/*
 * Copyright (C) 2004 Andrew Mihal
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef __ANNEAL_H__
#define __ANNEAL_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ext/slist>
#include <vector>
#include <boost/random.hpp>
#include <math.h>

#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"

using std::pair;
using std::vector;
using __gnu_cxx::slist;

using boost::uniform_01;
using boost::variate_generator;
using boost::mt19937;

using vigra::LineIterator;
using vigra::Point2D;
using vigra::Rect2D;

namespace enblend {

template<class UniformRandomNumberGenerator, class IntType>
struct uniform_smallint_integer
{
public:
    typedef UniformRandomNumberGenerator base_type;
    typedef IntType result_type;

    uniform_smallint_integer(base_type & rng, IntType min, IntType max) : _rng(&rng)
    { set(min, max); }

    void set(result_type min, result_type max);

    result_type min() const { return _min; }
    result_type max() const { return _max; }
    base_type& base() const { return *_rng; }

    result_type operator()() {
        // we must not use the low bits here, because LCGs get very bad then
        return (IntType)((((*_rng)() - _baseMin) / _factor) % _range) + _min;
    }

private:
    typedef typename base_type::result_type base_result;
    base_type * _rng;
    IntType _min, _max;
    base_result _range;
    base_result _baseMin;
    //int _factor;
    base_result _factor;
};

template<class UniformRandomNumberGenerator, class IntType>
void uniform_smallint_integer<UniformRandomNumberGenerator, IntType>::
set(result_type min, result_type max) 
{
    // skip all this if min and max are unchanged.
    if (min == _min && _max == max) return;
    _min = min;
    _max = max;
    assert(min < max);

    _range = static_cast<base_result>(_max-_min)+1;
    // FIXME there was a bug here where the variable would get shadowed
    _factor = 1;
  
    // LCGs get bad when only taking the low bits.
    // (probably put this logic into a partial template specialization)
    // Check how many low bits we can ignore before we get too much
    // quantization error.
    _baseMin = _rng->min();
    base_result r_base = _rng->max() - _baseMin;
    if(r_base == std::numeric_limits<base_result>::max()) {
        _factor = 2;
        r_base /= 2;
    }
    r_base += 1;
    if(r_base % _range == 0) {
        // No quantization effects, good
        _factor = r_base / _range;
    } else {
        // carefully avoid overflow; pessimizing heree
        for( ; r_base/_range/32 >= _range; _factor *= 2)
            r_base /= 2;
    }
    //cout << "min=" << _min << " max=" << _max << " range=" << _range << " factor=" << _factor << endl;
}

template <typename CostImage>
class AnnealConfiguration {
public:

    AnnealConfiguration(const CostImage* const d, slist<pair<bool, Point2D> > *v) :
            costImage(d),
            totalCost(0.0),
            maxStepSize(20),
            w(d->width()),
            h(d->height()),
            cacheMoveIndex(0),
            cacheMovePoint(0,0),
            cacheSegmentCost(0.0),
            cacheNextSegmentCost(0.0),
            cacheDeltaCost(0.0),
            moveablePointRNG(Twister, 0, 10),
            stepRNG(Twister, 0, 10) {

        unsigned int pointIndex = 0;
        slist<pair<bool, Point2D> >::iterator lastPoint = v->previous(v->end());
        for (slist<pair<bool, Point2D> >::iterator point = v->begin();
                point != v->end(); ++point, ++pointIndex) {

            originalPoints.push_back(point->second);
            currentPoints.push_back(point->second);

            if (point->first) {
                moveablePointIndices.push_back(pointIndex);
            }

            if (point->first || lastPoint->first) {
                // Segment between lastPoint and point is moveable.
                double segmentCost = evalSegmentCost(point->second, lastPoint->second);
                totalCost += segmentCost;
                segmentCosts.push_back(segmentCost);
            }
            else {
                segmentCosts.push_back(0.0);
            }

            lastPoint = point;
        }

        moveablePointRNG.set(0, moveablePointIndices.size() - 1);
    }

    double getCost() const {
        return totalCost;
    }

    int numMoveableVertices() const {
        return moveablePointIndices.size();
    }

    vector<Point2D>& getCurrentPoints() {
        return currentPoints;
    }

    void setMaxStepSize(unsigned int maxStepSize) {
        this->maxStepSize = maxStepSize;
    }

    // Propose a move, return the percent improvement
    double step();

    void commit() {
        currentPoints[cacheMoveIndex] = cacheMovePoint;
        segmentCosts[cacheMoveIndex] = cacheSegmentCost;
        segmentCosts[cacheMoveIndex+1] = cacheNextSegmentCost;
        totalCost += cacheDeltaCost;
    }

protected:

    const CostImage *costImage;

    // Cost of the entire snake
    double totalCost;

    // original locations of points
    vector<Point2D> originalPoints;

    // Current positions of all points
    vector<Point2D> currentPoints;

    // Indices of movable points in currentPoints and originalPoints
    vector<unsigned int> moveablePointIndices;

    // Moveable segment costs. Index is for segment behind corresponding currentPoint.
    vector<double> segmentCosts;

    int maxStepSize;

    int w;
    int h;

    // Calculate cost of segment between pointA and pointB
    double evalSegmentCost(const Point2D &pointA, const Point2D &pointB) const;

    // step function caches a proposed move here.
    unsigned int cacheMoveIndex;
    Point2D cacheMovePoint;
    double cacheSegmentCost;
    double cacheNextSegmentCost;
    double cacheDeltaCost;

    // Variate generator that picks a random moveable vertex.
    uniform_smallint_integer<mt19937, unsigned int> moveablePointRNG;
    uniform_smallint_integer<mt19937, int> stepRNG;

};

template <typename CostImage>
double AnnealConfiguration<CostImage>::evalSegmentCost(const Point2D& pointA, const Point2D &pointB) const {
    typedef typename CostImage::const_traverser CostIterator;
    typedef typename CostImage::PixelType CostType;

    double cost = 0.0;
    int lineLength = 0;

    if (pointA != pointB) {
        LineIterator<CostIterator> lineStart(costImage->upperLeft() + pointA,
                                             costImage->upperLeft() + pointB);
        LineIterator<CostIterator> lineEnd(costImage->upperLeft() + pointB,
                                           costImage->upperLeft() + pointB);

        do {
            CostType pointCost = *lineStart;
            if (pointCost == NumericTraits<CostType>::max()) cost += 20 * pointCost;
            else cost += 2 * pointCost;
            ++lineLength;
            ++lineStart;
        } while (lineStart != lineEnd);
    }

    if (lineLength < 30) return cost;
    if (lineLength < 50) return cost + (1 << (lineLength-30));
    else return cost + (1 << 20);
    //return cost + (1 << (std::min(20, lineLength-20)));
};

template <typename CostImage>
double AnnealConfiguration<CostImage>::step() {

    // Choose a moveable point randomly.
    cacheMoveIndex = moveablePointIndices[moveablePointRNG()];
    cacheMovePoint = currentPoints[cacheMoveIndex];

    // Limit the size of the neighborhood if we are close to the image edge.
    int minDeltaX = std::max((1 - cacheMovePoint.x), -maxStepSize);
    int maxDeltaX = std::min((w - 2 - cacheMovePoint.x), maxStepSize);
    int minDeltaY = std::max((1 - cacheMovePoint.y), -maxStepSize);
    int maxDeltaY = std::min((h - 2 - cacheMovePoint.y), maxStepSize);

    stepRNG.set(minDeltaX, maxDeltaX);
    int deltaX = stepRNG();

    stepRNG.set(minDeltaY, maxDeltaY);
    int deltaY = stepRNG();

    cacheMovePoint += Diff2D(deltaX, deltaY);

    unsigned int previousPointIndex = (cacheMoveIndex == 0) ? (currentPoints.size() - 1) : (cacheMoveIndex - 1);
    unsigned int nextPointIndex = (cacheMoveIndex + 1) % currentPoints.size();

    cacheSegmentCost = evalSegmentCost(cacheMovePoint, currentPoints[previousPointIndex]);
    cacheNextSegmentCost = evalSegmentCost(currentPoints[nextPointIndex], cacheMovePoint);

    double deltaCost = (cacheSegmentCost - segmentCosts[cacheMoveIndex])
                     + (cacheNextSegmentCost - segmentCosts[nextPointIndex]);

    // Measure distance between newLocation and originalPoints[pointIndex];
    Point2D currentLocation = currentPoints[cacheMoveIndex];
    Point2D originalLocation = originalPoints[cacheMoveIndex];

    double distanceDeltaCost = (cacheMovePoint - originalLocation).magnitude()
                             - (currentLocation - originalLocation).magnitude();

    cacheDeltaCost = deltaCost + 1.25 * distanceDeltaCost;

    return (cacheDeltaCost / totalCost);
};

template <typename CostImage>
class GDAConfiguration {
public:
    typedef typename CostImage::PixelType CostImagePixelType;
    typedef typename CostImage::const_traverser CostIterator;

    GDAConfiguration(/*const*/ CostImage* const d, slist<pair<bool, Point2D> > *v) : costImage(d) {

        costImageBounds = Rect2D(d->size());
        kMax = 1;

        // Copy original point locations into originalPoints and mfEstimates
        slist<pair<bool, Point2D> >::iterator lastPoint = v->previous(v->end());
        for (slist<pair<bool, Point2D> >::iterator currentPoint = v->begin(); currentPoint != v->end(); ) {
            convergedPoints.push_back(false);
            originalPoints.push_back(currentPoint->second);
            mfEstimates.push_back(currentPoint->second);

            if (!currentPoint->first) {
                // Point is not moveable.
                vector<Point2D> *stateSpace = new vector<Point2D>();
                stateSpace->push_back(currentPoint->second);
                pointStateSpaces.push_back(stateSpace);

                vector<double> *stateProbabilities = new vector<double>();
                stateProbabilities->push_back(1.0);
                pointStateProbabilities.push_back(stateProbabilities);

                lastPoint = currentPoint;
                ++currentPoint;
            }
            else {
                Point2D lastPoint2D = lastPoint->second;
                Point2D currentPoint2D = currentPoint->second;
                lastPoint = currentPoint;
                ++currentPoint;
                Point2D nextPoint2D = (currentPoint == v->end()) ? v->begin()->second : currentPoint->second;

                // Determine state space of currentPoint along normal vector
                Diff2D normal(lastPoint2D.y - nextPoint2D.y, nextPoint2D.x - lastPoint2D.x);
                normal *= std::min(costImageBounds.width(), costImageBounds.height()) / (3 * normal.magnitude());

                LineIterator<Diff2D> lineBegin(Diff2D(currentPoint2D) + normal, Diff2D(currentPoint2D) - normal);
                Diff2D lineEnd = Diff2D(currentPoint2D) - normal;

                vector<Point2D> *stateSpace = new vector<Point2D>();
                pointStateSpaces.push_back(stateSpace);
                while (*lineBegin != lineEnd) {
                    if (costImageBounds.contains(Point2D(*lineBegin))) {
                        if ((*costImage)[*lineBegin] != NumericTraits<CostImagePixelType>::max()) {
                            stateSpace->push_back(Point2D(*lineBegin));
                            //(*costImage)[*lineBegin] = 150;
                        }
                    }
                    ++lineBegin;
                    if (*lineBegin == lineEnd) break;
                    ++lineBegin;
                    if (*lineBegin == lineEnd) break;
                    ++lineBegin;
                }

                vector<double> *stateProbabilities = new vector<double>(stateSpace->size(), 1.0 / stateSpace->size());
                pointStateProbabilities.push_back(stateProbabilities);

                kMax = std::max(kMax, stateSpace->size());
            }
        }

        tau = 0.85;

        // Maximum cost change possible by any single annealing move
        deltaEMax = 500.0;
        deltaEMin = 1.0;
        double epsilon = 1.0 / (kMax * kMax);

        tInitial = ceil(deltaEMax / log((kMax - 1 + (kMax * kMax * epsilon)) / (kMax - 1 - (kMax * kMax * epsilon))));
        tFinal = deltaEMin / log((kMax - (kMax * epsilon) - 1) / (kMax * epsilon));
        cout << "tInitial=" << tInitial << endl;
        cout << "tFinal=" << tFinal << endl;
    }

    ~GDAConfiguration() {
        for (vector<vector<Point2D>* >::iterator i = pointStateSpaces.begin();
                i != pointStateSpaces.end();
                ++i) {
            delete *i;
        }

        for (vector<vector<double>* >::iterator i = pointStateProbabilities.begin();
                i != pointStateProbabilities.end();
                ++i) {
            delete *i;
        }
    }

    void run() {
        tCurrent = tInitial;
        while (tCurrent > tFinal) {
            double epsilon = 1.0 / kMax;
            eta = (int)ceil(log(epsilon) / log(((kMax - 2.0) / (2.0 * kMax) * exp(-tCurrent / deltaEMax)) + 0.5));
            cout << "tCurrent=" << tCurrent << " eta=" << eta << endl;
            for (unsigned int i = 0; i < eta; i++) {
                //cout << "i = " << i << endl;
                iterate();
            }
            tCurrent *= tau;
        }

        int convergeCount = 0;
        for (unsigned int i = 0; i < convergedPoints.size(); i++) {
            if (convergedPoints[i]) convergeCount++;
        }
        cout << "Total points=" << convergedPoints.size() << " converged=" << convergeCount << endl;
    }

    vector<Point2D> & getCurrentPoints() { return mfEstimates; }

    double originalCost() { return 0.0; }
    double currentCost() { return 0.0; }

protected:

    void iterate() {
        double E[kMax];
        double pi[kMax];

        unsigned int lastIndex = originalPoints.size() - 1;
        for (unsigned int index = 0; index < originalPoints.size(); ++index) {
            unsigned int nextIndex = (index + 1) % originalPoints.size();
            //cout << "lastIndex=" << lastIndex << " index=" << index << " nextIndex=" << nextIndex << " size=" << originalPoints.size() << endl;

            Point2D originalPoint = originalPoints[index];
            //Point2D currentPointEstimate = mfEstimates[index];
            Point2D lastPointEstimate = mfEstimates[lastIndex];
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            //cout << "op=" << originalPoint << " cp=" << currentPointEstimate << " lpe=" << lastPointEstimate << " npe=" << nextPointEstimate << endl;
            //bool lastPointInCostImage = costImageBounds.contains(lastPointEstimate);
            //bool nextPointInCostImage = costImageBounds.contains(nextPointEstimate);

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            unsigned int localK = stateSpace->size();

            lastIndex = index;
            if (localK == 1 || convergedPoints[index]) continue;

            // Calculate E values
            for (unsigned int i = 0; i < localK; ++i) {
                Point2D currentPoint = (*stateSpace)[i];
                E[i] = costImageCost(lastPointEstimate, currentPoint)
                        + costImageCost(currentPoint, nextPointEstimate)
                        + (currentPoint - originalPoint).magnitude();
            }

            for (unsigned int j = 0; j < localK; ++j) pi[j] = 0.0;

            // Calculate new stateProbabilities
            for (unsigned int j = 0; j < localK; ++j) {
                for (unsigned int i = j; i < localK; ++i) {
                    double piT = (*stateProbabilities)[i] + (*stateProbabilities)[j];
                    double An = 1.0 / (1.0 + exp((E[j] - E[i]) / tCurrent));
                    pi[j] += An * piT;
                    if (i != j) pi[i] += (1.0 - An) * piT;
                }
            }

            // Normalize
            for (unsigned int i = 0; i < localK; ++i) {
                (*stateProbabilities)[i] = pi[i] / localK;
            }

        }

        kMax = 1;
        // Make new mean field estimates.
        for (unsigned int index = 0; index < pointStateSpaces.size(); ++index) {
            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            double estimateX = 0.0;
            double estimateY = 0.0;
            unsigned int localK = stateSpace->size();
            if (localK == 1 || convergedPoints[index]) continue;
            for (unsigned int k = 0; k < localK; ++k) {
                double weight = (*stateProbabilities)[k];
                if (weight > 0.99) convergedPoints[index] = true;
                Point2D state = (*stateSpace)[k];
                //cout << "state=" << state << " weight=" << weight << endl;
                estimateX += weight * (double)state.x;
                estimateY += weight * (double)state.y;
            }
            mfEstimates[index] = Point2D((int)round(estimateX), (int)round(estimateY));
            //if (convergedPoints[index]) {
            //    cout << "point " << originalPoints[index] << " has converged to " << mfEstimates[index] << endl;
            //    for (unsigned int k = 0; k < localK; ++k) {
            //        cout << "    " << (*stateSpace)[k] << " = " << (*stateProbabilities)[k] << endl;
            //    }
            //}

            //cout << "op=" << originalPoints[index] << " localK=" << localK << " mfe=(" << estimateX << ", " << estimateY << ")" << "=" << mfEstimates[index] << endl;

            // Remove improbable solutions from the search space
            for (unsigned int k = 0; k < stateProbabilities->size(); ) {
                double weight = (*stateProbabilities)[k];
                if (weight < 0.0001) {
                    (*stateProbabilities)[k] = (*stateProbabilities)[stateProbabilities->size() - 1];
                    (*stateSpace)[k] = (*stateSpace)[stateSpace->size() - 1];
                    stateProbabilities->pop_back();
                    stateSpace->pop_back();
                } else {
                    ++k;
                }
            }

            kMax = std::max(kMax, stateProbabilities->size());

            // FIXME ensure new mfEstimate is inside costImage
        }

    }

    double costImageCost(const Point2D &start, const Point2D &end) {
        //cout << "costImageCost(" << start << ", " << end << ")" << endl;
        double cost = 0.0;
        int lineLength = 0;

        if (start != end) {
            LineIterator<CostIterator> lineStart(costImage->upperLeft() + start,
                                                 costImage->upperLeft() + end);
            LineIterator<CostIterator> lineEnd(costImage->upperLeft() + end,
                                               costImage->upperLeft() + end);

            do {
                CostImagePixelType pointCost = *lineStart;
                if (pointCost == NumericTraits<CostImagePixelType>::max()) cost += 20 * pointCost;
                else cost += 2 * pointCost;
                ++lineLength;
                ++lineStart;
            } while (lineStart != lineEnd);
        }

        if (lineLength < 30) return cost;
        if (lineLength < 50) return cost + (1 << (lineLength-30));
        else return cost + (1 << 20);
    }

    /*const*/ CostImage *costImage;
    Rect2D costImageBounds;

    vector<Point2D> originalPoints;
    vector<Point2D> mfEstimates;
    vector<vector<Point2D>* > pointStateSpaces;
    vector<vector<double>* > pointStateProbabilities;
    vector<bool> convergedPoints;

    // Initial Temperature
    double tInitial;

    // Final Temperature
    double tFinal;

    // Current Temperature
    double tCurrent;

    // Iterations per Temperature
    unsigned int eta;

    // Cooling constant
    double tau;

    double deltaEMax;
    double deltaEMin;
    unsigned int kMax;

};

template <typename CostImage>
void annealSnake(/*const*/ CostImage* const ci, slist<pair<bool, Point2D> > *snake) {

    GDAConfiguration<CostImage> cfg(ci, snake);
    cfg.run();

    slist<pair<bool, Point2D> >::iterator snakePoint = snake->begin();
    vector<Point2D>::iterator annealedPoint = cfg.getCurrentPoints().begin();
    for (; snakePoint != snake->end(); ++snakePoint, ++annealedPoint) {
        snakePoint->second = *annealedPoint;
    }

    //AnnealConfiguration<CostImage> cfg(ci, snake);

    //uniform_01<mt19937> thermal(Twister);

    //unsigned int stepsPerIteration = 256 * cfg.numMoveableVertices();

    //cout << "initial cost=" << cfg.getCost() << endl;
    //cout << "stepsPerIteration=" << stepsPerIteration << endl;

    //double initialTemperature = 1.0;
    //double temperatureCutoff = 0.01;
    //double temperatureDampening = 0.85;
    //double temperature = initialTemperature;

    //while (temperature > temperatureCutoff) {
    //    unsigned int stepSize = 3 + (unsigned int)(temperature * 20 /*.20 * std::min(ci->width(), ci->height())*/);
    //    cfg.setMaxStepSize(stepSize);
    //    unsigned int acceptedSteps = 0;
    //    double acceptCoefficient = 0.02 + (0.20 * temperature);
    //    for (unsigned int i = 0; i < stepsPerIteration; i++) {
    //        double percentDifference = cfg.step();
    //        double acceptProbability = acceptCoefficient * exp(-5.0 * percentDifference / temperature);
    //        if ((percentDifference < 0.0) || (thermal() < acceptProbability)) {
    //            cfg.commit();
    //            acceptedSteps++;
    //        }
    //    }
    //    cout << "temp=" << temperature
    //         << " stepSize=" << stepSize
    //         << " accepted steps=" << acceptedSteps
    //         << " accept ratio=" << ((double)acceptedSteps / (double)stepsPerIteration)
    //         << " cost=" << cfg.getCost()
    //         << endl;
    //    temperature *= temperatureDampening;
    //}

    //cout << "final cost=" << cfg.getCost() << endl;

    //slist<pair<bool, Point2D> >::iterator snakePoint = snake->begin();
    //vector<Point2D>::iterator annealedPoint = cfg.getCurrentPoints().begin();
    //for (; snakePoint != snake->end(); ++snakePoint, ++annealedPoint) {
    //    snakePoint->second = *annealedPoint;
    //}

};

} // namespace enblend

#endif /* __ANNEAL_H__ */
