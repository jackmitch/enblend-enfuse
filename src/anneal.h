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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>
#include <ext/slist>
#include <algorithm>
#include <vector>
#include <math.h>

#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"

using std::for_each;
using std::pair;
using std::vector;
using __gnu_cxx::slist;

using boost::lambda::bind;
using boost::lambda::_1;
using boost::lambda::delete_ptr;

using vigra::LineIterator;
using vigra::Point2D;
using vigra::Rect2D;

namespace enblend {

template <typename CostImage>
class GDAConfiguration {
public:
    typedef typename CostImage::PixelType CostImagePixelType;
    typedef typename CostImage::const_traverser CostIterator;

    GDAConfiguration(const CostImage* const d, slist<pair<bool, Point2D> > *v) : costImage(d) {

        kMax = 1;

        // Copy original point locations into originalPoints and mfEstimates
        slist<pair<bool, Point2D> >::iterator lastPoint = v->previous(v->end());
        for (slist<pair<bool, Point2D> >::iterator currentPoint = v->begin(); currentPoint != v->end(); ) {

            originalPoints.push_back(currentPoint->second);
            mfEstimates.push_back(currentPoint->second);

            vector<Point2D> *stateSpace = new vector<Point2D>();
            pointStateSpaces.push_back(stateSpace);

            vector<double> *stateProbabilities = new vector<double>();
            pointStateProbabilities.push_back(stateProbabilities);

            if (!currentPoint->first) {
                // Point is not moveable.
                stateSpace->push_back(currentPoint->second);
                stateProbabilities->push_back(1.0);

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
                normal *= std::min(costImage->width(), costImage->height()) / (3 * normal.magnitude());

                Diff2D lineEnd = Diff2D(currentPoint2D) - normal;
                LineIterator<Diff2D> lineBegin(Diff2D(currentPoint2D) + normal, lineEnd);

                while (lineBegin != lineEnd) {
                    if (costImage->isInside(*lineBegin)) {
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

                unsigned int localK = stateSpace->size();
                for (unsigned int i = 0; i < localK; ++i) stateProbabilities->push_back(1.0 / localK);

                kMax = std::max(kMax, localK);
            }

            convergedPoints.push_back(stateSpace->size() < 2);
        }

        tau = 0.75;
        deltaEMax = 300.0;
        deltaEMin = 5.0;
        double epsilon = 1.0 / (kMax * kMax);
        tInitial = ceil(deltaEMax / log((kMax - 1 + (kMax * kMax * epsilon)) / (kMax - 1 - (kMax * kMax * epsilon))));
        tFinal = deltaEMin / log((kMax - (kMax * epsilon) - 1) / (kMax * epsilon));
    }

    ~GDAConfiguration() {
        for_each(pointStateSpaces.begin(), pointStateSpaces.end(), bind(delete_ptr(), _1));
        for_each(pointStateProbabilities.begin(), pointStateProbabilities.end(), bind(delete_ptr(), _1));
    }

    void run() {
        tCurrent = tInitial;
        int numIterations = (int)ceil(log(tFinal/tInitial)/log(tau));
        while (tCurrent > tFinal) {
            double epsilon = 1.0 / kMax;
            unsigned int eta = (unsigned int)ceil(log(epsilon) / log(((kMax - 2.0) / (2.0 * kMax) * exp(-tCurrent / deltaEMax)) + 0.5));
            cout << "tCurrent=" << tCurrent << " eta=" << eta << endl;
            for (unsigned int i = 0; i < eta; i++) {
                iterate();
            }
            tCurrent *= tau;
        }
    }

    vector<Point2D> & getCurrentPoints() { return mfEstimates; }

    double currentCost() {
        double cost = 0.0;
        for (unsigned int index = 0; index < originalPoints.size(); ++index) {
            unsigned int nextIndex = (index + 1) % originalPoints.size();
            Point2D originalPoint = originalPoints[index];
            Point2D currentPointEstimate = mfEstimates[index];
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            if (costImage->isInside(currentPointEstimate) && costImage->isInside(nextPointEstimate)) {
                cost += costImageCost(currentPointEstimate, nextPointEstimate);
            }
            cost += (currentPointEstimate - originalPoint).magnitude() / 2.0;
        }
        return cost;
    }

protected:

    void iterate() {
        double E[kMax];
        double pi[kMax];

        unsigned int lastIndex = originalPoints.size() - 1;
        for (unsigned int index = 0; index < originalPoints.size(); ++index) {
            // Skip updating points that have already converged.
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            unsigned int localK = stateSpace->size();

            unsigned int nextIndex = (index + 1) % originalPoints.size();
            Point2D originalPoint = originalPoints[index];
            Point2D lastPointEstimate = mfEstimates[lastIndex];
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            lastIndex = index;

            // Calculate E values
            for (unsigned int i = 0; i < localK; ++i) {
                Point2D currentPoint = (*stateSpace)[i];
                E[i] = costImageCost(lastPointEstimate, currentPoint)
                        + costImageCost(currentPoint, nextPointEstimate)
                        + ((currentPoint - originalPoint).magnitude() / 2.0);
                pi[i] = 0.0;
            }

            // Calculate new stateProbabilities
            for (unsigned int j = 0; j < localK; ++j) {
                double piTj = (*stateProbabilities)[j];
                pi[j] += piTj;
                for (unsigned int i = (j+1); i < localK; ++i) {
                    double piT = (*stateProbabilities)[i] + piTj;
                    double An = 1.0 + exp((E[j] - E[i]) / tCurrent);
                    pi[j] += piT / An;
                    pi[i] += (1.0 - (1.0 / An)) * piT;
                }
                (*stateProbabilities)[j] = pi[j] / localK;
            }
        }

        kMax = 1;
        // Make new mean field estimates.
        for (unsigned int index = 0; index < pointStateSpaces.size(); ++index) {
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            unsigned int localK = stateSpace->size();
            double estimateX = 0.0;
            double estimateY = 0.0;

            for (unsigned int k = 0; k < localK; ++k) {
                double weight = (*stateProbabilities)[k];
                if (weight > 0.99) convergedPoints[index] = true;
                Point2D state = (*stateSpace)[k];
                estimateX += weight * (double)state.x;
                estimateY += weight * (double)state.y;
            }
            mfEstimates[index] = Point2D((int)round(estimateX), (int)round(estimateY));

            // Remove improbable solutions from the search space
            for (unsigned int k = 0; k < stateSpace->size(); ) {
                double weight = (*stateProbabilities)[k];
                if (weight < 0.0001) {
                    // Replace this state with last state
                    (*stateProbabilities)[k] = (*stateProbabilities)[stateProbabilities->size() - 1];
                    (*stateSpace)[k] = (*stateSpace)[stateSpace->size() - 1];

                    // Delete last state
                    stateProbabilities->pop_back();
                    stateSpace->pop_back();
                } else {
                    ++k;
                }
            }

            localK = stateSpace->size();
            if (localK < 2) convergedPoints[index] = true;
            kMax = std::max(kMax, stateProbabilities->size());

            // FIXME ensure new mfEstimate is inside costImage?
        }

    }

    int costImageCost(const Point2D &start, const Point2D &end) {
        int cost = 0;
        int lineLength = 0;

        if (start != end) {
            CostIterator endIterator = costImage->upperLeft() + end;
            LineIterator<CostIterator> lineStart(costImage->upperLeft() + start, endIterator);

            do {
                cost += *lineStart;
                ++lineLength;
                ++lineStart;
            } while (lineStart != endIterator);
        }

        if (lineLength < 8) cost += NumericTraits<CostImagePixelType>::max() * (8 - lineLength);

        return cost;
    }

    const CostImage *costImage;

    // Original point locations
    vector<Point2D> originalPoints;

    // Mean-field estimates of current point locations
    vector<Point2D> mfEstimates;

    // State spaces of each point
    vector<vector<Point2D>* > pointStateSpaces;

    // Probability vectors for each state space
    vector<vector<double>* > pointStateProbabilities;

    // Flags indicate which points have converged
    vector<bool> convergedPoints;

    // Initial Temperature
    double tInitial;

    // Final Temperature
    double tFinal;

    // Current Temperature
    double tCurrent;

    // Cooling constant
    double tau;

    // Maximum cost change possible by any single annealing move
    double deltaEMax;

    // Minimum cost change possible by any single annealing move
    double deltaEMin;

    // Largest state space over all points
    unsigned int kMax;

};

template <typename CostImage>
void annealSnake(const CostImage* const ci, slist<pair<bool, Point2D> > *snake) {

    GDAConfiguration<CostImage> cfg(ci, snake);
    //cout << "original cost = " << cfg.currentCost() << endl;
    cfg.run();
    //cout << "final cost = " << cfg.currentCost() << endl;

    slist<pair<bool, Point2D> >::iterator snakePoint = snake->begin();
    vector<Point2D>::iterator annealedPoint = cfg.getCurrentPoints().begin();
    for (; snakePoint != snake->end(); ++snakePoint, ++annealedPoint) {
        snakePoint->second = *annealedPoint;
    }

};

} // namespace enblend

#endif /* __ANNEAL_H__ */
