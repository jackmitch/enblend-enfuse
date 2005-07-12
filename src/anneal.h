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

#include <algorithm>
#include <vector>
#include <boost/random.hpp>
#include <math.h>

#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"

using std::pair;
using std::vector;

using boost::uniform_01;
//using boost::uniform_int;
//using boost::detail::uniform_smallint_integer;
using boost::variate_generator;
using boost::mt19937;

using vigra::LineIterator;
using vigra::Point2D;

namespace enblend {

template<class UniformRandomNumberGenerator, class IntType>
struct uniform_smallint_integer
{
public:
  typedef UniformRandomNumberGenerator base_type;
  typedef IntType result_type;

  uniform_smallint_integer(base_type & rng, IntType min, IntType max)
    : _rng(&rng)
  { set(min, max); }

  void set(result_type min, result_type max);
  
  result_type min() const { return _min; }
  result_type max() const { return _max; }
  base_type& base() const { return *_rng; }
  
  result_type operator()()
  {
    // we must not use the low bits here, because LCGs get very bad then
    //return (((*_rng)() - _rng->min()) / _factor) % _range + _min;
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

//FIXME turn this into a valid functor for random_shuffle
//template <class IntType=int>
//class STLCompatibleUniformSmallintFunctor {
//public:
//    inline IntType operator()(IntType maxNonInclusive) const {
//        uniform_smallint<IntType> us(0, maxNonInclusive - 1);
//        return us(Twister);
//    }
//};

template <typename DistanceCostImage, typename StitchCostImage>
class AnnealConfiguration {
public:
    typedef DistanceCostImage DistanceCostImageType;
    typedef StitchCostImage StitchCostImageType;

    AnnealConfiguration(const DistanceCostImage *d,
            const StitchCostImage *s,
            const vector<pair<bool, Point2D> > &v)
    : snake(v), origSnake(v), cost(0.0), segmentCost(v.size(), 0.0), moveablePointIndices(),
      dci(d), sci(s), maxStepSize(10),
      w(d->width()), h(d->height()),
      movingPrevIndex(0), movingIndex(0), movingPoint(),
      newSegmentCost(0.0), newPrevSegmentCost(0.0), newCost(0.0),
      moveablePointRNG(Twister, 0, 10), stepRNG(Twister, 0, 10) {

        // Calculate the initial snake cost
        // Populate the moveablePointIndices vector
        double c = 0.0;
        for (unsigned int i = 0; i < snake.size(); i++) {
            if (!snake[i].first) moveablePointIndices.push_back(i);
            double d = evalSegmentCost(i);
            segmentCost[i] = d;
            c += d;
        }
        cost = c;

        moveablePointRNG.set(0, moveablePointIndices.size() - 1);
        //cout << "moveablePoints=" << moveablePointIndices.size() << endl;
        //// Create the moveablePointIndex generators
        //range = new uniform_smallint<unsigned int>(0, moveablePointIndices.size() - 1);
        //rng = new variate_generator<mt19937&, uniform_smallint<unsigned int> > rng(Twister, range);
    }

    //AnnealConfiguration(const AnnealConfiguration &c)
    //: cost(c.cost), snake(c.snake), segmentCost(c.segmentCost),
    //  moveablePointIndices(c.moveablePointIndices), dci(c.dci), sci(c.sci) { }

    ~AnnealConfiguration() { }

    //AnnealConfiguration& operator=(const AnnealConfiguration& c) {
    //    cost = c.cost;
    //    snake.clear();
    //    snake.insert(snake.end(), c.snake.begin(), c.snake.end());
    //    segmentCost.clear();
    //    segmentCost.insert(segmentCost.end(), c.segmentCost.begin(), c.segmentCost.end());
    //    moveablePointIndices.clear();
    //    moveablePointIndices.insert(moveablePointIndices.end(), c.moveablePointIndices.begin(), c.moveablePointIndices.end());
    //    dci = c.dci;
    //    sci = c.sci;
    //    return *this;
    //}

    const vector<pair<bool, Point2D> >& getSnake() const {
        return snake;
    }

    double getCost() const {
        return cost;
    }

    void setMaxStepSize(unsigned int maxStepSize) {
        this->maxStepSize = maxStepSize;
    }

    // Propose a move, return the proposed new cost
    double step();

    void commit() {
        cost = newCost;
        segmentCost[movingPrevIndex] = newPrevSegmentCost;
        segmentCost[movingIndex] = newSegmentCost;
        snake[movingIndex].second = movingPoint;
    }

    void print() const;

protected:

    // Calculate cost of segment between pointIndex and pointIndex+1
    double evalSegmentCost(const unsigned int pointIndex) const;
    double evalSegmentCost(const unsigned int pointIndex, const Point2D &pointB) const;
    double evalSegmentCost(const unsigned int movingPointIndex, const Point2D &pointA, const unsigned int pointIndex) const;
    double evalSegmentCost(const unsigned int movingPointIndex, const Point2D &pointA, const Point2D &pointB) const;

    // The snake itself. First is moveable flag, second is point
    vector<pair<bool, Point2D> > snake;
    // The original snake point locations.
    vector<pair<bool, Point2D> > origSnake;

    // Total snake cost and cost of each segment.
    // segmentCost includes non-moveable segments.
    double cost;
    vector<double> segmentCost;

    // snake indices that are moveable.
    vector<unsigned int> moveablePointIndices;

    // Cost function inputs
    const DistanceCostImage *dci;
    const StitchCostImage *sci;

    int maxStepSize;
    int w;
    int h;

    // step function caches a proposed move here.
    unsigned int movingPrevIndex;
    unsigned int movingIndex;
    Point2D movingPoint;
    double newSegmentCost;
    double newPrevSegmentCost;
    double newCost;

    // Variate generator that picks a random moveable vertex.
    uniform_smallint_integer<mt19937, unsigned int> moveablePointRNG;
    uniform_smallint_integer<mt19937, int> stepRNG;

};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::evalSegmentCost(const unsigned int pointIndex) const {
    unsigned int nextPointIndex = (pointIndex + 1) % snake.size();

    // If segment is between two frozen points, it has zero cost.
    if (snake[pointIndex].first && snake[nextPointIndex].first) return 0.0;

    const Point2D &pointA = snake[pointIndex].second;
    const Point2D &pointB = snake[nextPointIndex].second;
    return evalSegmentCost(pointIndex, pointA, pointB);
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::evalSegmentCost(const unsigned int pointIndex, const Point2D &pointB) const {
    return evalSegmentCost(pointIndex, snake[pointIndex].second, pointB);
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::evalSegmentCost(const unsigned int movingPointIndex, const Point2D &pointA, const unsigned int pointIndex) const {
    return evalSegmentCost(movingPointIndex, pointA, snake[pointIndex].second);
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::evalSegmentCost(const unsigned int movingPointIndex, const Point2D &pointA, const Point2D &pointB) const {
    typedef typename DistanceCostImage::const_traverser DistanceIterator;
    typedef typename StitchCostImage::const_traverser StitchIterator;

    // Total cost of snake according to dci
    double distanceCost = 0.0;
    // Total cost of snake according to sci
    double stitchCost = 0.0;

    LineIterator<StitchIterator> lineS(sci->upperLeft() + pointA,
            sci->upperLeft() + pointB);
    LineIterator<StitchIterator> lineSEnd(sci->upperLeft() + pointB,
            sci->upperLeft() + pointB);

    unsigned int lineLength = 0;
    do {
        stitchCost += 256.0 * log(0.2 * *lineS + 1.0) / log(52.0);
        ++lineS;
        ++lineLength;
    } while (lineS != lineSEnd);

    // Measure how far away the pointA is from its original location.
    const Point2D &movingPointOriginalLocation = origSnake[movingPointIndex].second;
    Diff2D delta = pointA - movingPointOriginalLocation;
    double distance = delta.magnitude();

    // Measure how far away from a feature the original location of pointA was
    double originalFeatureDistance = sqrt((*dci)[movingPointOriginalLocation]);

    // It is costly for points to stray too far from their original locations,
    // depending on how near the center they are.
    if (distance > (originalFeatureDistance / 2)) {
        distanceCost = 2.0 * distance;
    }

    //cout << "pointA=" << pointA << " originalLoc=" << movingPointOriginalLocation << " distanceCost=" << distanceCost << " origFeatureDistance=" << originalFeatureDistance << endl;

    // Weighting
    double weightedCost = stitchCost + distanceCost;
    
    // Stop clumping by penalizing short lines.
    if (lineLength < 5) weightedCost *= 2.0;

    // Also penalize long lines
    if (lineLength > 30) weightedCost *= 2.0;

    return weightedCost;
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::step() {

    // Choose a moveable point randomly.
    movingIndex = moveablePointIndices[moveablePointRNG()];
    movingPrevIndex = (movingIndex + snake.size() - 1) % snake.size();
    movingPoint = snake[movingIndex].second;

    unsigned int nextPointIndex = (movingIndex + 1) % snake.size();

    // Limit the size of the neighborhood if we are close to the image edge.
    int minDeltaX = std::max((1 - movingPoint.x), -maxStepSize);
    int maxDeltaX = std::min((w - 2 - movingPoint.x), maxStepSize);
    int minDeltaY = std::max((1 - movingPoint.y), -maxStepSize);
    int maxDeltaY = std::min((h - 2 - movingPoint.y), maxStepSize);

    stepRNG.set(minDeltaX, maxDeltaX);
    int deltaX = stepRNG();

    stepRNG.set(minDeltaY, maxDeltaY);
    int deltaY = stepRNG();
    movingPoint += Diff2D(deltaX, deltaY);

    newPrevSegmentCost = evalSegmentCost(movingPrevIndex, movingPoint);
    newSegmentCost = evalSegmentCost(movingIndex, movingPoint, nextPointIndex); 

    newCost = cost - segmentCost[movingPrevIndex] - segmentCost[movingIndex]
            + newPrevSegmentCost + newSegmentCost;

    return newCost;
/*
    Point2D originalPoint = snake[index].second;
    double bestCostSoFar = NumericTraits<double>::max();
    Point2D bestMoveSoFar = snake[index].second;
    double prevSegmentNewCost = segmentCost[prev];
    double nextSegmentNewCost = segmentCost[index]; 

    //cout << "step considering moving point (" << originalPoint.x << ", " << originalPoint.y << ")" << endl;
    //cout << "original cost=" << cost << endl;
    //cout << "original segmentA=" << segmentCost[prev] << endl;
    //cout << "original segmentB=" << segmentCost[index] << endl;

    // Limit the size of the neighborhood if we are close to the image edge.
    int minDeltaX = std::max((1 - originalPoint.x), -maxStepSize);
    int maxDeltaX = std::min((dci->width() - 2 - originalPoint.x), maxStepSize);
    int minDeltaY = std::max((1 - originalPoint.y), -maxStepSize);
    int maxDeltaY = std::min((dci->height() - 2 - originalPoint.y), maxStepSize);

    // Find the best move in the neighborhood.
    // visit neighborhood moves randomly - this leads to a truly random choice if there are
    // several possibilities with the same cost.
    // FIXME this sucks use random_shuffle
    // FIXME consider minStepSize
    vector<int> permutedXMoves;
    for (int deltaX = minDeltaX; deltaX <= maxDeltaX; deltaX++) {
        vector<int>::iterator pos = permutedXMoves.begin();
        uniform_int<unsigned int> insertRange(0, permutedXMoves.size());
        variate_generator<mt19937&, uniform_int<unsigned int> > insertRNG(Twister, insertRange);
        pos += insertRNG();
        permutedXMoves.insert(pos, deltaX);
    }

    for (unsigned int deltaXIndex = 0; deltaXIndex < permutedXMoves.size(); deltaXIndex++) {
        int deltaX = permutedXMoves[deltaXIndex];

        vector<int> permutedYMoves;
        // FIXME this sucks use random_shuffle
        // FIXME consider minStepSize
        for (int deltaY = minDeltaY; deltaY <= maxDeltaY; deltaY++) {
            vector<int>::iterator pos = permutedYMoves.begin();
            uniform_int<unsigned int> insertRange(0, permutedYMoves.size());
            variate_generator<mt19937&, uniform_int<unsigned int> > insertRNG(Twister, insertRange);
            pos += insertRNG();
            permutedYMoves.insert(pos, deltaY);
        }

        for (unsigned int deltaYIndex = 0; deltaYIndex < permutedYMoves.size(); deltaYIndex++) {
            int deltaY = permutedYMoves[deltaYIndex];

            // This is not really a move - skip it.
            if (deltaX == 0 && deltaY == 0) continue;

            // Check for minStepSize (chessboard metric)
            if ((abs(deltaX) + abs(deltaY)) < minStepSize) continue;

            snake[index].second = originalPoint + Diff2D(deltaX, deltaY);
            double psCost = evalSegmentCost(prev);
            double nsCost = evalSegmentCost(index);
            double trialCost = psCost + nsCost;

            //cout << "move to (" << deltaX << ", " << deltaY << ") has cost " << trialCost << endl;

            if (trialCost < bestCostSoFar) {
                bestMoveSoFar = snake[index].second;
                bestCostSoFar = trialCost;
                prevSegmentNewCost = psCost;
                nextSegmentNewCost = nsCost;
            }
        }
    }

    // Update configuration state.
    snake[index].second = bestMoveSoFar;
*/
    //cost = cost - segmentCost[prev] - segmentCost[index]
    //        + prevSegmentNewCost + nextSegmentNewCost;
    //segmentCost[prev] = prevSegmentNewCost;
    //segmentCost[index] = nextSegmentNewCost;

    //cout << "moving point to (" << bestMoveSoFar.x << ", " << bestMoveSoFar.y << ")" << endl;
    //cout << "new cost=" << cost << endl;
    //cout << "new segmentA=" << segmentCost[prev] << endl;
    //cout << "new segmentB=" << segmentCost[index] << endl;
    //return index;
};

template <typename DistanceCostImage, typename StitchCostImage>
void AnnealConfiguration<DistanceCostImage, StitchCostImage>::print() const {
    for (unsigned int i = 0; i < snake.size(); i++) {
        pair<bool, Point2D> p = snake[i];
        if (!p.first) {
            cout << " (" << p.second.x << ", " << p.second.y << ")";
        }
    }
    return;
};

template <typename DistanceCostImage, typename StitchCostImage>
void annealSnake(const DistanceCostImage* const dci,
        const StitchCostImage* const sci,
        vector<pair<bool, Point2D> > *snake) {

    typedef AnnealConfiguration<DistanceCostImage, StitchCostImage> AnnealConfigurationType;
    AnnealConfigurationType cfg(dci, sci, *snake);

    uniform_01<mt19937> thermal(Twister);

    unsigned int stepsPerIteration = 256 * snake->size();
    //unsigned int stepsPerIteration = 10;
    cout << "initial cost=" << cfg.getCost() << endl;
    cout << "stepsPerIteration=" << stepsPerIteration << endl;

    double initialTemperature = 1.0;
    double temperatureCutoff = 0.0005;
    double temperatureDampening = 0.85;
    double temperature = initialTemperature;

    while (temperature > temperatureCutoff) {
        unsigned int stepSize = 5 + (unsigned int)(temperature * 70);
        cfg.setMaxStepSize(stepSize);
        unsigned int acceptedSteps = 0;
        double acceptCoefficient = 0.05 + (0.15 * temperature);
        for (unsigned int i = 0; i < stepsPerIteration; i++) {
            double prevCost = cfg.getCost();
            double newCost = cfg.step();
            double percentDifference = (newCost - prevCost) / prevCost;
            double acceptProbability = acceptCoefficient * exp(-5.0 * percentDifference / temperature);
            if ((newCost < prevCost) || (thermal() < acceptProbability)) {
                cfg.commit();
                acceptedSteps++;
                //cout << "i=" << i << " prevCost=" << prevCost << " newCost=" << cfg.getCost() << endl;
            }
        }
        cout << "temp=" << temperature
             << " stepSize=" << stepSize
             << " accepted steps=" << acceptedSteps
             << " accept ratio=" << ((double)acceptedSteps / (double)stepsPerIteration)
             << " cost=" << cfg.getCost()
             << endl;
        temperature *= temperatureDampening;
    }

    cout << "final cost=" << cfg.getCost() << endl;

    //vector<int> trials(snake->size(), 0);

/*
    int stepsPerTemperature = 30000;

    double initialTemperature = 0.30;
    double coolingFactor = 0.90;
    double finalTemperature = 0.05;
    int stepSize = 15;
    int initialMinStepSize = 13;
    double minStepSizeCoolingFactor = 0.80;

    uniform_01<mt19937> rng(Twister);

    // The overall best solution found
    // This starts out as the initial configuration.
    AnnealConfigurationType bestConfig(dci, sci, *snake);
    double bestCost = bestConfig.getCost();

//for (unsigned int i = 0; i < 25; i++) {
//    bestConfig.step(stepSize);
//    cout << "prevCost=" << bestCost << " newCost=" << bestConfig.getCost() << endl;
//    bestCost = bestConfig.getCost();
//}

    // The hill-climbing configuration
    // This starts out as the initial configuration.
    AnnealConfigurationType currentConfig(bestConfig);
    double currentCost = currentConfig.getCost();

    double currentTemperature = initialTemperature;
    double minStepSize = initialMinStepSize;
    while (currentTemperature > finalTemperature) {

        for (int i = 0; i < stepsPerTemperature; i++) {
            AnnealConfigurationType testConfig(currentConfig);
            int indexMoved = testConfig.step(stepSize, (int)floor(minStepSize));
            //trials[indexMoved]++;
            double testCost = testConfig.getCost();

            if (testCost <= bestCost) {
                bestConfig = testConfig;
                bestCost = testCost;
            }

            double percentDifference = (testCost - currentCost) / currentCost;
            double acceptProbability = exp(-1.0 * percentDifference / currentTemperature);
            bool accept = (rng() < acceptProbability);

            if ((testCost < currentCost) || accept) {
                currentConfig = testConfig;
                currentCost = testCost;
            }
        }

        cout << "temp=" << currentTemperature << " minStepSize=" << (int)floor(minStepSize) << " bestCost=" << bestCost << " currentCost=" << currentCost << endl;
        currentTemperature *= coolingFactor;
        minStepSize *= minStepSizeCoolingFactor;
    }
*/
    // Copy result to input snake.
    snake->clear();
    snake->insert(snake->end(), cfg.getSnake().begin(), cfg.getSnake().end());

    //cout << "move statistics:" << endl;
    //for (unsigned int i = 0; i < snake->size(); i++) {
    //    cout << "index " << i << " frozen=" << (*snake)[i].first << " moved=" << trials[i] << endl;
    //}

};

} // namespace enblend

#endif /* __ANNEAL_H__ */
