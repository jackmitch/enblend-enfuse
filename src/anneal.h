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

#include <vector>
#include <boost/random.hpp>
#include <math.h>

#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"

using std::pair;
using std::vector;

using boost::uniform_01;
using boost::uniform_int;
using boost::variate_generator;
using boost::mt19937;

using vigra::LineIterator;
using vigra::Point2D;

namespace enblend {

template <typename DistanceCostImage, typename StitchCostImage>
class AnnealConfiguration {
public:
    typedef DistanceCostImage DistanceCostImageType;
    typedef StitchCostImage StitchCostImageType;

    AnnealConfiguration(const DistanceCostImage *d,
            const StitchCostImage *s,
            const vector<pair<bool, Point2D> > &v)
    : cost(0.0), snake(v), segmentCost(v.size(), 0.0), moveablePointIndices(),
      dci(d), sci(s) {
        initSnakeCost();
    }

    AnnealConfiguration(const AnnealConfiguration &c)
    : cost(c.cost), snake(c.snake), segmentCost(c.segmentCost),
      moveablePointIndices(c.moveablePointIndices), dci(c.dci), sci(c.sci) { }

    ~AnnealConfiguration() { }

    AnnealConfiguration& operator=(const AnnealConfiguration& c) {
        cost = c.cost;
        snake.clear();
        snake.insert(snake.end(), c.snake.begin(), c.snake.end());
        segmentCost.clear();
        segmentCost.insert(segmentCost.end(), c.segmentCost.begin(), c.segmentCost.end());
        moveablePointIndices.clear();
        moveablePointIndices.insert(moveablePointIndices.end(), c.moveablePointIndices.begin(), c.moveablePointIndices.end());
        dci = c.dci;
        sci = c.sci;
        return *this;
    }

    const vector<pair<bool, Point2D> >& getSnake() const {
        return snake;
    }

    double getCost() const {
        return cost;
    }

    int step(int stepSize, int minStepSize);
    void print() const;

protected:

    // Calculate cost of entire snake
    void initSnakeCost();

    // Calculate cost of segment between pointIndex and pointIndex+1
    double evalSegmentCost(int pointIndex) const;

    // Total snake cost and cost of each segment.
    double cost;
    vector<pair<bool, Point2D> > snake;
    vector<double> segmentCost;
    vector<unsigned int> moveablePointIndices;
    const DistanceCostImage *dci;
    const StitchCostImage *sci;
};

template <typename DistanceCostImage, typename StitchCostImage>
void AnnealConfiguration<DistanceCostImage, StitchCostImage>::initSnakeCost() {
    double c = 0.0;
    for (unsigned int i = 0; i < snake.size(); i++) {
        if (!snake[i].first) moveablePointIndices.push_back(i);
        double d = evalSegmentCost(i);
        segmentCost[i] = d;
        c += d;
    }
    cost = c;
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::evalSegmentCost(int pointIndex) const {
    typedef typename DistanceCostImage::const_traverser DistanceIterator;
    typedef typename StitchCostImage::const_traverser StitchIterator;

    int nextPointIndex = (pointIndex + 1) % snake.size();
    pair<bool, Point2D> pointA = snake[pointIndex];
    pair<bool, Point2D> pointB = snake[nextPointIndex];

    // If segment is between two frozen points, it has zero cost.
    if (pointA.first && pointB.first) return 0.0;

    // Total cost of snake according to dci
    double distanceCost = 0.0;
    // Total cost of snake according to sci
    double stitchCost = 0.0;

    LineIterator<DistanceIterator> lineD(dci->upperLeft() + pointA.second,
            dci->upperLeft() + pointB.second);
    LineIterator<DistanceIterator> lineDEnd(dci->upperLeft() + pointB.second,
            dci->upperLeft() + pointB.second);
    LineIterator<StitchIterator> lineS(sci->upperLeft() + pointA.second,
            sci->upperLeft() + pointB.second);

        distanceCost += *lineD;
    int lineLength = 0;
    do {
        //distanceCost += *lineD;
        if (*lineS > (NumericTraits<typename StitchCostImage::value_type>::max() / 16))
            stitchCost += exp((6.0 * (*lineS / NumericTraits<typename StitchCostImage::value_type>::max())) - 1.0);
        ++lineD;
        ++lineS;
        ++lineLength;
    } while (lineD != lineDEnd);

    // Only lines exceeding a certain length are subject to penalty.
    if (lineLength < 3) lineLength = 15;
    else if (lineLength < 15) lineLength = 0;

    // Weighting
    return (1.0 * stitchCost) + ((double)lineLength / 5.0) /*+ ((double)distanceCost / (1<<19))*/ ;
};

template <typename DistanceCostImage, typename StitchCostImage>
int AnnealConfiguration<DistanceCostImage, StitchCostImage>::step(int maxStepSize, int minStepSize) {

    // Choose a moveable point randomly.
    uniform_int<unsigned int> range(0, moveablePointIndices.size() - 1);
    variate_generator<mt19937, uniform_int<unsigned int> > rng(Twister, range);
    unsigned int index = moveablePointIndices[rng()];

//vector<int> permutedMPIs;
//for (unsigned int rng = 0; rng < moveablePointIndices.size(); rng++) {
//    unsigned int index = moveablePointIndices[rng];
//    vector<int>::iterator pos = permutedMPIs.begin();
//    uniform_int<unsigned int> insertRange(0, permutedMPIs.size());
//    variate_generator<mt19937&, uniform_int<unsigned int> > insertRNG(Twister, insertRange);
//    pos += insertRNG();
//    permutedMPIs.insert(pos, index);
//}
//
//for (unsigned int rng = 0; rng < permutedMPIs.size(); rng++) {
//unsigned int index = permutedMPIs[rng];

    unsigned int prev = (index + snake.size() - 1) % snake.size();
    //unsigned int next = (i + 1) % snake.size();

    //double prevSegmentPrevCost = segmentCost[prev];
    //double nextSegmentPrevCost = segmentCost[next];

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
    cost = cost - segmentCost[prev] - segmentCost[index]
            + prevSegmentNewCost + nextSegmentNewCost;
    segmentCost[prev] = prevSegmentNewCost;
    segmentCost[index] = nextSegmentNewCost;

    //cout << "moving point to (" << bestMoveSoFar.x << ", " << bestMoveSoFar.y << ")" << endl;
    //cout << "new cost=" << cost << endl;
    //cout << "new segmentA=" << segmentCost[prev] << endl;
    //cout << "new segmentB=" << segmentCost[index] << endl;
//}
    return index;
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

    //vector<int> trials(snake->size(), 0);

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

    // Copy result to input snake.
    snake->clear();
    snake->insert(snake->end(), bestConfig.getSnake().begin(), bestConfig.getSnake().end());

    //cout << "move statistics:" << endl;
    //for (unsigned int i = 0; i < snake->size(); i++) {
    //    cout << "index " << i << " frozen=" << (*snake)[i].first << " moved=" << trials[i] << endl;
    //}
};

} // namespace enblend

#endif /* __ANNEAL_H__ */
