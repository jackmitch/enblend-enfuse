/*
 * Copyright (C) 2004-2009 Andrew Mihal
 * Copyright (C) 2009-2017 Christoph Spiel
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

#ifndef ANNEAL_H_INCLUDED
#define ANNEAL_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <limits>
#include <vector>

#ifdef _WIN32
#include <cmath>
#else
#include <math.h>
#endif

#include <vigra/diff2d.hxx>
#include <vigra/iteratoradapter.hxx>

#include "masktypedefs.h"
#include "muopt.h"
#include "opencl.h"
#include "opencl_anneal.h"
#include "openmp_lock.h"
#include "timer.h"


#ifdef OPENCL
namespace GPU
{
    extern std::unique_ptr<ocl::CalculateStateProbabilities> StateProbabilities;
}
#endif


namespace enblend {

inline static vigra::Diff2D
normal_vector(const vigra::Point2D& a_previous_point,
              const vigra::Point2D& a_current_point,
              const vigra::Point2D& a_next_point)
{
    // vp: vector from a_previous_point to a_current_point
    const vigra::Diff2D vp(a_current_point.x - a_previous_point.x, a_current_point.y - a_previous_point.y);
    // vn: vector from a_current_point to a_next_point
    const vigra::Diff2D vn(a_next_point.x - a_current_point.x, a_next_point.y - a_current_point.y);
    // np: normal to vp
    const vigra::Diff2D np(-vp.y, vp.x);
    // nn: normal to vn
    const vigra::Diff2D nn(-vn.y, vn.x);

    // Answer normal vector at a_current_point;
    // normal points to the left of vp and vn.
    return np + nn;
}


template <typename CostImage, typename VisualizeImage>
class GDAConfiguration
{
public:
    typedef typename CostImage::PixelType CostImagePixelType;
    typedef typename vigra::NumericTraits<CostImagePixelType>::Promote CostImagePromoteType;

    GDAConfiguration(const CostImage* const d, Segment* v, VisualizeImage* const vi) :
        costImage(d), visualizeStateSpaceImage(vi) {
        kMax = 1;
        distanceWeight = 1.0;
        mismatchWeight = 1.0;

        const int costImageShortDimension = std::min(costImage->width(), costImage->height());
        // Determine state space of currentPoint
        const int stateSpaceWidth = costImageShortDimension / 3;

        vigra::Point2D previousPoint = v->back().second;
        for (Segment::iterator current = v->begin(); current != v->end();) {
            bool currentMoveable = current->first;
            vigra::Point2D currentPoint = current->second;
            ++current;
            vigra::Point2D nextPoint = current == v->end() ? v->begin()->second : current->second;

            mfEstimates.push_back(currentPoint);

            std::vector<vigra::Point2D>* stateSpace = new std::vector<vigra::Point2D>();
            pointStateSpaces.push_back(stateSpace);

            std::vector<int>* stateDistances = new std::vector<int>();
            pointStateDistances.push_back(stateDistances);

            vigra::Diff2D normal = normal_vector(previousPoint, currentPoint, nextPoint);
            const double normal_magnitude = normal.magnitude();

            if (currentMoveable && normal_magnitude > std::numeric_limits<double>::epsilon()) {
                normal *= stateSpaceWidth / normal_magnitude;

                vigra::Diff2D leftPoint = currentPoint + normal;
                vigra::Diff2D rightPoint = currentPoint - normal;

                // Choose a reasonable number of state points between these extremes
                const int lineLength = std::max(std::abs(rightPoint.x - leftPoint.x),
                                                std::abs(rightPoint.y - leftPoint.y));
                const int spaceBetweenPoints =
                    static_cast<int>(ceil(lineLength / static_cast<double>(AnnealPara.kmax)));

                vigra::LineIterator<vigra::Diff2D> linePoint(currentPoint, leftPoint);
                for (int i = 0; i < (lineLength + 1) / 2; ++i, ++linePoint) {
                    // Stop searching along the line if we leave the
                    // cost image or enter a max-cost region.
                    if (!costImage->isInside(*linePoint)) {
                        break;
                    } else if ((*costImage)[*linePoint] == vigra::NumericTraits<CostImagePixelType>::max()) {
                        break;
                    } else if (i % spaceBetweenPoints == 0) {
                        stateSpace->push_back(vigra::Point2D(*linePoint));
                        stateDistances->push_back(std::max(std::abs(linePoint->x - currentPoint.x),
                                                           std::abs(linePoint->y - currentPoint.y)) / 2);
                        if (visualizeStateSpaceImage) {
                            (*visualizeStateSpaceImage)[*linePoint] = VISUALIZE_STATE_SPACE;
                        }
                    }
                }
                linePoint = vigra::LineIterator<vigra::Diff2D>(currentPoint, rightPoint);
                ++linePoint;
                for (int i = 1; i < 1 + lineLength / 2; ++i, ++linePoint) {
                    // Stop searching along the line if we leave the
                    // cost image or enter a max-cost region.
                    if (!costImage->isInside(*linePoint)) {
                        break;
                    } else if ((*costImage)[*linePoint] == vigra::NumericTraits<CostImagePixelType>::max()) {
                        break;
                    } else if (i % spaceBetweenPoints == 0) {
                        stateSpace->push_back(vigra::Point2D(*linePoint));
                        stateDistances->push_back(std::max(std::abs(linePoint->x - currentPoint.x),
                                                           std::abs(linePoint->y - currentPoint.y)) / 2);
                        if (visualizeStateSpaceImage) {
                            (*visualizeStateSpaceImage)[*linePoint] = VISUALIZE_STATE_SPACE;
                        }
                    }
                }
            }

            if (stateSpace->size() == 0) {
                stateSpace->push_back(currentPoint);
                stateDistances->push_back(0);
                if (visualizeStateSpaceImage && costImage->isInside(currentPoint)) {
                    (*visualizeStateSpaceImage)[currentPoint] = VISUALIZE_STATE_SPACE_INSIDE;
                }
            }

            const unsigned int localK = stateSpace->size();
            if (localK > AnnealPara.kmax) {
                std::cerr << command
                     << ": local k = " << localK << " > k_max = " << AnnealPara.kmax
                     << std::endl;
                exit(1);
            }

            kMax = std::max(kMax, localK);

            pointStateProbabilities.push_back(new std::vector<double>(localK, 1.0 / localK));

            convergedPoints.push_back(localK < 2);

            previousPoint = currentPoint;
        }

        tau = AnnealPara.tau;
        deltaEMax = AnnealPara.deltaEMax;
        deltaEMin = AnnealPara.deltaEMin;
        const double kmax = static_cast<double>(kMax);
        if (kmax > 2.5) {
            tInitial = ceil(deltaEMax / log(kmax / (kmax - 2.0)));
            tFinal = deltaEMin / log(kmax * kmax - kmax - 1.0);
        } else {
            tInitial = deltaEMax;
            tFinal = deltaEMin;
        }
    }

    virtual ~GDAConfiguration() {
        std::for_each(pointStateSpaces.begin(), pointStateSpaces.end(),
                      [](std::vector<vigra::Point2D>* x) {delete x;});
        std::for_each(pointStateProbabilities.begin(), pointStateProbabilities.end(),
                      [](std::vector<double>* x) {delete x;});
        std::for_each(pointStateDistances.begin(), pointStateDistances.end(),
                      [](std::vector<int>* x) {delete x;});
    }

    void run() {
        int progressIndicator = 1;
        int numIterations = static_cast<int>(ceil(log(tFinal / tInitial) / log(tau)));
        int iterationCount = 0;
        int iterationsPerTick = (numIterations + 3) / 4;

        tCurrent = tInitial;

        while (kMax > 1 && tCurrent > tFinal) {
            const double epsilon = 1.0 / kMax;
            const unsigned int eta =
                static_cast<unsigned int>(ceil(log(epsilon)
                                               / log(((kMax - 2.0) / (2.0 * kMax) * exp(-tCurrent / deltaEMax))
                                                     + 0.5)));

            if (Verbose >= VERBOSE_GDA_MESSAGES) {
                const std::ios_base::fmtflags ioFlags(std::cerr.flags());
                std::cerr << "\n"
                     << command
                     << ": info: t = " << std::scientific << std::setprecision(3) << tCurrent
                     << ", eta = " << std::setw(4) << eta
                     << ", k_max = " << std::setw(3) << kMax;
                std::cerr.flush();
                std::cerr.flags(ioFlags);
            }

            for (unsigned int i = 0; i < eta; i++) {
                iterate();
            }

            tCurrent *= tau;

            if (Verbose >= VERBOSE_GDA_MESSAGES) {
                int numConvergedPoints = 0;
                for (unsigned int i = 0; i < convergedPoints.size(); i++) {
                    if (convergedPoints[i]) {
                        numConvergedPoints++;
                    }
                }
                std::cerr << ", " << numConvergedPoints
                     << " of " << convergedPoints.size()
                     << " points converged";
                std::cerr.flush();
            }
            else if (Verbose >= VERBOSE_MASK_MESSAGES && iterationCount % iterationsPerTick == 0) {
                std::cerr << " " << progressIndicator << "/4";
                progressIndicator++;
                std::cerr.flush();
            }

            iterationCount++;
        }

        if (visualizeStateSpaceImage) {
            // Remaining unconverged state space points
            for (unsigned int i = 0; i < pointStateSpaces.size(); ++i) {
                std::vector<vigra::Point2D>* stateSpace = pointStateSpaces[i];
                for (unsigned int j = 0; j < stateSpace->size(); ++j) {
                    vigra::Point2D point = (*stateSpace)[j];
                    if (visualizeStateSpaceImage->isInside(point)) {
                        (*visualizeStateSpaceImage)[point] = VISUALIZE_STATE_SPACE_UNCONVERGED;
                    }
                }
            }
        }

        if (Verbose >= VERBOSE_GDA_MESSAGES) {
            std::cerr << std::endl;
            for (unsigned int i = 0; i < convergedPoints.size(); i++) {
                if (!convergedPoints[i]) {
                    std::cerr << command
                         << ": info: unconverged point: "
                         << std::endl;
                    std::vector<vigra::Point2D>* stateSpace = pointStateSpaces[i];
                    std::vector<double>* stateProbabilities = pointStateProbabilities[i];
                    const unsigned int localK = stateSpace->size();
                    for (unsigned int state = 0; state < localK; ++state) {
                        std::cerr << command
                             << ": info: state " << (*stateSpace)[state]
                             << ", weight = " << (*stateProbabilities)[state]
                             << std::endl;
                    }
                    std::cerr << command
                         << ": info: mfEstimate = " << mfEstimates[i]
                         << std::endl;
                }
            }
        }
    }

    const std::vector<vigra::Point2D>& getCurrentPoints() const {
        return mfEstimates;
    }

    void setOptimizerWeights(double aDistanceWeight, double aMismatchWeight) {
        // IMPLEMENTATION NOTE: We normalize to 2.0, because we
        // want to reproduce the results of Enblend-3.2.  Up to
        // Enblend-3.2 the distanceWeight and mismatchWeight were both
        // fixed at 1.0.
        const double sum = 0.5 * (aDistanceWeight + aMismatchWeight);
        assert(aDistanceWeight >= 0.0);
        assert(aMismatchWeight >= 0.0);
        assert(sum > 0.0);
        distanceWeight = aDistanceWeight / sum;
        mismatchWeight = aMismatchWeight / sum;
    }

protected:
    virtual void calculateStateProbabilities() {
        const int mf_size = static_cast<int>(mfEstimates.size());

#ifdef OPENMP
#pragma omp parallel
#endif
        {
            double* E = new double[kMax];
            double* Pi = new double[kMax];

#ifdef OPENMP
#pragma omp for nowait schedule(guided)
#endif
            for (int index = 0; index < mf_size; ++index) {
                // Skip updating points that have already converged.
                convergedPointsLock.set();
                if (convergedPoints[index]) {
                    convergedPointsLock.unset();
                    continue;
                }
                convergedPointsLock.unset();

                const std::vector<vigra::Point2D>* stateSpace = pointStateSpaces[index];
                std::vector<double>* stateProbabilities = pointStateProbabilities[index];
                const std::vector<int>* stateDistances = pointStateDistances[index];
                const unsigned int localK = stateSpace->size();

                const int lastIndex = index == 0 ? mf_size - 1 : index - 1;
                const unsigned int nextIndex = (index + 1) % mf_size;
                const vigra::Point2D lastPointEstimate = mfEstimates[lastIndex];
                const bool lastPointInCostImage = costImage->isInside(lastPointEstimate);
                const vigra::Point2D nextPointEstimate = mfEstimates[nextIndex];
                const bool nextPointInCostImage = costImage->isInside(nextPointEstimate);

                // Calculate E values.
                for (unsigned i = 0U; i < localK; ++i) {
                    const vigra::Point2D currentPoint = (*stateSpace)[i];
                    const int distanceCost = (*stateDistances)[i];
                    int mismatchCost = 0;
                    if (lastPointInCostImage) {
                        mismatchCost += costImageCost(lastPointEstimate, currentPoint);
                    }
                    if (nextPointInCostImage) {
                        mismatchCost += costImageCost(currentPoint, nextPointEstimate);
                    }

                    const double cost =
                        distanceWeight * static_cast<double>(distanceCost) +
                        mismatchWeight * static_cast<double>(mismatchCost);
                    E[i] = cost / tCurrent;
                    Pi[i] = 0.0;
                }

                timer::WallClock wall_clock;
                wall_clock.start();
                // Calculate new stateProbabilities
                // An = 1 / (1 + exp((E[j] - E[i]) / tCurrent))
                // pi[j]' = 1/K * sum_(0)_(k-1) An(i,j) * (pi[i] + pi[j])
                for (unsigned j = 0U; j < localK; ++j) {
                    const double piTj = (*stateProbabilities)[j];
                    Pi[j] += piTj;
                    const double ej = E[j];
                    for (unsigned i = j + 1U; i < localK; ++i) {
                        const double piT = (*stateProbabilities)[i] + piTj;
                        double piTAn = piT / (1.0 + std::exp(ej - E[i]));
                        if (EXPECT_RESULT(std::isnan(piTAn), false)) {
                            // exp term is infinity or zero.
                            piTAn = ej > E[i] ? 0.0 : piT;
                        }
                        Pi[j] += piTAn;
                        Pi[i] += piT - piTAn;
                    }
                    (*stateProbabilities)[j] = Pi[j] / localK;
                }
                wall_clock.stop();
                if (parameter::as_boolean("time-state-probabilities", false))
                {
                    ocl::StowFormatFlags _;

                    std::cerr <<
                        "\n" <<
                        command << ": timing: wall-clock runtime of `Calculate New State Probabilities' (CPU): " <<
                        std::setprecision(3) << 1e6 * wall_clock.value() << " µs\n" <<
                        std::endl;
                }
            }

            delete [] E;
            delete [] Pi;
        } // omp parallel
    }

    void iterate() {
        calculateStateProbabilities();

        kMax = 1;
        size_t kmax_local = 1;

#ifdef OPENMP
#pragma omp parallel firstprivate(kmax_local)
#endif
        {
#ifdef OPENMP
#pragma omp for nowait schedule(guided)
#endif
            for (int index = 0; index < static_cast<int>(pointStateSpaces.size()); ++index) {
                convergedPointsLock.set();
                if (convergedPoints[index]) {
                    convergedPointsLock.unset();
                    continue;
                }
                convergedPointsLock.unset();

                std::vector<vigra::Point2D>* stateSpace = pointStateSpaces[index];
                std::vector<double>* stateProbabilities = pointStateProbabilities[index];
                std::vector<int>* stateDistances = pointStateDistances[index];
                unsigned int localK = stateSpace->size();
                double estimateX = 0.0;
                double estimateY = 0.0;

                // Make new mean field estimates.
                double totalWeight = 0.0;
                bool hasHighWeightState = false;
                for (unsigned int k = 0; k < localK; ++k) {
                    const double weight = (*stateProbabilities)[k];
                    totalWeight += weight;
                    if (weight > 0.99) {
                        hasHighWeightState = true;
                    }
                    const vigra::Point2D state = (*stateSpace)[k];
                    estimateX += weight * static_cast<double>(state.x);
                    estimateY += weight * static_cast<double>(state.y);
                }
                estimateX /= totalWeight;
                estimateY /= totalWeight;

                vigra::Point2D newEstimate(vigra::NumericTraits<int>::fromRealPromote(estimateX),
                                           vigra::NumericTraits<int>::fromRealPromote(estimateY));

                // Sanity check
                if (!costImage->isInside(newEstimate)) {
                    cerrLock.set();
                    std::cerr << command
                              << ": warning: new mean field estimate outside cost image"
                              << std::endl;
                    for (unsigned int state = 0; state < localK; ++state) {
                        std::cerr << command
                                  << ": note: state " << (*stateSpace)[state]
                                  << " weight = "
                                  << (*stateProbabilities)[state]
                                  << std::endl;
                    }
                    std::cerr << command
                              << ": note: new estimate = " << newEstimate
                              << std::endl;
                    cerrLock.unset();

                    // Skip this point from now on.
                    convergedPointsLock.set();
                    convergedPoints[index] = true;
                    convergedPointsLock.unset();
                    continue;
                }

                mfEstimates[index] = newEstimate;

                // Remove improbable solutions from the search space
                double totalWeights = 0.0;
                const double cutoffWeight = hasHighWeightState ? 0.50 : 0.00001;
                for (unsigned int k = 0; k < stateSpace->size(); ) {
                    const double weight = (*stateProbabilities)[k];
                    if (weight < cutoffWeight) {
                        // Replace this state with last state
                        (*stateProbabilities)[k] = (*stateProbabilities)[stateProbabilities->size() - 1];
                        (*stateSpace)[k] = (*stateSpace)[stateSpace->size() - 1];
                        (*stateDistances)[k] = (*stateDistances)[stateDistances->size() - 1];

                        // Delete last state
                        stateProbabilities->pop_back();
                        stateSpace->pop_back();
                        stateDistances->pop_back();
                    } else {
                        totalWeights += weight;
                        ++k;
                    }
                }

                // Renormalize
                for (unsigned int k = 0; k < stateSpace->size(); ++k) {
                    (*stateProbabilities)[k] /= totalWeights;
                }

                localK = stateSpace->size();
                if (localK < 2) {
                    convergedPointsLock.set();
                    convergedPoints[index] = true;
                    convergedPointsLock.unset();
                }

                kmax_local = std::max(kmax_local, stateProbabilities->size());
            }

            kMaxLock.set();
            kMax = std::max(kMax, static_cast<unsigned int>(kmax_local));
            kMaxLock.unset();
        } // omp parallel
    }

    int costImageCost(const vigra::Point2D& start_point, const vigra::Point2D& end_point) const {
        typedef typename CostImage::ConstIterator CostIterator;

        const int shortLineThreshold = 8; // We penalize lines below this limit.
        const CostIterator end(costImage->upperLeft() + end_point);

        vigra::LineIterator<CostIterator> lineEnd(end, end);
        vigra::LineIterator<CostIterator> line(costImage->upperLeft() + start_point, end);
        int cost = 0;

        while (line != lineEnd) {
            cost += *line;
            ++line;
        }

        const int lineLength =
            std::max(std::abs(end_point.x - start_point.x), std::abs(end_point.y - start_point.y));

        if (lineLength < shortLineThreshold) {
            cost += vigra::NumericTraits<CostImagePixelType>::max() * (shortLineThreshold - lineLength);
        }

        return cost;
    }

    bool segmentIntersect(const vigra::Point2D& l1a, const vigra::Point2D& l1b,
                          const vigra::Point2D& l2a, const vigra::Point2D& l2b) const {
        int denom = (l2b.y - l2a.y) * (l1b.x - l1a.x) - (l2b.x - l2a.x) * (l1b.y - l1a.y);
        if (denom == 0) {
            return false;       // lines are parallel or coincident
        }

        int uaNum = (l2b.x - l2a.x) * (l1a.y - l2a.y) - (l2b.y - l2a.y) * (l1a.x - l2a.x);
        int ubNum = (l1b.x - l1a.x) * (l1a.y - l2a.y) - (l1b.y - l1a.y) * (l1a.x - l2a.x);
        if (denom < 0) {
            uaNum *= -1;
            ubNum *= -1;
            denom *= -1;
        }
        if (uaNum > 0 && uaNum < denom && ubNum > 0 && ubNum < denom) {
            return true;
        }
        return false;
    }

    const CostImage *costImage;
    VisualizeImage *visualizeStateSpaceImage;

    // Mean-field estimates of current point locations
    std::vector<vigra::Point2D> mfEstimates;

    // State spaces of each point
    std::vector<std::vector<vigra::Point2D>*> pointStateSpaces;

    // Probability vectors for each state space
    std::vector<std::vector<double>*> pointStateProbabilities;

    std::vector<std::vector<int>*> pointStateDistances;

    // Flags indicate which points have converged
    std::vector<bool> convergedPoints;
    omp::lock convergedPointsLock;

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
    omp::lock kMaxLock;

    // Weight factors for the distance of a point from the initial
    // seam line and the total mismatch accumulated along the seam
    // line segment.
    double distanceWeight;;
    double mismatchWeight;

    omp::lock cerrLock;
}; // class GDAConfiguration


#ifdef OPENCL

template <typename CostImage, typename VisualizeImage>
class GDAConfigurationGPU : public GDAConfiguration<CostImage, VisualizeImage>
{
    typedef GDAConfiguration<CostImage, VisualizeImage> super;

public:
    GDAConfigurationGPU(const CostImage* const d, Segment* v, VisualizeImage* const vi) : super(d, v, vi)
    {
#ifdef DEBUG
        std::cerr << command << ": info: choose OpenCL accelaration for Anneal Snake" << std::endl;
#endif
    }

protected:
    void calculateStateProbabilities() override // GPU version
    {
        const int mf_size = static_cast<int>(super::mfEstimates.size());

        const size_t maximum_probability_vector_size =
            (*std::max_element(super::pointStateProbabilities.begin(),
                               super::pointStateProbabilities.end(),
                               [] (const std::vector<double>* x, const std::vector<double>* y)
        {return x->size() < y->size();}))->size();

        // Method GPU::StateProbabilities->setup() allocates space for
        // `E' and `Pi' for us.  In particular it will use the GPU's
        // favorite global memory.
        float* E;
        float* Pi;
        GPU::StateProbabilities->setup(maximum_probability_vector_size, static_cast<size_t>(super::kMax),
                                       E, Pi);

        for (int index = 0; index < mf_size; ++index)
        {
            if (super::convergedPoints[index])
            {
                continue;
            }

            const std::vector<vigra::Point2D>* stateSpace = super::pointStateSpaces[index];
            std::vector<double>* stateProbabilities = super::pointStateProbabilities[index];
            const std::vector<int>* stateDistances = super::pointStateDistances[index];
            const int localK = static_cast<int>(stateSpace->size());

            const int lastIndex = (index == 0 ? mf_size : index) - 1;
            const int nextIndex = (index + 1) % mf_size;
            const vigra::Point2D lastPointEstimate = super::mfEstimates[lastIndex];
            const bool lastPointInCostImage = super::costImage->isInside(lastPointEstimate);
            const vigra::Point2D nextPointEstimate = super::mfEstimates[nextIndex];
            const bool nextPointInCostImage = super::costImage->isInside(nextPointEstimate);

            // Calculate E values.
            for (int i = 0; i < localK; ++i)
            {
                const vigra::Point2D currentPoint = (*stateSpace)[i];
                const int distanceCost = (*stateDistances)[i];
                int mismatchCost = 0;
                if (lastPointInCostImage)
                {
                    mismatchCost += super::costImageCost(lastPointEstimate, currentPoint);
                }
                if (nextPointInCostImage)
                {
                    mismatchCost += super::costImageCost(currentPoint, nextPointEstimate);
                }

                const double cost =
                    super::distanceWeight * static_cast<double>(distanceCost) +
                    super::mismatchWeight * static_cast<double>(mismatchCost);
                E[i] = static_cast<float>(cost / super::tCurrent);
                Pi[i] = 0.0f;
            } // for i

            timer::WallClock wall_clock;
            wall_clock.start();
            GPU::StateProbabilities->run(localK, stateProbabilities, super::kMax, E, Pi);
            wall_clock.stop();

            if (parameter::as_boolean("time-state-probabilities", false))
            {
                ocl::StowFormatFlags _;

                std::cerr <<
                    "\n" <<
                    command << ": timing: wall-clock runtime of `Calculate New State Probabilities' (GPU): " <<
                    std::setprecision(3) << 1e6 * wall_clock.value() << " µs\n" <<
                    std::endl;
            }
        } // for index

        GPU::StateProbabilities->teardown();
    }
}; // class GDAConfigurationGPU

#endif // OPENCL


template <typename CostImage, typename VisualizeImage>
void annealSnake(const CostImage* const ci,
                 const std::pair<double, double>& optimizerWeights,
                 Segment* snake,
                 VisualizeImage* const vi)
{
    timer::WallClock wall_clock;

    wall_clock.start();

#ifdef OPENCL
    const bool enable_kernel = parameter::as_boolean("gpu-kernel-anneal", true);

    std::unique_ptr<GDAConfiguration<CostImage, VisualizeImage> >
        cfg((GPUContext && enable_kernel) ?
            new GDAConfigurationGPU<CostImage, VisualizeImage>(ci, snake, vi) :
            new GDAConfiguration<CostImage, VisualizeImage>(ci, snake, vi));

    if (GPUContext && enable_kernel)
    {
        GPU::StateProbabilities->wait();     // Ensure that kernel was built.
    }
#else
    std::unique_ptr<GDAConfiguration<CostImage, VisualizeImage> >
        cfg(new GDAConfiguration<CostImage, VisualizeImage>(ci, snake, vi));
#endif

    cfg->setOptimizerWeights(optimizerWeights.first, optimizerWeights.second);
    cfg->run();

    std::vector<vigra::Point2D>::const_iterator annealedPoint = cfg->getCurrentPoints().begin();
    for (Segment::iterator snakePoint = snake->begin();
         snakePoint != snake->end();
         ++snakePoint, ++annealedPoint)
    {
        snakePoint->second = *annealedPoint;
    }

    wall_clock.stop();
    if (parameter::as_boolean("time-anneal-snake", false))
    {
        ocl::StowFormatFlags _;

        std::cerr <<
            "\n" <<
            command << ": timing: wall-clock runtime of `Anneal Snake': " <<
            std::setprecision(3) << 1000.0 * wall_clock.value() << " ms\n" <<
            std::endl;
    }
}

} // namespace enblend

#endif // ANNEAL_H_INCLUDED


// Local Variables:
// mode: c++
// End:
