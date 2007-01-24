/*
 * Copyright (C) 2004-2007 Andrew Mihal
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
//#include <list>
#ifdef _WIN32
#include <slist>
#else
#include <ext/slist>
#endif
#include <algorithm>
#include <vector>

#ifdef _WIN32
#include <cmath>
#else
#include <math.h>
#endif

#include "gpu.h"
#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"
#include "vigra_ext/XMIWrapper.h"

using std::for_each;
using std::pair;
using std::vector;
//using std::list;
#ifdef _WIN32
using std::slist;
#else
using __gnu_cxx::slist;
#endif

using boost::lambda::bind;
using boost::lambda::_1;
using boost::lambda::delete_ptr;

using vigra::LineIterator;
using vigra::Point2D;
using vigra::Rect2D;

using vigra_ext::copyPaintedSetToImage;

namespace enblend {

template <typename CostImage>
void drawDottedLine(CostImage & i, vector<Point2D> & l, typename CostImage::PixelType p) {
    typedef typename CostImage::PixelType CostImagePixelType;

    miPixel pixels[2];
    pixels[0] = p;
    pixels[1] = p;
    miGC *pGC = miNewGC(2, pixels);
    miPaintedSet *paintedSet = miNewPaintedSet();
    miPoint *mip = new miPoint[l.size()];

    int index = 0;
    for (vector<Point2D>::iterator points = l.begin(); points != l.end(); ++points, ++index) {
        mip[index].x = (*points).x;
        mip[index].y = (*points).y;
    }

    miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, index, mip);
    copyPaintedSetToImage(destImageRange(i), paintedSet, Diff2D(0,0));
    miClearPaintedSet(paintedSet);

    p = (p > (NumericTraits<CostImagePixelType>::max() / 2))
            ? NumericTraits<CostImagePixelType>::zero()
            : NumericTraits<CostImagePixelType>::max();
    pixels[0] = p;
    pixels[1] = p;
    miSetGCPixels(pGC, 2, pixels);
    
    miDrawPoints(paintedSet, pGC, MI_COORD_MODE_ORIGIN, index, mip);
    copyPaintedSetToImage(destImageRange(i), paintedSet, Diff2D(0,0));

    miDeleteGC(pGC);
    miDeletePaintedSet(paintedSet);
    delete[] mip;
}

template <typename CostImage, typename VisualizeImage>
class GDAConfiguration {
public:
    typedef typename CostImage::PixelType CostImagePixelType;
    typedef typename CostImage::const_traverser CostIterator;

    GDAConfiguration(const CostImage* const d, slist<pair<bool, Point2D> > *v, VisualizeImage* const vi)
            : costImage(d),
              visualizeStateSpaceImage(vi),
              E(NULL), Pi(NULL), EF(NULL), PiF(NULL) {

        kMax = 1;

        int costImageShortDimension = std::min(costImage->width(), costImage->height());
        // Determine state space of currentPoint
        int stateSpaceWidth = costImageShortDimension / 3;

        slist<pair<bool, Point2D> >::iterator last = v->previous(v->end());
        Point2D previousPoint = last->second;
        for (slist<pair<bool, Point2D> >::iterator current = v->begin(); current != v->end(); ) {

            bool currentMoveable = current->first;
            Point2D currentPoint = current->second;
            ++current;
            Point2D nextPoint = (current == v->end()) ? v->begin()->second : current->second;

            mfEstimates.push_back(currentPoint);

            vector<Point2D> *stateSpace = new vector<Point2D>();
            pointStateSpaces.push_back(stateSpace);

            vector<int> *stateDistances = new vector<int>();
            pointStateDistances.push_back(stateDistances);

            if (currentMoveable) {
                // vp = vector from previousPoint to currentPoint
                Diff2D vp(currentPoint.x - previousPoint.x, currentPoint.y - previousPoint.y);
                // vn = vector from currentPoint to nextPoint
                Diff2D vn(nextPoint.x - currentPoint.x, nextPoint.y - currentPoint.y);
                // np = normal to vp
                Diff2D np(-vp.y, vp.x);
                // nn = normal to vn
                Diff2D nn(-vn.y, vn.x);
                
                // normal = normal vector at currentPoint
                // normal points to the left of vp and vn.
                Diff2D normal = np + nn;
                normal *= (stateSpaceWidth / normal.magnitude());

                Diff2D leftPoint = currentPoint + normal;
                Diff2D rightPoint = currentPoint - normal;

                // Choose a reasonable number of state points between these extremes
                int lineLength = std::max(std::abs(rightPoint.x - leftPoint.x),
                                          std::abs(rightPoint.y - leftPoint.y));
                int spaceBetweenPoints = static_cast<int>(ceil(lineLength / (double)GDAKmax));

                LineIterator<Diff2D> linePoint(currentPoint, leftPoint);
                for (int i = 0; i < (lineLength+1)/2; ++i, ++linePoint) {
                    if ((*costImage)[*linePoint] == NumericTraits<CostImagePixelType>::max()) break;
                    else if (((i % spaceBetweenPoints) == 0) && costImage->isInside(*linePoint)) {
                        stateSpace->push_back(Point2D(*linePoint));
                        stateDistances->push_back(std::max(std::abs(linePoint->x - currentPoint.x),
                                                           std::abs(linePoint->y - currentPoint.y)) / 2);
                        if (visualizeStateSpaceImage) (*visualizeStateSpaceImage)[*linePoint].setBlue(255);
                    }
                }
                linePoint = LineIterator<Diff2D>(currentPoint, rightPoint);
                ++linePoint;
                for (int i=(lineLength+1)/2; i < lineLength; ++i, ++linePoint) {
                    if ((*costImage)[*linePoint] == NumericTraits<CostImagePixelType>::max()) break;
                    else if (((i % spaceBetweenPoints) == 0) && costImage->isInside(*linePoint)) {
                        stateSpace->push_back(Point2D(*linePoint));
                        stateDistances->push_back(std::max(std::abs(linePoint->x - currentPoint.x),
                                                           std::abs(linePoint->y - currentPoint.y)) / 2);
                        if (visualizeStateSpaceImage) (*visualizeStateSpaceImage)[*linePoint].setBlue(255);
                    }
                }

                //LineIterator<Diff2D> linePoint(leftPoint, rightPoint);
                //for (int i = 0; i < lineLength; ++i, ++linePoint) {
                //    if (((i % spaceBetweenPoints) == 0)
                //            && costImage->isInside(*linePoint)
                //            && ((*costImage)[*linePoint] != NumericTraits<CostImagePixelType>::max())) {
                //        stateSpace->push_back(Point2D(*linePoint));
                //        stateDistances->push_back(std::max(std::abs(linePoint->x - currentPoint.x),
                //                                           std::abs(linePoint->y - currentPoint.y)) / 2);

                //        if (visualizeStateSpaceImage) (*visualizeStateSpaceImage)[*linePoint].setBlue(255);
                //    }
                //}

            }

            if (stateSpace->size() == 0) {
                stateSpace->push_back(currentPoint);
                stateDistances->push_back(0);
                if (visualizeStateSpaceImage && costImage->isInside(currentPoint)) {
                    (*visualizeStateSpaceImage)[currentPoint].setBlue(200);
                }
            }

            unsigned int localK = stateSpace->size();

            if (localK > GDAKmax) {
                cerr << "enblend: localK=" << localK << " > GDAKmax=" << GDAKmax << endl;
                exit(1);
            }

            kMax = std::max(kMax, localK);

            pointStateProbabilities.push_back(new vector<double>(localK, 1.0 / localK));

            convergedPoints.push_back(localK < 2);

            previousPoint = currentPoint;
        }

        if (UseGPU) {
            EF = new float[kMax * mfEstimates.size()];
            PiF = new float[kMax * mfEstimates.size()];
        } else {
            E = new int[kMax];
            Pi = new double[kMax];
        }

        tau = 0.75;
        deltaEMax = 7000.0;
        deltaEMin = 5.0;
        double epsilon = 1.0 / (kMax * kMax);
        tInitial = ceil(deltaEMax / log((kMax - 1 + (kMax * kMax * epsilon)) / (kMax - 1 - (kMax * kMax * epsilon))));
        tFinal = deltaEMin / log((kMax - (kMax * epsilon) - 1) / (kMax * epsilon));

    }

    ~GDAConfiguration() {
        for_each(pointStateSpaces.begin(), pointStateSpaces.end(), bind(delete_ptr(),_1));
        for_each(pointStateProbabilities.begin(), pointStateProbabilities.end(), bind(delete_ptr(),_1));
        for_each(pointStateDistances.begin(), pointStateDistances.end(), bind(delete_ptr(),_1));
        delete[] E;
        delete[] Pi;
        delete[] EF;
        delete[] PiF;
    }

    void run() {
        int progressIndicator = 1;
        int numIterations = (int)ceil(log(tFinal/tInitial)/log(tau));
        int iterationCount = 0;
        int iterationsPerTick = (numIterations+3) / 4;

        if (UseGPU) configureGPUTextures(kMax, pointStateSpaces.size());

        tCurrent = tInitial;

        while (tCurrent > tFinal) {
            double epsilon = 1.0 / kMax;
            unsigned int eta = (unsigned int)ceil(log(epsilon)
                             / log(((kMax - 2.0) / (2.0 * kMax) * exp(-tCurrent / deltaEMax)) + 0.5));

            if (Verbose > VERBOSE_GDA_MESSAGES) {
                cout << endl << "tCurrent=" << tCurrent << " eta=" << eta << " kMax=" << kMax;
                cout.flush();
            }

            for (unsigned int i = 0; i < eta; i++) iterate();

            tCurrent *= tau;

            if (Verbose > VERBOSE_GDA_MESSAGES) {
                int numConvergedPoints = 0;
                for (unsigned int i = 0; i < convergedPoints.size(); i++) {
                    if (convergedPoints[i]) numConvergedPoints++;
                }
                cout << " converged=" << numConvergedPoints << "/" << convergedPoints.size();
                cout.flush();
            }
            else if ((Verbose > VERBOSE_MASK_MESSAGES) && (iterationCount % iterationsPerTick) == 0) {
                cout << " " << progressIndicator++ << "/4";
                cout.flush();
            }

            iterationCount++;

        }

        if (UseGPU) clearGPUTextures();

        if (visualizeStateSpaceImage) {
            // Remaining unconverged state space points
            for (unsigned int i = 0; i < pointStateSpaces.size(); ++i) {
                vector<Point2D> *stateSpace = pointStateSpaces[i];
                for (unsigned int j = 0; j < stateSpace->size(); ++j) {
                    Point2D point = (*stateSpace)[j];
                    if (visualizeStateSpaceImage->isInside(point)) (*visualizeStateSpaceImage)[point].setGreen(255);
                }
            }
            // Optimized contour
            //drawDottedLine(visualizeStateSpaceImage, mfEstimates, 225);
        }

        if (Verbose > VERBOSE_GDA_MESSAGES) {
            cout << endl;
            for (unsigned int i = 0; i < convergedPoints.size(); i++) {
                if (!convergedPoints[i]) {
                    cout << "Unconverged point:" << endl;
                    vector<Point2D> *stateSpace = pointStateSpaces[i];
                    vector<double> *stateProbabilities = pointStateProbabilities[i];
                    unsigned int localK = stateSpace->size();
                    for (unsigned int state = 0; state < localK; ++state) {
                        cout << "    state " << (*stateSpace)[state]
                             << " weight=" << (*stateProbabilities)[state] << endl;
                    }
                    cout << "    mfEstimate=" << mfEstimates[i] << endl;
                }
            }
        }

    }

    vector<Point2D> & getCurrentPoints() { return mfEstimates; }

protected:

    //virtual void calculateStateProbabilities() {
    inline void calculateStateProbabilities() {

        unsigned int lastIndex = mfEstimates.size() - 1;
        for (unsigned int index = 0; index < mfEstimates.size(); ++index) {
            // Skip updating points that have already converged.
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            vector<int> *stateDistances = pointStateDistances[index];
            unsigned int localK = stateSpace->size();

            unsigned int nextIndex = (index + 1) % mfEstimates.size();
            Point2D lastPointEstimate = mfEstimates[lastIndex];
            bool lastPointInCostImage = costImage->isInside(lastPointEstimate);
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            bool nextPointInCostImage = costImage->isInside(nextPointEstimate);
            lastIndex = index;

            // Calculate E values.
            // exp_a scaling factor is part of the Schraudolph approximation.
            // for all e_i, e_j, in E: -700 < e_j-e_i < 700
            double exp_a = 1512775 / tCurrent; // = (1048576 / M_LN2) / tCurrent;
            for (unsigned int i = 0; i < localK; ++i) {
                Point2D currentPoint = (*stateSpace)[i];
                E[i] = (*stateDistances)[i];
                if (lastPointInCostImage) E[i] += costImageCost(lastPointEstimate, currentPoint);
                if (nextPointInCostImage) E[i] += costImageCost(currentPoint, nextPointEstimate);
                E[i] = NumericTraits<int>::fromRealPromote(E[i] * exp_a);
                Pi[i] = 0.0;
            }

            // Calculate new stateProbabilities
            // An = 1 / (1 + exp((E[j] - E[i]) / tCurrent))
            // I am using an approximation of the exp function from:
            // Nicol N. Schraudolph. A Fast, Compact Approximation of the Exponential Function.
            // Neural Computation, vol. 11, pages 853--862, 1999.
            union {
                double d;
                #ifdef WORDS_BIGENDIAN
                    struct { int i, j; } n;
                #else
                    struct { int j, i; } n;
                #endif
            } eco;
            eco.n.j = 0;

            // An = 1 / (1 + exp( (E[j] - E[i]) / T )
            // pi[j]' = 1/K * sum_(0)_(k-1) An(i,j) * (pi[i] + pi[j])
            for (unsigned int j = 0; j < localK; ++j) {
                double piTj = (*stateProbabilities)[j];
                Pi[j] += piTj;
                int ej = E[j];
                for (unsigned int i = (j+1); i < localK; ++i) {
                    double piT = (*stateProbabilities)[i] + piTj;
                    eco.n.i = (ej - E[i]) + (1072693248 - 60801);
                    // FIXME eco.n.i is overflowing into NaN range!
                    double piTAn = piT / (1 + eco.d);
                    //if (isnan(piTAn)) {
                    //    cout << endl << "piTAn=" << piTAn << " piT=" << piT << " denom=" << (1+eco.d) << endl;
                    //    cout << "eco.n.i=" << eco.n.i << " ej=" << ej << " ei=" << E[i] << " ej-ei=" << (ej - E[i]) << endl;
                    //    printf("%08x         ej-ei\n", (ej - E[i]));
                    //    printf("%08x         adj\n", 1072693248 - 60801);
                    //    printf("%08x%08x\n", eco.n.i, eco.n.j);
                    //}
                    if (isnan(piTAn)) {
                        // exp term is infinity or zero.
                        if (ej > E[i]) piTAn = 0.0;
                        else piTAn = piT;
                    }
                    Pi[j] += piTAn;
                    Pi[i] += piT - piTAn;
                }
                double result = Pi[j] / localK;
                //if (isnan(result)) {
                //    cerr << endl << "nan! Pi[j]=" << Pi[j] << " localK=" << localK << endl;
                //    for (unsigned int n = 0; n < localK; ++n) {
                //        cerr << "    sp[" << n << "] = " << (*stateProbabilities)[n]
                //             << "\te[" << n << "] = " << E[n] << endl;
                //    }
                //}
                (*stateProbabilities)[j] = result;
            }
        }

    }

    //virtual void calculateStateProbabilitiesGPU() {
    inline void calculateStateProbabilitiesGPU() {

        unsigned int unconvergedPoints = 0;
        unsigned int lastIndex = mfEstimates.size() - 1;
        for (unsigned int index = 0; index < mfEstimates.size(); ++index) {
            // Skip updating points that have already converged.
            if (convergedPoints[index]) continue;

            unsigned int rowIndex = unconvergedPoints / 4;
            unsigned int vectorIndex = unconvergedPoints % 4;
            float *EFbase = &(EF[(rowIndex * kMax * 4) + vectorIndex]);
            float *PiFbase = &(PiF[(rowIndex * kMax * 4) + vectorIndex]);

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            vector<int> *stateDistances = pointStateDistances[index];
            unsigned int localK = stateSpace->size();

            unsigned int nextIndex = (index + 1) % mfEstimates.size();
            Point2D lastPointEstimate = mfEstimates[lastIndex];
            bool lastPointInCostImage = costImage->isInside(lastPointEstimate);
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            bool nextPointInCostImage = costImage->isInside(nextPointEstimate);
            lastIndex = index;

            // Calculate E values.
            for (unsigned int i = 0; i < localK; ++i) {
                Point2D currentPoint = (*stateSpace)[i];
                EFbase[4*i] = (*stateDistances)[i];
                if (lastPointInCostImage) EFbase[4*i] += costImageCost(lastPointEstimate, currentPoint);
                if (nextPointInCostImage) EFbase[4*i] += costImageCost(currentPoint, nextPointEstimate);
                PiFbase[4*i] = static_cast<float>((*stateProbabilities)[i]);
                if (isnan(PiFbase[4*i]) || isnan((*stateProbabilities)[i])) {
                    union {
                        double d;
                        #ifdef WORDS_BIGENDIAN
                            struct { int i, j; } n;
                        #else
                            struct { int j, i; } n;
                        #endif
                    } eco;
                    cout << "gpu incoming pi is nan: PiFbase=";
                    printf("%08x", PiFbase[4*i]);
                    cout << " spi=";
                    eco.d = (*stateProbabilities)[i];
                    printf("%08x%08x\n", eco.n.i, eco.n.j);
                }
                if (isnan(EFbase[4*i])) {
                    union {
                        double d;
                        float f;
                        #ifdef WORDS_BIGENDIAN
                            struct { int i, j; } n;
                        #else
                            struct { int j, i; } n;
                        #endif
                    } eco;
                    eco.n.i = 0;
                    eco.n.j = 0;
                    eco.f = EFbase[4*i];
                    cout << "gpu incoming EF is nan: EFbase=";
                    printf("%08x%08x\n", eco.n.i, eco.n.j);
                }
            }
            for (unsigned int i = localK; i < kMax; ++i) {
                PiFbase[4*i] = 0.0f;
            }

            unconvergedPoints++;
        }

        // Calculate all of the new PiF values on the GPU in parallel
        gpuGDAKernel(kMax, unconvergedPoints, tCurrent, EF, PiF, PiF);

        // Write the results back to pointStateProbabilities
        unconvergedPoints = 0;
        for (unsigned int index = 0; index < mfEstimates.size(); ++index) {
            // Skip updating points that have already converged.
            if (convergedPoints[index]) continue;

            unsigned int rowIndex = unconvergedPoints / 4;
            unsigned int vectorIndex = unconvergedPoints % 4;
            float *PiFbase = &(PiF[(rowIndex * kMax * 4) + vectorIndex]);

            vector<double> *stateProbabilities = pointStateProbabilities[index];
            unsigned int localK = stateProbabilities->size();

            for (unsigned int i = 0; i < localK; ++i) {
                (*stateProbabilities)[i] = static_cast<double>(PiFbase[4*i]);
                if (isnan(PiFbase[4*i]) || isnan((*stateProbabilities)[i])) {
                    union {
                        double d;
                        float f;
                        #ifdef WORDS_BIGENDIAN
                            struct { int i, j; } n;
                        #else
                            struct { int j, i; } n;
                        #endif
                    } eco;
                    cout << "gpu outgoing pi is nan: PiFbase=";
                    eco.n.i = eco.n.j = 0;
                    eco.f = PiFbase[4*i];
                    printf("%08x%08x", eco.n.i, eco.n.j);
                    cout << " spi=";
                    eco.d = (*stateProbabilities)[i];
                    printf("%08x%08x\n", eco.n.i, eco.n.j);
                }
            }

            unconvergedPoints++;
        }

    }

    //virtual void iterate() {
    void iterate() {

        if (UseGPU) {
            calculateStateProbabilitiesGPU();

            //// Copy GPU-calculated results
            //vector<vector<double> > gpuStateProbabilities;
            //for (unsigned int index = 0; index < pointStateProbabilities.size(); ++index) {
            //    gpuStateProbabilities.push_back(vector<double>(*(pointStateProbabilities[index])));
            //}

            //// Do regular CPU computations
            //calculateStateProbabilities();

            //// Compare
            //for (unsigned int index = 0; index < pointStateProbabilities.size(); ++index) {
            //    vector<double> &gpuProbs = gpuStateProbabilities[index];
            //    vector<double> &cpuProbs = *(pointStateProbabilities[index]);
            //    cout << "index " << index << endl;
            //    for (unsigned int k = 0; k < gpuProbs.size(); ++k) {
            //        double diff = std::abs(gpuProbs[k] - cpuProbs[k]);
            //        if (diff > 0.001) {
            //            cout << gpuProbs[k] << ", " << cpuProbs[k] << " abs=" << std::abs(gpuProbs[k] - cpuProbs[k]) << endl;
            //        }
            //    }
            //}

        } else {
            calculateStateProbabilities();
        }

        kMax = 1;
        for (unsigned int index = 0; index < pointStateSpaces.size(); ++index) {
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            vector<int> *stateDistances = pointStateDistances[index];
            unsigned int localK = stateSpace->size();
            double estimateX = 0.0;
            double estimateY = 0.0;

            // Make new mean field estimates.
            double totalWeight = 0.0;
            bool hasHighWeightState = false;
            for (unsigned int k = 0; k < localK; ++k) {
                double weight = (*stateProbabilities)[k];
                totalWeight += weight;
                if (weight > 0.99) hasHighWeightState = true;
                Point2D state = (*stateSpace)[k];
                estimateX += weight * (double)state.x;
                estimateY += weight * (double)state.y;
            }
            estimateX /= totalWeight;
            estimateY /= totalWeight;

            Point2D newEstimate(NumericTraits<int>::fromRealPromote(estimateX),
                                NumericTraits<int>::fromRealPromote(estimateY));

            // Sanity check
            if (!costImage->isInside(newEstimate)) {
                union {
                    double d;
                    #ifdef WORDS_BIGENDIAN
                        struct { int i, j; } n;
                    #else
                        struct { int j, i; } n;
                    #endif
                } eco;
                cout << endl << "enblend: optimizer warning: new mean field estimate is outside cost image." << endl;
                for (unsigned int state = 0; state < localK; ++state) {
                    cout << "    state " << (*stateSpace)[state]
                         << " weight=";
                    cout << (*stateProbabilities)[state] << " = ";
                    eco.d = (*stateProbabilities)[state];
                    printf("%08x%08x\n", eco.n.i, eco.n.j);
                }
                cout << "    mfEstimate=" << newEstimate << endl;
                // Skip this point from now on.
                convergedPoints[index] = true;
                continue;
            }

            mfEstimates[index] = newEstimate;

            // Remove improbable solutions from the search space
            double totalWeights = 0.0;
            double cutoffWeight = hasHighWeightState ? 0.50 : 0.00001;
            for (unsigned int k = 0; k < stateSpace->size(); ) {
                double weight = (*stateProbabilities)[k];
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
            if (localK < 2) convergedPoints[index] = true;
            kMax = std::max(kMax, stateProbabilities->size());

        }

    }

    //virtual int costImageCost(const Point2D &start, const Point2D &end) {
    inline int costImageCost(const Point2D &start, const Point2D &end) {
        //if (!(costImage->isInside(start) && costImage->isInside(end))) {
        //    cerr << "start and end points are not inside image: start=" << start << " end=" << end << endl;
        //    exit(-1);
        //}

        int cost = 0;

        int lineLength = std::max(std::abs(end.x - start.x), std::abs(end.y - start.y));

        if (lineLength > 0) {
            LineIterator<CostIterator> lineStart(costImage->upperLeft() + start, costImage->upperLeft() + end);
            for (int i = 0; i < lineLength; ++i) {
                cost += *lineStart;
                ++lineStart;
            }
        }

        if (lineLength < 8) cost += NumericTraits<CostImagePixelType>::max() * (8 - lineLength);

        return cost;
    }

    bool segmentIntersect(const Point2D & l1a, const Point2D & l1b, const Point2D & l2a, const Point2D & l2b) {
        int denom = (l2b.y - l2a.y)*(l1b.x - l1a.x) - (l2b.x - l2a.x)*(l1b.y - l1a.y);
        if (denom == 0) return false; // lines are parallel or coincident
        int uaNum = (l2b.x - l2a.x)*(l1a.y - l2a.y) - (l2b.y - l2a.y)*(l1a.x - l2a.x);
        int ubNum = (l1b.x - l1a.x)*(l1a.y - l2a.y) - (l1b.y - l1a.y)*(l1a.x - l2a.x);
        if (denom < 0) { uaNum *= -1; ubNum *= -1; denom *= -1; }
        if (uaNum > 0 && uaNum < denom && ubNum > 0 && ubNum < denom) return true;
        return false;
    }

    const CostImage *costImage;
    VisualizeImage *visualizeStateSpaceImage;

    // Mean-field estimates of current point locations
    vector<Point2D> mfEstimates;

    // State spaces of each point
    vector<vector<Point2D>* > pointStateSpaces;

    // Probability vectors for each state space
    vector<vector<double>* > pointStateProbabilities;

    vector<vector<int>* > pointStateDistances;

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

    // Data arrays for CPU probability calculations
    int *E;
    double *Pi;

    // Data arrays for GPU probability calculations
    float *EF;
    float *PiF;

};

template <typename CostImage, typename VisualizeImage>
void annealSnake(const CostImage* const ci, slist<pair<bool, Point2D> > *snake, VisualizeImage* const vi) {

    GDAConfiguration<CostImage, VisualizeImage> cfg(ci, snake, vi);

    cfg.run();

    slist<pair<bool, Point2D> >::iterator snakePoint = snake->begin();
    vector<Point2D>::iterator annealedPoint = cfg.getCurrentPoints().begin();
    for (; snakePoint != snake->end(); ++snakePoint, ++annealedPoint) {
        snakePoint->second = *annealedPoint;
    }

};

} // namespace enblend

#endif /* __ANNEAL_H__ */
