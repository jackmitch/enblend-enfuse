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
#include <list>
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

#include "vigra/diff2d.hxx"
#include "vigra/iteratoradapter.hxx"
#include "vigra_ext/XMIWrapper.h"

using std::for_each;
using std::pair;
using std::vector;
using std::list;
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

/*
template <typename CostImage>
void drawDottedLine(CostImage & i, list<Point2D> & l, typename CostImage::PixelType p) {
    typedef typename CostImage::PixelType CostImagePixelType;

    miPixel pixels[2];
    pixels[0] = p;
    pixels[1] = p;
    miGC *pGC = miNewGC(2, pixels);
    miPaintedSet *paintedSet = miNewPaintedSet();
    miPoint *mip = new miPoint[l.size()];

    int index = 0;
    for (list<Point2D>::iterator points = l.begin(); points != l.end(); ++points, ++index) {
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

template<typename CostImage>
void drawSegments(CostImage & i, list<Point2D> & l1, list<Point2D> & l2, typename CostImage::PixelType p) {
    typedef typename CostImage::PixelType CostImagePixelType;

    miPixel pixels[2];
    pixels[0] = p;
    pixels[1] = p;
    miGC *pGC = miNewGC(2, pixels);
    miPaintedSet *paintedSet = miNewPaintedSet();
    miPoint mip[2];

    for (list<Point2D>::iterator l1Point = l1.begin(), l2Point = l2.begin(); l1Point != l1.end(); ++l1Point, ++l2Point) {
        mip[0].x = (*l1Point).x;
        mip[0].y = (*l1Point).y;
        mip[1].x = (*l2Point).x;
        mip[1].y = (*l2Point).y;
        miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, 2, mip);
    }

    copyPaintedSetToImage(destImageRange(i), paintedSet, Diff2D(0,0));

    miDeleteGC(pGC);
    miDeletePaintedSet(paintedSet);
}
*/

template <typename CostImage>
class GDAConfiguration {
public:
    typedef typename CostImage::PixelType CostImagePixelType;
    typedef typename CostImage::const_traverser CostIterator;

    GDAConfiguration(const CostImage* const d, list<pair<bool, Point2D> > *v) : costImage(d) {

        kMax = 1;

        CostImage visualizeStateSpaceImage(*costImage);

        int costImageShortDimension = std::min(costImage->width(), costImage->height());
        // Determine state space of currentPoint
        int stateSpaceWidth = costImageShortDimension / 3;

/*
        list<bool> selfIntersects;
        list<Point2D> leftPoints;
        list<Point2D> middlePoints;
        list<Point2D> rightPoints;

        // FIXME assuming v is an open contour bounded by one nonmoveable vertex on each end.
        list<pair<bool, Point2D> >::iterator last = --(v->end());
        Point2D previousPoint = last->second;
        for (list<pair<bool, Point2D> >::iterator current = v->begin(); current != last; ) {

            bool currentMoveable = current->first;
            Point2D currentPoint = current->second;
            ++current;
            Point2D nextPoint = current->second;

            //cout << "previousPoint=" << previousPoint << " currentPoint=" << currentPoint << " nextPoint=" << nextPoint << endl;

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

                middlePoints.push_back(currentPoint);
                leftPoints.push_back(currentPoint + normal);
                rightPoints.push_back(currentPoint - normal);
                selfIntersects.push_back(false);

                //cout << "leftPoint=" << (currentPoint+normal) << " middlePoint=" << currentPoint << " rightPoint=" << (currentPoint - normal) << endl;
            }

            previousPoint = currentPoint;
        }

        drawDottedLine(visualizeStateSpaceImage, leftPoints, 50);
        drawDottedLine(visualizeStateSpaceImage, rightPoints, 75);
        drawDottedLine(visualizeStateSpaceImage, middlePoints, 100);

        list<Point2D>::iterator lastMiddle = middlePoints.begin();
        list<Point2D>::iterator lastLeft = leftPoints.begin();
        list<Point2D>::iterator lastRight = rightPoints.begin();
        list<Point2D>::iterator lastSelfIntersects = selfIntersects.begin();
        list<Point2D>::iterator currentSelfIntersects = ++(selfIntersects.begin());
        list<Point2D>::iterator currentMiddle = ++(middlePoints.begin());
        list<Point2D>::iterator currentLeft = ++(leftPoints.begin());
        list<Point2D>::iterator currentRight = ++(rightPoints.begin());
        for (; currentMiddle != middlePoints.end();
                ++lastMiddle, ++lastLeft, ++lastRight, ++currentMiddle, ++currentLeft, ++currentRight, ++lastSelfIntersects, ++currentSelfIntersects) {

            list<Point2D>::iterator backwardsSelfIntersects = lastSelfIntersects;
            list<Point2D>::iterator backwardsLeft = lastLeft;
            list<Point2D>::iterator backwardsMiddle = lastMiddle;
            //bool hasIntersect = false;
            while(segmentIntersect(*backwardsMiddle, *backwardsLeft, *currentMiddle, *currentLeft)) {
                //hasIntersect = true;
                *backwardsSelfIntersects = true;
                if (backwardsLeft == leftPoints.begin()) break;
                --backwardsLeft;
                --backwardsMiddle;
                --backwardsSelfIntersects;
            }
            //if (hasIntersect) {
            //    Point2D pointToUse = *backwardsLeft;
            //    for (++backwardsLeft, ++backwardsMiddle; backwardsLeft != currentLeft; ++backwardsLeft, ++backwardsMiddle) {
            //        *backwardsLeft = pointToUse;
            //        //int priorDistance = (*backwardsLeft - *backwardsMiddle).squaredMagnitude();
            //        //int newDistance = (*currentLeft - *backwardsMiddle).squaredMagnitude();
            //        //if (newDistance > priorDistance)
            //        //    *backwardsLeft = *currentLeft;
            //    }
            //    *currentLeft = pointToUse;
            //}

*/

/*
            if (segmentIntersect(*lastMiddle, *lastLeft, *currentMiddle, *currentLeft)) {
                //cout << "left intersect " << *lastMiddle << "->" << *lastLeft << "  " << *currentMiddle << "->" << *currentLeft << endl;

                //int currentDistance = (*currentLeft - *currentMiddle).squaredMagnitude();
                //int lastDistance = (*lastLeft - *currentMiddle).squaredMagnitude();
                //if (currentDistance < lastDistance) {
                    *currentLeft = *lastLeft;
                //}
            }

            //if (lastRight->x > 200 && lastRight->x < 220) cout << "right test " << *lastMiddle << "->" << *lastRight << "  " << *currentMiddle << "->" << *currentRight << endl;
            if (segmentIntersect(*lastMiddle, *lastRight, *currentMiddle, *currentRight)) {
                //cout << "right intersect " << *lastMiddle << "->" << *lastRight << "  " << *currentMiddle << "->" << *currentRight << endl;
                int currentDistance = (*currentRight - *currentMiddle).squaredMagnitude();
                int lastDistance = (*lastRight - *currentMiddle).squaredMagnitude();
                if (currentDistance < lastDistance) {
                    *currentRight = *lastRight;
                }
            }
*/
/*
        }
*/
/*
        list<Point2D>::iterator lastNonSelfIntersect = leftPoints.begin();
        currentSelfIntersects = selfIntersects.begin();
        currentLeft = leftPoints.begin();
        bool 
        for (; currentSelfIntersects != selfIntersects.end(); ++currentSelfIntersects, ++currentLeft) {
            if (!*currentSelfIntersects) {

        drawDottedLine(visualizeStateSpaceImage, leftPoints, 150);
        drawDottedLine(visualizeStateSpaceImage, rightPoints, 175);
        drawDottedLine(visualizeStateSpaceImage, middlePoints, 200);

        drawSegments(visualizeStateSpaceImage, middlePoints, leftPoints, 220);
        //miPixel pixels[2];
        //pixels[0] = NumericTraits<CostImagePixelType>::max();
        //pixels[1] = NumericTraits<CostImagePixelType>::max();
        //miGC *pGC = miNewGC(2, pixels);
        //miSetGCAttrib(pGC, MI_GC_LINE_WIDTH, stateSpaceWidth);
        //miSetGCMiterLimit(pGC, 1.2);
        //miPaintedSet *paintedSet = miNewPaintedSet();

        //int goodPoints = 0;
        //miPoint *leftParallelPoints = new miPoint[v->size()];
        //miPoint *rightParallelPoints = new miPoint[v->size()];
        //miPoint *centerPoints = new miPoint[v->size()];
*/


/*
        // Copy original point locations into originalPoints and mfEstimates
        //list<pair<bool, Point2D> >::iterator lastPoint = v->previous(v->end());
        list<pair<bool, Point2D> >::iterator lastPoint = v->end();
        --lastPoint;
        for (list<pair<bool, Point2D> >::iterator currentPoint = v->begin(); currentPoint != v->end(); ) {

        Point2D lastNormalPoint;
        Point2D lastInvNormalPoint;
        bool hasLastNormals = false;
*/
        list<pair<bool, Point2D> >::iterator lastPoint = v->end();
        --lastPoint;
        for (list<pair<bool, Point2D> >::iterator currentPoint = v->begin(); currentPoint != v->end(); ) {

            originalPoints.push_back(currentPoint->second);
            mfEstimates.push_back(currentPoint->second);

            vector<Point2D> *stateSpace = new vector<Point2D>();
            pointStateSpaces.push_back(stateSpace);

            vector<double> *stateProbabilities = new vector<double>();
            pointStateProbabilities.push_back(stateProbabilities);

            vector<int> *stateDistances = new vector<int>();
            pointStateDistances.push_back(stateDistances);

            if (!currentPoint->first) {
                // Point is not moveable.
                stateSpace->push_back(currentPoint->second);
                stateProbabilities->push_back(1.0);
                stateDistances->push_back(0);

                lastPoint = currentPoint;
                ++currentPoint;
                hasLastNormals = false;
            }
            else {
                Point2D lastPoint2D = lastPoint->second;
                Point2D currentPoint2D = currentPoint->second;
                lastPoint = currentPoint;
                ++currentPoint;
                Point2D nextPoint2D = (currentPoint == v->end()) ? v->begin()->second : currentPoint->second;

                // vp = vector from lastPoint to currentPoint
                Diff2D vp(currentPoint2D.x - lastPoint2D.x, currentPoint2D.y - lastPoint2D.y);
                // vn = vector from currentPoint to nextPoint
                Diff2D vn(nextPoint2D.x - currentPoint2D.x, nextPoint2D.y - currentPoint2D.y);
                // np = normal to vp
                Diff2D np(-vp.y, vp.x);
                // nn = normal to vn
                Diff2D nn(-vn.y, vn.x);

                // normal = normal vector at currentPoint
                // normal points to the left of vp and vn.
                Diff2D normal = np + nn;

                normal *= (stateSpaceWidth / normal.magnitude());
                Diff2D invNormal = -normal;

                if (

                hasLastNormals = true;
                lastNormalPoint = normal;
                lastInvNormalPoint = invNormal;
#if 0
                // Make polyline segment containing three points
                // First point: point halfway between lastPoint2D and currentPoint2D
                // Second point: currentPoint2D
                // Last point: point halfway between currentPoint2D and nextPoint2D

                // vp = vector from lastPoint to currentPoint
                Diff2D vp(currentPoint2D.x - lastPoint2D.x, currentPoint2D.y - lastPoint2D.y);
                Point2D previousHalf = lastPoint2D + (vp / 2);

                // vn = vector from currentPoint to nextPoint
                Diff2D vn(nextPoint2D.x - currentPoint2D.x, nextPoint2D.y - currentPoint2D.y);
                Point2D nextHalf = currentPoint2D + (vn / 2);

                miPoint points[3];
                points[0].x = previousHalf.x;
                points[0].y = previousHalf.y;
                points[1].x = currentPoint2D.x;
                points[1].y = currentPoint2D.y;
                points[2].x = nextHalf.x;
                points[2].y = nextHalf.y;

                // Rasterize polyline segment
                miClearPaintedSet(paintedSet);
                --(pixels[0]);
                --(pixels[1]);
                miSetGCPixels(pGC, 2, pixels);
                miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, 3, points);
if (pixels[1] == 208)
                copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));
                visualizeStateSpaceImage[currentPoint2D] = 255;
#endif

#if 0
                // Determine state space of currentPoint
                double stateSpaceWidth = costImageShortDimension / 3.0;

                // vp = vector from lastPoint to currentPoint
                Diff2D vp(currentPoint2D.x - lastPoint2D.x, currentPoint2D.y - lastPoint2D.y);

                // vn = vector from currentPoint to nextPoint
                Diff2D vn(nextPoint2D.x - currentPoint2D.x, nextPoint2D.y - currentPoint2D.y);

                // np = normal to vp
                Diff2D np(-vp.y, vp.x);

                // nn = normal to vn
                Diff2D nn(-vn.y, vn.x);

                // normal = normal vector at currentPoint
                // normal points to the left of vp and vn.
                Diff2D normal = np + nn;

                Diff2D normalScaled = normal * (stateSpaceWidth / normal.magnitude());
                Diff2D invNormalScaled = normal * (-1 * stateSpaceWidth / normal.magnitude());
                //// cosTheta = cosine of angle between vp and vn
                //double cosTheta = -1 * ((vp.x * vn.x) + (vp.y * vn.y)) / (vp.magnitude() * vn.magnitude());

                //// Acute junction distance = distance from currentPoint to closest point in wide polyline
                //// on the acute side
                //double ajd = stateSpaceWidth / sqrt((1 - cosTheta) / 2);
                //cout << "vp=" << vp << " vn=" << vn << " cosTheta=" << cosTheta << " ajd=" << ajd << endl;

                //// Check to see if normal is on the acute or obtuse side of the intersection between vp and vn.
                //// Take dot product of vn and np.
                //// If positive -> vn is to the left of vp -> normal is on acute side
                //// If negative -> vn is right of vp -> normal is on obtuse side
                //int vnDOTnp = (vn.x * np.x) + (vn.y * np.y);

                //Diff2D normalScaled = normal;
                //Diff2D invNormalScaled = -normal;

                //if (vnDOTnp > 0) {
                //    // Normal is acute
                //    // scale normal by ajd
                //    normalScaled *= ajd / normal.magnitude();
                //    // scale invNormal by stateSpaceWidth
                //    invNormalScaled *= stateSpaceWidth / normal.magnitude();
                //} else {
                //    // Normal is obtuse
                //    // scale normal by stateSpaceWidth
                //    normalScaled *= stateSpaceWidth / normal.magnitude();
                //    // scale invNormal by ajd
                //    invNormalScaled *= ajd / normal.magnitude();
                //}

                //cout << "currentPoint=" << currentPoint2D << " normalScaled=" << normalScaled << " invNormalScaled=" << invNormalScaled << endl;

                //// Find furthest good point along normalScaled
                //Diff2D normalFarPoint(currentPoint2D);
                //{
                //    Diff2D lineStart = Diff2D(currentPoint2D);
                //    Diff2D lineEnd = Diff2D(currentPoint2D) + normalScaled;
                //    LineIterator<Diff2D> lineIterator(lineStart, lineEnd);
                //    while (lineIterator != lineEnd) {
                //        if (costImage->isInside(*lineIterator) &&
                //            ((*costImage)[*lineIterator] != NumericTraits<CostImagePixelType>::max())) {
                //            normalFarPoint = *lineIterator;
                //            //visualizeStateSpaceImage[*lineIterator] = 125;
                //        }
                //        ++lineIterator;
                //    }
                //}
                ////visualizeStateSpaceImage[normalFarPoint] = 150;

                //// Find furthest good point along invNormalScaled
                //Diff2D invNormalFarPoint(currentPoint2D);
                //{
                //    Diff2D lineStart = Diff2D(currentPoint2D);
                //    Diff2D lineEnd = Diff2D(currentPoint2D) + invNormalScaled;
                //    LineIterator<Diff2D> lineIterator(lineStart, lineEnd);
                //    while (lineIterator != lineEnd) {
                //        if (costImage->isInside(*lineIterator) &&
                //            ((*costImage)[*lineIterator] != NumericTraits<CostImagePixelType>::max())) {
                //            invNormalFarPoint = *lineIterator;
                //            //visualizeStateSpaceImage[*lineIterator] = 175;
                //        }
                //        ++lineIterator;
                //    }
                //}
                ////visualizeStateSpaceImage[invNormalFarPoint] = 200;
                ////visualizeStateSpaceImage[currentPoint2D] = 255;

                Diff2D normalFarPoint = Diff2D(currentPoint2D) + normalScaled;
                Diff2D invNormalFarPoint = Diff2D(currentPoint2D) + invNormalScaled;
                leftParallelPoints[goodPoints].x = normalFarPoint.x;
                leftParallelPoints[goodPoints].y = normalFarPoint.y;
                rightParallelPoints[goodPoints].x = invNormalFarPoint.x;
                rightParallelPoints[goodPoints].y = invNormalFarPoint.y;
                centerPoints[goodPoints].x = currentPoint2D.x;
                centerPoints[goodPoints].y = currentPoint2D.y;
                ++goodPoints;
#endif

// Old state space enumeration
#if 1
                // Determine state space of currentPoint along normal vector
                // Normal points to the left of the vector
                Diff2D normal(lastPoint2D.y - nextPoint2D.y, nextPoint2D.x - lastPoint2D.x);
                normal *= costImageShortDimension / (3 * normal.magnitude());

                Diff2D lineStart = Diff2D(currentPoint2D) + normal;
                Diff2D lineEnd = Diff2D(currentPoint2D) - normal;
                LineIterator<Diff2D> lineBegin(lineStart, lineEnd);

                // Choose a reasonable number of state points along this line.
                int numberOfStatePoints = 32;
                int lineLength = std::max(std::abs(lineEnd.x - lineStart.x), std::abs(lineEnd.y - lineStart.y));
                int spaceBetweenPoints = lineLength / numberOfStatePoints;
                //while (lineLength > numberOfStatePoints) {
                //    ++spaceBetweenPoints;
                //    lineLength >>= 1;
                //}

                int pointNumber = 0;
                while (lineBegin != lineEnd) {
                    if ((pointNumber % spaceBetweenPoints) == 0) {
                        if (costImage->isInside(*lineBegin)) {
                            if ((*costImage)[*lineBegin] != NumericTraits<CostImagePixelType>::max()) {
                                Point2D linePoint(*lineBegin);
                                stateSpace->push_back(linePoint);
                                stateDistances->push_back((int)((currentPoint2D - linePoint).magnitude()) / 4);
                                visualizeStateSpaceImage[*lineBegin] = 100;
                            }
                        }
                    }
                    ++pointNumber;
                    ++lineBegin;
                }
#endif

                unsigned int localK = stateSpace->size();
                kMax = std::max(kMax, localK);

                for (unsigned int i = 0; i < localK; ++i) stateProbabilities->push_back(1.0 / localK);

                //if (localK > (unsigned int)numberOfStatePoints) {
                //    //int lineLength = std::max(std::abs(lineEnd.x - lineStart.x), std::abs(lineEnd.y - lineStart.y));
                //    cout << lineStart << " -> " << lineEnd << " = " << lineLength << " space=" << spaceBetweenPoints << " pointNumber=" << pointNumber << " localK=" << localK << endl;
                //}

            }

            convergedPoints.push_back(stateSpace->size() < 2);
        }
*/

#if 0
        pixels[0] = 0;
        pixels[1] = 1;
        miClearPaintedSet(paintedSet);
        miSetGCPixels(pGC, 2, pixels);
        miPixel pixels[2];
        pixels[0] = 200; //NumericTraits<CostImagePixelType>::max();
        pixels[1] = 200; //NumericTraits<CostImagePixelType>::max();
        miGC *pGC = miNewGC(2, pixels);
        miSetGCAttrib(pGC, MI_GC_LINE_WIDTH, stateSpaceWidth);
        miSetGCMiterLimit(pGC, 1.2);
        miPaintedSet *paintedSet = miNewPaintedSet();

        miPoint *points = new miPoint[v->size()];
        int goodPoints = 0;

        for (list<pair<bool, Point2D> >::iterator currentPoint = v->begin(); currentPoint != v->end(); ++currentPoint) {
            if (currentPoint->first) {
                points[goodPoints].x = currentPoint->second.x;
                points[goodPoints].y = currentPoint->second.y;
                ++goodPoints;
            }
        }

        miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, points);

        delete[] points;

        copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));

        miDeleteGC(pGC);
        miDeletePaintedSet(paintedSet);
#endif

#if 0
        miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, centerPoints);
        copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));
        miClearPaintedSet(paintedSet);
        pixels[0] = 200;
        pixels[1] = 200;
        miSetGCPixels(pGC, 2, pixels);
        miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, leftParallelPoints);
        copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));
        miClearPaintedSet(paintedSet);
        pixels[0] = 175;
        pixels[1] = 175;
        miSetGCPixels(pGC, 2, pixels);
        miDrawLines(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, rightParallelPoints);
        copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));
        miClearPaintedSet(paintedSet);
        pixels[0] = 150;
        pixels[1] = 150;
        miSetGCPixels(pGC, 2, pixels);
        miDrawPoints(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, leftParallelPoints);
        miDrawPoints(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, rightParallelPoints);
        miDrawPoints(paintedSet, pGC, MI_COORD_MODE_ORIGIN, goodPoints, centerPoints);
        copyPaintedSetToImage(destImageRange(visualizeStateSpaceImage), paintedSet, Diff2D(0,0));

        delete[] leftParallelPoints;
        delete[] rightParallelPoints;
        delete[] centerPoints;
        miDeleteGC(pGC);
        miDeletePaintedSet(paintedSet);
#endif

        ImageExportInfo visInfo("enblend_anneal_state_space.tif");
        exportImage(srcImageRange(visualizeStateSpaceImage), visInfo);

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
        for_each(pointStateDistances.begin(), pointStateDistances.end(), bind(delete_ptr(), _1));
    }

    void run() {
        tCurrent = tInitial;
        int numIterations = (int)ceil(log(tFinal/tInitial)/log(tau));
        while (tCurrent > tFinal) {
            double epsilon = 1.0 / kMax;
            unsigned int eta = (unsigned int)ceil(log(epsilon) / log(((kMax - 2.0) / (2.0 * kMax) * exp(-tCurrent / deltaEMax)) + 0.5));
            cout << "tCurrent=" << tCurrent << " eta=" << eta << " kMax=" << kMax;
            for (unsigned int i = 0; i < eta; i++) {
                iterate();
            }
            tCurrent *= tau;

            int numConvergedPoints = 0;
            for (unsigned int i = 0; i < convergedPoints.size(); i++) {
                if (convergedPoints[i]) numConvergedPoints++;
            }
            cout << " converged=" << numConvergedPoints << "/" << convergedPoints.size() << endl;
        }

        for (unsigned int i = 0; i < convergedPoints.size(); i++) {
            if (!convergedPoints[i]) {
                cout << "Unconverged original point=" << originalPoints[i] << endl;
                vector<Point2D> *stateSpace = pointStateSpaces[i];
                vector<double> *stateProbabilities = pointStateProbabilities[i];
                unsigned int localK = stateSpace->size();
                for (unsigned int state = 0; state < localK; ++state) {
                    cout << "    state " << (*stateSpace)[state] << " weight=" << (*stateProbabilities)[state] << endl;
                }
                cout << "    mfEstimate=" << mfEstimates[i] << endl;
            }
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
            cost += (currentPointEstimate - originalPoint).magnitude() / 4;
        }
        return cost;
    }

protected:

    void iterate() {
        int *E = new int[kMax];
        double *pi = new double[kMax];

        unsigned int lastIndex = originalPoints.size() - 1;
        for (unsigned int index = 0; index < originalPoints.size(); ++index) {
            // Skip updating points that have already converged.
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            vector<int> *stateDistances = pointStateDistances[index];
            unsigned int localK = stateSpace->size();

            unsigned int nextIndex = (index + 1) % originalPoints.size();
            Point2D originalPoint = originalPoints[index];
            Point2D lastPointEstimate = mfEstimates[lastIndex];
            bool lastPointInCostImage = costImage->isInside(lastPointEstimate);
            Point2D nextPointEstimate = mfEstimates[nextIndex];
            bool nextPointInCostImage = costImage->isInside(nextPointEstimate);
            lastIndex = index;

            // Calculate E values.
            // exp_a scaling factor is part of the Schraudolph approximation.
            double exp_a = (1048576 / M_LN2) / tCurrent;
            for (unsigned int i = 0; i < localK; ++i) {
                Point2D currentPoint = (*stateSpace)[i];
                E[i] = (*stateDistances)[i];
                if (lastPointInCostImage) E[i] += costImageCost(lastPointEstimate, currentPoint);
                if (nextPointInCostImage) E[i] += costImageCost(currentPoint, nextPointEstimate);
                E[i] = (int)(E[i] * exp_a);
                pi[i] = 0.0;
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

            for (unsigned int j = 0; j < localK; ++j) {
                double piTj = (*stateProbabilities)[j];
                pi[j] += piTj;
                for (unsigned int i = (j+1); i < localK; ++i) {
                    double piT = (*stateProbabilities)[i] + piTj;
                    eco.n.i = E[j] - E[i] + 1072693248 - 60801;
                    double piTAn = piT / (1 + eco.d);
                    pi[j] += piTAn;
                    pi[i] += piT - piTAn;
                }
                (*stateProbabilities)[j] = pi[j] / localK;
            }
        }

        delete[] E;
        delete[] pi;

        kMax = 1;
        // Make new mean field estimates.
        for (unsigned int index = 0; index < pointStateSpaces.size(); ++index) {
            if (convergedPoints[index]) continue;

            vector<Point2D> *stateSpace = pointStateSpaces[index];
            vector<double> *stateProbabilities = pointStateProbabilities[index];
            vector<int> *stateDistances = pointStateDistances[index];
            unsigned int localK = stateSpace->size();
            double estimateX = 0.0;
            double estimateY = 0.0;

            double totalWeight = 0.0;
            for (unsigned int k = 0; k < localK; ++k) {
                double weight = (*stateProbabilities)[k];
                totalWeight += weight;
                if (weight > 0.99) convergedPoints[index] = true;
                Point2D state = (*stateSpace)[k];
                estimateX += weight * (double)state.x;
                estimateY += weight * (double)state.y;
            }
            estimateX /= totalWeight;
            estimateY /= totalWeight;
			Point2D newEstimate(NumericTraits<int>::fromRealPromote(estimateX), NumericTraits<int>::fromRealPromote(estimateY));
            mfEstimates[index] = newEstimate;

            // Sanity check
            if (!costImage->isInside(newEstimate)) {
                cerr << "new mf estimate for point " << originalPoints[index] << " is outside cost image: " << newEstimate << endl;
                exit(-1);
            }

            // Remove improbable solutions from the search space
            for (unsigned int k = 0; k < stateSpace->size(); ) {
                double weight = (*stateProbabilities)[k];
                if (weight < 0.0001) {
                    // Replace this state with last state
                    (*stateProbabilities)[k] = (*stateProbabilities)[stateProbabilities->size() - 1];
                    (*stateSpace)[k] = (*stateSpace)[stateSpace->size() - 1];
                    (*stateDistances)[k] = (*stateDistances)[stateDistances->size() - 1];

                    // Delete last state
                    stateProbabilities->pop_back();
                    stateSpace->pop_back();
                    stateDistances->pop_back();
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

    // Original point locations
    vector<Point2D> originalPoints;

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

};

template <typename CostImage>
void annealSnake(const CostImage* const ci, slist<pair<bool, Point2D> > *snake) {

    list<pair<bool, Point2D> > listSnake;
    for (slist<pair<bool, Point2D> >::iterator p = snake->begin(); p != snake->end(); ++p) {
        listSnake.push_back(*p);
    }

    GDAConfiguration<CostImage> cfg(ci, &listSnake);
    //cout << "original cost = " << cfg.currentCost() << endl;
    //cfg.run();
    //cout << "final cost = " << cfg.currentCost() << endl;

    //slist<pair<bool, Point2D> >::iterator snakePoint = snake->begin();
    //vector<Point2D>::iterator annealedPoint = cfg.getCurrentPoints().begin();
    //for (; snakePoint != snake->end(); ++snakePoint, ++annealedPoint) {
    //    snakePoint->second = *annealedPoint;
    //}

};

} // namespace enblend

#endif /* __ANNEAL_H__ */
