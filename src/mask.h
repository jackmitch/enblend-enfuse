/*
 * Copyright (C) 2004-2005 Andrew Mihal
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
#ifndef __MASK_H__
#define __MASK_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <ext/slist>

#include "common.h"
//#include "anneal.h"
#include "nearest.h"
#include "path.h"

#include "vigra/contourcirculator.hxx"
#include "vigra/error.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra_ext/impexalpha.hxx"
#include "vigra_ext/XMIWrapper.h"

using std::make_pair;
using std::vector;
using __gnu_cxx::slist;

using vigra::combineThreeImages;
using vigra::combineTwoImages;
using vigra::CrackContourCirculator;
using vigra::exportImage;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::initImageIf;
using vigra::NumericTraits;
using vigra::Point2D;
using vigra::RGBToGrayAccessor;
using vigra::transformImage;
using vigra::transformImageIf;

using vigra::functor::Arg1;
using vigra::functor::Arg2;
using vigra::functor::Arg3;
using vigra::functor::ifThenElse;
using vigra::functor::Param;

using vigra_ext::copyPaintedSetToImage;

namespace enblend {

/** Calculate a blending mask between whiteImage and blackImage.
 */
template <typename ImageType, typename AlphaType, typename MaskType>
MaskType *createMask(ImageType *white,
        ImageType *black,
        AlphaType *whiteAlpha,
        AlphaType *blackAlpha,
        EnblendROI &uBB,
        EnblendROI &iBB,
        bool wraparound,
        EnblendROI &mBB) {

    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

    //// Read mask from a file instead of calculating it.
    //// Be sure to still calculate mBB below.
    //MaskType *fileMask = new MaskType(uBB.size());
    //ImageImportInfo fileMaskInfo("enblend_mask.tif");
    //importImage(fileMaskInfo, destImage(*fileMask));
    //MaskType *mask = fileMask;

    // uBB rounded up to multiple of 8 pixels in each direction
    Diff2D stride8_size(1 + ((uBB.size().x + 7) >> 3), 1 + ((uBB.size().y + 7) >> 3));

    // range of stride8 pixels that intersect uBB
    EnblendROI stride8_initBB;
    stride8_initBB.setCorners(Diff2D(0, 0), Diff2D(1 + (uBB.size().x >> 3), 1 + (uBB.size().y >> 3)));

    // Stride 8 mask
    MaskType *maskInit = new MaskType(stride8_size);
    // Mask initializer pixel values:
    // 0 = outside both black and white image, or inside both images.
    // 1 = inside white image only.
    // 255 = inside black image only.
    //MaskType *maskInit = new MaskType(uBB.size());
    // mem xsection = BImage*ubb

    // Set maskInit = 1 at all pixels where whiteImage contributes.
    //initImageIf(destImageRange(*maskInit),
    //        maskIter(whiteAlpha->upperLeft() + uBB.getUL()),
    //        NumericTraits<MaskPixelType>::one());
    initImageIf(stride8_initBB.apply(destImageRange(*maskInit)),
            stride(8, 8, maskIter(whiteAlpha->upperLeft() + uBB.getUL())),
            NumericTraits<MaskPixelType>::one());

    // maskInit = maskInit + 255 at all pixels where blackImage contributes.
    // if whiteImage also contributes, this wraps around to zero.
    //transformImageIf(srcImageRange(*maskInit),
    //        maskIter(blackAlpha->upperLeft() + uBB.getUL()),
    //        destImage(*maskInit),
    //        (Param(NumericTraits<MaskPixelType>::one()) * Arg1()) + Param(NumericTraits<MaskPixelType>::max()));
    transformImageIf(stride8_initBB.apply(srcImageRange(*maskInit)),
            stride(8, 8, maskIter(blackAlpha->upperLeft() + uBB.getUL())),
            stride8_initBB.apply(destImage(*maskInit)),
            (Param(NumericTraits<MaskPixelType>::one()) * Arg1()) + Param(NumericTraits<MaskPixelType>::max()));

    // Mask transform replaces 0 areas with either 1 or 255.
    //MaskType *maskTransform = new MaskType(uBB.size());
    //MaskType *maskTransform = new MaskType(stride8_size);
    MaskType *mask = new MaskType(stride8_size + Diff2D(2,2));
    // mem xsection = 2*BImage*ubb
    // ignore 1-pixel border around maskInit and maskTransform
    nearestFeatureTransform(wraparound,
            srcImageRange(*maskInit),
            destIter(mask->upperLeft() + Diff2D(1,1)));
    // mem xsection = 2*BImage*ubb + 2*UIImage*ubb

    delete maskInit;
    // mem xsection = BImage*ubb

    //MaskType *mask = new MaskType(uBB.size());
    // Add 1-pixel border all around so we can use the crackcontourcirculator
    //MaskType *mask = new MaskType(stride8_size + Diff2D(2,2));
    // mem xsection = BImage*ubb + MaskType*ubb

    // Dump maskTransform into mask
    // maskTransform = 1, then mask = max value (white image)
    // maskTransform != 1, then mask = zero - (black image)
    //transformImage(srcImageRange(*maskTransform),
    //        destIter(mask->upperLeft() + Diff2D(1,1)),
    //        ifThenElse(Arg1() == Param(NumericTraits<MaskPixelType>::one()),
    //                Param(NumericTraits<MaskPixelType>::max()),
    //                Param(NumericTraits<MaskPixelType>::zero())));

    //delete maskTransform;
    // mem xsection = MaskType*ubb

    // 0 = uninitialized border region
    // 1 = white image
    // 255 = black image
    // Vectorize white regions in mask
    int distance = 4;
    Point2D borderUL(1,1);
    Point2D borderLR(mask->width()-1, mask->height()-1);
    vector<slist<pair<bool, Point2D> > *> snakes;
    MaskIteratorType my = mask->upperLeft() + Diff2D(1,1);
    MaskIteratorType mend = mask->lowerRight() + Diff2D(-1, -1);
    for (int y = 1; my.y < mend.y; ++y, ++my.y) {
        MaskIteratorType mx = my;
        MaskPixelType lastColor = NumericTraits<MaskPixelType>::max();

        for (int x = 1; mx.x < mend.x; ++x, ++mx.x) {
            if ((*mx == NumericTraits<MaskPixelType>::one()) && (lastColor == NumericTraits<MaskPixelType>::max())) {

                // Found the corner of a previously unvisited white region.
                // Create a snake to hold the border of this region.
                vector<Point2D> excessPoints;
                slist<pair<bool, Point2D> > *snake = new slist<pair<bool, Point2D> >();
                snakes.push_back(snake);

                // Walk around border of white region.
                CrackContourCirculator<MaskIteratorType> crack(mx);
                CrackContourCirculator<MaskIteratorType> crackEnd(crack);
                bool lastPointFrozen = false;
                int distanceLastPoint = 0;
                do {
                    Point2D currentPoint = *crack + Diff2D(x,y);
                    crack++;
                    Point2D nextPoint = *crack + Diff2D(x,y);

                    // See if currentPoint lies on border.
                    if ((currentPoint.x == borderUL.x) || (currentPoint.x == borderLR.x)
                            || (currentPoint.y == borderUL.y) || (currentPoint.y == borderLR.y)) {

                        // See if currentPoint is in a corner.
                        if ((currentPoint.x == borderUL.x && currentPoint.y == borderUL.y)
                                || (currentPoint.x == borderUL.x && currentPoint.y == borderLR.y)
                                || (currentPoint.x == borderLR.x && currentPoint.y == borderUL.y)
                                || (currentPoint.x == borderLR.x && currentPoint.y == borderLR.y)) {
                            snake->push_front(make_pair(false, currentPoint));
                            distanceLastPoint = 0;
                        }
                        else if (!lastPointFrozen
                                || ((nextPoint.x != borderUL.x) && (nextPoint.x != borderLR.x)
                                        && (nextPoint.y != borderUL.y) && (nextPoint.y != borderLR.y))) {
                            snake->push_front(make_pair(false, currentPoint));
                            distanceLastPoint = 0;
                        }
                        lastPointFrozen = true;
                    }
                    else {
                        // Current point is not frozen.
                        // FIXME need to fix painting routine below if pixels are to be skipped.
                        if ((distanceLastPoint % distance) == 0) {
                            snake->push_front(make_pair(true, currentPoint));
                            distanceLastPoint = 0;
                        } else {
                            excessPoints.push_back(currentPoint);
                        }
                        lastPointFrozen = false;
                    }
                    distanceLastPoint++;
                } while (crack != crackEnd);

                // Paint the border so this region will not be found again
                for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
                         vertexIterator != snake->end(); ++vertexIterator) {
                    (*mask)[vertexIterator->second] = NumericTraits<MaskPixelType>::zero();
                }
                for (vector<Point2D>::iterator vertexIterator = excessPoints.begin();
                        vertexIterator != excessPoints.end(); ++vertexIterator) {
                    (*mask)[*vertexIterator] = NumericTraits<MaskPixelType>::zero();
                }
            }

            lastColor = *mx;
        }
    }
    
    // Fill mask with union region
    mask->init(NumericTraits<MaskPixelType>::zero());
    combineTwoImages(stride(8, 8, srcIterRange(whiteAlpha->upperLeft() + uBB.getUL(), whiteAlpha->upperLeft() + uBB.getLR())),
            stride(8, 8, maskIter(blackAlpha->upperLeft() + uBB.getUL())),
            destIter(mask->upperLeft() + Diff2D(1,1)),
            ifThenElse(Arg1() || Arg2(), Param(NumericTraits<MaskPixelType>::max()), Param(NumericTraits<MaskPixelType>::zero())));

    // Mark movable snake vertices (vertices inside union region)
    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;
        for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
                vertexIterator != snake->end(); ++vertexIterator) {
            if (vertexIterator->first &&
                    ((*mask)(vertexIterator->second.x, vertexIterator->second.y) == NumericTraits<MaskPixelType>::zero())) {
                vertexIterator->first = false;
            }
        }
    }

    //ImageExportInfo smallMaskInfo("enblend_small_mask.tif");
    //exportImage(srcImageRange(*mask), smallMaskInfo);
    delete mask;

    // Convert snake vertices to root-relative vertices
    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;
        for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
                vertexIterator != snake->end(); ++vertexIterator) {
            vertexIterator->second = uBB.getUL() + (8 * (vertexIterator->second + Diff2D(-1,-1)));
        }
    }

    // Find extent of moveable snake vertices, and vertices bordering moveable vertices
    int leftExtent = NumericTraits<int>::max();
    int rightExtent = NumericTraits<int>::min();
    int topExtent = NumericTraits<int>::max();
    int bottomExtent = NumericTraits<int>::min();
    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;
        slist<pair<bool, Point2D> >::iterator lastVertex = snake->previous(snake->end());
        //--lastVertex;
        for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
                vertexIterator != snake->end(); ++vertexIterator) {
            if (lastVertex->first || vertexIterator->first) {
                leftExtent = std::min(leftExtent, lastVertex->second.x);
                rightExtent = std::max(rightExtent, lastVertex->second.x);
                topExtent = std::min(topExtent, lastVertex->second.y);
                bottomExtent = std::max(bottomExtent, lastVertex->second.y);
                leftExtent = std::min(leftExtent, vertexIterator->second.x);
                rightExtent = std::max(rightExtent, vertexIterator->second.x);
                topExtent = std::min(topExtent, vertexIterator->second.y);
                bottomExtent = std::max(bottomExtent, vertexIterator->second.y);
            }
            lastVertex = vertexIterator;
        }
    }

    // Vertex bounding box
    EnblendROI vBB;
    vBB.setCorners(Point2D(leftExtent, topExtent), Point2D(rightExtent+1, bottomExtent+1));

    // Make sure that vertex bounding box is bigger than iBB by one pixel in each direction
    EnblendROI iBBPlus;
    iBBPlus.setCorners(iBB.getUL() + Diff2D(-1,-1), iBB.getLR() + Diff2D(1,1));
    vBB.unite(iBBPlus, vBB);

    // Vertex-Union bounding box: portion of uBB inside vBB.
    EnblendROI uvBB;
    vBB.intersect(uBB, uvBB);

    // Offset between vBB and uvBB
    Diff2D uvBBOffset = uvBB.getUL() - vBB.getUL();

    // Push ul corner of vBB so that there is an even number of pixels between vBB and uvBB.
    // This is necessary for striding by two over vBB.
    if (uvBBOffset.x % 2) vBB.setUpperLeft(vBB.getUL() + Diff2D(-1,0));
    if (uvBBOffset.y % 2) vBB.setUpperLeft(vBB.getUL() + Diff2D(0,-1));
    uvBBOffset = uvBB.getUL() - vBB.getUL();
    Diff2D uvBBOffsetHalf = uvBBOffset / 2;

    // Create stitch mismatch image
    typedef UInt8 MismatchImagePixelType;
    EnblendNumericTraits<MismatchImagePixelType>::ImageType mismatchImage((vBB.size() + Diff2D(1,1)) / 2,
            NumericTraits<MismatchImagePixelType>::zero());
    //cout << "mismatch image w=" << mismatchImage.width() << " h=" << mismatchImage.height() << endl;
    //cout << "uvBBOffset =" << uvBBOffset << " half=" << uvBBOffsetHalf << endl;
    combineTwoImages(stride(2, 2, uvBB.apply(srcImageRange(*white, RGBToGrayAccessor<ImagePixelType>()))),
                     stride(2, 2, uvBB.apply(srcImage(*black, RGBToGrayAccessor<ImagePixelType>()))),
                     destIter(mismatchImage.upperLeft() + uvBBOffsetHalf),
                     abs(Arg1() - Arg2()));
    transformImage(srcImageRange(mismatchImage), destImage(mismatchImage),
                     ifThenElse(Arg1() > Param(10), Arg1(), Param(NumericTraits<MismatchImagePixelType>::one())));
    // Areas where only one image contribute have maximum cost
    combineThreeImages(stride(2, 2, uvBB.apply(srcImageRange(*whiteAlpha))),
                     stride(2, 2, uvBB.apply(srcImage(*blackAlpha))),
                     srcIter(mismatchImage.upperLeft() + uvBBOffsetHalf),
                     destIter(mismatchImage.upperLeft() + uvBBOffsetHalf),
                     ifThenElse(Arg1() ^ Arg2(), Param(NumericTraits<MismatchImagePixelType>::max()), Arg3()));

    // Anneal snakes over mismatch image
    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;
    };

    // Use Dijkstra to route between moveable snake vertices over mismatchImage.
    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;
        slist<pair<bool, Point2D> >::iterator lastVertex = snake->previous(snake->end());
        //--lastVertex;
        //bool lastVertexInVBB = vBB.includes(lastVertex->second);

        for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
            vertexIterator != snake->end(); ++vertexIterator) {
            Point2D point = vertexIterator->second;
            //bool pointInVBB = vBB.includes(point);

            if (/*lastVertexInVBB && pointInVBB &&*/ (lastVertex->first || vertexIterator->first)) {
                // FIXME do dijkstra between these points
                // route from point to lastPoint
                // backtrack when solution is found from lastPoint to point
                // insert new vertices before vertexIterator

                // Move point relative to vBB stride 2
                //cout << "point=" << point;
                point = (point - vBB.getUL()) / 2;
                //cout << " moved=" << point << endl;
                //mismatchImage[point] = vertexIterator->first ? 128 : 200;

                // Last point relative to vBB stride 2
                Point2D lastPoint = (lastVertex->second - vBB.getUL()) / 2;
                //mismatchImage[lastPoint] = lastVertex->first ? 128 : 200;

                int radius = 25;
                EnblendROI pointSurround;
                pointSurround.setCorners(point - Diff2D(radius,radius), point + Diff2D(radius,radius));
                //cout << "pointSurroundFirst=" << pointSurround << endl;
                EnblendROI lastPointSurround;
                lastPointSurround.setCorners(lastPoint - Diff2D(radius,radius), lastPoint + Diff2D(radius,radius));
                //cout << "lastPointSurroundFirst=" << lastPointSurround << endl;
                pointSurround.unite(lastPointSurround, pointSurround);
                //cout << "union=" << pointSurround << endl;
                EnblendROI withinVBB;
                withinVBB.setCorners(Diff2D(0,0), (vBB.size() + Diff2D(1,1)) / 2);
                pointSurround.intersect(withinVBB, pointSurround);
                //cout << "intersect=" << pointSurround << endl;

                //cout << "vBB=" << vBB << endl;
                //cout << "pointSurround=" << pointSurround << endl;
                //cout << "point=" << point << " moved=" << (point-pointSurround.getUL()) << endl;
                //cout << "lastPoint=" << lastPoint << " moved=" << (lastPoint-pointSurround.getUL()) << endl;

                vector<Point2D> *shortPath = minCostPath(pointSurround.apply(srcImageRange(mismatchImage)),
                        point - pointSurround.getUL(),
                        lastPoint - pointSurround.getUL());

                for (vector<Point2D>::iterator shortPathPoint = shortPath->begin();
                        shortPathPoint != shortPath->end();
                        ++shortPathPoint) {
                    mismatchImage[*shortPathPoint + pointSurround.getUL()] = 50;
                    snake->insert_after(lastVertex, make_pair(false, ((*shortPathPoint + pointSurround.getUL()) * 2) + vBB.getUL()));
                }

                delete shortPath;

                mismatchImage[point] = vertexIterator->first ? 128 : 200;
                mismatchImage[lastPoint] = lastVertex->first ? 128 : 200;

            }
            //else if (pointInVBB) {
            //    point = (point - vBB.getUL()) / 2;
            //    mismatchImage[point] = 90;
            //}

            lastVertex = vertexIterator;
            //lastVertexInVBB = pointInVBB;
        }
    }

    ImageExportInfo mismatchInfo("enblend_mismatch.tif");
    exportImage(srcImageRange(mismatchImage), mismatchInfo);

    // Fill snakes on uBB-sized mask
    miPixel pixels[2];
    pixels[0] = NumericTraits<MaskPixelType>::max();
    pixels[1] = NumericTraits<MaskPixelType>::max();
    miGC *pGC = miNewGC(2, pixels);
    miPaintedSet *paintedSet = miNewPaintedSet();

    mask = new MaskType(uBB.size());

    for (vector<slist<pair<bool, Point2D> > *>::iterator snakeIterator = snakes.begin();
            snakeIterator != snakes.end(); ++snakeIterator) {
        slist<pair<bool, Point2D> > *snake = *snakeIterator;

        miPoint *points = new miPoint[snake->size()];

        int i = 0;
        for (slist<pair<bool, Point2D> >::iterator vertexIterator = snake->begin();
                vertexIterator != snake->end(); ++vertexIterator, ++i) {
            Point2D vertex = vertexIterator->second - uBB.getUL();
            points[i].x = vertex.x;
            points[i].y = vertex.y;
        }

        miFillPolygon(paintedSet, pGC, MI_SHAPE_GENERAL, MI_COORD_MODE_ORIGIN, snake->size(), points);

        delete[] points;
    }

    copyPaintedSetToImage(destImageRange(*mask), paintedSet, Diff2D(0,0));

    miDeleteGC(pGC);
    miDeletePaintedSet(paintedSet);

    // Find the bounding box of the mask transition line and put it in mBB.
    MaskIteratorType firstMulticolorColumn = mask->lowerRight();
    MaskIteratorType lastMulticolorColumn = mask->upperLeft();
    MaskIteratorType firstMulticolorRow = mask->lowerRight();
    MaskIteratorType lastMulticolorRow = mask->upperLeft();

    MaskIteratorType myPrev = mask->upperLeft();
    my = mask->upperLeft() + Diff2D(0,1);
    mend = mask->lowerRight();
    for (; my.y < mend.y; ++my.y, ++myPrev.y) {
        MaskIteratorType mxLeft = my;
        MaskIteratorType mx = my + Diff2D(1,0);
        MaskIteratorType mxUpLeft = myPrev;
        MaskIteratorType mxUp = myPrev + Diff2D(1,0);

        if (*mxUpLeft != *mxLeft) {
            // Transition line is between mxUpLeft and mxLeft.
            if (firstMulticolorRow.y > mxUpLeft.y) firstMulticolorRow = mxUpLeft;
            if (lastMulticolorRow.y < mxLeft.y) lastMulticolorRow = mxLeft;
        }

        for (; mx.x < mend.x; ++mx.x, ++mxLeft.x, ++mxUp.x) {
            if (*mxLeft != *mx || *mxUp != *mx) {
                // Transition line is between mxLeft and mx and between mx and mxUp
                if (firstMulticolorColumn.x > mxLeft.x) firstMulticolorColumn = mxLeft;
                if (lastMulticolorColumn.x < mx.x) lastMulticolorColumn = mx;
                if (firstMulticolorRow.y > mxUp.y) firstMulticolorRow = mxUp;
                if (lastMulticolorRow.y < mx.y) lastMulticolorRow = mx;
            }
        }
    }

    // Check that mBB is well-defined.
    if ((firstMulticolorColumn.x >= lastMulticolorColumn.x)
            || (firstMulticolorRow.y >= lastMulticolorRow.y)) {
        // No transition pixels were found in the mask at all.
        // This means that one image has no contribution.
        if (*(mask->upperLeft()) == NumericTraits<MaskPixelType>::zero()) {
            // If the mask is entirely black, then inspectOverlap should have caught this.
            // It should have said that the white image is redundant.
            vigra_fail("Mask is entirely black, but white image was not identified as redundant.");
        }
        else {
            // If the mask is entirely white, then the black image would have been identified
            // as redundant if black and white were swapped.
            // Set mBB to the full size of the mask.
            mBB.setCorners(uBB.getUL(), uBB.getLR());
            // Explain why the black image disappears completely.
            cerr << "enblend: the previous images are completely overlapped by the current images"
                 << endl;
        }
    }
    else {
        // Move mBB lower right corner out one pixel, per VIGRA convention.
        ++lastMulticolorColumn.x;
        ++lastMulticolorRow.y;

        // mBB is defined relative to the inputUnion origin.
        mBB.setCorners(
                uBB.getUL() + Diff2D(firstMulticolorColumn.x - mask->upperLeft().x,
                                     firstMulticolorRow.y - mask->upperLeft().y),
                uBB.getUL() + Diff2D(lastMulticolorColumn.x - mask->upperLeft().x,
                                     lastMulticolorRow.y - mask->upperLeft().y));
    }

    if (Verbose > VERBOSE_ROIBB_SIZE_MESSAGES) {
        cout << "Mask transition line bounding box: ("
             << mBB.getUL().x
             << ", "
             << mBB.getUL().y
             << ") -> ("
             << mBB.getLR().x
             << ", "
             << mBB.getLR().y
             << ")" << endl;
    }

    return mask;

}

} // namespace enblend

#endif /* __MASK_H__ */
