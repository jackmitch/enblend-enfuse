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
#ifndef __MASK_H__
#define __MASK_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <vector>
#include <queue>
#include <ext/hash_set>

#include "common.h"
#include "anneal.h"
#include "nearest.h"

#include "vigra/contourcirculator.hxx"
#include "vigra/error.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/imageiteratoradapter.hxx"
#include "vigra/impex.hxx"
#include "vigra/impexalpha.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra/transformimage.hxx"
#include "vgl/vgl_polygon_scan_iterator.h"

using __gnu_cxx::hash_set;

using std::make_pair;
using std::priority_queue;
using std::vector;

using vigra::BasicImage;
using vigra::BCFImage;
using vigra::BImage;
using vigra::BlueAccessor;
using vigra::BRGBCFImage;
using vigra::combineThreeImages;
using vigra::combineTwoImages;
using vigra::CrackContourCirculator;
using vigra::Diff2D;
using vigra::exportImage;
using vigra::FindMinMax;
using vigra::GreenAccessor;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::initImageBorder;
using vigra::initImageIf;
using vigra::inspectImage;
using vigra::linearIntensityTransform;
using vigra::LineIterator;
using vigra::NumericTraits;
using vigra::Point2D;
using vigra::RedAccessor;
using vigra::RGBToGrayAccessor;
using vigra::RGBValue;
using vigra::transformImage;
using vigra::transformImageIf;
using vigra::UIImage;

using vigra::functor::Arg1;
using vigra::functor::Arg2;
using vigra::functor::Arg3;
using vigra::functor::ifThenElse;
using vigra::functor::Param;

// Hash function to support hash_set<Diff2D>
namespace __gnu_cxx {
template <>
class hash<Diff2D> {
public:
    size_t operator()(const Diff2D& x) const {
        return (size_t)(x.x ^ x.y);
    }
};
}

namespace enblend {

/** Calculate a blending mask between whiteImage and blackImage.
 */
template <typename ImageType, typename AlphaType, typename MaskType>
MaskType *createMask(const pair<const ImageType*, const AlphaType*> whitePair, /*AlphaType * const whiteAlpha,*/
        const pair<const ImageType*, const AlphaType*> blackPair, /*const AlphaType * const blackAlpha,*/
        const EnblendROI &iBB,
        const EnblendROI &uBB,
        const bool wraparound,
        EnblendROI &mBB) {

    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

    const ImageType *whiteImage = whitePair.first;
    const AlphaType *whiteAlpha = whitePair.second;
    const ImageType *blackImage = blackPair.first;
    const AlphaType *blackAlpha = blackPair.second;

    //// Read mask from a file instead of calculating it.
    //MaskType *fileMask = new MaskType(uBB.size());
    //ImageImportInfo fileMaskInfo("enblend_mask.tif");
    //importImage(fileMaskInfo, destImage(*fileMask));
    //maskBounds(srcImageRange(*fileMask), uBB, mBB);
    //return fileMask;

    // Mask initializer pixel values:
    // 0 = outside both black and white image, or inside both images.
    // 1 = inside white image only.
    // 255 = inside black image only.
    #ifdef ENBLEND_CACHE_IMAGES
    BCFImage *maskInit = new BCFImage(uBB.size());
    #else
    BImage *maskInit = new BImage(uBB.size());
    #endif
    // mem xsection = BImage*ubb

    // Set maskInit = 1 at all pixels where whiteImage contributes.
    initImageIf(destImageRange(*maskInit),
            maskIter(whiteAlpha->upperLeft() + uBB.getUL()),
            1);

    // maskInit = maskInit + 255 at all pixels where blackImage contributes.
    // if whiteImage also contributes, this wraps around to zero.
    transformImageIf(srcImageRange(*maskInit),
            maskIter(blackAlpha->upperLeft() + uBB.getUL()),
            destImage(*maskInit),
            linearIntensityTransform(1, 255));

    // Mask transform replaces 0 areas with either 1 or 255.
    #ifdef ENBLEND_CACHE_IMAGES
    BCFImage *maskTransform = new BCFImage(uBB.size());
    // FIXME should this be a cf image?
    UIImage *maskDistance = new UIImage(uBB.size());
    #else
    BImage *maskTransform = new BImage(uBB.size());
    UIImage *maskDistance = new UIImage(uBB.size());
    #endif
    // mem xsection = 2*BImage*ubb
    nearestFeatureTransform(wraparound,
            srcImageRange(*maskInit),
            destImage(*maskTransform),
            destImage(*maskDistance));
    // mem xsection = 2*BImage*ubb + 2*UIImage*ubb

    delete maskInit;
    // mem xsection = BImage*ubb + UIImage*ubb

    MaskType *mask = new MaskType(uBB.size());
    // mem xsection = BImage*ubb + MaskType*ubb + UIImage*ubb

    // Dump maskTransform into mask
    // FIXME replace this with white flood fill of snakes polygons.
    // maskTransform = 1, then mask = max value (white image)
    // maskTransform != 1, then mask = zero - (black image)
    transformImage(srcImageRange(*maskTransform),
            destImage(*mask),
            ifThenElse(Arg1() == Param(1),
                    Param(NumericTraits<MaskPixelType>::max()),
                    Param(NumericTraits<MaskPixelType>::zero())));

    delete maskTransform;
    // mem xsection = MaskType*ubb + UIImage*uBB

    EnblendROI iBB_uBB;
    iBB_uBB.setCorners(iBB.getUL() - uBB.getUL(), iBB.getLR() - uBB.getUL());
    cout << "iBB_uBB= (" << iBB_uBB.getUL().x << ", " << iBB_uBB.getUL().y << ")"
         << " -> "
         << "(" << iBB_uBB.getLR().x << ", " << iBB_uBB.getLR().y << ")" << endl;

    // Map distances to costs. Max distance = cost 0. Min distance (feature) = cost 1<<20
    UIImage *distanceCostImage = new UIImage(iBB.size().x + 2, iBB.size().y + 2, NumericTraits<unsigned int>::max());
    //// FIXME nearestFeatureTransform should return max distance
    //// then we can remove inspectImage call
    //FindMinMax<UIImage::value_type> minmax;
    //inspectImage(iBB_uBB.apply(srcImageRange(*maskDistance)), minmax);
    //// FIXME combine these two calls using one functor
    //transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
    //        destIter(distanceCostImage->upperLeft() + Diff2D(1,1)),
    //        linearRangeMapping(minmax.min, minmax.max, 1<<20, 0));
    //transformImage(srcImageRange(*distanceCostImage),
    //               destImage(*distanceCostImage),
    //               ifThenElse(Arg1() >= Param((unsigned int)(1<<20)),
    //                          Param(NumericTraits<unsigned int>::max()),
    //                          Arg1()));
    // FIXME just have nearestFeatureTransform put the data here to start with.
    copyImage(iBB_uBB.apply(srcImageRange(*maskDistance)), destIter(distanceCostImage->upperLeft() + Diff2D(1,1)));
    delete maskDistance;

    // Map stitch mismatches to costs.
    BImage *stitchCostImage = new BImage(iBB.size().x + 2, iBB.size().y + 2, NumericTraits<short>::max());
    combineTwoImages(iBB.apply(srcImageRange(*whiteImage, RGBToGrayAccessor<typename ImageType::value_type>())),
            iBB.apply(srcImage(*blackImage, RGBToGrayAccessor<typename ImageType::value_type>())),
            destIter(stitchCostImage->upperLeft() + Diff2D(1,1)),
            abs(Arg1() - Arg2()));

    //// Step N: min distance (feature) = 1<<20 mapped to infinite cost

    //// Debug: describe cost function as RGB image for visualization.
    //BRGBCFImage *maskDistanceB = new BRGBCFImage(uBB.size());

    //transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
    //        iBB_uBB.apply(destImage(*maskDistanceB, BlueAccessor<typename BRGBCFImage::value_type>())),
    //        linearRangeMapping(0, 1<<20, 0, 200));

    //// Step 2: map pixel discrepancies to costs.
    //combineThreeImages(iBB.apply(srcImageRange(*whiteImage, RGBToGrayAccessor<typename ImageType::value_type>())),
    //        iBB.apply(srcImage(*blackImage, RGBToGrayAccessor<typename ImageType::value_type>())),
    //        iBB_uBB.apply(srcImage(*maskDistance)),
    //        iBB_uBB.apply(destImage(*maskDistance)),
    //        Arg3() + (Param(1<<25) * (abs(Arg1() - Arg2()) / Param(255))));

    //// Map costs to 0-255 pixel range.
    //transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
    //        iBB_uBB.apply(destImage(*maskDistanceB, GreenAccessor<typename BRGBCFImage::value_type>())),
    //        linearRangeMapping(0, 1<<25, 0, 255));

    //// Debug: infinite cost points mapped to 255 red pixel value.
    //combineTwoImages(iBB_uBB.apply(srcImageRange(*maskDistance)),
    //                 iBB_uBB.apply(srcImage(*maskDistanceB)),
    //                 iBB_uBB.apply(destImage(*maskDistanceB)),
    //                 ifThenElse(Arg1() == Param(NumericTraits<unsigned int>::max()),
    //                            Param(RGBValue<unsigned char>(255,0,0)),
    //                            Arg2()));

    //// Solve cost function
    //dijkstra(iBB_uBB.apply(srcImageRange(*maskDistance)),
    //        iBB_uBB.apply(destImage(*maskDistanceB, GreenAccessor<typename BRGBCFImage::value_type>())),
    //        entryExitPoints);

    // Mask crack contains iBB surrounded by a 1-pixel black border.
    MaskType *maskCrack = new MaskType(iBB.size() + Diff2D(2, 2));
    copyImage(iBB_uBB.apply(srcImageRange(*mask)), destIter(maskCrack->upperLeft() + Diff2D(1, 1)));

    vector<vector<pair<bool, Point2D> > *> snakes;
    MaskIteratorType my = maskCrack->upperLeft() + Diff2D(1,1);
    MaskIteratorType mbegin = my;
    MaskIteratorType mend = maskCrack->lowerRight() + Diff2D(-1, -1);
    for (; my.y != mend.y; ++my.y) {
        MaskIteratorType mx = my;
        MaskPixelType lastColor = NumericTraits<MaskPixelType>::zero();
        for (; mx.x != mend.x; ++mx.x) {
            if ((*mx == NumericTraits<MaskPixelType>::max())
                    && (lastColor == NumericTraits<MaskPixelType>::zero())) {
                // Previously un-visited white region.
                // Create a snake
                vector<pair<bool, Point2D> > *snake = new vector<pair<bool, Point2D> >();
                snakes.push_back(snake);

                // Iterate over copy of image, modifying original
                MaskType *maskCrackCopy = new MaskType(*maskCrack);
                MaskIteratorType ci = maskCrackCopy->upperLeft() + (mx - maskCrack->upperLeft());
                CrackContourCirculator<MaskIteratorType> crack(ci);
                CrackContourCirculator<MaskIteratorType> crackEnd(crack);
                int distanceLastPoint = 0;
                bool lastPointFrozen = false;
                do {
                    CrackContourCirculator<MaskIteratorType> crackThis = crack;
                    CrackContourCirculator<MaskIteratorType> crackNext = ++crack;
                    MaskIteratorType crackPoint = mx + *crackThis;
                    MaskIteratorType crackNextPoint = mx + *crackNext;

                    // Paint crackPoint.
                    *crackPoint = NumericTraits<MaskPixelType>::max() / 2;

                    //cout << "crackThis=(" << (*crackThis).x << ", " << (*crackThis).y << ") "
                    //     << "crackNext=(" << (*crackNext).x << ", " << (*crackNext).y << ") ";

                    // Check if current point is frozen.
                    if ((crackPoint.x == mbegin.x)
                            || (crackPoint.x == mend.x)
                            || (crackPoint.y == mbegin.y)
                            || (crackPoint.y == mend.y)) {
                        if ((crackPoint.x == mbegin.x && crackPoint.y == mbegin.y)
                                || (crackPoint.x == mend.x && crackPoint.y == mbegin.y)
                                || (crackPoint.x == mbegin.x && crackPoint.y == mend.y)
                                || (crackPoint.x == mend.x && crackPoint.y == mend.y)) {
                            // Crack point is frozen on a corner.
                            // Make a frozen point here.
                            snake->push_back(make_pair(true, *crackThis + (mx - maskCrack->upperLeft())));
                            distanceLastPoint = 0;
                            //cout << "frozen reason 1" << endl;
                        }
                        else {
                            // It is frozen on an edge.
                            // This is a frozen point if the next point is not frozen.
                            // or if the last point is not frozen.
                            if (!lastPointFrozen
                                    || ((crackNextPoint.x != mbegin.x)
                                        && (crackNextPoint.x != mend.x)
                                        && (crackNextPoint.y != mbegin.y)
                                        && (crackNextPoint.y != mend.y))) {
                                snake->push_back(make_pair(true, *crackThis + (mx - maskCrack->upperLeft())));
                                distanceLastPoint = 0;
                                //cout << "frozen reason 2" << endl;
                            } else {
                                //cout << "frozen but not marked" << endl;
                            }
                        }
                        lastPointFrozen = true;
                    }
                    else {
                        // Current point is not frozen.
                        if ((distanceLastPoint % 10) == 0) {
                            // Make non-frozen point if we are far enough away from last non-frozen point.
                            snake->push_back(make_pair(false, *crackThis + (mx - maskCrack->upperLeft())));
                            distanceLastPoint = 0;
                            //cout << "unfrozen" << endl;
                        } else {
                            //cout << "unfrozen not marked" << endl;
                        }
                        lastPointFrozen = false;
                    }
                    distanceLastPoint++;
                } while (crack != crackEnd);

                delete maskCrackCopy;
            }
            lastColor = *mx;
        }
    }


    BRGBCFImage *maskCrackVisualize = new BRGBCFImage(iBB.size() + Diff2D(2, 2));
    //transformImage(srcImageRange(*maskCrack),
    //        destImage(*maskCrackVisualize, RedAccessor<typename BRGBCFImage::value_type>()),
    //        linearRangeMapping(NumericTraits<MaskPixelType>::min(),
    //                           NumericTraits<MaskPixelType>::max(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::min(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::max()));
    //transformImage(srcImageRange(*maskCrack),
    //        destImage(*maskCrackVisualize, GreenAccessor<typename BRGBCFImage::value_type>()),
    //        linearRangeMapping(NumericTraits<MaskPixelType>::min(),
    //                           NumericTraits<MaskPixelType>::max(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::min(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::max()));
    //transformImage(srcImageRange(*maskCrack),
    //        destImage(*maskCrackVisualize, BlueAccessor<typename BRGBCFImage::value_type>()),
    //        linearRangeMapping(NumericTraits<MaskPixelType>::min(),
    //                           NumericTraits<MaskPixelType>::max(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::min(),
    //                           NumericTraits<typename BRGBCFImage::value_type::value_type>::max()));
    transformImage(srcImageRange(*stitchCostImage),
            destImage(*maskCrackVisualize, GreenAccessor<BRGBCFImage::value_type>()),
            Arg1());
    combineTwoImages(srcImageRange(*distanceCostImage),
            srcImage(*maskCrackVisualize, RedAccessor<BRGBCFImage::value_type>()),
            destImage(*maskCrackVisualize, RedAccessor<BRGBCFImage::value_type>()),
            //ifThenElse(Arg1() == Param(NumericTraits<unsigned int>::max()),
            ifThenElse(Arg1() == Param(NumericTraits<unsigned int>::zero()),
                    Param(0xff), Arg2()));

    // Visualize snakes
    for (unsigned int i = 0; i < snakes.size(); i++) {
        vector<pair<bool, Point2D> > *snake = snakes[i];
        annealSnake<UIImage, BImage>(distanceCostImage, stitchCostImage, snake);
        for (unsigned int j = 0; j < snake->size(); j++) {
            unsigned int next = (j+1) % snake->size();
            pair<bool, Point2D> pointA = (*snake)[j];
            pair<bool, Point2D> pointB = (*snake)[next];
            LineIterator<typename BRGBCFImage::traverser> line(maskCrackVisualize->upperLeft() + pointA.second,
                    maskCrackVisualize->upperLeft() + pointB.second);
            LineIterator<typename BRGBCFImage::traverser> end(maskCrackVisualize->upperLeft() + pointB.second,
                    maskCrackVisualize->upperLeft() + pointB.second);
            for (; line != end; ++line) {
                *line = RGBValue<unsigned char>(0, 0, 0xff);
            }
        }
        for (unsigned int j = 0; j < snake->size(); j++) {
            pair<bool, Point2D> point = (*snake)[j];
            (*maskCrackVisualize)[point.second] = point.first ? RGBValue<unsigned char>(0xff, 0, 0) : RGBValue<unsigned char>(0, 0xff, 0);
        }
    }

    ImageExportInfo maskCrackInfo("enblend_mask_crack.tif");
    exportImage(srcImageRange(*maskCrackVisualize), maskCrackInfo);

    delete distanceCostImage;
    delete stitchCostImage;
    delete maskCrackVisualize;
    delete maskCrack;

    initImage(iBB_uBB.apply(destImageRange(*mask)), NumericTraits<MaskPixelType>::zero());
    vgl_polygon snakesPoly;
    for (unsigned int i = 0; i < snakes.size(); i++) {
        // Flood fill each snake white.
        vgl_polygon_sheet sheet;
        vector<pair<bool, Point2D> > *snake = snakes[i];
        for (unsigned int j = 0; j < snake->size(); j++) {
            sheet.push_back((*snake)[j].second);
        }
        snakesPoly.push_back(sheet);
    }
    //cout << "snakesPoly.size()=" << snakesPoly.size() << endl;
    //for (unsigned int i = 0; i < snakesPoly.size(); i++) {
    //    cout << "snakesPoly[" << i << "].size()=" << snakesPoly[i].size() << endl;
    //}
    // Iterator that points to 0,0 of snake coord space = mask->upperLeft() + iBB_uBB.getUL();
    vgl_polygon_scan_iterator<MaskIteratorType> fill(mask->upperLeft() + iBB_uBB.getUL() + Diff2D(-1,-1), snakesPoly);
    vgl_polygon_scan_iterator<MaskIteratorType> fillEnd(mask->upperLeft() + iBB_uBB.getUL() + Diff2D(-1,-1), snakesPoly);
    do {
        *fill = NumericTraits<MaskPixelType>::max();
    } while (++fill != fillEnd);
    //for (unsigned int i = 0; i < 100; i++) {
    //    *fill = NumericTraits<MaskPixelType>::max();
    //    ++fill;
    //}

    // Delete all snakes.
    for (unsigned int i = 0; i < snakes.size(); i++) {
        delete snakes[i];
    }

    //// Debug: combine cost visualization with mask.
    //BRGBCFImage *visualization = new BRGBCFImage(uBB.size());
    //transformImage(srcImageRange(*mask), destImage(*visualization),
    //        ifThenElse(Arg1() == Param(0),
    //                   Param(RGBValue<unsigned char>(0)),
    //                   Param(RGBValue<unsigned char>(255))));
    //initImageBorder(destIterRange(visualization->upperLeft() + iBB_uBB.getUL() + Diff2D(-1, -1),
    //                             visualization->upperLeft() + iBB_uBB.getLR() + Diff2D(1, 1)),
    //                1, RGBValue<unsigned char>(0,0,0));
    //copyImage(iBB_uBB.apply(srcImageRange(*maskDistanceB)),
    //        iBB_uBB.apply(destImage(*visualization)));

    //Diff2D point = *entryExitPointsCopy.begin() + iBB_uBB.getUL();
    //CrackContourCirculator<typename BRGBCFImage::traverser> crack(visualization->upperLeft() + point);
    //CrackContourCirculator<typename BRGBCFImage::traverser> crackend(crack);
    //vector<Diff2D> snakePoints;
    //unsigned int everyTenthPoint = 0;
    //do {
    //    //hash_set<Diff2D>::iterator topLocation = entryExitPointsCopy.find(*crack + point - iBB_uBB.getUL());
    //    Diff2D crack_iBB_uBB = *crack + point - iBB_uBB.getUL();
    //    if (topLocation == entryExitPointsCopy.end()) {
    //        if ((everyTenthPoint % 10) == 0) snakePoints.push_back(*crack);
    //        everyTenthPoint++;
    //    } else {
    //        // Special unmovable point
    //        cout << "found special unmovable point." << endl;
    //        snakePoints.push_back(*crack);
    //        everyTenthPoint = 1;
    //    }
    //    //FIXME check if point is on iBB boundary. - make it unmovable.
    //} while (++crack != crackend);
    //for (unsigned int i = 0; i < snakePoints.size(); i++) {
    //    (*visualization)[point + snakePoints[i]] = RGBValue<unsigned char>(0xff, 0, 0);
    //}

    // Debug: mark entry/exit points in yellow
    //for (hash_set<Diff2D>::iterator i = entryExitPointsCopy.begin(); i != entryExitPointsCopy.end(); ++i) {
    //    Diff2D point = *i + iBB_uBB.getUL();
    //    cout << "entry/exit point (" << point.x << ", " << point.y << ")" << endl;
    //    (*visualization)[point] = RGBValue<unsigned char>(0xff, 0xff, 0x0);
    //}

    //// Debug: write visualization out as tif
    //ImageExportInfo visualizationInfo("enblend_distance.tif");
    //visualizationInfo.setPosition(uBB.getUL());
    //exportImage(srcImageRange(*visualization), visualizationInfo);
    ////delete maskDistanceB;
    //delete visualization;

    
    //delete entryExitPoints;

    // Calculate mask bounds.
    maskBounds(srcImageRange(*mask), uBB, mBB);

    return mask;
};

// Find the (white) points where the transition line enters/leaves the iBB.
// FIXME only white points are allowed. - these will become anchor points of the snake.
template <class MaskImageIterator, class MaskAccessor>
hash_set<Diff2D> *findTransitionLineEntryExitPoints(
        MaskImageIterator mask_upperleft,
        MaskImageIterator mask_lowerright,
        MaskAccessor ma) {
    // FIXME what if transition line is a closed curve
    // FIXME what if transition line is a combination of closed curves and lines
    typedef typename MaskAccessor::value_type MaskPixelType;

    hash_set<Diff2D> *points = new hash_set<Diff2D>();

    // First try top row.
    MaskImageIterator mbegin = mask_upperleft;
    MaskImageIterator mx = mbegin;
    Diff2D dx;
    MaskPixelType lastColor = ma(mx);
    MaskImageIterator mend = mask_lowerright;
    for (; mx.x != mend.x; ++mx.x, ++dx.x) {
        if (ma(mx) != lastColor) {
            points->insert(dx);
            lastColor = ma(mx);
        }
    }
    // Right column.
    for (--mx.x, --dx.x; mx.y != mend.y; ++mx.y, ++dx.y) {
        if (ma(mx) != lastColor) {
            points->insert(dx);
            lastColor = ma(mx);
        }
    }
    // Bottom row.
    for (--mx.y, --dx.y; mx.x >= mbegin.x; --mx.x, --dx.x) {
        if (ma(mx) != lastColor) {
            points->insert(dx);
            lastColor = ma(mx);
        }
    }
    // Left column.
    for (++mx.x, ++dx.x; mx.y >= mbegin.y; --mx.y, --dx.y) {
        if (ma(mx) != lastColor) {
            points->insert(dx);
            lastColor = ma(mx);
        }
    }

    return points;
};

// Version with argument object factories
template <class MaskImageIterator, class MaskAccessor>
inline hash_set<Diff2D> *findTransitionLineEntryExitPoints(
        triple<MaskImageIterator, MaskImageIterator, MaskAccessor> mask) {
    return findTransitionLineEntryExitPoints(mask.first, mask.second, mask.third);
};

template <typename Point, typename Image>
class dijkstra_compare {
public:
    dijkstra_compare(Image *i) : image(i) {}
    bool operator()(const Point &a, const Point &b) {
        //cout << "comparing a=(" << a.x << ", " << a.y << ")=" << (*image)[a]
        //     << " b=(" << b.x << ", " << b.y << ")=" << (*image)[b] << endl;
        // want the priority queue sorted in ascending order.
        return ((*image)[a] > (*image)[b]);
    }
protected:
    Image *image;
};

// Find max cost path between each pair of points.
template <class CostImageIterator, class CostAccessor,
          class DestPathImageIterator, class DestPathImageAccessor>
void dijkstra(CostImageIterator cost_upperleft,
        CostImageIterator cost_lowerright,
        CostAccessor ca,
        DestPathImageIterator dest_upperleft,
        DestPathImageAccessor da,
        hash_set<Diff2D> *entryExitPoints) {

    typedef typename CostAccessor::value_type InputCostType;
    typedef float WorkingCostType;
    typedef BasicImage<WorkingCostType> WorkingCostImageType;
    typedef priority_queue<Diff2D, vector<Diff2D>, dijkstra_compare<Diff2D, WorkingCostImageType> > PQ;

    int w = cost_lowerright.x - cost_upperleft.x;
    int h = cost_lowerright.y - cost_upperleft.y;

    while (!entryExitPoints->empty()) {
        Diff2D entryExitPoint = *(entryExitPoints->begin());
        entryExitPoints->erase(entryExitPoints->begin());

        // 4-bit direction encoding {up, down, left, right}
        //  A  8  9
        //  2  0  1
        //  6  4  5
        const unsigned char neighborArray[] = {0xA, 8, 9, 1, 5, 4, 6, 2};
        const unsigned char neighborArrayInverse[] = {5, 4, 6, 2, 0xA, 8, 9, 1};
        BImage *pathNextHop = new BImage(w, h);
        WorkingCostImageType *costSoFar = new WorkingCostImageType(w, h, NumericTraits<WorkingCostType>::max());
        PQ *pq = new PQ(dijkstra_compare<Diff2D, WorkingCostImageType>(costSoFar));

        // Set costSoFar to cost function at entry point
        (*costSoFar)[entryExitPoint] = ca(cost_upperleft + entryExitPoint);
        pq->push(entryExitPoint);

        cout << "entry point (" << entryExitPoint.x << ", " << entryExitPoint.y << ")" << endl;

        while (!pq->empty()) {
            Diff2D top = pq->top();
            pq->pop();

            WorkingCostType costToTop = (*costSoFar)[top];

            //cout << "visiting top (" << top.x << ", " << top.y << ")"
            //     << " costToTop=" << costToTop << endl;
            //if (pq->empty()) cout << "pq is now empty." << endl;

            // Check if top is an exit point
            hash_set<Diff2D>::iterator topLocation = entryExitPoints->find(top);
            if (topLocation == entryExitPoints->end()) {
                // else for each 8-neighbor of top with costSoFar==0 do relax.
                for (int i = 0; i < 8; i++) {
                    // Get the neighbor
                    unsigned char neighborDirection = neighborArray[i];
                    Diff2D neighborPoint = top;
                    if (neighborDirection & 0x8) --neighborPoint.y;
                    if (neighborDirection & 0x4) ++neighborPoint.y;
                    if (neighborDirection & 0x2) --neighborPoint.x;
                    if (neighborDirection & 0x1) ++neighborPoint.x;

                    // Make sure neighbor is in valid region
                    if (neighborPoint.y < 0) continue;
                    if (neighborPoint.y >= h) continue;
                    if (neighborPoint.x < 0) continue;
                    if (neighborPoint.x >= w) continue;

                    // See if the neighbor has already been visited.
                    // If neighbor has maximal cost, it has not been visited.
                    // If so skip it.
                    WorkingCostType neighborPreviousCost = (*costSoFar)[neighborPoint];
                    if (neighborPreviousCost != NumericTraits<WorkingCostType>::max()) continue;

                    WorkingCostType neighborCost = ca(cost_upperleft + neighborPoint);
                    if (neighborCost == NumericTraits<InputCostType>::max()) continue;

                    //cout << "visiting neighbor (" << neighborPoint.x << ", " << neighborPoint.y << ")" << endl;
                    //cout << "neighborCost=" << neighborCost << " neighborPreviousCost=" << neighborPreviousCost << endl;
                    // Calculate new cost to neighbor (with saturating arithmetic)
                    WorkingCostType newNeighborCost = neighborCost + costToTop;
                    //if (newNeighborCost < neighborCost) { // wraparound occured.
                    //    newNeighborCost = NumericTraits<CostType>::max();
                    //}
                    if (newNeighborCost < neighborPreviousCost) {
                        // We have found the shortest path to neighbor.
                        (*costSoFar)[neighborPoint] = newNeighborCost;
                        (*pathNextHop)[neighborPoint] = neighborArrayInverse[i];
                        pq->push(neighborPoint);
                        //cout << "pushed neighbor new top=(" << pq->top().x << ", " << pq->top().y << ")" << endl;
                    }
                }

                //     if cost of neighbor > 0
                //         calculate new costSoFar for neighbor
                //         set pathNextHop for neighbor to top
                //         push neighbor on pq
                //     else
                //         neighbor is out of bounds.
            }
            else {
                // If yes then follow back to beginning using pathNextHop
                //     remove top from points
                cout << "exit point (" << top.x << ", " << top.y << ")" << endl;
                entryExitPoints->erase(topLocation);
                //     do {
                //         draw top in dest
                //         get pathNextHop for top
                //         move top to pathNextHop
                //     } (while pathNextHop != 0)
                //     break;
                unsigned char nextHop = 0;
                do {
                    da.set(NumericTraits<typename DestPathImageAccessor::value_type>::max(),
                            dest_upperleft + top);
                    nextHop = (*pathNextHop)[top];
                    if (nextHop & 0x8) --top.y;
                    if (nextHop & 0x4) ++top.y;
                    if (nextHop & 0x2) --top.x;
                    if (nextHop & 0x1) ++top.x;
                } while (nextHop != 0);
                break;
            }
        }

        delete pathNextHop;
        delete costSoFar;
        delete pq;
    }

};

// Version with argument object factories
template <class CostImageIterator, class CostAccessor,
          class DestPathImageIterator, class DestPathImageAccessor>
inline void dijkstra(triple<CostImageIterator, CostImageIterator, CostAccessor> cost,
        pair<DestPathImageIterator, DestPathImageAccessor> dest,
        hash_set<Diff2D> *entryExitPoints) {

    dijkstra(cost.first, cost.second, cost.third,
            dest.first, dest.second,
            entryExitPoints);

};

} // namespace enblend

#endif /* __MASK_H__ */
