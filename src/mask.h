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
#include "nearest.h"

#include "vigra/error.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/impexalpha.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra/transformimage.hxx"

using __gnu_cxx::hash_set;

using std::priority_queue;
using std::vector;

using vigra::BasicImage;
using vigra::BCFImage;
using vigra::BImage;
using vigra::BlueAccessor;
using vigra::BRGBCFImage;
using vigra::Diff2D;
using vigra::exportImage;
using vigra::FindMinMax;
using vigra::GreenAccessor;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::initImageIf;
using vigra::inspectImage;
using vigra::linearIntensityTransform;
using vigra::NumericTraits;
using vigra::transformImage;
using vigra::transformImageIf;
using vigra::UIImage;

using vigra::functor::Arg1;
using vigra::functor::Arg2;
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
template <typename AlphaType, typename MaskType>
MaskType *createMask(const AlphaType * const whiteAlpha,
        const AlphaType * const blackAlpha,
        const EnblendROI &iBB,
        const EnblendROI &uBB,
        const bool wraparound,
        EnblendROI &mBB) {

    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

    //// Read mask from a file instead of calculating it.
    //MaskType *fileMask = new MaskType(uBB.size());
    //ImageImportInfo fileMaskInfo("enblend_mask.tif");
    //importImage(fileMaskInfo, destImage(*fileMask));
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
    // maskTransform = 1, then mask = max value (white image)
    // maskTransform != 1, then mask = zero - (black image)
    transformImage(srcImageRange(*maskTransform),
            destImage(*mask),
            ifThenElse(Arg1() == Param(1),
                    Param(NumericTraits<MaskPixelType>::max()),
                    Param(NumericTraits<MaskPixelType>::zero())));

    delete maskTransform;
    // mem xsection = MaskType*ubb + UIImage*uBB

    // Calculate mask bounds.
    maskBounds(srcImageRange(*mask), uBB, mBB);

    EnblendROI iBB_uBB;
    iBB_uBB.setCorners(iBB.getUL() - uBB.getUL(), iBB.getLR() - uBB.getUL());
    cout << "iBB_uBB= (" << iBB_uBB.getUL().x << ", " << iBB_uBB.getUL().y << ")"
         << " -> "
         << "(" << iBB_uBB.getLR().x << ", " << iBB_uBB.getLR().y << ")" << endl;

    // Find points
    hash_set<Diff2D> *entryExitPoints =
            findTransitionLineEntryExitPoints(iBB_uBB.apply(srcImageRange(*mask)));

    // Debug: write mask cost function to output image file
    FindMinMax<UIImage::value_type> minmax;
    inspectImage(srcImageRange(*maskDistance), minmax);
    transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
            iBB_uBB.apply(destImage(*maskDistance)),
            linearRangeMapping(minmax.min, minmax.max, 1<<20, 0));
    BRGBCFImage *maskDistanceB = new BRGBCFImage(uBB.size());
    transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
            iBB_uBB.apply(destImage(*maskDistanceB, GreenAccessor<typename BRGBCFImage::value_type>())),
            linearRangeMapping(0, 1<<20, 0, 200));
    // Mark feature points as infinite cost (green=255)
    combineTwoImages(iBB_uBB.apply(srcImageRange(*maskDistance)),
                     iBB_uBB.apply(srcImage(*maskDistanceB, GreenAccessor<typename BRGBCFImage::value_type>())),
                     iBB_uBB.apply(destImage(*maskDistanceB, GreenAccessor<typename BRGBCFImage::value_type>())),
                     ifThenElse(Arg1() == Param((unsigned int)(1<<20)), Param((unsigned char)255), Arg2()));
    transformImage(iBB_uBB.apply(srcImageRange(*maskDistance)),
                     iBB_uBB.apply(destImage(*maskDistance)),
                    ifThenElse(Arg1() == Param((unsigned int)(1<<20)), Param(NumericTraits<unsigned int>::max()), Arg1()));
    // Mark entry/exit points in red
    for (hash_set<Diff2D>::iterator i = entryExitPoints->begin(); i != entryExitPoints->end(); ++i) {
        Diff2D point = *i + iBB_uBB.getUL();
        cout << "entry/exit point (" << point.x << ", " << point.y << ")" << endl;
        ((*maskDistanceB)[point]).setRed(0xFF);
    }

    dijkstra(iBB_uBB.apply(srcImageRange(*maskDistance)),
            iBB_uBB.apply(destImage(*maskDistanceB, BlueAccessor<typename BRGBCFImage::value_type>())),
            entryExitPoints);

    ImageExportInfo maskDistanceUSInfo("enblend_distance.tif");
    maskDistanceUSInfo.setPosition(uBB.getUL());
    exportImage(srcImageRange(*maskDistanceB), maskDistanceUSInfo);
    delete maskDistanceB;

    delete entryExitPoints;
    delete maskDistance;

    return mask;
};

// Find the points where the transition line enters/leaves the iBB.
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

    typedef typename CostAccessor::value_type CostType;
    typedef BasicImage<CostType> CostImageType;
    typedef priority_queue<Diff2D, vector<Diff2D>, dijkstra_compare<Diff2D, CostImageType> > PQ;

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
        CostImageType *costSoFar = new CostImageType(w, h, NumericTraits<CostType>::max());
        PQ *pq = new PQ(dijkstra_compare<Diff2D, CostImageType>(costSoFar));

        // Set costSoFar to cost function at entry point
        (*costSoFar)[entryExitPoint] = ca(cost_upperleft + entryExitPoint);
        pq->push(entryExitPoint);

        cout << "entry point (" << entryExitPoint.x << ", " << entryExitPoint.y << ")" << endl;

        while (!pq->empty()) {
            Diff2D top = pq->top();
            pq->pop();

            CostType costToTop = (*costSoFar)[top];

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
                    CostType neighborPreviousCost = (*costSoFar)[neighborPoint];
                    if (neighborPreviousCost != NumericTraits<CostType>::max()) continue;

                    CostType neighborCost = ca(cost_upperleft + neighborPoint);
                    if (neighborCost == NumericTraits<CostType>::max()) continue;

                    //cout << "visiting neighbor (" << neighborPoint.x << ", " << neighborPoint.y << ")" << endl;
                    //cout << "neighborCost=" << neighborCost << " neighborPreviousCost=" << neighborPreviousCost << endl;
                    // Calculate new cost to neighbor (with saturating arithmetic)
                    CostType newNeighborCost = neighborCost + costToTop;
                    if (newNeighborCost < neighborCost) { // wraparound occured.
                        newNeighborCost = NumericTraits<CostType>::max();
                    }
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
