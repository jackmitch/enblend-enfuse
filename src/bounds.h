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
#ifndef __BOUNDS_H__
#define __BOUNDS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "common.h"
#include "pyramid.h"

using std::cerr;
using std::cout;
using std::endl;
using std::min;

namespace enblend {

/** Characterize the overlap between src1 and src2.
 *  The images may overlap completely (CompleteOverlap),
 *  partially (PartialOverlap),
 *  or not at all (NoOverlap).
 */
template <typename SrcImageIterator, typename SrcAccessor>
Overlap inspectOverlap(
        SrcImageIterator src1_upperleft,
        SrcImageIterator src1_lowerright,
        SrcAccessor s1a,
        SrcImageIterator src2_upperleft,
        SrcAccessor s2a) {

    SrcImageIterator s1y = src1_upperleft;
    SrcImageIterator s2y = src2_upperleft;
    SrcImageIterator send = src1_lowerright;

    bool foundOverlap = false;
    bool foundDistinctS2 = false;

    for (; s1y.y != send.y; ++s1y.y, ++s2y.y) {

        SrcImageIterator s1x = s1y;
        SrcImageIterator s2x = s2y;

        for (; s1x.x != send.x; ++s1x.x, ++s2x.x) {
            if (s1a(s1x) && s2a(s2x)) {
                foundOverlap = true;
            } else if (s2a(s2x)) {
                foundDistinctS2 = true;
            }
            if (foundOverlap && foundDistinctS2) {
                // If we have found a pixel where there is overlap,
                // and also a pixel where src2 alone contributes,
                // then we know it's PartialOverlap and we can quit.
                return PartialOverlap;
            }
        }
    }

    if (foundOverlap) {
        return CompleteOverlap;
    } else {
        return NoOverlap;
    }

};

// Argument object factory version.
template <typename SrcImageIterator, typename SrcAccessor>
Overlap inspectOverlap(
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src1,
        pair<SrcImageIterator, SrcAccessor> src2) {
    return inspectOverlap(src1.first, src1.second, src1.third,
                          src2.first, src2.second);
};

/** Determine the region-of-interest and number of blending levels to use,
 *  given the current mask-bounding-box and intersection-bounding-box.
 *  We also need to know if the image is a 360-degree pano so we can check
 *  for the case that the ROI wraps around the left and right edges.
 */
template <typename PyramidPixelType>
unsigned int roiBounds(const EnblendROI &inputUnion,
        const EnblendROI &iBB,
        const EnblendROI &mBB,
        const EnblendROI &uBB,
        EnblendROI &roiBB,
        bool wraparoundForMask) {

    unsigned int levels = 1;

    if (ExactLevels == 0) {
        // Estimate the number of blending levels to use based on the size of the iBB.
        // Assume the transition line runs approximately down the center of the iBB.
        // Choose a number of levels that makes the mask spread out to the edges
        // of the iBB.
        // Calculate short dimension of iBB.
        unsigned int shortDimension = min(iBB.size().x, iBB.size().y);
        while (levels < 30) {
            unsigned int extent = filterHalfWidth<PyramidPixelType>(levels + 1);
            if ((2 * extent) > shortDimension) {
                // levels + 1 is too many levels.
                break;
            }
            levels++;
        }

        if (levels == 1) {
            cerr << "enblend: overlap region is too small to make "
                 << "more than one pyramid level."
                 << endl;
        }

    } else {
        levels = ExactLevels;
    }

    unsigned int extent = filterHalfWidth<PyramidPixelType>(levels);
    Diff2D extentDiff(extent, extent);
    roiBB.setCorners(mBB.getUL() - extentDiff, mBB.getLR() + extentDiff);

    if (wraparoundForMask &&
            (roiBB.getUL().x < 0 || roiBB.getLR().x > uBB.getLR().x)) {
        // If the ROI goes off either edge of the uBB,
        // and the uBB is the full size of the output image,
        // and the wraparound flag is specified,
        // then make roiBB the full width of uBB.
        roiBB.setCorners(Diff2D(0, roiBB.getUL().y),
                Diff2D(uBB.getLR().x, roiBB.getLR().y));
    }

    // ROI must not be bigger than uBB.
    uBB.intersect(roiBB, roiBB);

    // Verify the number of levels based on the size of the ROI.
    unsigned int roiShortDimension = min(roiBB.size().x, roiBB.size().y);
    unsigned int allowableLevels;
    for (allowableLevels = 1; allowableLevels < levels; allowableLevels++) {
        if (roiShortDimension <= 8) {
            // ROI dimensions preclude using more levels than allowableLevels.
            break;
        }
        roiShortDimension = (roiShortDimension + 1) >> 1;
    }

    if (allowableLevels < ExactLevels) {
        cerr << "enblend: image geometry precludes using more than "
             << allowableLevels
             << " levels." << endl;
    }

    if (Verbose > VERBOSE_NUMLEVELS_MESSAGES) {
        cout << "Using " << allowableLevels << " blending levels" << endl;
    }
    if (Verbose > VERBOSE_ROIBB_SIZE_MESSAGES) {
        cout << "Region of Interest bounding box: ("
             << roiBB.getUL().x
             << ", "
             << roiBB.getUL().y
             << ") -> ("
             << roiBB.getLR().x
             << ", "
             << roiBB.getLR().y
             << ")" << endl;
    }

    return allowableLevels;
};

/** Find the bounding box of the mask transition line and put it in mBB.
 */
template <typename MaskImageIterator, typename MaskAccessor>
void maskBounds(
        MaskImageIterator mask_upperleft,
        MaskImageIterator mask_lowerright,
        MaskAccessor ma,
        const EnblendROI &uBB,
        EnblendROI &mBB) {

    MaskImageIterator firstMulticolorColumn = mask_lowerright;
    MaskImageIterator lastMulticolorColumn = mask_upperleft;
    MaskImageIterator firstMulticolorRow = mask_lowerright;
    MaskImageIterator lastMulticolorRow = mask_upperleft;

    MaskImageIterator myPrev = mask_upperleft;
    MaskImageIterator my = mask_upperleft + Diff2D(0,1);
    MaskImageIterator mend = mask_lowerright;
    for (; my.y != mend.y; ++my.y, ++myPrev.y) {
        MaskImageIterator mxLeft = my;
        MaskImageIterator mx = my + Diff2D(1,0);
        MaskImageIterator mxUpLeft = myPrev;
        MaskImageIterator mxUp = myPrev + Diff2D(1,0);

        if (ma(mxUpLeft) != ma(mxLeft)) {
            // Transition line is between mxUpLeft and mxLeft.
            if (firstMulticolorRow.y > mxUpLeft.y) firstMulticolorRow = mxUpLeft;
            if (lastMulticolorRow.y < mxLeft.y) lastMulticolorRow = mxLeft;
        }

        for (; mx.x != mend.x; ++mx.x, ++mxLeft.x, ++mxUp.x) {
            if (ma(mxLeft) != ma(mx) || ma(mxUp) != ma(mx)) {
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
        vigra_fail("Mask transition line bounding box undefined.");
    }

    // Move mBB lower right corner out one pixel, per VIGRA convention.
    ++lastMulticolorColumn.x;
    ++lastMulticolorRow.y;

    // mBB is defined relative to the inputUnion origin.
    mBB.setCorners(
            uBB.getUL() + Diff2D(firstMulticolorColumn.x - mask_upperleft.x,
                                 firstMulticolorRow.y - mask_upperleft.y),
            uBB.getUL() + Diff2D(lastMulticolorColumn.x - mask_upperleft.x,
                                 lastMulticolorRow.y - mask_upperleft.y));

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

    return;
};

// Version with argument object factories.
template <typename MaskImageIterator, typename MaskAccessor>
inline void maskBounds(
        triple<MaskImageIterator, MaskImageIterator, MaskAccessor> mask,
        const EnblendROI &uBB,
        EnblendROI &mBB) {

    maskBounds(mask.first, mask.second, mask.third, uBB, mBB);

};

} // namespace enblend

#endif /* __BOUNDS_H__ */
