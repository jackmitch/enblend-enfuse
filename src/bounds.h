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
 *  given the current intersection-bounding-box.
 *  We also need to know if the image is a 360-degree pano so we can check
 *  for the case that the ROI wraps around the left and right edges.
 */
template <typename PyramidPixelType>
unsigned int roiBounds(const EnblendROI &inputUnion,
        const EnblendROI &iBB,
        const EnblendROI &uBB,
        EnblendROI &roiBB) {

    // Calculate how many levels we can create, and determine the size of
    // the region-of-interest on level 0 os that all of the pixels that
    // influence blending on the bottom level are included.
    unsigned int levels = 1;
    while (true) {

        // filterHalfWidth = how many pixels the transition line spreads
        // out into layer levels-1, given the precision of the mask.
        // filterHalfWidth = how many pixels an image feature spreads out into
        // layer levels-1, given the precision of the image.
        // Add these up to find the spread of the bounding box that is
        // necessary to include all of the level 0 pixels that will influence
        // blending on level levels-1.
        // Here I am assuming that the mask uses the same data type as a single
        // channel of the RGB pyramid.
        unsigned int extent = 2 * filterHalfWidth<PyramidPixelType>(levels - 1);

        // Set roiBB equal to iBB plus extent pixels in every direction.
        // I am conservatively estimating that the transition line touches
        // opposing corners of the iBB.
        // In some cases the transition line may touch adjacent corners of
        // the iBB. In those cases this code overestimates the size of the ROI.
        Diff2D extentDiff(extent, extent);
        roiBB.setCorners(iBB.getUL() - extentDiff, iBB.getLR() + extentDiff);

        // Check for wraparound condition.
        if (Wraparound
                && uBB.size().x == inputUnion.size().x
                && ((uBB.getUL().x > roiBB.getUL().x)
                    || (roiBB.getLR().x > uBB.getLR().x))) {
            // If the image is a 360 pano,
            // and the images in this step wrap around,
            // and the roi extends beyond the left or right edge of the image,
            // Then make the roi the full union bounding box size.
            roiBB = uBB;
        }
        else {
            // roiBB = uBB intersect roiBB
            // This guarantees that roiBB is no bigger than uBB.
            uBB.intersect(roiBB, roiBB);
        }

        if (levels == MaximumLevels) {
            // Hit the user-specified level limit.
            break;
        }
        else {
            // Calculate short dimension of ROI.
            unsigned int shortDimension = min(roiBB.size().x, roiBB.size().y);

            for (unsigned int i = 1; i < levels; i++) {
                shortDimension = (shortDimension + 1) >> 1;
            }

            if (shortDimension > 8) {
                // Try another level.
                levels++;
            }
            else {
                // quit - this is a good number of levels.
                break;
            }
        }
    }

    if (Verbose > VERBOSE_NUMLEVELS_MESSAGES) {
        cout << "Using " << levels << " blending levels" << endl;
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

    if (levels == 1 && levels != MaximumLevels) {
        cerr << "enblend: intersection of images is too small to make "
             << "more than one pyramid level."
             << endl;
    }

    return levels;
}

} // namespace enblend

#endif /* __BOUNDS_H__ */
