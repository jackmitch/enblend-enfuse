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

using std::cout;
using std::cerr;
using std::endl;
using std::min;

namespace enblend {

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

        if (levels == (unsigned int)MaximumLevels) {
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

    if (Verbose > 0) {
        cout << "Using " << levels << " blending levels" << endl;
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

    if (levels == 1 && levels != (unsigned int)MaximumLevels) {
        cerr << "enblend: intersection of images is too small to make "
             << "more than one pyramid level."
             << endl;
    }

    return levels;
}

} // namespace enblend

#endif /* __BOUNDS_H__ */
