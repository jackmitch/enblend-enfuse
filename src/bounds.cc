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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

// Union bounding box.
extern uint32 UBBFirstX;
extern uint32 UBBLastX;
extern uint32 UBBFirstY;
extern uint32 UBBLastY;

// Region of interest for pyramids.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

/** Caluclate the region of interest for pyramids and the number of
 *  levels we're going to make based on the blend mask.
 */
uint32 bounds(MaskPixel *mask) {
    // First calculate the bounding box of the transition line.
    int32 firstMulticolorRow = UBBLastY;
    int32 lastMulticolorRow = UBBFirstY;
    int32 firstMulticolorColumn = UBBLastX;
    int32 lastMulticolorColumn = UBBFirstX;

    // Scan rows.
    for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
        MaskPixel *maskP = &mask[y * OutputWidth + UBBFirstX];
        uint8 leftmostColor = maskP->r;
        for (uint32 x = UBBFirstX; x <= UBBLastX; x++) {
            if (maskP->r != leftmostColor) {
                firstMulticolorRow = min(firstMulticolorRow, (int32)y);
                lastMulticolorRow = max(lastMulticolorRow, (int32)y);
                break;
            }
            maskP++;
        }
    }

    // Scan columns.
    for (uint32 x = UBBFirstX; x <= UBBLastX; x++) {
        MaskPixel *maskP = &mask[UBBFirstY * OutputWidth + x];
        uint8 topmostColor = maskP->r;
        for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
            if (maskP->r != topmostColor) {
                firstMulticolorColumn = min(firstMulticolorColumn, (int32)x);
                lastMulticolorColumn = max(lastMulticolorColumn, (int32)x);
                break;
            }
            maskP += OutputWidth;
        }
    }

    if (Verbose > 0) {
        cout << "Transition line bounding box = ("
             << firstMulticolorColumn << ", "
             << firstMulticolorRow << ") -> ("
             << lastMulticolorColumn << ", "
             << lastMulticolorRow << ")"
             << endl;
    }

    // Caluclate how many levels we can create, and determine the size of
    // the region-of-interest on level 0 so that all of the pixels that
    // influence blending on the bottom level are included.
    uint32 maximumLevels = 1;
    while (true) {

        // filterHalfWidth = how many pixels the transition line spreads
        // out into layer maximumLevels-1, given the precision of the mask.
        // filterHalfWidth = how many pixels an image feature spreads
        // out into layer maximumLevels-1, given the precision of the image.
        // Add these up to find the spread of the bounding box that is
        // necessary to include all of the level 0 pixels that will influence
        // blending on level maximumLevels-1.
        // Since image precision = mask precision = 8 bpp, I'm just doubling
        // the value.
        int32 extent = 2 * filterHalfWidth(maximumLevels - 1, 255);

        ROIFirstX = max(firstMulticolorColumn - extent, (int32)UBBFirstX);
        ROILastX = min((int32)UBBLastX, lastMulticolorColumn + extent);
        ROIFirstY = max(firstMulticolorRow - extent, (int32)UBBFirstY);
        ROILastY = min((int32)UBBLastY, lastMulticolorRow + extent);

        //cout << "Level " << maximumLevels << " tl bb = ("
        //        << firstExtentColumn << ", "
        //        << firstExtentRow << ") -> ("
        //        << lastExtentColumn << ", "
        //        << lastExtentRow << ")"
        //        << endl;

        // caluclate short dimension.
        uint32 shortDimension = min(ROILastY - ROIFirstY + 1,
                ROILastX - ROIFirstX + 1);

        if ((shortDimension >> (maximumLevels - 1)) > 8) {
            // Try another level.
            maximumLevels++;
        }
        else {
            // quit - this is a good maximumLevels.
            break;
        }
    }

    if (Verbose > 0) {
        cout << "MaximumLevels = " << maximumLevels << endl;
        cout << "Region of Interest bounding box = ("
             << ROIFirstX << ", "
             << ROIFirstY << ") -> ("
             << ROILastX << ", "
             << ROILastY << ")"
             << endl;
    }

    if (maximumLevels == 1) {
        cerr << "enblend: intersection of images is too small to make "
             << "more than one pyramid level."
             << endl;
    }

    return maximumLevels;
}

/** Copy pixels from inside the UBB but outside the ROI from srcImage
 *  into dstImage. Ignore transparent pixels.
 */
void copyExcludedPixels(uint32 *dstImage, uint32 *srcImage) {

    for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
        for (uint32 x = UBBFirstX; x <= UBBLastX; x++) {
            if (y < ROIFirstY
                    || y > ROILastY
                    || x < ROIFirstX
                    || x > ROILastX) {
                uint32 index = (y * OutputWidth) + x;
                uint32 srcPixel = srcImage[index];
                if (TIFFGetA(srcPixel)) {
                    dstImage[index] = srcImage[index];
                }
            }
        }
    }

    return;
}

