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
 * along with Foobar; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <tiffio.h>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;
extern double StitchMismatchThreshold;

// Region of interest for this operation.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

/** Calculate a blending mask between outputBuf and inputTIFF.
 */
uint32 *createMask(uint32 *outputBuf, TIFF *inputTIFF) {

    // Allocate memory for the blending mask.
    uint32 *mask = (uint32*)_TIFFmalloc(
            OutputWidth * OutputHeight * sizeof(uint32));
    if (mask == NULL) {
        cerr << "enblend: malloc failed for mask" << endl;
        exit(1);
    }

    // Region of interest boundaries
    ROIFirstX = OutputWidth - 1;
    ROILastX = 0;
    ROIFirstY = OutputHeight - 1;
    ROILastY = 0;

    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFReadScanline(inputTIFF,
                &(mask[i * OutputWidth]),
                i,
                8);

        for (uint32 j = 0; j < OutputWidth; j++) {
            uint32 outPixel = outputBuf[(i * OutputWidth) + j];
            uint32 inPixel = mask[(i * OutputWidth) + j];

            if (outPixel != 0 || inPixel != 0) {
                // Pixel is in region of interest.
                ROIFirstX = min(j, ROIFirstX);
                ROILastX = max(j, ROILastX);
                ROIFirstY = min(i, ROIFirstY);
                ROILastY = max(i, ROILastY);
            }

            if (outPixel == 0 && inPixel == 0) {
                // Pixel is not in the union of the two images.
                // Make the mask pixel black.
                mask[(i * OutputWidth) + j] = 0x00000001;
            }
            else if (outPixel != 0 && inPixel == 0) {
                // Pixel is in out only.
                // Make the mask pixel blue.
                mask[(i * OutputWidth) + j] = 0xFF00FF00;
            }
            else if (outPixel == 0 && inPixel != 0) {
                // Pixel is in input only.
                // Make the mask pixel green.
                mask[(i * OutputWidth) + j] = 0xFFFF0000;
            }
            else {
                // Pixel is in the intersection of the two images.

                /* FIXME stitch mismatch avoidance is broken.
                // Calculate the stitch mismatch at this pixel.
                double rDiff = abs((int32)TIFFGetR(outPixel) - (int32)TIFFGetR(inPixel))
                        / 255.0;
                double bDiff = abs((int32)TIFFGetB(outPixel) - (int32)TIFFGetB(inPixel))
                        / 255.0;
                double gDiff = abs((int32)TIFFGetG(outPixel) - (int32)TIFFGetG(inPixel))
                        / 255.0;

                if (rDiff > StitchMismatchThreshold
                        || bDiff > StitchMismatchThreshold
                        || gDiff > StitchMismatchThreshold) {
                    // Make the mask pixel red.
                    mask[(i * OutputWidth) + j] = RED;
                } else {
                    // Make the mask pixel black.
                    mask[(i * OutputWidth) + j] = BLACK;
                }
                */

                // Make the mask pixel black.
                mask[(i * OutputWidth) + j] = 0xFF000001;

            }
        }
    }

    // Run the thinning transform on the mask inside the ROI.
    // This will replace the black pixels with either green or blue
    // based on how close each pixel is to a green or blue region.
    thinMask(mask);

    // Swap blue for white and green for black.
    // Only needs to be done within region of interest.
    // Outside ROI make mask transparent.
    for (uint32 i = 0; i < (OutputWidth * OutputHeight); i++) {
        //uint32 x = i % OutputWidth;
        //uint32 y = i / OutputWidth;
        //if (x < ROIFirstX || x > ROILastX
        //        || y < ROIFirstY || y > ROILastY) {
        //    mask[i] = TRANS;
        //}
        //else if (mask[i] == BLUE) {
        //    mask[i] = WHITE;
        //}
        //else if (mask[i] == GREEN) {
        //    mask[i] = TRANS;
        //}
        if (mask[i] & 0x0000FF00) {
            // White - pixel is near out
            mask[i] |= 0x00FFFFFF;
        } else if (mask[i] & 0x00FF0000) {
            // Black - pixel is near in
            mask[i] &= 0xFF000000;
        } else {
            // pixel was not remarked. Keep trans info, set to out - white.
            mask[i] |= 0x00FFFFFF;
        }
    }

    return mask;
}

