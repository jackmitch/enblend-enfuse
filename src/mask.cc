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

/** Calculate a blending mask between whiteImage and blackImage.
 *  Sets the union bounding box of whiteImage and blackImage.
 */
MaskPixel *createMask(uint32 *whiteImage, uint32 *blackImage) {

    // Allocate memory for the blending mask.
    MaskPixel *mask = (MaskPixel*)calloc(OutputWidth * OutputHeight,
            sizeof(MaskPixel));
    if (mask == NULL) {
        cerr << "enblend: malloc failed in mask for mask" << endl;
        exit(1);
    }

    // Region of interest boundaries
    UBBFirstX = OutputWidth - 1;
    UBBLastX = 0;
    UBBFirstY = OutputHeight - 1;
    UBBLastY = 0;

    uint32 *whitePixel = whiteImage;
    uint32 *blackPixel = blackImage;
    MaskPixel *maskPixel = mask;

    for (uint32 i = 0; i < (OutputWidth * OutputHeight); i++) {

        if (*whitePixel != 0 || *blackPixel != 0) {
            // Pixel is in region of interest.
            uint32 x = i % OutputWidth;
            uint32 y = i / OutputWidth;
            UBBFirstX = min(x, UBBFirstX);
            UBBLastX = max(x, UBBLastX);
            UBBFirstY = min(y, UBBFirstY);
            UBBLastY = max(y, UBBLastY);
        }

        if (*whitePixel == 0 && *blackPixel == 0) {
            // Pixel is not in the union of the two images.
            // Mark the pixel as thinnable.
            maskPixel->r = 1;
        }
        else if (*whitePixel == 0) {
            // Pixel is in blackImage but not whiteImage.
            // Make the pixel green.
            maskPixel->g = 255;
            maskPixel->a = 255;
        }
        else if (*blackPixel == 0) {
            // Pixel is in whiteImage but not blackImage.
            // Make the pixel blue.
            maskPixel->b = 255;
            maskPixel->a = 255;
        }
        else {
            // Pixel is in both images.
            // Mark the pixel as thinnable.
            maskPixel->r = 1;
            maskPixel->a = 255;
        }

        whitePixel++;
        blackPixel++;
        maskPixel++;

    }

    // Run the nearest feature transform on the mask inside the UBB.
    // This will replace the thinnable pixels with either green or blue
    // based on how close each pixel is to a green or blue region.
    nearestFeatureTransform(mask);

    // Remark all blue pixels as white. These pixels are closer to
    // whiteImage than blackImage.
    // Remark all green pixels as black. These pixels are closer to
    // blackImage than whiteImage.
    // Remark remaining thinnable pixels as white.
    // Do not change alpha channel - this stores the union of
    // whiteImage and blackImage.
    maskPixel = mask;
    for (uint32 i = 0; i < (OutputWidth * OutputHeight); i++) {
        if (maskPixel->g != 0) {
            maskPixel->r = 0;
            maskPixel->g = 0;
            maskPixel->b = 0;
        }
        else {
            maskPixel->r = 255;
            maskPixel->g = 255;
            maskPixel->b = 255;
        }

        maskPixel++;
    }

    return mask;
}

