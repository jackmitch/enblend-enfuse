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
FILE *createMask(FILE *whiteImageFile, FILE *blackImageFile) {

    uint32 ubbWidth = UBBLastX - UBBFirstX + 1;
    uint32 ubbHeight = UBBLastY - UBBFirstY + 1;

    // Allocate memory for the blending mask.
    MaskPixel *mask = (MaskPixel*)calloc(ubbWidth * ubbHeight,
            sizeof(MaskPixel));
    if (mask == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for mask)" << endl;
        exit(1);
    }

    uint32 *whiteLine = (uint32*)malloc(ubbWidth * sizeof(uint32));
    if (whiteLine == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for whiteLine)" << endl;
        exit(1);
    }

    uint32 *blackLine = (uint32*)malloc(ubbWidth * sizeof(uint32));
    if (blackLine == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for blackLine)" << endl;
        exit(1);
    }

    MaskPixel *maskPixel = mask;
    for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
        // Move to the pixel at (UBBFirstX, y)
        fseek(whiteImageFile, (y * OutputWidth + UBBFirstX) * sizeof(uint32), SEEK_SET);
        fseek(blackImageFile, (y * OutputWidth + UBBFirstX) * sizeof(uint32), SEEK_SET);

        // Read lines from white and black images.
        readFromTmpfile(whiteLine, sizeof(uint32), ubbWidth, whiteImageFile);
        readFromTmpfile(blackLine, sizeof(uint32), ubbWidth, blackImageFile);

        uint32 *whitePixel = whiteLine;
        uint32 *blackPixel = blackLine;
        for (uint32 x = 0; x < ubbWidth; x++) {

            bool isInWhiteImage = (TIFFGetA(*whitePixel) == 255);
            bool isInBlackImage = (TIFFGetA(*blackPixel) == 255);

            if (!isInWhiteImage && !isInBlackImage) {
                // Pixel is not in the union of the two images.
                // Mark the pixel as not a feature pixel.
                maskPixel->r = 1;
            }
            else if (isInBlackImage && !isInWhiteImage) {
                // Pixel is in blackImage but not whiteImage.
                // Make the pixel green.
                maskPixel->g = 255;
                maskPixel->a = 255;
            }
            else if (isInWhiteImage && !isInBlackImage) {
                // Pixel is in whiteImage but not blackImage.
                // Make the pixel blue.
                maskPixel->b = 255;
                maskPixel->a = 255;
            }
            else {
                // Pixel is in both images.
                // Mark the pixel as not a feature pixel.
                maskPixel->r = 1;
                maskPixel->a = 255;
            }

            whitePixel++;
            blackPixel++;
            maskPixel++;
        }
    }

    free(whiteLine);
    free(blackLine);

    // Run the nearest feature transform on the mask inside the UBB.
    // This will replace the not-a-feature pixels with either green or blue
    // based on how close each pixel is to a green or blue region.
    nearestFeatureTransform(mask);

    // Remark all blue pixels as white. These pixels are closer to
    // whiteImage than blackImage.
    // Remark all green pixels as black. These pixels are closer to
    // blackImage than whiteImage.
    // Remark remaining thinnable pixels as white.
    // Do not change alpha channel - this stores the outline of the union of
    // whiteImage and blackImage.
    maskPixel = mask;
    for (uint32 i = 0; i < (ubbWidth * ubbHeight); i++) {
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

    //saveMaskAsTIFF(mask);

    // Dump mask to file.
    return dumpToTmpfile(mask, sizeof(MaskPixel), ubbWidth * ubbHeight);
}

