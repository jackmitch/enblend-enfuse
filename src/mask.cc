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
#include <tiffio.h>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

/** Calculate a blending mask between outputBuf and inputTIFF.
 *
 *  Based on an algorithm in "Efficient Binary Image Thinning Using
 *  Neighborhood Maps" by Joseph M. Cychosz, from "Graphics Gems IV",
 *  ed. Paul Heckbert. Academic Press, Inc. 1994.
 */
uint32 *createMask(uint32 *outputBuf, TIFF *inputTIFF) {

    // Allocate memory for the blending mask.
    uint32 *mask = (uint32*)_TIFFmalloc(
            OutputWidth * OutputHeight * sizeof(uint32));
    if (mask == NULL) {
        cerr << "enblend: malloc failed for mask" << endl;
        exit(1);
    }

    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFReadScanline(inputTIFF,
                &(mask[i * OutputWidth]),
                i,
                8);

        for (uint32 j = 0; j < OutputWidth; j++) {
            uint32 outPixel = outputBuf[(i * OutputWidth) + j];
            uint32 inPixel = mask[(i * OutputWidth) + j];

            if (outPixel == 0 && inPixel == 0) {
                // Pixel is not in the union of the two images.
                // Make the mask pixel black.
                mask[(i * OutputWidth) + j] = 0xFF000000;
            }
            else if (outPixel != 0 && inPixel == 0) {
                // Pixel is in out only.
                // Make the mask pixel blue.
                mask[(i * OutputWidth) + j] = 0xFFFF0000;
            }
            else if (outPixel == 0 && inPixel != 0) {
                // Pixel is in input only.
                // Make the mask pixel green.
                mask[(i * OutputWidth) + j] = 0xFF00FF00;
            }
            else {
                // Pixel is in the intersection of the two images.
                // Make the mask pixel black.
                //TODO: red if different above threshold.
                mask[(i * OutputWidth) + j] = 0xFF000000;
            }
        }
    }

    return mask;
}

