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

// Region of interest for this operation.
extern uint32 UBBFirstX;
extern uint32 UBBLastX;
extern uint32 UBBFirstY;
extern uint32 UBBLastY;

			        	/* Direction masks:		*/
			        	/*   N	   S	 W     E	*/
static	uint16	directionMasks[]	= { 0200, 0002, 0040, 0010 };

/*	True if pixel neighbor map indicates the pixel is 8-simple and	*/
/*	not an end point and thus can be remarked.  The neighborhood	*/
/*	map is defined as an integer of bits abcdefghi with a non-zero	*/
/*	bit representing a non-zero pixel.  The bit assignment for the	*/
/*	neighborhood is:						*/
/*									*/
/*				a b c					*/
/*				d e f					*/
/*				g h i					*/

static	uint8   remarkTable[512] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

/** This function analyzes all of the thinnable pixels in the input mask
 *  and determines how far away each one is from a green or blue
 *  region. If the pixel is closer to a green region than any other region,
 *  the thinnable pixel is changed to green. Likewise for blue.
 *  This function only operates within the region-of-interest.
 *
 *  Based on an algorithm in "Efficient Binary Image Thinning Using
 *  Neighborhood Maps" by Joseph M. Cychosz, from "Graphics Gems IV",
 *  ed. Paul Heckbert. Academic Press, Inc. 1994.
 */
void thinMask(MaskPixel *mask) {

    uint32 passCount = 0;
    uint32 passIndex = 0;
    uint32 modifiedPixelCount = 1;
    uint32 thinnablePixelsRemaining = 0;

    uint16 pMap;
    uint16 qMap;
    uint16 *qbMapArray;
    uint16 directionMask;

    int32 directionIndexOffsets[] = {
            -OutputWidth, // NORTH
            OutputWidth,  // SOUTH
            -1,           // WEST
            1};           // EAST
    int32 directionIndexOffset;

    uint32 roiWidth = UBBLastX - UBBFirstX + 1;
    uint32 roiHeight = UBBLastY - UBBFirstY + 1;

    MaskPixel *firstPixel = &mask[UBBFirstY * OutputWidth + UBBFirstX];

    qbMapArray = (uint16*)malloc(roiWidth * sizeof(uint16));
    if (qbMapArray == NULL) {
        cerr << "enblend: out of memory for qbMapArray" << endl;
        exit(1);
    }

    while (modifiedPixelCount != 0) {
        passCount++;
        modifiedPixelCount = 0;

        for (passIndex = 0; passIndex < 4; passIndex++) {
            thinnablePixelsRemaining = 0;

            directionMask = directionMasks[passIndex];
            directionIndexOffset = directionIndexOffsets[passIndex];

            // Build initial previous scan buffer.
            pMap = 0776 | firstPixel->r;

            uint32 x;
            for (x = 0; x < (roiWidth-1); x++) {
                pMap = ((pMap << 1) & 0666) | 0110 | firstPixel[x+1].r;
                qbMapArray[x] = pMap;
            }
            // Right edge pixel is implicitly thinnable.
            qbMapArray[x] = ((pMap << 1) & 0666) | 0111;

            MaskPixel *currentPixel;

            // Scan image for pixel remarking candidates.
            uint32 y;
            for (y = 0; y < (roiHeight-1); y++) {
                currentPixel = &firstPixel[y * OutputWidth];

                // Calculate first pMap.
                qMap = qbMapArray[0];
                pMap = ((qMap << 2) & 0110) | 0666
                        | currentPixel[OutputWidth].r;

                // Process pixels across row.
                for (x = 0; x < (roiWidth-1); x++) {
                    qMap = qbMapArray[x];
                    pMap = ((pMap << 1) & 0666)
                            | ((qMap << 3) & 0110)
                            | currentPixel[OutputWidth + 1].r;

                    qbMapArray[x] = pMap;

                    if (((pMap & directionMask) == 0)
                            && remarkTable[pMap]) {
                        // Modify pixel.
                        modifiedPixelCount++;
                        // Remove thinnable bit.
                        currentPixel->r = 0;
                        // Set color mark bit.
                        currentPixel->g = currentPixel[directionIndexOffset].g;
                        currentPixel->b = currentPixel[directionIndexOffset].b;
                    }
                    else if (pMap & 0020) {
                        thinnablePixelsRemaining++;
                    }

                    currentPixel++;
                }

                // Right edge of row.
                qMap = qbMapArray[x];
                pMap = ((pMap << 1) & 0666)
                        | ((qMap << 3) & 0110)
                        | 1;

                qbMapArray[x] = pMap;

                if (((pMap & directionMask) == 0)
                        && remarkTable[pMap]) {
                    // Modify pixel.
                    modifiedPixelCount++;
                    // Remove thinnable bit.
                    currentPixel->r = 0;
                    // Set color mark bit.
                    currentPixel->g = currentPixel[directionIndexOffset].g;
                    currentPixel->b = currentPixel[directionIndexOffset].b;
                }
                else if (pMap & 0020) {
                    thinnablePixelsRemaining++;
                }
            }

            // Bottom scanline
            currentPixel = &firstPixel[y * OutputWidth];

            // Initial pMap.
            qMap = qbMapArray[0];
            pMap = ((qMap << 2) & 0110) | 0667;

            for (x = 0; x < roiWidth; x++) {
                qMap = qbMapArray[x];
                pMap = ((pMap << 1) & 0666) | ((qMap << 3) & 0110) | 0001;

                if (((pMap & directionMask) == 0)
                        && remarkTable[pMap]) {
                    // Modify pixel.
                    modifiedPixelCount++;
                    // Remove thinnable bit.
                    currentPixel->r = 0;
                    // Set color mark bit.
                    currentPixel->g = currentPixel[directionIndexOffset].g;
                    currentPixel->b = currentPixel[directionIndexOffset].b;
                }
                else if (pMap & 0020) {
                    thinnablePixelsRemaining++;
                }

                currentPixel++;
            }
        }

        if (Verbose > 0) {
            cout << "thin: pass " << passCount << ", "
                 << modifiedPixelCount
                 << " pixels thinned, "
                 << thinnablePixelsRemaining
                 << " pixels to go."
                 << endl;

        }

    }

    free(qbMapArray);

    return;

}
