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

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

// Region of interest for this operation.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

			        	/* Direction masks:		*/
			        	/*   N	   S	 W     E	*/
static	uint32	directionMasks[]	= { 0200, 0002, 0040, 0010 };

/*	True if pixel neighbor map indicates the pixel is 8-simple and	*/
/*	not an end point and thus can be remarked.  The neighborhood	*/
/*	map is defined as an integer of bits abcdefghi with a non-zero	*/
/*	bit representing a non-zero pixel.  The bit assignment for the	*/
/*	neighborhood is:						*/
/*									*/
/*				a b c					*/
/*				d e f					*/
/*				g h i					*/

static	unsigned char	remarkTable[512] = {
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

/** Convert x,y coordinates to an index in the 1-d data array.
 */
inline int coord(int x, int y) {
    return ((y * OutputWidth) + x);
}

/** This function analyzes all of the black pixels in the input mask
 *  and determines how far away each one is from a red, green, or blue
 *  region. If the pixel is closer to a red region than any other region,
 *  the black pixel is changed to red. Likewise for green and blue.
 *
 *  Based on an algorithm in "Efficient Binary Image Thinning Using
 *  Neighborhood Maps" by Joseph M. Cychosz, from "Graphics Gems IV",
 *  ed. Paul Heckbert. Academic Press, Inc. 1994.
 */
void thinMask(uint32 *mask) {

    int passCount = 0;
    int passIndex = 0;
    int modifiedPixelCount = 1;
    int blackPixelsRemaining = 0;

    uint32 pMap;
    uint32 qMap;
    uint32 *qbMapArray;
    uint32 directionMask;

    int32 directionIndexOffsets[] = {
            -OutputWidth, // NORTH
            OutputWidth,  // SOUTH
            -1,           // WEST
            1};           // EAST
    int32 directionIndexOffset;

    uint32 x;
    uint32 y;

    uint32 roiWidth = ROILastX - ROIFirstX + 1;

    qbMapArray = (uint32*)malloc(roiWidth * sizeof(uint32));
    if (qbMapArray == NULL) {
        cerr << "enblend: malloc failed for qbMapArray" << endl;
        exit(1);
    }

    while (modifiedPixelCount != 0) {
        passCount++;
        modifiedPixelCount = 0;

        for (passIndex = 0; passIndex < 4; passIndex++) {
            directionMask = directionMasks[passIndex];
            directionIndexOffset = directionIndexOffsets[passIndex];

            blackPixelsRemaining = 0;

            // Build initial previous scan buffer.
            pMap = (mask[coord(ROIFirstX,ROIFirstY)] == BLACK) ? 0777 : 0776;
            for (x = ROIFirstX; x < ROILastX; x++) {
                pMap = ((pMap << 1) & 0666)
                        | ((mask[coord(x+1,ROIFirstY)] == BLACK) ? 0111 : 0110);
                qbMapArray[x-ROIFirstX] = pMap;
            }
            // Right edge pixel
            qbMapArray[x-ROIFirstX] = ((pMap << 1) & 0666) | 0111;

            // Scan image for pixel remarking candidates.
            for (y = ROIFirstY; y < ROILastY; y++) {

                // Calculate first pMap.
                qMap = qbMapArray[0];
                pMap = ((qMap << 2) & 0110) | 0666
                        | ((mask[coord(ROIFirstX,y+1)] == BLACK) ? 0001 : 0000);

                // Process pixels across row.
                for (x = ROIFirstX; x <= ROILastX ; x++) {
                    qMap = qbMapArray[x-ROIFirstX];
                    pMap = ((pMap << 1) & 0666) | ((qMap << 3) & 0110);

                    if (x != OutputWidth - 1) {
                        pMap |= ((mask[coord(x+1,y+1)] == BLACK) ? 0001 : 0000);
                    } else {
                        // Right edge.
                        pMap |= 0001;
                    }

                    qbMapArray[x-ROIFirstX] = pMap;

                    if (((pMap & directionMask) == 0)
                            && remarkTable[pMap]) {
                        // Modify pixel.
                        modifiedPixelCount++;
                        mask[coord(x,y)] =
                                mask[coord(x,y) + directionIndexOffset];
                    }
                    else if ((pMap & 0020) != 0) {
                        blackPixelsRemaining++;
                    }
                }
            }

            // Process bottom scanline.
            // Initial pMap.
            qMap = qbMapArray[0];
            pMap = ((qMap << 2) & 0110) | 0667;

            for (x = ROIFirstX; x <= ROILastX; x++) {
                qMap = qbMapArray[x-ROIFirstX];
                pMap = ((pMap << 1) & 0666) | ((qMap << 3) & 0110) | 0001;

                if (((pMap & directionMask) == 0)
                        && remarkTable[pMap]) {
                    // Modify pixel.
                    modifiedPixelCount++;
                    mask[coord(x,y)] =
                            mask[coord(x,y) + directionIndexOffset];
                }
                else if ((pMap & 0020) != 0) {
                    blackPixelsRemaining++;
                }
            }
        }

        if (Verbose > 0) {
            cout << "thin: pass " << passCount << ", "
                 << modifiedPixelCount
                 << " pixels remarked, "
                 << blackPixelsRemaining
                 << " pixels to go."
                 << endl;

        }

    }

    free(qbMapArray);

    return;

}
