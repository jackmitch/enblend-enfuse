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
extern int MaximumLevels;
extern bool Wraparound;
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

/** Caluclate the union bounding box of these two images.
 *  Files are OutputWidth * OutputHeight * uint32.
 */
void ubbBounds(FILE *uint32File1, FILE *uint32File2) {
    UBBFirstX = OutputWidth - 1;
    UBBLastX = 0;
    UBBFirstY = OutputHeight - 1;
    UBBLastY = 0;

    uint32 *line1 = (uint32*)malloc(OutputWidth * sizeof(uint32));
    if (line1 == NULL) {
        cerr << endl
             << "enblend: out of memory (in ubbBounds for line1)" << endl;
        exit(1);
    }

    uint32 *line2 = (uint32*)malloc(OutputWidth * sizeof(uint32));
    if (line2 == NULL) {
        cerr << endl
             << "enblend: out of memory (in ubbBounds for line2)" << endl;
        exit(1);
    }

    // Make sure we're at the beginning of the files.
    rewind(uint32File1);
    rewind(uint32File2);

    for (uint32 y = 0; y < OutputHeight; y++) {
        readFromTmpfile(line1, sizeof(uint32), OutputWidth, uint32File1);
        readFromTmpfile(line2, sizeof(uint32), OutputWidth, uint32File2);

        uint32 *pixel1 = line1;
        uint32 *pixel2 = line2;
        for (uint32 x = 0; x < OutputWidth; x++) {
            if (*pixel1 != 0 || *pixel2 != 0) {
                // Pixel is in the union bounding box.
                UBBFirstX = min(x, UBBFirstX);
                UBBLastX = max(x, UBBLastX);
                UBBFirstY = min(y, UBBFirstY);
                UBBLastY = max(y, UBBLastY);
            }
            pixel1++;
            pixel2++;
        }
    }

    free(line1);
    free(line2);

    if (Verbose > 0) {
        cout << "Union bounding box = ("
             << UBBFirstX << ", " << UBBFirstY
             << ") -> ("
             << UBBLastX << ", " << UBBLastY
             << ")" << endl;
    }

    return;
}

/** Caluclate the region of interest for pyramids and the number of
 *  levels we're going to make based on the blend mask.
 */
uint32 roiBounds(FILE *maskFile) {

    uint32 ubbWidth = UBBLastX - UBBFirstX + 1;

    // First calculate the bounding box of the transition line.
    int32 firstMulticolorRow = UBBLastY;
    int32 lastMulticolorRow = UBBFirstY;
    int32 firstMulticolorColumn = UBBLastX;
    int32 lastMulticolorColumn = UBBFirstX;

    // Allocate memory for two rows.
    MaskPixel *firstRow = (MaskPixel*)malloc(ubbWidth * sizeof(MaskPixel));
    if (firstRow == NULL) {
        cerr << endl
             << "enblend: out of memory (in roiBounds for firstRow)" << endl;
        exit(1);
    }

    MaskPixel *row = (MaskPixel*)malloc(ubbWidth * sizeof(MaskPixel));
    if (row == NULL) {
        cerr << endl
             << "enblend: out of memory (in roiBounds for row)" << endl;
        exit(1);
    }

    // Make sure we're at the beginning of the mask file.
    rewind(maskFile);

    bool foundMulticolorRow = false;
    bool foundMulticolorColumn = false;
    for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
        readFromTmpfile(row, sizeof(MaskPixel), ubbWidth, maskFile);
        if (y == UBBFirstY) {
            memcpy(firstRow, row, ubbWidth * sizeof(MaskPixel));
        }

        MaskPixel *maskPFirst = firstRow;
        MaskPixel *maskP = row;
        uint8 leftmostColor = maskP->r;
        for (uint32 x = UBBFirstX; x <= UBBLastX; x++) {
            if (maskP->r != leftmostColor) {
                // Row y is a multicolor row.
                firstMulticolorRow = min(firstMulticolorRow, (int32)y);
                lastMulticolorRow = max(lastMulticolorRow, (int32)y);
                foundMulticolorRow = true;
            }
            if (maskP->r != maskPFirst->r) {
                // Column x is a multicolor column.
                firstMulticolorColumn = min(firstMulticolorColumn, (int32)x);
                lastMulticolorColumn = max(lastMulticolorColumn, (int32)x);
                foundMulticolorColumn = true;
            }
            maskP++;
            maskPFirst++;
        }
    }

    free(row);
    free(firstRow);

    if (!foundMulticolorRow && !foundMulticolorColumn) {
        // Sanity check for situation where there is no transition line
        // this means one image has no contribution.
        return 0;
    }
    else if (foundMulticolorRow && !foundMulticolorColumn) {
        // The transition line is vertical.
        firstMulticolorColumn = UBBFirstX;
        lastMulticolorColumn = UBBLastX;
    }
    else if (foundMulticolorColumn && !foundMulticolorRow) {
        // The transition line is horizontal.
        firstMulticolorRow = UBBFirstY;
        lastMulticolorRow = UBBLastY;
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
    uint32 levels = 1;
    while (true) {

        // filterHalfWidth = how many pixels the transition line spreads
        // out into layer levels-1, given the precision of the mask.
        // filterHalfWidth = how many pixels an image feature spreads
        // out into layer levels-1, given the precision of the image.
        // Add these up to find the spread of the bounding box that is
        // necessary to include all of the level 0 pixels that will influence
        // blending on level levels-1.
        // Since image precision = mask precision = 8 bpp, I'm just doubling
        // the value.
        int32 extent = 2 * filterHalfWidth(levels - 1, 255);

        if (Wraparound && ubbWidth == OutputWidth
                && ((firstMulticolorColumn - extent < (int32)UBBFirstX)
                        || (lastMulticolorColumn + extent > (int32)UBBLastX))) {
            // If the image is a 360 pano,
            // and the images in this blend step wrap around,
            // and the transition line is within extent pixels of the edge,
            // Then make the roi the full output width.
            ROIFirstX = UBBFirstX;
            ROILastX = UBBLastX;
        } else {
            ROIFirstX = max(firstMulticolorColumn - extent, (int32)UBBFirstX);
            ROILastX = min((int32)UBBLastX, lastMulticolorColumn + extent);
        }
        ROIFirstY = max(firstMulticolorRow - extent, (int32)UBBFirstY);
        ROILastY = min((int32)UBBLastY, lastMulticolorRow + extent);

        //cout << "Level " << levels << " tl bb = ("
        //        << firstExtentColumn << ", "
        //        << firstExtentRow << ") -> ("
        //        << lastExtentColumn << ", "
        //        << lastExtentRow << ")"
        //        << endl;

        // caluclate short dimension.
        uint32 shortDimension = min(ROILastY - ROIFirstY + 1,
                ROILastX - ROIFirstX + 1);

        if (levels == (uint32)MaximumLevels) {
            // Hit the user-specified level limit.
            break;
        }
        else if ((shortDimension >> (levels - 1)) > 8) {
            // Try another level.
            levels++;
        }
        else {
            // quit - this is a good number of levels.
            break;
        }
    }

    // The outside border of the ROI must be solid white or black in the mask,
    // or else there will be a boundary condition problem in expand.
    // Give ourselves a one pixel pad to ensure against this.
    ROIFirstX = (ROIFirstX == 0) ? 0 : max(ROIFirstX - 1, UBBFirstX);
    ROILastX = min(ROILastX + 1, UBBLastX);
    ROIFirstY = (ROIFirstY == 0) ? 0 : max(ROIFirstY - 1, UBBFirstY);
    ROILastY = min(ROILastY + 1, UBBLastY);

    if (Verbose > 0) {
        cout << "Using " << levels << " blending levels" << endl;
        cout << "Region of Interest bounding box = ("
             << ROIFirstX << ", "
             << ROIFirstY << ") -> ("
             << ROILastX << ", "
             << ROILastY << ")"
             << endl;
    }

    if (levels == 1 && levels != (uint32)MaximumLevels) {
        cerr << "enblend: intersection of images is too small to make "
             << "more than one pyramid level."
             << endl;
    }

    return levels;
}

/** Copy pixels from inside the UBB but outside the ROI from srcImage
 *  into dstImage. Ignore transparent pixels. Ignore pixels that aren't black
 *  in the mask.
 *  src and dst files are OutputWidth * OutputHeight * uint32.
 *  mask file is UBBWidth * UBBHeight * MaskPixel
 */
void copyExcludedPixels(FILE *dst, FILE *src, FILE *mask) {

    uint32 ubbWidth = UBBLastX - UBBFirstX + 1;

    uint32 *srcLine = (uint32*)malloc(ubbWidth * sizeof(uint32));
    if (srcLine == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for srcLine)" << endl;
        exit(1);
    }

    uint32 *dstLine = (uint32*)malloc(ubbWidth * sizeof(uint32));
    if (dstLine == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for dstLine)" << endl;
        exit(1);
    }

    MaskPixel *maskLine = (MaskPixel*)malloc(ubbWidth * sizeof(uint32));
    if (maskLine == NULL) {
        cerr << endl
             << "enblend: out of memory (in createMask for maskLine)" << endl;
        exit(1);
    }

    for (uint32 y = UBBFirstY; y <= UBBLastY; y++) {
        fseek(src, (y * OutputWidth + UBBFirstX) * sizeof(uint32), SEEK_SET);
        fseek(dst, (y * OutputWidth + UBBFirstX) * sizeof(uint32), SEEK_SET);
        fseek(mask, (y - UBBFirstY) * ubbWidth * sizeof(MaskPixel), SEEK_SET);
        readFromTmpfile(srcLine, sizeof(uint32), ubbWidth, src);
        readFromTmpfile(dstLine, sizeof(uint32), ubbWidth, dst);
        readFromTmpfile(maskLine, sizeof(MaskPixel), ubbWidth, mask);
        
        uint32 *srcPixel = srcLine;
        uint32 *dstPixel = dstLine;
        MaskPixel *maskPixel = maskLine;
        for (uint32 x = UBBFirstX; x <= UBBLastX; x++) {
            if (y < ROIFirstY
                    || y > ROILastY
                    || x < ROIFirstX
                    || x > ROILastX) {
                if (maskPixel->a == 255 && maskPixel->r == 0) {
                    *dstPixel = *srcPixel;
                }
            }
            srcPixel++;
            dstPixel++;
            maskPixel++;
        }

        // Write the modified line back to dst file.
        fseek(dst, -(ubbWidth * sizeof(uint32)), SEEK_CUR);
        writeToTmpfile(dstLine, sizeof(uint32), ubbWidth, dst);
    }

    free(srcLine);
    free(dstLine);
    free(maskLine);

    return;
}

/** Copy ROI-sized area of LPPixels into output-sized tmpfile.
 *  Ignore pixels that are transparent in UBB-sized mask tmpfile.
 */
void copyROIToOutputWithMask(LPPixel *roi, FILE *uint32File, FILE *maskFile) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 ubbWidth = UBBLastX - UBBFirstX + 1;

    uint32 *outputLine = (uint32*)malloc(roiWidth * sizeof(uint32));
    if (outputLine == NULL) {
        cerr << endl
             << "enblend: out of memory in copy for outputLine" << endl;
        exit(1);
    }

    MaskPixel *maskLine = (MaskPixel*)malloc(roiWidth * sizeof(MaskPixel));
    if (maskLine == NULL) {
        cerr << endl
             << "enblend: out of memory in copy for maskLine" << endl;
        exit(1);
    }

    for (uint32 y = ROIFirstY; y <= ROILastY; y++) {
        uint32 ubbY = y - UBBFirstY;
        uint32 ubbX = ROIFirstX - UBBFirstX;

        // Read roiWidth-sized line from uint32File.
        fseek(uint32File, (y * OutputWidth + ROIFirstX) * sizeof(uint32), SEEK_SET);
        readFromTmpfile(outputLine, sizeof(uint32), roiWidth, uint32File);

        // Read roiWidth-sized line from maskFile.
        fseek(maskFile, (ubbY * ubbWidth + ubbX) * sizeof(MaskPixel), SEEK_SET);
        readFromTmpfile(maskLine, sizeof(MaskPixel), roiWidth, maskFile);

        MaskPixel *maskPixel = maskLine;
        uint32 *outputPixel = outputLine;
        for (uint32 x = 0; x < roiWidth; x++) {
            if (maskPixel->a == 255) {
                // Convert back to uint8 data from int16 data.
                roi->r = min(255, max(0, (int)roi->r));
                roi->g = min(255, max(0, (int)roi->g));
                roi->b = min(255, max(0, (int)roi->b));

                *outputPixel =
                        (roi->r & 0xFF)
                        | ((roi->g & 0xFF) << 8)
                        | ((roi->b & 0xFF) << 16)
                        | (0xFF << 24);

            }
            outputPixel++;
            maskPixel++;
            roi++;
        }

        // Write back results.
        fseek(uint32File, -(roiWidth * sizeof(uint32)), SEEK_CUR);
        writeToTmpfile(outputLine, sizeof(uint32), roiWidth, uint32File);

    }

    free(outputLine);
    free(maskLine);

    return;
}
