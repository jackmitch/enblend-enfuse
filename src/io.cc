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
#include <list>
#ifdef __GW32C__
#include <map>
#endif
#include <vector>

#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#include <unistd.h>

#include "enblend.h"

using namespace std;

// Global values from command line parameters.
extern int Verbose;
extern int MaximumLevels;
extern bool OneAtATime;
extern uint32 OutputWidth;
extern uint32 OutputHeight;
extern double StitchMismatchThreshold;
extern uint16 PlanarConfig;
extern uint16 Photometric;

// Union bounding box.
extern uint32 UBBFirstX;
extern uint32 UBBLastX;
extern uint32 UBBFirstY;
extern uint32 UBBLastY;

// The region of interest for pyramids.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

#ifdef __GW32C__
map<FILE*, char*> FileToFilenameMap;
#endif

/** Close a temporary file, delete if necessary on GnuWin32 systems. */
void closeTmpfile(FILE *f) {

    fclose(f);

#ifdef __GW32C__
    // Delete the file on GnuWin32.
    if (unlink(FileToFilenameMap[f]) < 0) {
        cerr << "enblend: error deleting temporary file." << endl;
        perror("reason");
        exit(1);
    }
    free(FileToFilenameMap[f]);
    FileToFilenameMap.erase(f);
#endif

    return;
}

/** Create a temporary file. */
FILE *createTmpfile() {

    char tmpFilename[] = ".enblend_tmpXXXXXX";

#ifdef HAVE_MKSTEMP
    int tmpFD = mkstemp(tmpFilename);
    if (tmpFD < 0) {
        cerr << "enblend: error creating temporary file." << endl;
        exit(1);
    }

    FILE *f = fdopen(tmpFD, "wb+");
#else
    char *tmpReturn = mktemp(tmpFilename);
    if (tmpReturn == NULL) {
        cerr << "enblend: error creating temporary file." << endl;
        exit(1);
    }

    FILE *f = fopen(tmpFilename, "wb+");
#endif

    if (f == NULL) {
        cerr << "enblend: error opening temporary file." << endl;
        perror("reason");
        exit(1);
    }

#ifdef __GW32C__
    int tmpFilenameLen = strlen(tmpFilename) + 1;
    char *tmpFilenameCopy = (char*)malloc(tmpFilenameLen * sizeof(char));
    if (tmpFilenameCopy == NULL) {
        cerr << endl
             << "enblend: out of memory (in createTmpfile for tmpFilenameCopy)"
             << endl;
        exit(1);
    }
    strncpy(tmpFilenameCopy, tmpFilename, tmpFilenameLen);
    FileToFilenameMap[f] = tmpFilenameCopy;
#else
    // Arrange to have the file deleted on close.
    if (unlink(tmpFilename) < 0) {
        cerr << "enblend: error deleting temporary file." << endl;
        perror("reason");
        exit(1);
    }
#endif

    return f;
}

/** Read data from a temporary file.
 *  Puts size * nmemb bytes into ptr from the current file position.
 *  File position is moved forward the same number of bytes.
 */
void readFromTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t itemsRead = fread(ptr, size, nmemb, stream);
    if (itemsRead < nmemb) {
        cerr << "enblend: error reading from temporary file." << endl;
        perror("reason");
        exit(1);
    }
}

/** Write data to a temporary file.
 *  Writes size*nmemb bytes from ptr.
 *  File position is moved forward the same number of bytes.
 */
void writeToTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t itemsWritten = fwrite(ptr, size, nmemb, stream);
    if (itemsWritten < nmemb) {
        cerr << "enblend: error writing to temporary file." << endl;
        perror("reason");
        exit(1);
    }
}

/** Write data to a temporary file.
 *  Writes size*nmemb bytes of data from ptr.
 *  Keeps the file open and rewinds it to the beginning.
 *  Frees ptr.
 */
FILE *dumpToTmpfile(void *ptr, size_t size, size_t nmemb) {
    FILE *f = createTmpfile();

    writeToTmpfile(ptr, size, nmemb, f);
    rewind(f);

    free(ptr);

    return f;
}

/** Dump a pyramid to a temporary file.
 *  Keeps the file open and rewinds it to the beginning.
 *  Frees the space used by the pyramid.
 */
FILE *dumpPyramidToTmpfile(vector<LPPixel*> &v) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Open output temporary file.
    FILE *pyramidFile = createTmpfile();
    if (pyramidFile == NULL) {
        cerr << "enblend: error opening temporary file." << endl;
        perror("reason");
        exit(1);
    }

    for (uint32 l = 0; l < v.size(); l++) {
        uint32 lPixels = (roiWidth >> l) * (roiHeight >> l);
        writeToTmpfile(v[l], sizeof(LPPixel), lPixels, pyramidFile);
        free(v[l]);
    }

    rewind(pyramidFile);

    return pyramidFile;
}

/** Save a mask as a TIFF called mask.tif.
 */
void saveMaskAsTIFF(MaskPixel *mask, char *filename) {

    uint32 *image = (uint32*)calloc(OutputWidth * OutputHeight, sizeof(uint32));
    if (image == NULL) {
        cerr << endl
             << "enblend: out of memory (in saveMaskAsTIFF for image)" << endl;
        exit(1);
    }

    TIFF *outputTIFF = TIFFOpen(filename, "w");
    if (outputTIFF == NULL) {
        cerr << "enblend: error opening TIFF file \""
             << filename << "\"" << endl;
        exit(1);
    }

    TIFFSetField(outputTIFF, TIFFTAG_ORIENTATION, 1);
    TIFFSetField(outputTIFF, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(outputTIFF, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(outputTIFF, TIFFTAG_IMAGEWIDTH, OutputWidth);
    TIFFSetField(outputTIFF, TIFFTAG_IMAGELENGTH, OutputHeight);
    TIFFSetField(outputTIFF, TIFFTAG_PHOTOMETRIC, Photometric);
    TIFFSetField(outputTIFF, TIFFTAG_PLANARCONFIG, PlanarConfig);

    MaskPixel *maskP = mask;
    uint32 *imageP = image;
    for (uint32 y = 0; y < OutputHeight; y++) {
        for (uint32 x = 0; x < OutputWidth; x++) {
            if (y >= UBBFirstY && y <= UBBLastY
                     && x >= UBBFirstX && x <= UBBLastX) {
                *imageP = 
                        (maskP->r & 0xFF)
                        | ((maskP->g & 0xFF) << 8)
                        | ((maskP->b & 0xFF) << 16)
                        | ((maskP->a & 0xFF) << 24);
                maskP++;
            }
            else *imageP = 0;
            imageP++;
        }
    }
    
    for (uint32 scan = 0; scan < OutputHeight; scan++) {
        TIFFWriteScanline(outputTIFF,
                &(image[scan * OutputWidth]),
                scan,
                8);
    }

    free(image);

    TIFFClose(outputTIFF);

}

/** Save a pyramid as a multilayered tiff.
 */
void savePyramidAsTIFF(vector<LPPixel*> &p, char *filename) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Make a copy of p.
    vector<LPPixel*> pCopy;
    for (uint32 i = 0; i < p.size(); i++) {
        uint32 levelPixels = (roiWidth >> i) * (roiHeight >> i);
        LPPixel *level = (LPPixel*)malloc(levelPixels * sizeof(LPPixel));
        if (level == NULL) {
            cerr << endl
                 << "enblend: out of memory (in savePyramidAsTIFF for level)"
                 << endl;
            exit(1);
        }
        pCopy.push_back(level);
        memcpy(level, p[i], levelPixels * sizeof(LPPixel));
    }

    // Output image.
    uint32 *image = (uint32*)malloc(OutputWidth * OutputHeight * sizeof(uint32));
    if (image == NULL) {
        cerr << endl
             << "enblend: out of memory (in savePyramidAsTIFF for image)"
             << endl;
        exit(1);
    }

    for (uint32 i = 0; i < p.size(); i++) {

        char filenameBuf[512];
        snprintf(filenameBuf, 512, "%s%04u.tif", filename, (unsigned int)i);
        cout << filenameBuf << endl;
        TIFF *outputTIFF = TIFFOpen(filenameBuf, "w");
        if (outputTIFF == NULL) {
            cerr << "enblend: error opening TIFF file \""
                 << filenameBuf << "\"" << endl;
            exit(1);
        }


        vector<LPPixel*> pCopySubset;

        // clear levels up to i.
        for (uint32 j = 0; j < i; j++) {
            uint32 jPixels = (roiWidth >> j) * (roiHeight >> j);
            memset(pCopy[j], 0, jPixels * sizeof(LPPixel));
            pCopySubset.push_back(pCopy[j]);
        }
        pCopySubset.push_back(pCopy[i]);

        // Collapse from i up.
        collapsePyramid(pCopySubset);

        // Write pCopySubset into a layer of the tiff.
        //TIFFSetDirectory(outputTIFF, i);
        TIFFSetField(outputTIFF, TIFFTAG_ORIENTATION, 1);
        TIFFSetField(outputTIFF, TIFFTAG_SAMPLESPERPIXEL, 4);
        TIFFSetField(outputTIFF, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(outputTIFF, TIFFTAG_IMAGEWIDTH, OutputWidth);
        TIFFSetField(outputTIFF, TIFFTAG_IMAGELENGTH, OutputHeight);
        TIFFSetField(outputTIFF, TIFFTAG_PHOTOMETRIC, Photometric);
        TIFFSetField(outputTIFF, TIFFTAG_PLANARCONFIG, PlanarConfig);

        memset(image, 0, OutputWidth * OutputHeight * sizeof(uint32));

        LPPixel *pixel = pCopy[0];
        for (uint32 y = ROIFirstY; y <= ROILastY; y++) {
            for (uint32 x = ROIFirstX; x <= ROILastX; x++) {
                pixel->r = min(255, max(0, (int)abs(pixel->r)));
                pixel->g = min(255, max(0, (int)abs(pixel->g)));
                pixel->b = min(255, max(0, (int)abs(pixel->b)));
                image[y * OutputWidth + x] =
                        (pixel->r & 0xFF)
                        | ((pixel->g & 0xFF) << 8)
                        | ((pixel->b & 0xFF) << 16)
                        | (0xFF << 24);
                pixel++;
            }
        }

        for (uint32 line = 0; line < OutputHeight; line++) {
            TIFFWriteScanline(outputTIFF,
                    &(image[line * OutputWidth]),
                    line,
                    8);
        }

        //TIFFWriteDirectory(outputTIFF);
        TIFFClose(outputTIFF);

    }

    //TIFFClose(outputTIFF);

    for (uint32 i = 0; i < p.size(); i++) {
        free(pCopy[i]);
    }

    free(image);

    return;
}
