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
#include <vector>

#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>

#include "enblend.h"

using namespace std;

// Global values from command line parameters.
int Verbose = 0;
int MaximumLevels = 0;
bool OneAtATime = false;
uint32 OutputWidth = 0;
uint32 OutputHeight = 0;
double StitchMismatchThreshold = 0.4;
uint16 PlanarConfig;
uint16 Photometric;

// Union bounding box.
uint32 UBBFirstX;
uint32 UBBLastX;
uint32 UBBFirstY;
uint32 UBBLastY;

// The region of interest for pyramids.
uint32 ROIFirstX;
uint32 ROILastX;
uint32 ROIFirstY;
uint32 ROILastY;

/** Print the usage information and quit. */
void printUsageAndExit() {
    cout << "==== enblend, version " << VERSION << " ====" << endl;
    cout << "Usage: enblend [options] -o OUTPUT INPUTS" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << " -o filename       Write output to file" << endl;
    //TODO stitch mismatch avoidance is work-in-progress.
    //cout << " -t float          Stitch mismatch threshold, [0.0, 1.0]" << endl;
    cout << " -l number         Maximum number of levels to use" << endl;
    cout << " -s                Blend images one at a time" << endl;
    cout << " -v                Verbose" << endl;
    cout << " -h                Print this help message" << endl;
    exit(1);
}

int main(int argc, char** argv) {

    // The name of the output file.
    char *outputFileName = NULL;

    // List of input files.
    list<char*> inputFileNameList;
    list<char*>::iterator listIterator;

    // Parse command line.
    int c;
    extern char *optarg;
    extern int optind;
    while ((c = getopt(argc, argv, "sl:o:t:vh")) != -1) {
        switch (c) {
            case 'h': {
                printUsageAndExit();
                break;
            }
            case 'v': {
                Verbose++;
                break;
            }
            case 't': {
                StitchMismatchThreshold = strtod(optarg, NULL);
                if (StitchMismatchThreshold < 0.0
                        || StitchMismatchThreshold > 1.0) {
                    cerr << "enblend: threshold must be between "
                         << "0.0 and 1.0 inclusive."
                         << endl;
                    printUsageAndExit();
                }
                break;
            }
            case 'o': {
                if (outputFileName != NULL) {
                    cerr << "enblend: more than one output file specified."
                         << endl;
                    printUsageAndExit();
                    break;
                }
                int len = strlen(optarg) + 1;
                outputFileName = (char*)malloc(len * sizeof(char));
                if (outputFileName == NULL) {
                    cerr << "enblend: out of memory for outputFileName" << endl;
                    exit(1);
                }
                strncpy(outputFileName, optarg, len);
                break;
            }
            case 'l': {
                MaximumLevels = atoi(optarg);
                if (MaximumLevels < 1) {
                    cerr << "enblend: maximum levels must be 1 or more."
                         << endl;
                    printUsageAndExit();
                }
                break;
            }
            case 's': {
                OneAtATime = true;
                break;
            }
            default: {
                printUsageAndExit();
                break;
            }
        }
    }

    // Make sure mandatory output file name parameter given.
    if (outputFileName == NULL) {
        cerr << "enblend: no output file specified." << endl;
        printUsageAndExit();
    }

    // Remaining parameters are input files.
    if (optind < argc) {
        while (optind < argc) {
            inputFileNameList.push_back(argv[optind++]);
        }
    } else {
        cerr << "enblend: no input files specified." << endl;
        printUsageAndExit();
    }

    // Create output TIFF object.
    TIFF *outputTIFF = TIFFOpen(outputFileName, "w");
    if (outputTIFF == NULL) {
        cerr << "enblend: error opening output TIFF file \""
             << outputFileName
             << "\"" << endl;
        exit(1);
    }

    // Set basic parameters for output TIFF.
    TIFFSetField(outputTIFF, TIFFTAG_ORIENTATION, 1);
    TIFFSetField(outputTIFF, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(outputTIFF, TIFFTAG_BITSPERSAMPLE, 8);

    // Give output tiff the same parameters as first input tiff.
    // Check that all input tiffs are the same size.
    listIterator = inputFileNameList.begin();
    while (listIterator != inputFileNameList.end()) {
        TIFF *inputTIFF = TIFFOpen(*listIterator, "r");

        if (inputTIFF == NULL) {
            // Error opening tiff.
            cerr << "enblend: error opening input TIFF file \""
                 << *listIterator
                 << "\"" << endl;
            exit(1);
        }

        if (listIterator == inputFileNameList.begin()) {
            // The first input tiff
            TIFFGetField(inputTIFF, TIFFTAG_PLANARCONFIG, &PlanarConfig);
            TIFFGetField(inputTIFF, TIFFTAG_PHOTOMETRIC, &Photometric);
            TIFFGetField(inputTIFF, TIFFTAG_IMAGEWIDTH, &OutputWidth);
            TIFFGetField(inputTIFF, TIFFTAG_IMAGELENGTH, &OutputHeight);

            TIFFSetField(outputTIFF, TIFFTAG_IMAGEWIDTH, OutputWidth);
            TIFFSetField(outputTIFF, TIFFTAG_IMAGELENGTH, OutputHeight);
            TIFFSetField(outputTIFF, TIFFTAG_PHOTOMETRIC, Photometric);
            TIFFSetField(outputTIFF, TIFFTAG_PLANARCONFIG, PlanarConfig);

            // We already forced this to be a 4-channel-per pixel (above), so 
            // make the fourth channel be alpha if RGB image (see tiff spec)
            // If not specified, Photoshop doesn't open up with transparent area
            if (Photometric == PHOTOMETRIC_RGB) {
                uint16 v[1];				
                v[0] = EXTRASAMPLE_ASSOCALPHA;
                TIFFSetField(outputTIFF, TIFFTAG_EXTRASAMPLES, 1, v);
            }			

            if (Verbose > 0) {
                cout << "output size = "
                     << OutputWidth
                     << "x"
                     << OutputHeight << endl;
            }
        }
        else {
            uint32 otherWidth;
            uint32 otherHeight;
            TIFFGetField(inputTIFF, TIFFTAG_IMAGEWIDTH, &otherWidth);
            TIFFGetField(inputTIFF, TIFFTAG_IMAGELENGTH, &otherHeight);

            if (otherWidth != OutputWidth || otherHeight != OutputHeight) {
                cerr << "enblend: input TIFF \""
                     << *listIterator
                     << "\" size "
                     << otherWidth << "x" << otherHeight
                     << " is not the same size as the first TIFF \""
                     << inputFileNameList.front()
                     << "\" size "
                     << OutputWidth << "x" << OutputHeight
                     << endl;
                exit(1);
            }
        }
            
        TIFFClose(inputTIFF);
        listIterator++;
    }

    // Create the initial white image.
    uint32 *whiteImage = assemble(inputFileNameList, false);
    // mem = 1 fullsize uint32

    // Main blending loop
    while (!inputFileNameList.empty()) {

        // Write whiteImage to disk as output file.
        // free whiteImage.
        // mem = 0

        // Create the black image.
        // write to disk as temp file.
        // free blackImage.
        uint32 *blackImage = assemble(inputFileNameList, OneAtATime);
        // mem = 0

        // Create the blend mask.
        // load blackImage and whiteImage from disk one row at a time.
        MaskPixel *mask = createMask(whiteImage, blackImage);
        // maybe an extra-conservative nnt algo?
        // mem = 1 ubb MaskPixel + 2 ubb uint32 distance matrices
        // mem = 1 ubb MaskPixel

        // Calculate the ROI bounds and number of levels.
        uint32 levels = bounds(mask);
        if (MaximumLevels > 0) {
            levels = min(levels, (uint32)MaximumLevels);
        }
        // write mask to disk as temp file.
        // free mask.
        // mem = 0

        // Copy parts of blackImage outside of ROI into whiteImage.
        // work on blackImage and whiteImage from files one row at a time.
        copyExcludedPixels(whiteImage, blackImage);
        // mem = 0.

        // Build Laplacian pyramid from blackImage
        // load blackImage from disk one row at a time.
        // dump result to disk.
        vector<LPPixel*> *blackLP = laplacianPyramid(blackImage, levels);
        // mem = 4/3 roi LPPixel
        // mem = 0
        //savePyramid(*blackLP, "black");

        // Free allocated memory.
        // done with blackImage temp file.
        _TIFFfree(blackImage);
        // mem = 0

        // Build Gaussian pyramid from mask.
        // load mask from disk one row at a time.
        // dump result to disk.
        vector<LPPixel*> *maskGP = gaussianPyramid(mask, levels);
        // mem = 4/3 roi LPPixel
        // mem = 0
        savePyramid(*maskGP, "mask");

        // Build Laplacian pyramid from whiteImage
        // load whiteImage from disk one row at a time.
        vector<LPPixel*> *whiteLP = laplacianPyramid(whiteImage, levels);
        // mem = 4/3 roi LPPixel
        //savePyramid(*whiteLP, "white");

        // Blend pyramids
        // store only whiteLP in memory.
        // load blackLP and maskGP from disk one row at a time.
        blend(*whiteLP, *blackLP, *maskGP);

        // Free allocated memory.
        // done with maskGP and blackLP temp files.
        for (uint32 i = 0; i < levels; i++) {
            free((*maskGP)[i]);
            free((*blackLP)[i]);
        }
        delete maskGP;
        delete blackLP;
        // mem = 1 fullsize uint32 + fullsize MaskPixel + 4/3 roi LPPixel
        // should mem = 4/3 roi LPPixel

        // Collapse result back into whiteImage.
        collapsePyramid(*whiteLP, whiteImage, mask);
        // copy whiteLP[0] into whiteImage on disk.
        // can be done one row at a time, don't neet to load whiteImage or maskImage entirely.

        free(mask);
        // mem = 1 fullsize uint32 + 4/3 roi LPPixel
        for (uint32 i = 0; i < levels; i++) {
            free((*whiteLP)[i]);
        }
        delete whiteLP;
        // mem = 1 fullsize uint32

    }

    // don't need to do any of this - whiteImage is already complete.

    if (Verbose > 0) {
        cout << "Writing output file." << endl;
    }

    // dump output scanlines.
    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFWriteScanline(outputTIFF,
                &(whiteImage[i * OutputWidth]),
                i,
                8);
    }

    // close outputTIFF.
    TIFFClose(outputTIFF);

    _TIFFfree(whiteImage);

    return 0;
}
