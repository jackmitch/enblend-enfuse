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
#include <list>
#include <vector>

#include <getopt.h>
#include <string.h>
#include <tiffio.h>

#include "enblend.h"

using namespace std;

// Global values from command line parameters.
int Verbose = 0;
uint32 OutputWidth = 0;
uint32 OutputHeight = 0;
double StitchMismatchThreshold = 0.4;

// The region of interest for the current operation.
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
    //FIXME stitch mismatch avoidance is broken.
    //cout << " -t float          Stitch mismatch threshold, [0.0, 1.0]" << endl;
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
    while ((c = getopt(argc, argv, "o:t:vh")) != -1) {
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
                if (StitchMismatchThreshold < 0.0 || StitchMismatchThreshold > 1.0) {
                    cerr << "enblend: threshold must be between 0.0 and 1.0 inclusive."
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
                    cerr << "enblend: malloc failed for outputFileName" << endl;
                    exit(1);
                }
                strncpy(outputFileName, optarg, len);
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
            uint16 planarConfig;
            uint16 photometric;

            TIFFGetField(inputTIFF, TIFFTAG_PLANARCONFIG, &planarConfig);
            TIFFGetField(inputTIFF, TIFFTAG_PHOTOMETRIC, &photometric);
            TIFFGetField(inputTIFF, TIFFTAG_IMAGEWIDTH, &OutputWidth);
            TIFFGetField(inputTIFF, TIFFTAG_IMAGELENGTH, &OutputHeight);

            TIFFSetField(outputTIFF, TIFFTAG_IMAGEWIDTH, OutputWidth);
            TIFFSetField(outputTIFF, TIFFTAG_IMAGELENGTH, OutputHeight);
            TIFFSetField(outputTIFF, TIFFTAG_PHOTOMETRIC, photometric);
            TIFFSetField(outputTIFF, TIFFTAG_PLANARCONFIG, planarConfig);

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

    // Allocate memory for the output image.
    uint32 *outputBuf = (uint32*)_TIFFmalloc(OutputWidth * OutputHeight * sizeof(uint32));
    if (outputBuf == NULL) {
        cerr << "enblend: malloc failed for outputBuf" << endl;
        exit(1);
    }

    // Remove the first TIFF file name.
    char *fileZeroName = inputFileNameList.front();
    inputFileNameList.pop_front();

    if (Verbose > 0) {
        cout << "Starting with first image \""
             << fileZeroName
             << "\"" << endl;
    }

    TIFF *inputTIFF = TIFFOpen(fileZeroName, "r");
    if (inputTIFF == NULL) {
        // Error opening tiff.
        cerr << "enblend: error opening input TIFF file \""
             << fileZeroName
             << "\"" << endl;
        exit(1);
    }

    // The first TIFF gets copied directly into the output.
    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFReadScanline(inputTIFF,
                &(outputBuf[i * OutputWidth]),
                i,
                8);
    }

    // Done with the first TIFF.
    TIFFClose(inputTIFF);

    //TODO: Sort the images in the list.

    // Main blending loop
    while (!inputFileNameList.empty()) {
        char *inputFileName = inputFileNameList.front();
        inputFileNameList.pop_front();

        if (Verbose > 0) {
            cout << "Next image: \"" << inputFileName << "\"" << endl;
        }

        inputTIFF = TIFFOpen(inputFileName, "r");
        if (inputTIFF == NULL) {
            // Error opening tiff.
            cerr << "enblend: error opening input TIFF file \""
                 << inputFileName
                 << "\"" << endl;
            exit(1);
        }

        // Create blend mask. Black pixels are the inputTIFF zone,
        // White pixels are the outputBuf zone
        uint32 *mask = createMask(outputBuf, inputTIFF);

        // Count max number of levels we can make from ROI size.
        int32 shortDimension = min(ROILastX - ROIFirstX + 1,
                ROILastY - ROIFirstY + 1);
        if (shortDimension < 4) {
            cerr << "enblend: union of images is too small to make even one pyramid."
                 << endl;
            exit(1);
        }
        int32 maximumLevels = 1;
        while (shortDimension > 8) {
            shortDimension = shortDimension >> 1;
            maximumLevels++;
        }

        // Build Gaussian pyramid from mask.
        vector<uint32*> *maskPyramid = gaussianPyramid(mask, maximumLevels);

        // TODO: find a good level to stop at where the blending zone does not
        // overrun the overlap zone too much.

        // Build Laplacian pyramid from outputBuf
        vector<uint32*> *outputLP = laplacianPyramid(outputBuf, maximumLevels);

        // Build Laplacian pyramid from inputTIFF
        vector<uint32*> *inputLP = laplacianPyramid(inputTIFF, maximumLevels);

        // blend.
        blend(maskPyramid, inputLP, outputLP);

        // collapse result into outputBuf
        collapsePyramid(outputLP, outputBuf);

        TIFFClose(inputTIFF);

        //_TIFFfree(mask);

        for (int32 i = 0; i < maximumLevels; i++) {
            uint32 *g = (*maskPyramid)[i];
            free(g);
            g = (*outputLP)[i];
            free(g);
            g = (*inputLP)[i];
            free(g);
        }
        delete maskPyramid;
        delete outputLP;
        delete inputLP;

        //TODO: sort list again.
    }

    // dump output scanlines.
    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFWriteScanline(outputTIFF,
                &(outputBuf[i * OutputWidth]),
                i,
                8);
    }

    // close outputTIFF.
    TIFFClose(outputTIFF);

    return 0;
}
