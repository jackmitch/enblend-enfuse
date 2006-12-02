/*
 * Copyright (C) 2004-2005 Andrew Mihal
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

#ifdef _WIN32
#include <win32helpers\win32config.h>
#endif // _WIN32

// Defines lrint for fast fromRealPromotes
#include "float_cast.h"

#ifdef _WIN32
// Make sure we bring in windows.h the right way
#define _STLP_VERBOSE_AUTO_LINK
#define _USE_MATH_DEFINES
#define NOMINMAX
#define VC_EXTRALEAN 
#include <windows.h>
#undef DIFFERENCE
#endif  //  _WIN32

#ifdef __GW32C__
#undef malloc
#define BOOST_NO_STDC_NAMESPACE 1
#endif

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

#ifndef _WIN32
#include <getopt.h>
#else
extern "C" int getopt(int nargc, char** nargv, char* ostr);
#endif

extern "C" char *optarg;
extern "C" int optind;

#ifndef _WIN32
#include <fenv.h>
#endif

#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <tiffconf.h>

#ifdef _WIN32
#include <io.h>
#endif

#include <boost/random/mersenne_twister.hpp>
#include <lcms.h>

// Globals

// Random number generator for dithering
boost::mt19937 Twister;

// Global values from command line parameters.
int Verbose = 1;
unsigned int ExactLevels = 0;
bool OneAtATime = true;
bool Wraparound = false;
double StitchMismatchThreshold = 0.4;
bool GimpAssociatedAlphaHack = false;
bool UseCIECAM = false;
bool UseLZW = false;
bool OutputSizeGiven = false;
int OutputWidthCmdLine = 0;
int OutputHeightCmdLine = 0;
bool Checkpoint = false;
bool UseGPU = false;

// Objects for ICC profiles
cmsHPROFILE InputProfile = NULL;
cmsHPROFILE XYZProfile = NULL;
cmsHTRANSFORM InputToXYZTransform = NULL;
cmsHTRANSFORM XYZToInputTransform = NULL;
cmsViewingConditions ViewingConditions;
LCMSHANDLE CIECAMTransform = NULL;

#include "common.h"
#include "enblend.h"

#include "vigra/impex.hxx"
#include "vigra/sized_int.hxx"

#include <tiffio.h>
using std::cerr;
using std::cout;
using std::endl;
using std::list;

using vigra::CachedFileImageDirector;
using vigra::Diff2D;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::Rect2D;
using vigra::StdException;

using enblend::enblendMain;

#ifdef _WIN32
#define strdup _strdup
#endif

/** Print the usage information and quit. */
void printUsageAndExit() {
    cout << "==== enblend, version " << VERSION << " ====" << endl;
    cout << "Usage: enblend [options] -o OUTPUT INPUTS" << endl;
    cout << endl;
    cout << "Common options:" << endl;
    cout << " -a                Pre-assemble non-overlapping images" << endl;
    cout << " -h                Print this help message" << endl;
    cout << " -l number         Number of levels to use (1 to 29)" << endl;
    cout << " -o filename       Write output to file" << endl;
    cout << " -v                Verbose" << endl;
    cout << " -w                Blend across -180/+180 boundary" << endl;
    cout << " -z                Use LZW compression" << endl;
    cout << " -x                Checkpoint partial results" << endl;

    cout << endl << "Extended options:" << endl;
    cout << " -b kilobytes      Image cache block size (default=2MiB)" << endl;
    cout << " -c                Use CIECAM02 to blend colors" << endl;
    cout << " -g                Associated alpha hack for Gimp (ver. < 2) and Cinepaint" << endl;
    cout << " -f WIDTHxHEIGHT   Manually set the size of the output image." << endl
         << "                   Useful for cropped and shifted input TIFF images," << endl
         << "                   such as those produced by Nona." << endl;
    cout << " -m megabytes      Use this much memory before going to disk (default=1GiB)" << endl;

    //TODO stitch mismatch avoidance is work-in-progress.
    //cout << " -t float          Stitch mismatch threshold, [0.0, 1.0]" << endl;

    // deprecated
    //cout << " -s                Blend images one at a time, in the order given" << endl;
    exit(1);
}

/** Make sure all cached file images get destroyed,
 *  and hence the temporary files deleted,
 *  if we are killed.
 */
void sigint_handler(int sig) {
    CachedFileImageDirector::v().~CachedFileImageDirector();
    #if !defined(__GW32C__) && !defined(_WIN32)
    signal(SIGINT, SIG_DFL);
    kill(getpid(), SIGINT);
    #else
    exit(0);
    #endif
}

int main(int argc, char** argv) {

#ifdef _WIN32
    // Make sure the FPU is set to rounding mode so that the lrint
    // functions in float_cast.h will work properly.
    // See changes in vigra numerictraits.hxx
    _controlfp( _RC_NEAR, _MCW_RC );
#else
    fesetround(FE_TONEAREST);
#endif
    
    signal(SIGINT, sigint_handler);

    // Make sure libtiff is compiled with TIF_PLATFORM_CONSOLE
    // to avoid interactive warning dialogs.
    //TIFFSetWarningHandler(NULL);
    //TIFFSetErrorHandler(NULL);

    // The name of the output file.
    char *outputFileName = NULL;

    // List of input files.
    list<char*> inputFileNameList;
    list<char*>::iterator inputFileNameIterator;

    // Parse command line.
    int c;
    while ((c = getopt(argc, argv, "ab:cf:ghl:m:o:st:vwxyz")) != -1) {
        switch (c) {
            case 'a': {
                OneAtATime = false;
                break;
            }
            case 'b': {
                int kilobytes = atoi(optarg);
                if (kilobytes < 1) {
                    cerr << "enblend: cache block size must be 1 or more."
                         << endl;
                    printUsageAndExit();
                }
                CachedFileImageDirector::v().setBlockSize(
                        (long long)kilobytes << 10);
                break;
            }
            case 'c': {
                UseCIECAM = true;
                break;
            }
            case 'f': {
                OutputSizeGiven = true;
                int nP = sscanf(optarg, "%dx%d", &OutputWidthCmdLine, &OutputHeightCmdLine);
                if (nP != 2) {
                    cerr << "enblend: the -f option requires a parameter of the form WIDTHxHEIGHT"
                         << endl;
                    printUsageAndExit();
                }
                break;
            }
            case 'g': {
                GimpAssociatedAlphaHack = true;
                break;
            }
            case 'h': {
                printUsageAndExit();
                break;
            }
            case 'l': {
                int levels = atoi(optarg);
                if (levels < 1 || levels > 29) {
                    cerr << "enblend: levels must in the range 1 to 29."
                         << endl;
                    printUsageAndExit();
                }
                ExactLevels = (unsigned int)levels;
                break;
            }
            case 'm': {
                int megabytes = atoi(optarg);
                if (megabytes < 1) {
                    cerr << "enblend: memory limit must be 1 or more."
                         << endl;
                    printUsageAndExit();
                }
                CachedFileImageDirector::v().setAllocation(
                        (long long)megabytes << 20);
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
                try {
                    outputFileName = new char[len];
                } catch (std::bad_alloc& e) {
                    cerr << endl << "enblend: out of memory"
                         << endl << e.what()
                         << endl;
                    exit(1);
                }
                strncpy(outputFileName, optarg, len);
                break;
            }
            case 's': {
                // Deprecated sequential blending flag.
                OneAtATime = true;
                cerr << "enblend: the -s flag is deprecated." << endl;
                break;
            }
            case 't': {
                printUsageAndExit();
                //StitchMismatchThreshold = strtod(optarg, NULL);
                //if (StitchMismatchThreshold < 0.0
                //        || StitchMismatchThreshold > 1.0) {
                //    cerr << "enblend: threshold must be between "
                //         << "0.0 and 1.0 inclusive."
                //         << endl;
                //    printUsageAndExit();
                //}
                break;
            }
            case 'v': {
                Verbose++;
                break;
            }
            case 'w': {
                Wraparound = true;
                break;
            }
            case 'x': {
                Checkpoint = true;
                break;
            }
            case 'y': {
                UseGPU = true;
                break;
            }
            case 'z': {
                UseLZW = true;
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
#ifdef _WIN32
            // There has got to be an easier way...
            char drive[_MAX_DRIVE];
            char dir[_MAX_DIR];
            char fname[_MAX_FNAME];
            char ext[_MAX_EXT];
            char newFile[_MAX_PATH];

            _splitpath(argv[optind], drive, dir, NULL, NULL);

            struct _finddata_t finddata;
            intptr_t findhandle;
            int stop = 0;

            findhandle = _findfirst(argv[optind], &finddata);
            if (findhandle != -1) {
                do {
                    _splitpath(finddata.name, NULL, NULL, fname, ext);
                    _makepath(newFile, drive, dir, fname, ext);
                    
                    // TODO (jbeda): This will leak -- the right way to 
                    // fix this is to make this a list of std::string.  
                    // I'll look into this after we get things working
                    // on Win32
                    inputFileNameList.push_back(strdup(newFile));
                } while (_findnext(findhandle, &finddata) == 0);
                _findclose(findhandle);
            }

            optind++;
#else
            inputFileNameList.push_back(argv[optind++]);
#endif
        }
    } else {
        cerr << "enblend: no input files specified." << endl;
        printUsageAndExit();
    }

    // Check that more than one input file was given.
    if (inputFileNameList.size() <= 1) {
        cerr << "enblend: only one input file given. "
             << "Enblend needs two or more overlapping input images in order "
             << "to do blending calculations. The output will be the same as "
             << "the input."
             << endl;
    }

    // List of info structures for each input image.
    list<ImageImportInfo*> imageInfoList;
    list<ImageImportInfo*>::iterator imageInfoIterator;

    bool isColor = false;
    const char *pixelType = NULL;
    ImageImportInfo::ICCProfile iccProfile;
    Rect2D inputUnion;

    // Check that all input images have the same parameters.
    inputFileNameIterator = inputFileNameList.begin();
    while (inputFileNameIterator != inputFileNameList.end()) {

        ImageImportInfo *inputInfo = NULL;
        try {
            inputInfo = new ImageImportInfo(*inputFileNameIterator);
        } catch (StdException& e) {
            cerr << endl << "enblend: error opening input file \""
                 << *inputFileNameIterator << "\":"
                 << endl << e.what()
                 << endl;
            exit(1);
        }

        // Save this image info in the list.
        imageInfoList.push_back(inputInfo);

        if (Verbose > VERBOSE_INPUT_IMAGE_INFO_MESSAGES) {
            cout << "Input image \""
                 << *inputFileNameIterator
                 << "\" ";

            if (inputInfo->isColor()) cout << "RGB ";

            if (!inputInfo->getICCProfile().empty()) cout << "ICC ";

            cout << inputInfo->getPixelType() << " "
                 << "position="
                 << inputInfo->getPosition().x
                 << "x"
                 << inputInfo->getPosition().y
                 << " "
                 << "size="
                 << inputInfo->width()
                 << "x"
                 << inputInfo->height()
                 << endl;
        }

        if (inputInfo->numExtraBands() < 1) {
            // Complain about lack of alpha channel.
            cerr << "enblend: Input image \""
                 << *inputFileNameIterator << "\" does not have an alpha "
                 << "channel. This is required to determine which pixels "
                 << "contribute to the final image."
                 << endl;
            exit(1);
        }

        // Get input image's position and size.
        Rect2D imageROI(Point2D(inputInfo->getPosition()),
                Size2D(inputInfo->width(), inputInfo->height()));

        if (inputFileNameIterator == inputFileNameList.begin()) {
            // The first input image.
            inputUnion = imageROI;
            isColor = inputInfo->isColor();
            pixelType = inputInfo->getPixelType();
            iccProfile = inputInfo->getICCProfile();
            if (!iccProfile.empty()) {
                InputProfile = cmsOpenProfileFromMem(iccProfile.data(), iccProfile.size());
                if (InputProfile == NULL) {
                    cerr << endl << "enblend: error parsing ICC profile data from file\""
                         << *inputFileNameIterator
                         << "\"" << endl;
                    exit(1);
                }
            }
        }
        else {
            // second and later images.
            inputUnion |= imageROI;

            if (isColor != inputInfo->isColor()) {
                cerr << "enblend: Input image \""
                     << *inputFileNameIterator << "\" is "
                     << (inputInfo->isColor() ? "color" : "grayscale")
                     << " but previous images are "
                     << (isColor ? "color" : "grayscale")
                     << "." << endl;
                exit(1);
            }
            if (strcmp(pixelType, inputInfo->getPixelType())) {
                cerr << "enblend: Input image \""
                     << *inputFileNameIterator << "\" has pixel type "
                     << inputInfo->getPixelType()
                     << " but previous images have pixel type "
                     << pixelType
                     << "." << endl;
                exit(1);
            }
            if (!std::equal(iccProfile.begin(), iccProfile.end(), inputInfo->getICCProfile().begin())) {
                ImageImportInfo::ICCProfile mismatchProfile = inputInfo->getICCProfile();
                cmsHPROFILE newProfile = NULL;
                if (!mismatchProfile.empty()) {
                    newProfile = cmsOpenProfileFromMem(mismatchProfile.data(), mismatchProfile.size());
                    if (newProfile == NULL) {
                        cerr << endl << "enblend: error parsing ICC profile data from file\""
                             << *inputFileNameIterator
                             << "\"" << endl;
                        exit(1);
                    }
                }

                cerr << endl << "enblend: Input image \""
                     << *inputFileNameIterator
                     << "\" has ";
                if (newProfile) {
                    cerr << " ICC profile \""
                         << cmsTakeProductName(newProfile)
                         << " "
                         << cmsTakeProductDesc(newProfile)
                         << "\"";
                } else {
                    cerr << " no ICC profile";
                }
                cerr << " but previous images have ";
                if (InputProfile) {
                    cerr << " ICC profile \""
                         << cmsTakeProductName(InputProfile)
                         << " "
                         << cmsTakeProductDesc(InputProfile)
                         << "\"." << endl;
                } else {
                    cerr << " no ICC profile." << endl;
                }
                cerr << "enblend: Blending images with different color spaces may have unexpected results."
                     << endl;

            }
        }

        inputFileNameIterator++;
    }

    // Make sure that inputUnion is at least as big as given by the -f paramater.
    if (OutputSizeGiven) {
        inputUnion |= Rect2D(Size2D(OutputWidthCmdLine, OutputHeightCmdLine));
    }

    // Create the Info for the output file.
    ImageExportInfo outputImageInfo(outputFileName);
    if (UseLZW) outputImageInfo.setCompression("LZW");

    // Pixel type of the output image is the same as the input images.
    outputImageInfo.setPixelType(pixelType);

    // Set the output image ICC profile
    outputImageInfo.setICCProfile(iccProfile);

    if (UseCIECAM) {
        if (InputProfile == NULL) {
            cerr << "enblend: Input images do not have ICC profiles. Assuming sRGB." << endl;
            InputProfile = cmsCreate_sRGBProfile();
        }
        XYZProfile = cmsCreateXYZProfile();

        InputToXYZTransform = cmsCreateTransform(InputProfile, TYPE_RGB_DBL,
                                                 XYZProfile, TYPE_XYZ_DBL,
                                                 INTENT_PERCEPTUAL, cmsFLAGS_NOTPRECALC);
        if (InputToXYZTransform == NULL) {
            cerr << "enblend: Error building color transform from \""
                 << cmsTakeProductName(InputProfile)
                 << " "
                 << cmsTakeProductDesc(InputProfile)
                 << "\" to XYZ." << endl;
            exit(1);
        }

        XYZToInputTransform = cmsCreateTransform(XYZProfile, TYPE_XYZ_DBL,
                                                 InputProfile, TYPE_RGB_DBL,
                                                 INTENT_PERCEPTUAL, cmsFLAGS_NOTPRECALC);
        if (XYZToInputTransform == NULL) {
            cerr << "enblend: Error building color transform from XYZ to \""
                 << cmsTakeProductName(InputProfile)
                 << " "
                 << cmsTakeProductDesc(InputProfile)
                 << "\"." << endl;
            exit(1);
        }

        // P2 Viewing Conditions: D50, 500 lumens
        ViewingConditions.whitePoint.X = 96.42;
        ViewingConditions.whitePoint.Y = 100.0;
        ViewingConditions.whitePoint.Z = 82.49;
        ViewingConditions.Yb = 20.0;
        ViewingConditions.La = 31.83;
        ViewingConditions.surround = AVG_SURROUND;
        ViewingConditions.D_value = 1.0;

        CIECAMTransform = cmsCIECAM02Init(&ViewingConditions);
        if (!CIECAMTransform) {
            cerr << endl << "enblend: Error initializing CIECAM02 transform." << endl;
            exit(1);
        }
    }

    // The size of the output image.
    if (Verbose > VERBOSE_INPUT_UNION_SIZE_MESSAGES) {
        cout << "Output image size: " << inputUnion << endl;
    }

    // Set the output image position and resolution.
    outputImageInfo.setXResolution(300.0);
    outputImageInfo.setYResolution(300.0);
    outputImageInfo.setPosition(inputUnion.upperLeft());

    // Sanity check on the output image file.
    try {
        // This seems to be a reasonable way to check if
        // the output file is going to work after blending
        // is done.
        encoder(outputImageInfo);
    } catch (StdException & e) {
        cerr << endl << "enblend: error opening output file \""
             << outputFileName
             << "\":"
             << endl << e.what()
             << endl;
        exit(1);
    }

    // Invoke templatized blender.
    try {
        if (isColor) {
            if (strcmp(pixelType,      "UINT8" ) == 0) enblendMain<RGBValue<UInt8 > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT8"  ) == 0) enblendMain<RGBValue<Int8  > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT16") == 0) enblendMain<RGBValue<UInt16> >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT16" ) == 0) enblendMain<RGBValue<Int16 > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT32") == 0) enblendMain<RGBValue<UInt32> >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT32" ) == 0) enblendMain<RGBValue<Int32 > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT64") == 0) enblendMain<RGBValue<UInt64> >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT64" ) == 0) enblendMain<RGBValue<Int64 > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "FLOAT" ) == 0) enblendMain<RGBValue<float > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "DOUBLE") == 0) enblendMain<RGBValue<double> >(imageInfoList, outputImageInfo, inputUnion);
            else {
                cerr << "enblend: images with pixel type \""
                     << pixelType
                     << "\" are not supported."
                     << endl;
                exit(1);
            }
        } else {
            //if (strcmp(pixelType,      "UINT8" ) == 0) enblendMain<UInt8 >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT8"  ) == 0) enblendMain<Int8  >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT16") == 0) enblendMain<UInt16>(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT16" ) == 0) enblendMain<Int16 >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT32") == 0) enblendMain<UInt32>(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT32" ) == 0) enblendMain<Int32 >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT64") == 0) enblendMain<UInt64>(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT64" ) == 0) enblendMain<Int64 >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "FLOAT" ) == 0) enblendMain<float >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "DOUBLE") == 0) enblendMain<double>(imageInfoList, outputImageInfo, inputUnion);
            //else {
                cerr << "enblend: images with pixel type \""
                     << pixelType
                     << "\" are not supported."
                     << endl;
                exit(1);
            //}
        }

        // delete entries in imageInfoList, in case
        // enblend loop returned early.
        imageInfoIterator = imageInfoList.begin();
        while (imageInfoIterator != imageInfoList.end()) {
            delete *imageInfoIterator++;
        }

        delete[] outputFileName;

    } catch (std::bad_alloc& e) {
        cerr << endl << "enblend: out of memory"
             << endl << e.what()
             << endl;
        exit(1);
    } catch (StdException& e) {
        cerr << endl << "enblend: an exception occured"
             << endl << e.what()
             << endl;
        exit(1);
    }

    if (CIECAMTransform) cmsCIECAM02Done(CIECAMTransform);
    if (InputToXYZTransform) cmsDeleteTransform(InputToXYZTransform);
    if (XYZToInputTransform) cmsDeleteTransform(XYZToInputTransform);
    if (XYZProfile) cmsCloseProfile(XYZProfile);
    if (InputProfile) cmsCloseProfile(InputProfile);

    // Success.
    return 0;
}
