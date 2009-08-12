/*
 * Copyright (C) 2004-2007 Andrew Mihal
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

#ifdef _MSC_VER
#include <win32helpers\win32config.h>
#define isnan _isnan
#endif // _MSC_VER

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
//extern "C" int getopt(int nargc, char** nargv, char* ostr);
extern "C" {
#include <win32helpers/getopt_long.h>
}
#endif

extern "C" char *optarg;
extern "C" int optind;

#ifndef _MSC_VER
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
std::string OutputFileName;
unsigned int ExactLevels = 0;
bool OneAtATime = true;
bool Wraparound = false;
bool GimpAssociatedAlphaHack = false;
bool UseCIECAM = false;
bool OutputSizeGiven = false;
int OutputWidthCmdLine = 0;
int OutputHeightCmdLine = 0;
int OutputOffsetXCmdLine = 0;
int OutputOffsetYCmdLine = 0;
bool Checkpoint = false;
bool UseGPU = false;
bool OptimizeMask = true;
bool CoarseMask = true;
std::string SaveMaskFileName;
std::string LoadMaskFileName;
std::string VisualizeMaskFileName;
unsigned int GDAKmax = 32;
unsigned int DijkstraRadius = 25;
unsigned int MaskVectorizeDistance = 0;
std::string OutputCompression;

// Globals related to catching SIGINT
#ifndef _WIN32
sigset_t SigintMask;
#endif

// Objects for ICC profiles
cmsHPROFILE InputProfile = NULL;
cmsHPROFILE XYZProfile = NULL;
cmsHTRANSFORM InputToXYZTransform = NULL;
cmsHTRANSFORM XYZToInputTransform = NULL;
cmsViewingConditions ViewingConditions;
LCMSHANDLE CIECAMTransform = NULL;

#include "common.h"
#include "enblend.h"
#ifdef HAVE_LIBGLEW
#include "gpu.h"
#endif

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


/** Print information on the current version and some configuration
 * details. */
void printVersionAndExit() {
    cout <<
        "enblend " << VERSION << "\n" <<
#ifdef ENBLEND_CACHE_IMAGES
        "Extra feature: image cache\n" <<
#endif
#ifdef HAVE_LIBGLEW
        "Extra feature: GPU acceleration\n" <<
#endif
        "\n" <<
        "Copyright (C) 2004-2008 Andrew Mihal.\n" <<
        "License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n" <<
        "This is free software: you are free to change and redistribute it.\n" <<
        "There is NO WARRANTY, to the extent permitted by law.\n" <<
        "\n" <<
        "Written by Andrew Mihal and others." <<
        endl;

    exit(0);
}


/** Print the usage information and quit. */
void printUsageAndExit(const bool error = true) {
    cout <<
        "Usage: enblend [options] -o OUTPUT INPUT...\n" <<
        "Blend INPUT images into a single OUTPUT image.\n" <<
        "\n" <<
        "Common options:\n" <<
        "  -V, --version          Output version information and exit\n" <<
        "  -a                     Pre-assemble non-overlapping images\n" <<
        "  -h, --help             Print this help message and exit\n" <<
        "  -l LEVELS              Number of blending levels to use (1 to 29)\n" <<
        "  -o FILENAME            Write output to FILENAME\n" <<
        "  -v, --verbose          Verbosely report progress; repeat to\n" <<
        "                            increase verbosity\n" <<
        "  -w                     Blend across -180/+180 degrees boundary\n" <<
        "  -x                     Checkpoint partial results\n" <<
        "  -z                     Use LZW compression (TIFF only).\n" <<
        "                           Kept for backward compatability with older scripts\n" <<
        "  --compression=COMP     Set compression of output image to COMP,\n" <<
        "                           where COMP is:\n" <<
        "                             NONE, PACKBITS, LZW, DEFLATE (for TIFF files)\n" <<
        "                             0-100 (for JPEG files)\n" <<
        "\n" <<
        "Extended options:\n" <<
        "  -b BLOCKSIZE           Image cache BLOCKSIZE in Kilobytes.  Default: " <<
        (CachedFileImageDirector::v().getBlockSize() / 1024LL) << "KB\n" <<
        "  -c                     Use CIECAM02 to blend colors\n" <<
        "  -g                     Associated-alpha hack for Gimp (before version 2)\n" <<
        "                           and Cinepaint\n" <<
#ifdef HAVE_LIBGLEW
        "  --gpu                  Use the graphics card to accelerate some computations\n" <<
#endif
        "  -f WIDTHxHEIGHT[+xXOFFSET+yYOFFSET]\n" <<
        "                         Manually set the size and position of the output\n" <<
        "                           image.  Useful for cropped and shifted input\n" <<
        "                           TIFF images, such as those produced by Nona.\n" <<
        "  -m CACHESIZE           Images CACHESIZE in Megabytes.  Default: " <<
        (CachedFileImageDirector::v().getAllocation() / 1048576LL) << "MB\n" <<
        "  --visualize=FILENAME   Save results of optimizer FILENAME\n" <<
        "\n" <<
        "Mask generation options:\n" <<
        "  --coarse-mask          Use an approximation to speedup mask generation.  Default\n" <<
        "  --fine-mask            Enable detailed mask generation.  Slow!  Use if\n" <<
        "                           overlap regions are very narrow.\n" <<
        "  --optimize             Turn on mask optimization.  This is the default.\n" <<
        "  --no-optimize          Turn off mask optimization.\n" <<
        "  --save-mask=FILENAME   Save the generated mask to FILENAME.\n" <<
        "  --load-mask=FILENAME   Use the mask in FILENAME instead of generating one.\n" <<
        endl;

    exit(error ? 1 : 0);
}

/** Make sure all cached file images get destroyed,
 *  and hence the temporary files deleted,
 *  if we are killed.
 */
void sigint_handler(int sig) {
    cout << endl << "Interrupted." << endl;
    // FIXME what if this occurs in a CFI atomic section?
    // This is no longer necessary, temp files are unlinked during creation.
    //CachedFileImageDirector::v().~CachedFileImageDirector();
#ifdef HAVE_LIBGLEW
    if (UseGPU) {
        // FIXME what if this occurs in a GL atomic section?
        //cout << "Cleaning up GPU state..." << endl;
        wrapupGPU();
    }
#endif
    #if !defined(__GW32C__) && !defined(_WIN32)
    struct sigaction action;
    action.sa_handler = SIG_DFL;
    sigemptyset(&(action.sa_mask));
    sigaction(SIGINT, &action, NULL);
    raise(SIGINT);
    #else
    exit(0);
    #endif
}


int process_options(int argc, char** argv) {
    enum OptionArgumentKind {
        NoArgument,
        StringArgument,
        FloatArgument,
        IntegerArgument
    };

    // NOTE: An OptionId is the index of a command line option in
    // "long_options".  Every change in "enum OptionId" must be
    // reflected in "long_options" and vice versa.
    enum OptionId {
        UseGpuId,                   //  0
        CoarseMaskId,               //  1
        FineMaskId,                 //  2
        OptimizeMaskId,             //  3
        NoOptimizeMaskId,           //  4
        SaveMaskId,                 //  5
        LoadMaskId,                 //  6
        VisualizeId,                //  7
        GdaKmaxId,                  //  8
        DijkstraRadiusId,           //  9
        MaskVectorizeDistanceId,    // 10
        CompressionId,              // 11
        VerboseId,                  // 12
        HelpId,                     // 13
        VersionId                   // 14
    };

    // NOTE: See note attached to "enum OptionId" above.
    static struct option long_options[] = {
        {"gpu", no_argument, 0, NoArgument},                                   //  0
        {"coarse-mask", no_argument, 0, NoArgument},                           //  1
        {"fine-mask", no_argument, 0, NoArgument},                             //  2
        {"optimize", no_argument, 0, NoArgument},                              //  3
        {"no-optimize", no_argument, 0, NoArgument},                           //  4
        {"save-mask", required_argument, 0, StringArgument},                   //  5
        {"load-mask", required_argument, 0, StringArgument},                   //  6
        {"visualize", required_argument, 0, StringArgument},                   //  7
        {"gda-kmax", required_argument, 0, IntegerArgument},                   //  8
        {"dijkstra-radius", required_argument, 0, IntegerArgument},            //  9
        {"mask-vectorize-distance", required_argument, 0, IntegerArgument},    // 10
        {"compression", required_argument, 0, StringArgument},                 // 11
        {"verbose", no_argument, 0, NoArgument},                               // 12
        {"help", no_argument, 0, NoArgument},                                  // 13
        {"version", no_argument, 0, NoArgument},                               // 14
        {0, 0, 0, 0}
    };

    // Parse command line.
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "Vab:cf:ghl:m:o:svwxz",
                            long_options, &option_index)) != -1) {
        switch (c) {
        case NoArgument: {
            if (long_options[option_index].flag != 0) break;
            switch (option_index) {
            case UseGpuId:
                UseGPU = true;
                break;
            case CoarseMaskId:
                CoarseMask = true;
                break;
            case FineMaskId:
                CoarseMask = false;
                break;
            case OptimizeMaskId:
                OptimizeMask = true;
                break;
            case NoOptimizeMaskId:
                OptimizeMask = false;
                break;
            case VerboseId:
                Verbose++;
                break;
            case HelpId:
                printUsageAndExit(false);
                break;          // never reached
            case VersionId:
                printVersionAndExit();
                break;          // never reached
            default:
                cerr << "enblend: internal error: unhandled \"NoArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case NoArgument"

        case StringArgument: {
            if (long_options[option_index].flag != 0) break;
            switch (option_index) {
            case SaveMaskId:
                SaveMaskFileName = optarg;
                break;
            case LoadMaskId:
                LoadMaskFileName = optarg;
                break;
            case VisualizeId:
                VisualizeMaskFileName = optarg;
                break;
            case CompressionId:
                OutputCompression = optarg;
                break;
            default:
                cerr << "enblend: internal error: unhandled \"StringArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case StringArgument"

        // case FloatArgument: {
        // ...
        // } // end of "case FloatArgument"

        case IntegerArgument: {
            if (long_options[option_index].flag != 0) break;
            unsigned int* optionUInt = NULL;
            switch (option_index) {
            case GdaKmaxId:
                optionUInt = &GDAKmax;
                break;
            case DijkstraRadiusId:
                optionUInt = &DijkstraRadius;
                break;
            case MaskVectorizeDistanceId:
                optionUInt = &MaskVectorizeDistance;
                break;
            default:
                cerr << "enblend: internal error: unhandled \"IntegerArgument\" option"
                     << endl;
                exit(1);
            }

            const int value = atoi(optarg);
            if (value < 1) {
                cerr << "enblend: " << long_options[option_index].name
                     << " must be 1 or more." << endl;
                printUsageAndExit();
            }
            *optionUInt = static_cast<unsigned int>(value);
            break;
        } // end of "case IntegerArgument"

        case 'V':
            printVersionAndExit();
            break;          // never reached
        case 'a':
            OneAtATime = false;
            break;
        case 'b': {
            const int kilobytes = atoi(optarg);
            if (kilobytes < 1) {
                cerr << "enblend: cache block size must be 1 or more." << endl;
                printUsageAndExit();
            }
            CachedFileImageDirector::v().setBlockSize(static_cast<long long>(kilobytes << 10));
            break;
        }
        case 'c':
            UseCIECAM = true;
            break;
        case 'f': {
            OutputSizeGiven = true;
            const int nP = sscanf(optarg, "%dx%d+%d+%d",
                                  &OutputWidthCmdLine, &OutputHeightCmdLine,
                                  &OutputOffsetXCmdLine, &OutputOffsetYCmdLine);
            if (nP == 4) {
                ; // full geometry string
            } else if (nP == 2) {
                OutputOffsetXCmdLine=0;
                OutputOffsetYCmdLine=0;
            } else {
                cerr << "enblend: the -f option requires a parameter "
                     << "of the form WIDTHxHEIGHT+X0+Y0 or WIDTHxHEIGHT" << endl;
                printUsageAndExit();
            }
            break;
        }
        case 'g':
            GimpAssociatedAlphaHack = true;
            break;
        case 'h':
            printUsageAndExit(false);
            break;
        case 'l': {
            const int levels = atoi(optarg);
            if (levels < 1 || levels > 29) {
                cerr << "enblend: levels must in the range 1 to 29." << endl;
                printUsageAndExit();
            }
            ExactLevels = static_cast<unsigned int>(levels);
            break;
        }
        case 'm': {
            const int megabytes = atoi(optarg);
            if (megabytes < 1) {
                cerr << "enblend: memory limit must be 1 or more." << endl;
                printUsageAndExit();
            }
            CachedFileImageDirector::v().setAllocation(static_cast<long long>(megabytes) << 20);
            break;
        }
        case 'o':
            if (!OutputFileName.empty()) {
                cerr << "enblend: more than one output file specified." << endl;
                printUsageAndExit();
                break;
            }
            OutputFileName = optarg;
            break;
        case 's':
            // Deprecated sequential blending flag.
            OneAtATime = true;
            cerr << "enblend: warning: Flag \"-s\" is deprecated." << endl;
            break;
        case 'v':
            Verbose++;
            break;
        case 'w':
            Wraparound = true;
            break;
        case 'x':
            Checkpoint = true;
            break;
        case 'z':
            OutputCompression = "LZW";
            break;

        default:
            printUsageAndExit();
            break;
        }
    }

    return optind;
}


int main(int argc, char** argv) {
#ifdef _MSC_VER
    // Make sure the FPU is set to rounding mode so that the lrint
    // functions in float_cast.h will work properly.
    // See changes in vigra numerictraits.hxx
    _controlfp(_RC_NEAR, _MCW_RC);
#else
    fesetround(FE_TONEAREST);
#endif

#ifndef _WIN32
    sigemptyset(&SigintMask);
    sigaddset(&SigintMask, SIGINT);

    struct sigaction action;
    action.sa_handler = sigint_handler;
    sigemptyset(&(action.sa_mask));
    sigaction(SIGINT, &action, NULL);
#else
    signal(SIGINT, sigint_handler);
#endif

    // Make sure libtiff is compiled with TIF_PLATFORM_CONSOLE
    // to avoid interactive warning dialogs.
    //TIFFSetWarningHandler(NULL);
    //TIFFSetErrorHandler(NULL);

    // List of input files.
    list<char*> inputFileNameList;
    list<char*>::iterator inputFileNameIterator;

    int optind;
    try {optind = process_options(argc, argv);}
    catch (StdException& e) {
        cerr << "enblend: error while processing command line options\n"
             << "enblend:     " << e.what()
             << endl;
        exit(1);
    }

    // Make sure mandatory output file name parameter given.
    if (OutputFileName.empty()) {
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

#ifdef HAVE_LIBGLEW
    if (UseGPU) {
      initGPU(&argc,argv);
    }
#endif

    //if (CachedFileImageDirector::v()->getManagedBlocks() < 4) {
    //    // Max simultaneous image access is in:
    //    // 4 in any of many calls to combineThreeImages
    //    // 4 gaussian pyramid init (src image layer, src alpha layer, dest pyramid image layer 0, dest pyramid alpha layer 0)
    //    // 4 in reduce (src image layer N, src alpha layer N, dest image layer N+1, dest alpha layer N+1)
    //    // FIXME complain or automatically adjust blocksize to get ManagedBlocks above 4?
    //}

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
    const char* pixelType = NULL;
    ImageImportInfo::ICCProfile iccProfile;
    Rect2D inputUnion;

    // Check that all input images have the same parameters.
    inputFileNameIterator = inputFileNameList.begin();
    int minDim = INT_MAX;
    while (inputFileNameIterator != inputFileNameList.end()) {

        ImageImportInfo* inputInfo = NULL;
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
            minDim = std::min(inputInfo->width(), inputInfo->height());
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
            if (!std::equal(iccProfile.begin(), iccProfile.end(),
                            inputInfo->getICCProfile().begin())) {
                ImageImportInfo::ICCProfile mismatchProfile = inputInfo->getICCProfile();
                cmsHPROFILE newProfile = NULL;
                if (!mismatchProfile.empty()) {
                    newProfile = cmsOpenProfileFromMem(mismatchProfile.data(),
                                                       mismatchProfile.size());
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
            if (inputInfo->width() < minDim) {
                minDim = inputInfo->width();
            }
            if (inputInfo->height() < minDim) {
                minDim = inputInfo->height();
            }
        }

        inputFileNameIterator++;
    }

    // Switch to fine mask, if the smallest coarse mask would be less
    // than 64 pixels wide or high.
    if (minDim / 8 < 64 && CoarseMask) {
        cout << "Input images to small for coarse mask, switching to fine mask." << endl;
        CoarseMask = false;
    }

    if (MaskVectorizeDistance == 0) {
        MaskVectorizeDistance = CoarseMask ? 4 : 20;
    }

    // Make sure that inputUnion is at least as big as given by the -f paramater.
    if (OutputSizeGiven) {
        inputUnion |= Rect2D(OutputOffsetXCmdLine, OutputOffsetYCmdLine,
                             OutputOffsetXCmdLine + OutputWidthCmdLine,
                             OutputOffsetYCmdLine + OutputHeightCmdLine);
    }

    // Create the Info for the output file.
    ImageExportInfo outputImageInfo(OutputFileName.c_str());
    if (!OutputCompression.empty()) {
        outputImageInfo.setCompression(OutputCompression.c_str());
    }

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
             << OutputFileName
             << "\":"
             << endl << e.what()
             << endl;
        exit(1);
    }

    // Sanity check on the LoadMaskFileName
    if (!LoadMaskFileName.empty()) try {
        ImageImportInfo maskInfo(LoadMaskFileName.c_str());
    } catch (StdException& e) {
        cerr << endl << "enblend: error opening load-mask input file \""
             << LoadMaskFileName << "\":"
             << endl << e.what()
             << endl;
        exit(1);
    }

    // Sanity check on the SaveMaskFileName
    if (!SaveMaskFileName.empty()) try {
        ImageExportInfo maskInfo(SaveMaskFileName.c_str());
        encoder(maskInfo);
    } catch (StdException& e) {
        cerr << endl << "enblend: error opening save-mask output file \""
             << SaveMaskFileName << "\":"
             << endl << e.what()
             << endl;
        exit(1);
    }

    // Sanity check on the VisualizeMaskFileName
    if (!VisualizeMaskFileName.empty()) try {
        ImageExportInfo maskInfo(VisualizeMaskFileName.c_str());
        encoder(maskInfo);
    } catch (StdException& e) {
        cerr << endl << "enblend: error opening visualize output file \""
             << VisualizeMaskFileName << "\":"
             << endl << e.what()
             << endl;
        exit(1);
    }

    if (!VisualizeMaskFileName.empty() && !OptimizeMask) {
        cerr << endl << "enblend: --visualize does nothing without --optimize."
             << endl;
    }

    // Invoke templatized blender.
    try {
        if (isColor) {
            if (strcmp(pixelType,      "UINT8" ) == 0) enblendMain<RGBValue<UInt8 > >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "FLOAT" ) == 0) enblendMain<RGBValue<float > >(imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (strcmp(pixelType, "INT8"  ) == 0) enblendMain<RGBValue<Int8  > >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "UINT16") == 0) enblendMain<RGBValue<UInt16> >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "INT16" ) == 0) enblendMain<RGBValue<Int16 > >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "UINT32") == 0) enblendMain<RGBValue<UInt32> >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "INT32" ) == 0) enblendMain<RGBValue<Int32 > >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT64") == 0) enblendMain<RGBValue<UInt64> >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT64" ) == 0) enblendMain<RGBValue<Int64 > >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "DOUBLE") == 0) enblendMain<RGBValue<double> >(imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                cerr << "enblend: images with pixel type \""
                     << pixelType
                     << "\" are not supported."
                     << endl;
                exit(1);
            }
        } else {
            if (strcmp(pixelType,      "UINT8" ) == 0) enblendMain<UInt8 >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "FLOAT" ) == 0) enblendMain<float >(imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (strcmp(pixelType, "INT8"  ) == 0) enblendMain<Int8  >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "UINT16") == 0) enblendMain<UInt16>(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "INT16" ) == 0) enblendMain<Int16 >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "UINT32") == 0) enblendMain<UInt32>(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "INT32" ) == 0) enblendMain<Int32 >(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "UINT64") == 0) enblendMain<UInt64>(imageInfoList, outputImageInfo, inputUnion);
            //else if (strcmp(pixelType, "INT64" ) == 0) enblendMain<Int64 >(imageInfoList, outputImageInfo, inputUnion);
            else if (strcmp(pixelType, "DOUBLE") == 0) enblendMain<double>(imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                cerr << "enblend: images with pixel type \""
                     << pixelType
                     << "\" are not supported."
                     << endl;
                exit(1);
            }
        }

        // delete entries in imageInfoList, in case
        // enblend loop returned early.
        imageInfoIterator = imageInfoList.begin();
        while (imageInfoIterator != imageInfoList.end()) {
            delete *imageInfoIterator++;
        }

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

#ifdef HAVE_LIBGLEW
    if (UseGPU) {
        wrapupGPU();
    }
#endif

    // Success.
    return 0;
}
