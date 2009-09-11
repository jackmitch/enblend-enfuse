/*
 * Copyright (C) 2004-2009 Andrew Mihal
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
#ifndef HAVE_CONFIG_H
#include <win32helpers\win32config.h>
#endif
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
#include <set>
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

#include "global.h"

// Globals
const std::string command("enfuse");

// Random number generator for dithering
boost::mt19937 Twister;

// Global values from command line parameters.
std::string OutputFileName(DEFAULT_OUTPUT_FILENAME);
int Verbose = DEFAULT_VERBOSITY;
unsigned int ExactLevels = 0;
bool OneAtATime = true;
boundary_t WrapAround = OpenBoundaries;
bool GimpAssociatedAlphaHack = false;
bool UseCIECAM = false;
bool OutputSizeGiven = false;
int OutputWidthCmdLine = 0;
int OutputHeightCmdLine = 0;
int OutputOffsetXCmdLine = 0;
int OutputOffsetYCmdLine = 0;
std::string OutputCompression;
std::string OutputPixelType;
double WExposure = 1.0;
double WContrast = 0.0;
double WSaturation = 0.2;
double WEntropy = 0.0;
double WMu = 0.5;
double WSigma = 0.2;
bool WSaturationIsDefault = true;
int ContrastWindowSize = 5;
std::string GrayscaleProjector;
struct EdgeFilterConfiguration {double edgeScale, lceScale, lceFactor;} FilterConfig = {0.0, 0.0, 0.0};
struct AlternativePercentage MinCurvature = {0.0, false};
int EntropyWindowSize = 3;
struct AlternativePercentage EntropyLowerCutoff = {0.0, true};
struct AlternativePercentage EntropyUpperCutoff = {100.0, true};
bool UseHardMask = false;
bool SaveMasks = false;
std::string SoftMaskTemplate("softmask-%n.tif");
std::string HardMaskTemplate("hardmask-%n.tif");
TiffResolution ImageResolution;

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
#include "enfuse.h"

#include "vigra/imageinfo.hxx"
#include "vigra/impex.hxx"
#include "vigra/sized_int.hxx"

#include <tiffio.h>

#ifdef DMALLOC
#include "dmalloc.h"            // must be last #include
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::hex;
using std::list;

using vigra::CachedFileImageDirector;
using vigra::Diff2D;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::Rect2D;
using vigra::StdException;

using enblend::enfuseMain;

#ifdef _WIN32
#define strdup _strdup
#endif

// Initialize data structures for precomputed entropy and logarithm.
template <typename InputPixelType, typename ResultPixelType>
size_t enblend::Histogram<InputPixelType, ResultPixelType>::precomputedSize = 0;
template <typename InputPixelType, typename ResultPixelType>
double* enblend::Histogram<InputPixelType, ResultPixelType>::precomputedLog = NULL;
template <typename InputPixelType, typename ResultPixelType>
double* enblend::Histogram<InputPixelType, ResultPixelType>::precomputedEntropy = NULL;


/** Print information on the current version and some configuration
 * details. */
void printVersionAndExit() {
    cout << "enfuse " << VERSION << "\n\n";

    if (Verbose >= VERBOSE_VERSION_REPORTING) {
        cout <<
            "Extra feature: dmalloc support: " <<
#ifdef DMALLOC
            "yes" <<
#else
            "no" <<
#endif
            "\n";

        cout <<
            "Extra feature: image cache: " <<
#ifdef ENBLEND_CACHE_IMAGES
            "yes" <<
#else
            "no" <<
#endif
            "\n";

#ifdef OPENMP
        const bool openmp_nested = omp_get_nested();
        omp_set_nested(true);
        const bool openmp_dynamic = omp_get_dynamic();
        omp_set_dynamic(true);
        cout <<
            "Extra feature: OpenMP: " <<
            "yes - version " << OPENMP_YEAR << '-' << OPENMP_MONTH << "\n" <<
            "                           - " << (omp_get_nested() ? "" : "no ") <<
            "support for nested parallelism\n" <<
            "                           - " << (omp_get_dynamic() ? "" : "no ") <<
            "support for dynamic adjustment of the number of threads\n" <<
            "                           - using " <<
            omp_get_num_procs() << " processor" << (omp_get_num_procs() >= 2 ? "s" : "") << " and up to " <<
            omp_get_max_threads() << " thread" << (omp_get_max_threads() >= 2 ? "s" : "") << "\n\n";
        omp_set_dynamic(openmp_dynamic);
        omp_set_nested(openmp_nested);
#else
        cout << "Extra feature: OpenMP: " << "no\n\n";
#endif

        cout <<
            "Supported image formats: " << vigra::impexListFormats() << "\n" <<
            "Supported file extensions: " << vigra::impexListExtensions() << "\n\n";
    }

    cout <<
        "Copyright (C) 2004-2009 Andrew Mihal.\n" <<
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
        "Usage: enfuse [options] [--output=IMAGE] INPUT...\n" <<
        "Fuse INPUT images into a single IMAGE.\n" <<
        "\n" <<
        "Common options:\n" <<
        "  -V, --version          output version information and exit\n" <<
        "  -h, --help             print this help message and exit\n" <<
        "  -l LEVELS              number of blending LEVELS to use (1 to 29)\n" <<
        "  -o, --output=FILE      write output to FILE; default: \"" << OutputFileName << "\"\n" <<
        "  -v, --verbose[=LEVEL]  verbosely report progress; repeat to\n" <<
        "                         increase verbosity or directly set to LEVEL\n" <<
        "  -w, --wrap[=MODE]      wrap around image boundary, where MODE is\n" <<
        "                         NONE, HORIZONTAL, VERTICAL, or BOTH; default: " <<
        enblend::stringOfWraparound(WrapAround) << ";\n" <<
        "  --compression=COMP     set compression of output image to COMP,\n" <<
        "                         where COMP is:\n" <<
        "                           NONE, PACKBITS, LZW, DEFLATE for TIFF files and\n" <<
        "                           0 to 100 for JPEG files\n" <<
        "\n" <<
        "Extended options:\n" <<
        "  -b BLOCKSIZE           image cache BLOCKSIZE in kilobytes; default: " <<
        (CachedFileImageDirector::v().getBlockSize() / 1024LL) << "KB\n" <<
        "  -c                     use CIECAM02 to blend colors\n" <<
        "  -d, --depth=DEPTH      set the number of bits per channel of the output\n" <<
        "                         image, where DEPTH is 8, 16, 32, r32, or r64\n" <<
        "  -g                     associated-alpha hack for Gimp (before version 2)\n" <<
        "                         and Cinepaint\n" <<
        "  -f WIDTHxHEIGHT[+xXOFFSET+yYOFFSET]\n" <<
        "                         manually set the size and position of the output\n" <<
        "                         image; useful for cropped and shifted input\n" <<
        "                         TIFF images, such as those produced by Nona\n" <<
        "  -m CACHESIZE           set image CACHESIZE in megabytes; default: " <<
        (CachedFileImageDirector::v().getAllocation() / 1048576LL) << "MB\n" <<
        "\n" <<
        "Fusion options:\n" <<
        "  --wExposure=WEIGHT     weight given to well-exposed pixels\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WExposure << "\n" <<
        "  --wSaturation=WEIGHT   weight given to highly-saturated pixels\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WSaturation << "\n" <<
        "  --wContrast=WEIGHT     weight given to high-contrast pixels\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WContrast << "\n" <<
        "  --wEntropy=WEIGHT      weight given to high entropy regions\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WEntropy << "\n" <<
        "  --wMu=MEAN             center also known as MEAN of Gaussian weighting\n" <<
        "                         function (0 <= MEAN <= 1); default: " << WMu << "\n" <<
        "  --wSigma=SIGMA         standard deviation of Gaussian weighting\n" <<
        "                         function (SIGMA > 0); default: " << WSigma << "\n" <<
        "  --SoftMask             average over all masks; this is the default\n" <<
        "  --HardMask             force hard blend masks and no averaging on finest\n" <<
        "                         scale; this is especially useful for focus\n" <<
        "                         stacks with thin and high contrast features,\n" <<
        "                         but leads to increased noise\n" <<
        "\n" <<
        "Expert options:\n" <<
        "  --ContrastWindowSize=SIZE\n" <<
        "                         set window SIZE for local-contrast analysis\n" <<
        "                         (SIZE >= 3); default: " << ContrastWindowSize  << "\n" <<
        "  --GrayProjector=OPERATOR\n" <<
        "                         apply gray-scale projection OPERATOR, where\n" <<
        "                         OPERATOR is one of \"average\", \"l-star\",\n" <<
        "                         \"lightness\", \"value\", \"luminance\", or\n" <<
        "                         \"channel-mixer:RED-WEIGHT:GREEN-WEIGHT:BLUE-WEIGHT\";\n" <<
        "                         default: \"" <<
        enblend::MultiGrayscaleAccessor<UInt8, NumericTraits<UInt8>::Promote>::defaultGrayscaleAccessorName() << "\"\n" <<
        "  --EdgeScale=EDGESCALE[:LCESCALE[:LCEFACTOR]]\n" <<
        "                         set scale on which to look for edges; positive\n" <<
        "                         LCESCALE switches on local contrast enhancement\n" <<
        "                         by LCEFACTOR (EDGESCALE, LCESCALE, LCEFACTOR >= 0);\n" <<
        "                         append \"%\" to LCESCALE for values relative to\n" <<
        "                         EDGESCALE; append \"%\" to LCEFACTOR for relative\n" <<
        "                         value; defaults: " <<
        FilterConfig.edgeScale << ":" << FilterConfig.lceScale << ":" << FilterConfig.lceFactor << "\n" <<
        "  --MinCurvature=CURVATURE\n" <<
        "                         minimum CURVATURE for an edge to qualify; append\n" <<
        "                         \"%\" for relative values; default: " << MinCurvature.str() << "\n" <<
        "  --EntropyWindowSize=SIZE\n" <<
        "                         set window SIZE for local entropy analysis\n" <<
        "                         (SIZE >= 3); default: " << EntropyWindowSize  << "\n" <<
        "  --EntropyCutoff=LOWERCUTOFF[:UPPERCUTOFF]\n" <<
        "                         LOWERCUTOFF is the value below of which pixels are\n" <<
        "                         treated as black and UPPERCUTOFF is the value above\n" <<
        "                         of which pixels are treated as white in the entropy\n" <<
        "                         weighting; append \"%\" signs for relative values;\n" <<
        "                         default: " <<
        EntropyLowerCutoff.str() << ":" << EntropyUpperCutoff.str() << "\n" <<
        "  --SaveMasks[=SOFT-TEMPLATE[:HARD-TEMPLATE]]\n" <<
        "                         save weight masks in SOFT-TEMPLATE and HARD-TEMPLATE;\n" <<
        "                         conversion chars: %i: mask index, %n: mask number,\n" <<
        "                         %p: full path, %d: dirname, %b: basename,\n" <<
        "                         %f: filename, %e: extension; lowercase characters\n" <<
        "                         refer to input images uppercase to output image;\n" <<
        "                         default: \"" << SoftMaskTemplate << "\":\"" << HardMaskTemplate << "\"\n" <<
        "  --debug                output mask images for debugging (deprecated);\n" <<
        "                         use \"--SaveMasks\" instead\n" <<
        "\n" <<
        "Report bugs at <" PACKAGE_BUGREPORT ">." <<
        endl;

    exit(error ? 1 : 0);
}

/** Make sure all cached file images get destroyed,
 *  and hence the temporary files deleted,
 *  if we are killed.
 */
void sigint_handler(int sig) {
    cerr << endl << "Interrupted." << endl;
    // FIXME what if this occurs in a CFI atomic section?
    // This is no longer necessary, temp files are unlinked during creation.
    //CachedFileImageDirector::v().~CachedFileImageDirector();
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


enum AllPossibleOptions {
    VersionOption, HelpOption, LevelsOption, OutputOption, VerboseOption,
    WrapAroundOption /* -w */, CompressionOption, LZWCompressionOption,
    BlockSizeOption, CIECAM02Option /* -c */,
    DepthOption, AssociatedAlphaOption /* -g */,
    SizeAndPositionOption /* -f */, CacheSizeOption,
    ExposureWeightOption, SaturationWeightOption,
    ContrastWeightOption, EntropyWeightOption,
    ExposureMuOption /* -wMu */, ExposureSigmaOption /* -wSigma */,
    SoftMaskOption, HardMaskOption,
    ContrastWindowSizeOption, GrayProjectorOption, EdgeScaleOption,
    MinCurvatureOption, EntropyWindowSizeOption, EntropyCutoffOption,
    DebugOption, SaveMasksOption
};

typedef std::set<enum AllPossibleOptions> OptionSetType;

bool contains(const OptionSetType& optionSet,
              enum AllPossibleOptions anOption)
{
    return optionSet.count(anOption) != 0;
}


/** Warn if options given at the command line have no effect. */
void warn_of_ineffective_options(const OptionSetType& optionSet)
{
    if (WExposure == 0.0 && contains(optionSet, ExposureWeightOption)) {
        if (contains(optionSet, ExposureMuOption)) {
            cerr << command <<
                ": warning: option \"--wMu\" has no effect as exposure weight\n" <<
                command <<
                ": warning:     is zero" <<
                endl;
        }
        if (contains(optionSet, ExposureSigmaOption)) {
            cerr << command <<
                ": warning: option \"--wSigma\" has no effect as exposure weight\n" <<
                command <<
                ": warning:     is zero" <<
                endl;
        }
    }

    if (WContrast == 0.0 && contains(optionSet, ContrastWindowSizeOption)) {
        cerr << command <<
            ": warning: option \"--ContrastWindowSize\" has no effect as contrast\n" <<
            command <<
            ": warning:     weight is zero" <<
            endl;
    }

    if (WExposure == 0.0 && WContrast == 0.0 && contains(optionSet, GrayProjectorOption)) {
        cerr << command <<
            ": warning: option \"--GrayProjector\" has no effect as exposure\n" <<
            command <<
            ": warning:     and contrast weight both are zero" <<
            endl;
    }

    if (WContrast == 0.0) {
        if (contains(optionSet, EdgeScaleOption)) {
            cerr << command <<
                ": warning: option \"--EdgeScale\" has no effect as contrast\n" <<
                command <<
                ": warning:     weight is zero" <<
                endl;
        }
        if (contains(optionSet, MinCurvatureOption)) {
            cerr << command <<
                ": warning: option \"--MinCurvature\" has no effect as contrast\n" <<
                command <<
                ": warning:     weight is zero" <<
                endl;
        }
    } else {
        if (FilterConfig.edgeScale > 0.0 &&
            contains(optionSet, ContrastWindowSizeOption) && MinCurvature.value <= 0.0) {
            cerr << command <<
                ": warning: option \"--ContrastWindowSize\" has no effect as\n" <<
                command <<
                ": warning:     EDGESCALE in \"--EdgeScale\" is positive and" <<
                command <<
                ": warning:     \"--MinCurvature\" is non-positive" <<
                endl;
        }
    }

    if (WEntropy == 0.0) {
        if (contains(optionSet, EntropyWindowSizeOption)) {
            cerr << command <<
                ": warning: option \"--EntropyWindowSize\" has no effect as\n" <<
                command <<
                ": warning:     entropy weight is zero" <<
                endl;
        }
        if (contains(optionSet, EntropyCutoffOption)) {
            cerr << command <<
                ": warning: option \"--EntropyCutoff\" has no effect as entropy\n" <<
                command <<
                ": warning:     weight is zero" <<
                endl;
        }
    }

    if (contains(optionSet, CompressionOption) &&
        !(enblend::getFileType(OutputFileName) == "TIFF" ||
          enblend::getFileType(OutputFileName) == "JPEG")) {
        cerr << command <<
            ": warning: compression is not supported with output\n" <<
            command <<
            ": warning:     file type \"" <<
            enblend::getFileType(OutputFileName) << "\"" <<
            endl;
    }

    if (contains(optionSet, AssociatedAlphaOption) &&
        enblend::getFileType(OutputFileName) != "TIFF") {
        cerr << command <<
            ": warning: option \"-g\" has no effect with output\n" <<
            command <<
            ": warning:     file type \"" <<
            enblend::getFileType(OutputFileName) << "\"" <<
            endl;
    }

#ifndef ENBLEND_CACHE_IMAGES
    if (contains(optionSet, CacheSizeOption)) {
        cerr << command <<
            ": warning: option \"-m\" has no effect in this version of " << command << ",\n" <<
            command <<
            ": warning:     because it was compiled without image cache" <<
            endl;
    }

    if (contains(optionSet, BlockSizeOption)) {
        cerr << command <<
            ": warning: option \"-b\" has no effect in this version of " << command << ",\n" <<
            command <<
            ": warning:     because it was compiled without image cache" <<
            endl;
    }
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
        CompressionId,           //  0
        WeightExposureId,        //  1
        WeightContrastId,        //  2
        WeightSaturationId,      //  3
        WeightMuId,              //  4
        WeightSigmaId,           //  5
        MinCurvatureId,          //  6
        EdgeScaleId,             //  7
        ContrastWindowSizeId,    //  8
        HardMaskId,              //  9
        GrayProjectorId,         // 10
        DebugId,                 // 11
        WeightEntropyId,         // 12
        EntropyWindowSizeId,     // 13
        EntropyCutoffId,         // 14
        SoftMaskId,              // 15
        VerboseId,               // 16
        HelpId,                  // 17
        VersionId,               // 18
        DepthId,                 // 19
        OutputId,                // 20
        SaveMasksId,             // 21
        WrapAroundId             // 22
    };

    // NOTE: See note attached to "enum OptionId" above.
    static struct option long_options[] = {
        {"compression", required_argument, 0, StringArgument},           //  0
        {"wExposure", required_argument, 0, FloatArgument},              //  1
        {"wContrast", required_argument, 0, FloatArgument},              //  2
        {"wSaturation", required_argument, 0, FloatArgument},            //  3
        {"wMu", required_argument, 0, FloatArgument},                    //  4
        {"wSigma", required_argument, 0, FloatArgument},                 //  5
        {"MinCurvature", required_argument, 0, StringArgument},          //  6
        {"EdgeScale", required_argument, 0, StringArgument},             //  7
        {"ContrastWindowSize", required_argument, 0, IntegerArgument},   //  8
        {"HardMask", no_argument, 0, NoArgument},                        //  9
        {"GrayProjector", required_argument, 0, StringArgument},         // 10
        {"debug", no_argument, 0, NoArgument},                           // 11
        {"wEntropy", required_argument, 0, FloatArgument},               // 12
        {"EntropyWindowSize", required_argument, 0, IntegerArgument},    // 13
        {"EntropyCutoff", required_argument, 0, StringArgument},         // 14
        {"SoftMask", no_argument, 0, NoArgument},                        // 15
        {"verbose", optional_argument, 0, IntegerArgument},              // 16
        {"help", no_argument, 0, NoArgument},                            // 17
        {"version", no_argument, 0, NoArgument},                         // 18
        {"depth", required_argument, 0, StringArgument},                 // 19
        {"output", required_argument, 0, StringArgument},                // 20
        {"SaveMasks", optional_argument, 0, StringArgument},             // 21
        {"wrap", optional_argument, 0, StringArgument},                  // 22
        {0, 0, 0, 0}
    };

    bool justPrintVersion = false;
    bool justPrintUsage = false;
    OptionSetType optionSet;

    // Parse command line.
    int option_index = 0;
    int c;
    opterr = 0;       // we have our own "unrecognized option" message
    while ((c = getopt_long(argc, argv, "Vb:cd:f:ghl:m:o:v::w::z",
                            long_options, &option_index)) != -1) {
        switch (c) {
        case NoArgument: {
            if (long_options[option_index].flag != 0) {
                break;
            }
            switch (option_index) {
            case HardMaskId:
                UseHardMask = true;
                optionSet.insert(HardMaskOption);
                break;
            case SoftMaskId:
                UseHardMask = false;
                optionSet.insert(SoftMaskOption);
                break;
            case DebugId:
                cerr << command
                     << ": warning: option \"--debug\" is deprecated; use \"--SaveMasks\" instead"
                     << endl;
                SaveMasks = true;
                optionSet.insert(DebugOption);
                break;
            case HelpId:
                justPrintUsage = true;
                optionSet.insert(HelpOption);
                break;
            case VersionId:
                justPrintVersion = true;
                optionSet.insert(VersionOption);
                break;
            default:
                cerr << command
                     << ": internal error: unhandled \"NoArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case NoArgument"

        case StringArgument: {
            if (long_options[option_index].flag != 0) {
                break;
            }
            switch (option_index) {
            case WrapAroundId:
                if (optarg != NULL && *optarg != 0) {
                    WrapAround = enblend::wraparoundOfString(optarg);
                    if (WrapAround == UnknownWrapAround) {
                        cerr << command
                             << ": unrecognized wrap-around mode \"" << optarg << "\"\n" << endl;
                        exit(1);
                    }
                } else {
                    WrapAround = HorizontalStrip;
                }
                optionSet.insert(WrapAroundOption);
                break;
            case MinCurvatureId: {
                char *tail;
                errno = 0;
                MinCurvature.value = strtod(optarg, &tail);
                if (errno == 0) {
                    if (*tail == 0) {
                        MinCurvature.isPercentage = false;
                    } else if (strcmp(tail, "%") == 0) {
                        MinCurvature.isPercentage = true;
                    } else {
                        cerr << command << ": unrecognized minimum gradient \""
                             << optarg << "\" specification." << endl;
                        exit(1);
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << optarg << "\" for minimum gradient: "
                         << enblend::errorMessage(errno) << endl;
                    exit(1);
                }
                optionSet.insert(MinCurvatureOption);
                break;
            }

            case EdgeScaleId: {
                char* s = new char[strlen(optarg) + 1];
                strcpy(s, optarg);
                char* save_ptr = NULL;
                char* token = enblend::strtoken_r(s, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                char* tail;

                if (token == NULL || *token == 0) {
                    cerr << command << ": no scale given to --EdgeScale.  "
                         << "EdgeScale is required." << endl;
                    exit(1);
                }
                errno = 0;
                FilterConfig.edgeScale = strtod(token, &tail);
                if (errno == 0) {
                    if (*tail != 0) {
                        cerr << command << ": could not decode \"" << tail
                             << "\" in edge scale specification \""
                             << token << "\" for EdgeScale." << endl;
                        exit(1);
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << token << "\" for EdgeScale: "
                         << enblend::errorMessage(errno) << endl;
                    exit(1);
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    FilterConfig.lceScale = strtod(token, &tail);
                    if (errno == 0) {
                        if (strcmp(tail, "%") == 0) {
                            FilterConfig.lceScale *= FilterConfig.edgeScale / 100.0;
                        } else if (*tail != 0) {
                            cerr << command << ": could not decode \"" << tail
                                 << "\" in specification \"" << token
                                 << "\" for LCE-scale." << endl;
                            exit(1);
                        }
                    } else {
                        cerr << command << ": illegal numeric format \""
                             << token << "\" for LCE-Scale: "
                             << enblend::errorMessage(errno) << endl;
                        exit(1);
                    }
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    FilterConfig.lceFactor = strtod(token, &tail);
                    if (errno == 0) {
                        if (strcmp(tail, "%") == 0) {
                            FilterConfig.lceFactor /= 100.0;
                        } else if (*tail != 0) {
                            cerr << command << ": could not decode \"" << tail
                                 << "\" in specification \"" << token
                                 << "\" for LCE-factor." << endl;
                            exit(1);
                        }
                    } else {
                        cerr << command << ": illegal numeric format \""
                             << token << "\" for LCE-factor: "
                             << enblend::errorMessage(errno) << endl;
                        exit(1);
                    }
                }

                if (save_ptr != NULL && *save_ptr != 0) {
                    cerr << command << ": warning: ignoring trailing garbage \""
                         << save_ptr << "\" in argument to --EdgeScale" << endl;
                }

                delete [] s;
                optionSet.insert(EdgeScaleOption);
                break;
            }

            case EntropyCutoffId: {
                char* s = new char[strlen(optarg) + 1];
                strcpy(s, optarg);
                char* save_ptr = NULL;
                char* token = enblend::strtoken_r(s, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                char* tail;

                if (token == NULL || *token == 0) {
                    cerr << command << ": no scale given to --EntropyCutoff.  "
                         << "LowerCutOff is required." << endl;
                    exit(1);
                }
                errno = 0;
                EntropyLowerCutoff.value = strtod(token, &tail);
                if (errno == 0) {
                    if (*tail == 0) {
                        EntropyLowerCutoff.isPercentage = false;
                    } else if (strcmp(tail, "%") == 0) {
                        EntropyLowerCutoff.isPercentage = true;
                    } else {
                        cerr << command << ": unrecognized entropy's lower cutoff \""
                             << tail << "\" in \"" << token << "\"" << endl;
                        exit(1);
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << token << "\" of entropy's lower cutoff: "
                         << enblend::errorMessage(errno) << endl;
                    exit(1);
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    EntropyUpperCutoff.value = strtod(token, &tail);
                    if (errno == 0) {
                        if (*tail == 0) {
                            EntropyUpperCutoff.isPercentage = false;
                        } else if (strcmp(tail, "%") == 0) {
                            EntropyUpperCutoff.isPercentage = true;
                        } else {
                            cerr << command << ": unrecognized entropy's upper cutoff \""
                                 << tail << "\" in \"" << token << "\"" << endl;
                            exit(1);
                        }
                    } else {
                        cerr << command << ": illegal numeric format \""
                             << token << "\" of entropy's upper cutoff: "
                             << enblend::errorMessage(errno) << endl;
                        exit(1);
                    }
                }

                if (save_ptr != NULL && *save_ptr != 0) {
                    cerr << command << ": warning: ignoring trailing garbage \""
                         << save_ptr << "\" in argument to --EntropyCutoff" << endl;
                }

                delete [] s;
                optionSet.insert(EntropyCutoffOption);
                break;
            }

            case CompressionId:
                OutputCompression = optarg;
                optionSet.insert(CompressionOption);
                break;

            case GrayProjectorId:
                GrayscaleProjector = optarg;
                optionSet.insert(GrayProjectorOption);
                break;

            case DepthId:
                OutputPixelType = enblend::outputPixelTypeOfString(optarg);
                optionSet.insert(DepthOption);
                break;

            case OutputId:
                if (contains(optionSet, OutputOption)) {
                    cerr << command
                         << ": warning: more than one output file specified"
                         << endl;
                }
                OutputFileName = optarg;
                optionSet.insert(OutputOption);
                break;

            case SaveMasksId:
                SaveMasks = true;
                optionSet.insert(SaveMasksOption);
                if (optarg != NULL && *optarg != 0) {
                    char* s = new char[strlen(optarg) + 1];
                    strcpy(s, optarg);
                    char* save_ptr = NULL;
                    char* token = enblend::strtoken_r(s, PATH_OPTION_DELIMITERS, &save_ptr);
                    SoftMaskTemplate = token;
                    token = enblend::strtoken_r(NULL, PATH_OPTION_DELIMITERS, &save_ptr);
                    if (token != NULL && *token != 0) {
                        HardMaskTemplate = token;
                    }
                    token = enblend::strtoken_r(NULL, PATH_OPTION_DELIMITERS, &save_ptr);
                    if (token != NULL && *token != 0) {
                        cerr << command
                             << ": warning: ignoring trailing garbage in --SaveMasks" << endl;
                    }
                    delete [] s;
                }
                break;

            default:
                cerr << command << ": internal error: unhandled \"StringArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case StringArgument"

        case FloatArgument: {
            if (long_options[option_index].flag != 0) {
                break;
            }
            double *optionDouble = NULL;
            switch (option_index) {
            case WeightExposureId:
                optionDouble = &WExposure;
                optionSet.insert(ExposureWeightOption);
                break;
            case WeightContrastId:
                optionDouble = &WContrast;
                optionSet.insert(ContrastWeightOption);
                break;
            case WeightSaturationId:
                optionDouble = &WSaturation;
                WSaturationIsDefault = false;
                optionSet.insert(SaturationWeightOption);
                break;
            case WeightMuId:
                optionDouble = &WMu;
                optionSet.insert(ExposureMuOption);
                break;
            case WeightSigmaId:
                optionDouble = &WSigma;
                optionSet.insert(ExposureSigmaOption);
                break;
            case WeightEntropyId:
                optionDouble = &WEntropy;
                optionSet.insert(EntropyWeightOption);
                break;
            default:
                cerr << command
                     << ": internal error: unhandled \"FloatArgument\" option"
                     << endl;
                exit(1);
            }

            char *lastChar = NULL;
            double value = strtod(optarg, &lastChar);
            if ((lastChar == optarg || value < 0.0 || value > 1.0) &&
                (option_index == WeightExposureId || option_index == WeightContrastId ||
                 option_index == WeightSaturationId || option_index == WeightEntropyId)) {
                cerr << command << ": " << long_options[option_index].name
                     << " must be in the range [0.0, 1.0]." << endl;
                exit(1);
            }

            *optionDouble = value;
            break;
        } // end of "case FloatArgument"

        case IntegerArgument: {
            if (long_options[option_index].flag != 0) {
                break;
            }

            switch (option_index) {
            case VerboseId:
                if (optarg != NULL && *optarg != 0) {
                    Verbose =
                        enblend::numberOfString(optarg,
                                                _1 >= 0,
                                                "verbosity level less than 0; will use 0",
                                                0);
                } else {
                    Verbose++;
                }
                optionSet.insert(VerboseOption);
                break;

            case ContrastWindowSizeId:
                ContrastWindowSize =
                    enblend::numberOfString(optarg,
                                            _1 >= 3,
                                            "contrast window size to small; will use size = 3",
                                            3);
                if (ContrastWindowSize % 2 != 1) {
                    cerr << command << ": warning: contrast window size \""
                         << ContrastWindowSize << "\" is even; increasing size to next odd number"
                         << endl;
                    ContrastWindowSize++;
                }
                optionSet.insert(ContrastWindowSizeOption);
                break;

            case EntropyWindowSizeId:
                EntropyWindowSize =
                    enblend::numberOfString(optarg,
                                            _1 >= 3,
                                            "entropy window size to small; will use size = 3",
                                            3);
                if (EntropyWindowSize % 2 != 1) {
                    cerr << command << ": warning: entropy window size \""
                         << EntropyWindowSize << "\" is even; increasing size to next odd number"
                         << endl;
                    EntropyWindowSize++;
                }
                optionSet.insert(EntropyWindowSizeOption);
                break;

            default:
                cerr << command << ": internal error: unhandled \"IntegerArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case IntegerArgument"

        case 'V':
            justPrintVersion = true;
            optionSet.insert(VersionOption);
            break;
        case 'b': {
            const int cache_block_size =
                enblend::numberOfString(optarg,
                                        _1 >= 1,
                                        "cache block size must be 1 KB or more; will use 1 KB",
                                        1);
            CachedFileImageDirector::v().setBlockSize(static_cast<long long>(cache_block_size) << 10);
            optionSet.insert(BlockSizeOption);
            break;
        }
        case 'c':
            UseCIECAM = true;
            optionSet.insert(CIECAM02Option);
            break;
        case 'd':
            OutputPixelType = enblend::outputPixelTypeOfString(optarg);
            optionSet.insert(DepthOption);
            break;
        case 'f': {
            OutputSizeGiven = true;
            const int nP = sscanf(optarg,
                                  "%dx%d+%d+%d",
                                  &OutputWidthCmdLine, &OutputHeightCmdLine,
                                  &OutputOffsetXCmdLine, &OutputOffsetYCmdLine);
            if (nP == 4) {
                ; // ok: full geometry string
            } else if (nP == 2) {
                OutputOffsetXCmdLine = 0;
                OutputOffsetYCmdLine = 0;
            } else {
                cerr << command << ": option \"-f\" requires a parameter\n"
                     << "Try \"enfuse --help\" for more information." << endl;
                exit(1);
            }
            optionSet.insert(SizeAndPositionOption);
            break;
        }
        case 'g':
            GimpAssociatedAlphaHack = true;
            optionSet.insert(AssociatedAlphaOption);
            break;
        case 'h':
            justPrintUsage = true;
            optionSet.insert(HelpOption);
            break;
        case 'l': {
            // We take care of "too many levels" in "bounds.h".
            ExactLevels =
                enblend::numberOfString(optarg,
                                        _1 >= 1U,
                                        "too few levels; will use one level",
                                        1U);
            optionSet.insert(LevelsOption);
            break;
        }
        case 'm': {
            const int cache_size =
                enblend::numberOfString(optarg,
                                        _1 >= 1,
                                        "cache memory limit less than 1 MB; will use 1 MB",
                                        1);
            CachedFileImageDirector::v().setAllocation(static_cast<long long>(cache_size) << 20);
            optionSet.insert(CacheSizeOption);
            break;
        }
        case 'o':
            if (contains(optionSet, OutputOption)) {
                cerr << command
                     << ": warning: more than one output file specified"
                     << endl;
            }
            OutputFileName = optarg;
            optionSet.insert(OutputOption);
            break;
        case 'v':
            if (optarg != NULL && *optarg != 0) {
                Verbose =
                    enblend::numberOfString(optarg,
                                            _1 >= 0,
                                            "verbosity level less than 0; will use 0",
                                            0);
            } else {
                Verbose++;
            }
            optionSet.insert(VerboseOption);
            break;
        case 'w':
            if (optarg != NULL && *optarg != 0) {
                WrapAround = enblend::wraparoundOfString(optarg);
                if (WrapAround == UnknownWrapAround) {
                    cerr << command
                         << ": unrecognized wrap-around mode \"" << optarg << "\"\n" << endl;
                    exit(1);
                }
            } else {
                WrapAround = HorizontalStrip;
            }
            optionSet.insert(WrapAroundOption);
            break;
        case 'z':
            cerr << command
                 << ": info: flag \"-z\" is deprecated; use \"--compression=LZW\" instead"
                 << endl;
            OutputCompression = "LZW";
            optionSet.insert(LZWCompressionOption);
            break;
        case '?':
            switch (optopt) {
                case 0: // unknown long option
                    cerr << command
                         << ": unknown option \""
                         << argv[optind - 1]
                         << "\"\n";
                    break;
                case 'b':           // FALLTHROUGH
                case 'd':           // FALLTHROUGH
                case 'f':           // FALLTHROUGH
                case 'l':           // FALLTHROUGH
                case 'm':           // FALLTHROUGH
                case 'o':
                    cerr << command
                         << ": option \"-"
                         << static_cast<char>(optopt)
                         << "\" requires an argument"
                         << endl;
                    break;
                default:
                    cerr << command << ": unknown option ";
                    if (isprint(optopt)) {
                        cerr << "\"-" << static_cast<char>(optopt) << "\"";
                    } else {
                        cerr << "character 0x" << hex << optopt;
                    }
                    cerr << endl;
            }
            cerr << "Try \"enfuse --help\" for more information." << endl;
            exit(1);

        default:
            cerr << command
                 << ": internal error: unhandled command line option"
                 << endl;
            exit(1);
        }
    }

    if (justPrintUsage) {
        printUsageAndExit(false);
        // never reached
    }

    if (justPrintVersion) {
        printVersionAndExit();
        // never reached
    }

    warn_of_ineffective_options(optionSet);

    return optind;
}


int main(int argc, char** argv)
{
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
    try {
        optind = process_options(argc, argv);
    } catch (StdException& e) {
        cerr << command << ": error while processing command line options\n"
             << command << ":     " << e.what()
             << endl;
        exit(1);
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
        cerr << command << ": no input files specified.\n";
        exit(1);
    }

    // List of info structures for each input image.
    list<ImageImportInfo*> imageInfoList;
    list<ImageImportInfo*>::iterator imageInfoIterator;

    bool isColor = false;
    std::string pixelType;
    TiffResolution resolution;
    ImageImportInfo::ICCProfile iccProfile;
    Rect2D inputUnion;

    // Check that all input images have the same parameters.
    inputFileNameIterator = inputFileNameList.begin();
    unsigned layer = 0;
    unsigned layers = 0;
    while (inputFileNameIterator != inputFileNameList.end()) {
        ImageImportInfo* inputInfo = NULL;
        try {
            std::string filename(*inputFileNameIterator);
            ImageImportInfo info(filename.c_str());
            if (layers == 0) { // OPTIMIZATION: call only once per file
                layers = info.numLayers();
            }
            if (layers >= 2) {
                filename = vigra::join_filename_layer(*inputFileNameIterator, layer);
            }
            ++layer;
            inputInfo = new ImageImportInfo(filename.c_str());
        } catch (StdException& e) {
            cerr << '\n' << command
                 << ": error opening input file \"" << *inputFileNameIterator << "\":\n"
                 << e.what() << endl;
            exit(1);
        }

        // Save this image info in the list.
        imageInfoList.push_back(inputInfo);

        if (Verbose >= VERBOSE_INPUT_IMAGE_INFO_MESSAGES) {
            cerr << command
                 << ": info: input image \""
                 << *inputFileNameIterator
                 << "\" "
                 << layer << '/' << layers << ' ';

            if (inputInfo->isColor()) {
                cerr << "RGB ";
            }

            if (!inputInfo->getICCProfile().empty()) {
                cerr << "ICC ";
            }

            cerr << inputInfo->getPixelType()
                 << " position="
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
            cerr << command
                 << ": warning: input image \""
                 << *inputFileNameIterator
                 << "\" does not have an alpha channel;\n"
                 << command
                 << ": warning: assuming all pixels should "
                 << "contribute to the final image"
                 << endl;
        }

        // Get input image's position and size.
        Rect2D imageROI(Point2D(inputInfo->getPosition()),
                        Size2D(inputInfo->width(), inputInfo->height()));

        if (inputFileNameIterator == inputFileNameList.begin()) {
            // First input image
            inputUnion = imageROI;
            isColor = inputInfo->isColor();
            pixelType = inputInfo->getPixelType();
            resolution = TiffResolution(inputInfo->getXResolution(),
                                        inputInfo->getYResolution());
            iccProfile = inputInfo->getICCProfile();
            if (!iccProfile.empty()) {
                InputProfile = cmsOpenProfileFromMem(iccProfile.data(), iccProfile.size());
                if (InputProfile == NULL) {
                    cerr << endl
                         << command << ": error parsing ICC profile data from file \""
                         << *inputFileNameIterator
                         << "\"" << endl;
                    exit(1);
                }
            }
        } else {
            // Second and later images
            inputUnion |= imageROI;

            if (isColor != inputInfo->isColor()) {
                cerr << command << ": input image \""
                     << *inputFileNameIterator << "\" is "
                     << (inputInfo->isColor() ? "color" : "grayscale") << "\n"
                     << command << ":   but previous images are "
                     << (isColor ? "color" : "grayscale")
                     << endl;
                exit(1);
            }
            if (pixelType != inputInfo->getPixelType()) {
                cerr << command << ": input image \""
                     << *inputFileNameIterator << "\" has pixel type "
                     << inputInfo->getPixelType() << ",\n"
                     << command << ":   but previous images have pixel type "
                     << pixelType
                     << endl;
                exit(1);
            }
            if (resolution !=
                TiffResolution(inputInfo->getXResolution(), inputInfo->getYResolution())) {
                cerr << command << ": info: input image \""
                     << *inputFileNameIterator << "\" has resolution "
                     << inputInfo->getXResolution() << " dpi x "
                     << inputInfo->getYResolution() << " dpi,\n"
                     << command << ": info:   but first image has resolution "
                     << resolution.x << " dpi x " << resolution.y << " dpi"
                     << endl;
            }
            if (!std::equal(iccProfile.begin(),
                            iccProfile.end(),
                            inputInfo->getICCProfile().begin())) {
                ImageImportInfo::ICCProfile mismatchProfile = inputInfo->getICCProfile();
                cmsHPROFILE newProfile = NULL;
                if (!mismatchProfile.empty()) {
                    newProfile = cmsOpenProfileFromMem(mismatchProfile.data(),
                                                       mismatchProfile.size());
                    if (newProfile == NULL) {
                        cerr << endl
                             << command << ": error parsing ICC profile data from file \""
                             << *inputFileNameIterator
                             << "\"" << endl;
                        exit(1);
                    }
                }

                cerr << endl << command << ": input image \""
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
                         << "\"" << endl;
                } else {
                    cerr << " no ICC profile" << endl;
                }
                cerr << command
                     << ": warning: blending images with different color spaces\n"
                     << command
                     << ": warning:     may have unexpected results"
                     << endl;
            }
        }

        if (layers == 1 || layer == layers)
        {
            layer = 0;
            layers = 0;
            inputFileNameIterator++;
        }
        else
        {
            // We are about to process the next layer in the _same_
            // image.  The imageInfoList already has been updated, but
            // inputFileNameList still lacks the filename.
            inputFileNameList.insert(inputFileNameIterator, *inputFileNameIterator);
        }
    }

    vigra_postcondition(imageInfoList.size() == inputFileNameList.size(),
                        "filename list and image info list are inconsistent");

    // Check that more than one input file was given.
    if (imageInfoList.size() <= 1) {
        cerr << command
             << ": warning: only one input image given.\n"
             << command
             << ": warning: Enfuse needs two or more overlapping input images in order to do\n"
             << command
             << ": warning: blending calculations.  The output will be the same as the input."
             << endl;
    }

    if (resolution == TiffResolution()) {
        cerr << command
             << ": warning: no usable resolution found in first image \""
             << *inputFileNameList.begin() << "\";\n"
             << command
             << ": warning:   will use " << DEFAULT_TIFF_RESOLUTION << " dpi"
             << endl;
        ImageResolution = TiffResolution(DEFAULT_TIFF_RESOLUTION,
                                         DEFAULT_TIFF_RESOLUTION);
    } else {
        ImageResolution = resolution;
    }

    // Make sure that inputUnion is at least as big as given by the -f paramater.
    if (OutputSizeGiven) {
        inputUnion |= Rect2D(OutputOffsetXCmdLine,
                             OutputOffsetYCmdLine,
                             OutputOffsetXCmdLine + OutputWidthCmdLine,
                             OutputOffsetYCmdLine + OutputHeightCmdLine);
    }

    // Create the Info for the output file.
    ImageExportInfo outputImageInfo(OutputFileName.c_str());
    if (!OutputCompression.empty()) {
        outputImageInfo.setCompression(OutputCompression.c_str());
    }

    // If not overridden by the command line, the pixel type of the
    // output image is the same as the input images'.  If the pixel
    // type is not supported by the output format, replace it with the
    // best match.
    {
        const std::string outputFileType = enblend::getFileType(OutputFileName);
        const std::string neededPixelType =
            OutputPixelType.empty() ? std::string(pixelType) : OutputPixelType;
        const std::string bestPixelType =
            enblend::bestPixelType(outputFileType, neededPixelType);
        if (neededPixelType != bestPixelType) {
            cerr << command
                 << ": warning: "
                 << (OutputPixelType.empty() ? "deduced" : "requested")
                 << " output pixel type is \""
                 << enblend::toLowercase(neededPixelType)
                 << "\", but image type \""
                 << enblend::toLowercase(outputFileType)
                 << "\"\n"
                 << command << ": warning:   supports \""
                 << enblend::toLowercase(bestPixelType)
                 << "\" at best;  will use \""
                 << enblend::toLowercase(bestPixelType)
                 << "\""
                 << endl;
        }
        outputImageInfo.setPixelType(bestPixelType.c_str());
        pixelType = enblend::maxPixelType(pixelType, bestPixelType);
    }

    // Set the output image ICC profile
    outputImageInfo.setICCProfile(iccProfile);

    if (UseCIECAM) {
        if (InputProfile == NULL) {
            cerr << command
                 << ": warning: input images do not have ICC profiles; assuming sRGB."
                 << endl;
            InputProfile = cmsCreate_sRGBProfile();
        }
        XYZProfile = cmsCreateXYZProfile();

        InputToXYZTransform = cmsCreateTransform(InputProfile, TYPE_RGB_DBL,
                                                 XYZProfile, TYPE_XYZ_DBL,
                                                 INTENT_PERCEPTUAL, cmsFLAGS_NOTPRECALC);
        if (InputToXYZTransform == NULL) {
            cerr << command << ": error building color transform from \""
                 << cmsTakeProductName(InputProfile)
                 << " "
                 << cmsTakeProductDesc(InputProfile)
                 << "\" to XYZ" << endl;
            exit(1);
        }

        XYZToInputTransform = cmsCreateTransform(XYZProfile, TYPE_XYZ_DBL,
                                                 InputProfile, TYPE_RGB_DBL,
                                                 INTENT_PERCEPTUAL, cmsFLAGS_NOTPRECALC);
        if (XYZToInputTransform == NULL) {
            cerr << command
                 << ": error building color transform from XYZ to \""
                 << cmsTakeProductName(InputProfile)
                 << " "
                 << cmsTakeProductDesc(InputProfile)
                 << "\"" << endl;
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
            cerr << endl
                 << command
                 << ": error initializing CIECAM02 transform"
                 << endl;
            exit(1);
        }
    }

    // The size of the output image.
    if (Verbose >= VERBOSE_INPUT_UNION_SIZE_MESSAGES) {
        cerr << command
             << ": info: output image size: "
             << inputUnion
             << endl;
    }

    // Set the output image position and resolution.
    outputImageInfo.setXResolution(ImageResolution.x);
    outputImageInfo.setYResolution(ImageResolution.y);
    outputImageInfo.setPosition(inputUnion.upperLeft());

    // Sanity check on the output image file.
    try {
        // This seems to be a reasonable way to check if
        // the output file is going to work after blending
        // is done.
        encoder(outputImageInfo);
    } catch (StdException & e) {
        cerr << endl
             << command
             << ": error opening output file \""
             << OutputFileName
             << "\";\n"
             << command
             << ": "
             << e.what()
             << endl;
        exit(1);
    }

    if (!OutputPixelType.empty()) {
        pixelType = enblend::maxPixelType(pixelType, OutputPixelType);
    }

    // Invoke templatized blender.
    try {
        if (isColor) {
            if      (pixelType == "UINT8")  enfuseMain<RGBValue<UInt8 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "INT8")   enfuseMain<RGBValue<Int8  > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT16") enfuseMain<RGBValue<UInt16> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enfuseMain<RGBValue<Int16 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enfuseMain<RGBValue<UInt32> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enfuseMain<RGBValue<Int32 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enfuseMain<RGBValue<float > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enfuseMain<RGBValue<double> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                cerr << command << ": RGB images with pixel type \""
                     << pixelType
                     << "\" are not supported"
                     << endl;
                exit(1);
            }
        } else {
            if (!WSaturationIsDefault && (WSaturation != 0.0)) {
                cerr << command
                     << ": warning: --WSaturation is not applicable to grayscale images;\n"
                     << command
                     << ": warning:   this parameter will have no effect"
                     << endl;
                WSaturation = 0.0;
            }
            if      (pixelType == "UINT8")  enfuseMain<UInt8 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "INT8")   enfuseMain<Int8  >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT16") enfuseMain<UInt16>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enfuseMain<Int16 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enfuseMain<UInt32>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enfuseMain<Int32 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enfuseMain<float >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enfuseMain<double>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                cerr << command
                     << ": black&white images with pixel type \""
                     << pixelType
                     << "\" are not supported"
                     << endl;
                exit(1);
            }
        }

        for (list<ImageImportInfo*>::iterator i = imageInfoList.begin();
             i != imageInfoList.end();
             ++i) {
            delete *i;
        }
    } catch (std::bad_alloc& e) {
        cerr << endl
             << command << ": out of memory\n"
             << command << ": " << e.what()
             << endl;
        exit(1);
    } catch (StdException& e) {
        cerr << endl
             << command << ": an exception occured\n"
             << command << ": " << e.what()
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
