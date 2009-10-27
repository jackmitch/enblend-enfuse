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

#include <getopt.h>
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
#include "signature.h"
#include "self_test.h"

// Globals
const std::string command("enfuse");

// Random number generator for dithering
boost::mt19937 Twister;

// Global values from command line parameters.
std::string OutputFileName(DEFAULT_OUTPUT_FILENAME);
int Verbose = 1;                //< src::default-verbosity-level 1
unsigned int ExactLevels = 0U;
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
double WExposure = 1.0;         //< src::default-weight-exposure 1.0
double WContrast = 0.0;         //< src::default-weight-contrast 0.0
double WSaturation = 0.2;       //< src::default-weight-saturation 0.2
double WEntropy = 0.0;          //< src::default-weight-entropy 0.0
double WMu = 0.5;               //< src::default-exposure-mu 0.5
double WSigma = 0.2;            //< src::default-exposure-sigma 0.2
bool WSaturationIsDefault = true;
int ContrastWindowSize = 5;
std::string GrayscaleProjector;
struct EdgeFilterConfiguration {double edgeScale, lceScale, lceFactor;} FilterConfig = {
    0.0,                        //< src::default-edge-scale 0.0
    0.0,                        //< src::default-lce-scale 0.0
    0.0                         //< src::default-lce-factor 0.0
};
struct AlternativePercentage MinCurvature = {0.0, false}; //< src::default-minimum-curvature 0
int EntropyWindowSize = 3;      //< src::default-entropy-window-size 3
struct AlternativePercentage EntropyLowerCutoff = {0.0, true}; //< src::default-entropy-lower-cutoff 0%
struct AlternativePercentage EntropyUpperCutoff = {100.0, true}; //< src::default-entropy-upper-cutoff 100%
bool UseHardMask = false;
bool SaveMasks = false;
std::string SoftMaskTemplate("softmask-%n.tif"); //< src::default-soft-mask-template softmask-%n.tif
std::string HardMaskTemplate("hardmask-%n.tif"); //< src::default-hard-mask-template hardmask-%n.tif

TiffResolution ImageResolution;
bool OutputIsValid = true;

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

Signature sig;

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


#define DUMP_GLOBAL_VARIABLES(...) dump_global_variables(__FILE__, __LINE__, ##__VA_ARGS__)
void dump_global_variables(const char* file, unsigned line,
                           std::ostream& out = std::cout)
{
    out <<
        "+ " << file << ":" << line << ": state of global variables\n" <<
        "+ Verbose = " << Verbose << ", option \"--verbose\"\n" <<
        "+ OutputFileName = <" << OutputFileName << ">\n" <<
        "+ ExactLevels = " << ExactLevels << "\n" <<
        "+ OneAtATime = " << enblend::stringOfBool(OneAtATime) << ", option \"-a\"\n" <<
        "+ WrapAround = " << enblend::stringOfWraparound(WrapAround) << ", option \"--wrap\"\n" <<
        "+ GimpAssociatedAlphaHack = " << enblend::stringOfBool(GimpAssociatedAlphaHack) <<
        ", option \"-g\"\n" <<
        "+ UseCIECAM = " << enblend::stringOfBool(UseCIECAM) << ", option \"-c\"\n" <<
        "+ OutputSizeGiven = " << enblend::stringOfBool(OutputSizeGiven) << ", option \"-f\"\n" <<
        "+     OutputWidthCmdLine = " << OutputWidthCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputHeightCmdLine = " << OutputHeightCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputOffsetXCmdLine = " << OutputOffsetXCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputOffsetYCmdLine = " << OutputOffsetYCmdLine << ", argument to option \"-f\"\n" <<
        "+ WExposure = " << WExposure << ", argument to option \"--exposure-weight\"\n" <<
        "+     WMu = " << WMu  << ", argument to option \"--exposure-mu\"\n" <<
        "+     WSigma = " << WSigma << ", argument to option \"--exposure-sigma\"\n" <<
        "+ WContrast = " << WContrast << ", argument to option \"--contrast-weight\"\n" <<
        "+ WSaturation = " << WSaturation << ", argument to option \"--saturation-weight\"\n" <<
        "+ WEntropy = " << WEntropy << ", argument to option \"--entropy-weight\"\n" <<
        "+ WSaturationIsDefault = " << enblend::stringOfBool(WSaturationIsDefault) << "\n" <<
        "+ ContrastWindowSize = " << ContrastWindowSize <<
        ", argument to option \"--contrast-window-size\"\n" <<
        "+ GrayscaleProjector = <" << GrayscaleProjector <<
        ">, argument to option \"--gray-projector\"\n" <<
        "+ FilterConfig = {\n" <<
        "+     edgeScale = " << FilterConfig.edgeScale << ",\n" <<
        "+     lceScale = " << FilterConfig.lceScale <<  ",\n" <<
        "+     lceFactor = " << FilterConfig.lceFactor <<  "\n" <<
        "+ }, arguments to option \"--contrast-edge-scale\"\n" <<
        "+ MinCurvature = {\n"
        "+     value = " << MinCurvature.value << ",\n" <<
        "+     isPercentage = " << enblend::stringOfBool(MinCurvature.isPercentage) << "\n" <<
        "+ }, arguments to option \"--contrast-min-curvature\"\n" <<
        "+ EntropyWindowSize = " << EntropyWindowSize <<
        ", argument to option \"--entropy-window-size\"\n" <<
        "+ EntropyLowerCutoff = {\n"
        "+     value = " << EntropyLowerCutoff.value << ",\n" <<
        "+     isPercentage = " << enblend::stringOfBool(EntropyLowerCutoff.isPercentage) << "\n" <<
        "+ }, first argument to option \"--entropy-cutoff\"\n" <<
        "+ EntropyUpperCutoff = {\n"
        "+     value = " << EntropyUpperCutoff.value << ",\n" <<
        "+     isPercentage = " << enblend::stringOfBool(EntropyUpperCutoff.isPercentage) << "\n" <<
        "+ }, second argument to option \"--entropy-cutoff\"\n" <<
        "+ UseHardMask = " << enblend::stringOfBool(UseHardMask) <<
        ", option \"--hard-mask\" or \"--soft-mask\"\n" <<
        "+ SaveMasks = " << enblend::stringOfBool(SaveMasks) << ", option \"--save-masks\"\n" <<
        "+     SoftMaskTemplate = <" << SoftMaskTemplate <<
        ">, first argument to option \"--save-masks\"\n" <<
        "+     HardMaskTemplate = <" << HardMaskTemplate <<
        ">, second argument to option \"--save-masks\"\n" <<
        "+ OutputCompression = <" << OutputCompression << ">, option \"--compression\"\n" <<
        "+ OutputPixelType = <" << OutputPixelType << ">, option \"--depth\"\n" <<
        "+ end of global variable dump\n";
}


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

#ifdef CACHE_IMAGES
        cout << "Extra feature: image cache: yes\n";
        {
#ifdef WIN32
            char lpPathBuffer[MAX_PATH];
            const DWORD dwRetVal = GetTempPath(MAX_PATH, lpPathBuffer);
            if (dwRetVal <= MAX_PATH && dwRetVal != 0) {
                cout << "  - cache file located in \"" << lpPathBuffer << "\"\n";
            }
#else
            const char* tmpdir = getenv("TMPDIR");
            cout << "  - environment variable TMPDIR ";
            if (tmpdir == NULL) {
                cout << "not set, cache file in default directory \"/tmp\"\n";
            } else {
                cout << "set, cache file located in \"" << tmpdir << "\"\n";
            }
#endif
        }
#else
        cout << "Extra feature: image cache: no\n";
#endif

#ifdef OPENMP
        const bool have_nested = have_openmp_nested();
        const bool have_dynamic = have_openmp_dynamic();
        cout <<
            "Extra feature: OpenMP: yes\n" <<
            "  - version " << OPENMP_YEAR << '-' << OPENMP_MONTH << "\n" <<
            "  - " << (have_nested ? "" : "no ") <<
            "support for nested parallelism;\n" <<
            "    nested parallelism " <<
            (have_nested && omp_get_nested() ? "enabled" : "disabled") << " by default\n" <<
            "  - " << (have_dynamic ? "" : "no ") <<
            "support for dynamic adjustment of the number of threads;\n" <<
            "    dynamic adjustment " <<
            (have_dynamic && omp_get_dynamic() ? "enabled" : "disabled") << " by default\n" <<
            "  - using " <<
            omp_get_num_procs() << " processor" << (omp_get_num_procs() >= 2 ? "s" : "") << " and up to " <<
            omp_get_max_threads() << " thread" << (omp_get_max_threads() >= 2 ? "s" : "") << "\n";
#else
        cout << "Extra feature: OpenMP: no\n";
#endif

        cout <<
            "\n" <<
            "Supported image formats: " << vigra::impexListFormats() << "\n" <<
            "Supported file extensions: " << vigra::impexListExtensions() << "\n\n";
    }

    if (Verbose >= VERBOSE_SIGNATURE_REPORTING) {
        cout.flush();
        std::wcout << sig.message() << L"\n\n";
        std::wcout.flush();
    }

    cout <<
        "Copyright (C) 2004-2009 Andrew Mihal.\n" <<
        "License GPLv2+: GNU GPL version 2 or later <http://www.gnu.org/licenses/gpl.html>\n" <<
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
        "  --compression=COMPRESSION\n" <<
        "                         set compression of output image to COMPRESSION,\n" <<
        "                         where COMPRESSION is:\n" <<
        "                         NONE, PACKBITS, LZW, DEFLATE for TIFF files and\n" <<
        "                         0 to 100 for JPEG files\n" <<
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
        "  --exposure-weight=WEIGHT\n" <<
        "                         weight given to well-exposed pixels\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WExposure << "\n" <<
        "  --saturation-weight=WEIGHT\n" <<
        "                         weight given to highly-saturated pixels\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WSaturation << "\n" <<
        "  --contrast-weight=WEIGHT\n" <<
        "                         weight given to pixels in high-contrast neighborhoods \n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WContrast << "\n" <<
        "  --entropy-weight=WEIGHT\n" <<
        "                         weight given to pixels in high entropy neighborhoods\n" <<
        "                         (0 <= WEIGHT <= 1); default: " << WEntropy << "\n" <<
        "  --exposure-mu=MEAN     center also known as MEAN of Gaussian weighting\n" <<
        "                         function (0 <= MEAN <= 1); default: " << WMu << "\n" <<
        "  --exposure-sigma=SIGMA standard deviation of Gaussian weighting\n" <<
        "                         function (SIGMA > 0); default: " << WSigma << "\n" <<
        "  --soft-mask            average over all masks; this is the default\n" <<
        "  --hard-mask            force hard blend masks and no averaging on finest\n" <<
        "                         scale; this is especially useful for focus\n" <<
        "                         stacks with thin and high contrast features,\n" <<
        "                         but leads to increased noise\n" <<
        "\n" <<
        "Expert options:\n" <<
        "  --contrast-window-size=SIZE\n" <<
        "                         set window SIZE for local-contrast analysis\n" <<
        "                         (SIZE >= 3); default: " << ContrastWindowSize  << "\n" <<
        "  --gray-projector=OPERATOR\n" <<
        "                         apply gray-scale projection OPERATOR in exposure\n" <<
        "                         or contrast weighing, where OPERATOR is one of\n" <<
        "                         \"average\", \"l-star\", \"lightness\", \"value\",\n" <<
        "                         \"luminance\", or\n" <<
        "                         \"channel-mixer:RED-WEIGHT:GREEN-WEIGHT:BLUE-WEIGHT\";\n" <<
        "                         default: \"" <<
        enblend::MultiGrayscaleAccessor<UInt8, NumericTraits<UInt8>::Promote>::defaultGrayscaleAccessorName() << "\"\n" <<
        "  --contrast-edge-scale=EDGESCALE[:LCESCALE[:LCEFACTOR]]\n" <<
        "                         set scale on which to look for edges; positive\n" <<
        "                         LCESCALE switches on local contrast enhancement\n" <<
        "                         by LCEFACTOR (EDGESCALE, LCESCALE, LCEFACTOR >= 0);\n" <<
        "                         append \"%\" to LCESCALE for values relative to\n" <<
        "                         EDGESCALE; append \"%\" to LCEFACTOR for relative\n" <<
        "                         value; defaults: " <<
        FilterConfig.edgeScale << ":" << FilterConfig.lceScale << ":" << FilterConfig.lceFactor << "\n" <<
        "  --contrast-min-curvature=CURVATURE\n" <<
        "                         minimum CURVATURE for an edge to qualify; append\n" <<
        "                         \"%\" for relative values; default: " << MinCurvature.str() << "\n" <<
        "  --entropy-window-size=SIZE\n" <<
        "                         set window SIZE for local entropy analysis\n" <<
        "                         (SIZE >= 3); default: " << EntropyWindowSize  << "\n" <<
        "  --entropy-cutoff=LOWERCUTOFF[:UPPERCUTOFF]\n" <<
        "                         LOWERCUTOFF is the value below of which pixels are\n" <<
        "                         treated as black and UPPERCUTOFF is the value above\n" <<
        "                         of which pixels are treated as white in the entropy\n" <<
        "                         weighting; append \"%\" signs for relative values;\n" <<
        "                         default: " <<
        EntropyLowerCutoff.str() << ":" << EntropyUpperCutoff.str() << "\n" <<
        "  --save-masks[=SOFT-TEMPLATE[:HARD-TEMPLATE]]\n" <<
        "                         save weight masks in SOFT-TEMPLATE and HARD-TEMPLATE;\n" <<
        "                         conversion chars: %i: mask index, %n: mask number,\n" <<
        "                         %p: full path, %d: dirname, %b: basename,\n" <<
        "                         %f: filename, %e: extension; lowercase characters\n" <<
        "                         refer to input images uppercase to the output image;\n" <<
        "                         default: \"" << SoftMaskTemplate << "\":\"" << HardMaskTemplate << "\"\n" <<
        "\n" <<
        "Report bugs at <" PACKAGE_BUGREPORT ">." <<
        endl;

    exit(error ? 1 : 0);
}


void cleanup_output(void)
{
#if DEBUG
    std::cout << "+ cleanup_output\n";
#endif

    if (!OutputIsValid) {
        std::cerr << command << ": info: remove invalid output image \"" << OutputFileName << "\"\n";
        errno = 0;
        if (unlink(OutputFileName.c_str()) != 0) {
            cerr << command <<
                ": warning: could not remove invalid output image \"" << OutputFileName << "\": " <<
                enblend::errorMessage(errno) << "\n";
        }
    }
}


/** Make sure all cached file images get destroyed, and hence the
 *  temporary files deleted, if we are killed.
 */
void sigint_handler(int sig)
{
    cerr << endl << command << ": interrupted" << endl;

    cleanup_output();

#if !defined(__GW32C__) && !defined(_WIN32)
    struct sigaction action;
    action.sa_handler = SIG_DFL;
    sigemptyset(&(action.sa_mask));
    sigaction(sig, &action, NULL);
#else
    signal(sig, SIG_DFL);
#endif
    raise(sig);
}


enum AllPossibleOptions {
    VersionOption, HelpOption, LevelsOption, OutputOption, VerboseOption,
    WrapAroundOption /* -w */, CompressionOption, LZWCompressionOption,
    BlockSizeOption, CIECAM02Option /* -c */,
    DepthOption, AssociatedAlphaOption /* -g */,
    SizeAndPositionOption /* -f */, CacheSizeOption,
    ExposureWeightOption, SaturationWeightOption,
    ContrastWeightOption, EntropyWeightOption,
    ExposureMuOption /* --contrast-mu */, ExposureSigmaOption /* --contrast-sigma */,
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
                ": warning: option \"--exposure-mu\" has no effect as exposure weight\n" <<
                command <<
                ": warning:     is zero" <<
                endl;
        }
        if (contains(optionSet, ExposureSigmaOption)) {
            cerr << command <<
                ": warning: option \"--exposure-sigma\" has no effect as exposure weight\n" <<
                command <<
                ": warning:     is zero" <<
                endl;
        }
    }

    if (WContrast == 0.0 && contains(optionSet, ContrastWindowSizeOption)) {
        cerr << command <<
            ": warning: option \"--contrast-window-size\" has no effect as contrast\n" <<
            command <<
            ": warning:     weight is zero" <<
            endl;
    }

    if (WExposure == 0.0 && WContrast == 0.0 && contains(optionSet, GrayProjectorOption)) {
        cerr << command <<
            ": warning: option \"--gray-projector\" has no effect as exposure\n" <<
            command <<
            ": warning:     and contrast weight both are zero" <<
            endl;
    }

    if (WContrast == 0.0) {
        if (contains(optionSet, EdgeScaleOption)) {
            cerr << command <<
                ": warning: option \"--contrast-edge-scale\" has no effect as contrast\n" <<
                command <<
                ": warning:     weight is zero" <<
                endl;
        }
        if (contains(optionSet, MinCurvatureOption)) {
            cerr << command <<
                ": warning: option \"--contrast-min-curvature\" has no effect as contrast\n" <<
                command <<
                ": warning:     weight is zero" <<
                endl;
        }
    } else {
        if (FilterConfig.edgeScale > 0.0 &&
            contains(optionSet, ContrastWindowSizeOption) && MinCurvature.value <= 0.0) {
            cerr << command <<
                ": warning: option \"--contrast-window-size\" has no effect as\n" <<
                command <<
                ": warning:     EDGESCALE in \"--contrast-edge-scale\" is positive and" <<
                command <<
                ": warning:     \"--contrast-min-curvature\" is non-positive" <<
                endl;
        }
    }

    if (WEntropy == 0.0) {
        if (contains(optionSet, EntropyWindowSizeOption)) {
            cerr << command <<
                ": warning: option \"--entropy-window-size\" has no effect as\n" <<
                command <<
                ": warning:     entropy weight is zero" <<
                endl;
        }
        if (contains(optionSet, EntropyCutoffOption)) {
            cerr << command <<
                ": warning: option \"--entropy-cutoff\" has no effect as entropy\n" <<
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

#ifndef CACHE_IMAGES
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


int process_options(int argc, char** argv)
{
    enum OptionId {
        OPTION_ID_OFFSET = 1023,    // Ids start at 1024
        CompressionId,
        WeightExposureId, Deprecated_WeightExposureId,
        WeightContrastId, Deprecated_WeightContrastId,
        WeightSaturationId, Deprecated_WeightSaturationId,
        WeightMuId, Deprecated_WeightMuId,
        WeightSigmaId, Deprecated_WeightSigmaId,
        MinCurvatureId, Deprecated_MinCurvatureId,
        EdgeScaleId, Deprecated_EdgeScaleId,
        ContrastWindowSizeId, Deprecated_ContrastWindowSizeId,
        HardMaskId, Deprecated_HardMaskId,
        GrayProjectorId, Deprecated_GrayProjectorId,
        WeightEntropyId, Deprecated_WeightEntropyId,
        EntropyWindowSizeId, Deprecated_EntropyWindowSizeId,
        EntropyCutoffId, Deprecated_EntropyCutoffId,
        SoftMaskId, Deprecated_SoftMaskId,
        VerboseId,
        HelpId,
        VersionId,
        DepthId,
        OutputId,
        SaveMasksId, Deprecated_SaveMasksId,
        WrapAroundId
    };

    static struct option long_options[] = {
        {"compression", required_argument, 0, CompressionId},
        {"exposure-weight", required_argument, 0, WeightExposureId},
        {"wExposure", required_argument, 0, Deprecated_WeightExposureId},
        {"contrast-weight", required_argument, 0, WeightContrastId},
        {"wContrast", required_argument, 0, Deprecated_WeightContrastId},
        {"saturation-weight", required_argument, 0, WeightSaturationId},
        {"wSaturation", required_argument, 0, Deprecated_WeightSaturationId},
        {"exposure-mu", required_argument, 0, WeightMuId},
        {"wMu", required_argument, 0, Deprecated_WeightMuId},
        {"exposure-sigma", required_argument, 0, WeightSigmaId},
        {"wSigma", required_argument, 0, Deprecated_WeightSigmaId},
        {"contrast-min-curvature", required_argument, 0, MinCurvatureId},
        {"MinCurvature", required_argument, 0, Deprecated_MinCurvatureId},
        {"contrast-edge-scale", required_argument, 0, EdgeScaleId},
        {"EdgeScale", required_argument, 0, Deprecated_EdgeScaleId},
        {"contrast-window-size", required_argument, 0, ContrastWindowSizeId},
        {"ContrastWindowSize", required_argument, 0, Deprecated_ContrastWindowSizeId},
        {"hard-mask", no_argument, 0, HardMaskId},
        {"HardMask", no_argument, 0, Deprecated_HardMaskId},
        {"gray-projector", required_argument, 0, GrayProjectorId},
        {"GrayProjector", required_argument, 0, Deprecated_GrayProjectorId},
        {"entropy-weight", required_argument, 0, WeightEntropyId},
        {"wEntropy", required_argument, 0, Deprecated_WeightEntropyId},
        {"entropy-window-size", required_argument, 0, EntropyWindowSizeId},
        {"EntropyWindowSize", required_argument, 0, Deprecated_EntropyWindowSizeId},
        {"entropy-cutoff", required_argument, 0, EntropyCutoffId},
        {"EntropyCutoff", required_argument, 0, Deprecated_EntropyCutoffId},
        {"soft-mask", no_argument, 0, SoftMaskId},
        {"SoftMask", no_argument, 0, Deprecated_SoftMaskId},
        {"verbose", optional_argument, 0, VerboseId},
        {"help", no_argument, 0, HelpId},
        {"version", no_argument, 0, VersionId},
        {"depth", required_argument, 0, DepthId},
        {"output", required_argument, 0, OutputId},
        {"save-mask", optional_argument, 0, SaveMasksId}, // singular form: not documented, not deprecated
        {"save-masks", optional_argument, 0, SaveMasksId},
        {"SaveMasks", optional_argument, 0, Deprecated_SaveMasksId},
        {"wrap", optional_argument, 0, WrapAroundId},
        {0, 0, 0, 0}
    };

    bool failed = false;
    bool justPrintVersion = false;
    bool justPrintUsage = false;
    OptionSetType optionSet;

    opterr = 0;       // we have our own "unrecognized option" message
    while (true) {
        int option_index;
        const int code = getopt_long(argc, argv, "Vb:cd:f:ghl:m:o:v::w::",
                                     long_options, &option_index);

        if (code == -1) {
            break;
        }

        switch (code) {
        case Deprecated_HardMaskId:
            PENALIZE_DEPRECATED_OPTION("--HardMask", "--hard-mask");
            // FALLTHROUGH
        case HardMaskId:
            UseHardMask = true;
            optionSet.insert(HardMaskOption);
            break;

        case Deprecated_SoftMaskId:
            PENALIZE_DEPRECATED_OPTION("--SoftMask", "--soft-mask");
            // FALLTHROUGH
        case SoftMaskId:
            UseHardMask = false;
            optionSet.insert(SoftMaskOption);
            break;

        case 'h': // FALLTHROUGH
        case HelpId:
            justPrintUsage = true;
            optionSet.insert(HelpOption);
            break;

        case 'V': // FALLTHROUGH
        case VersionId:
            justPrintVersion = true;
            optionSet.insert(VersionOption);
            break;

        case 'w': // FALLTHROUGH
        case WrapAroundId:
            if (optarg != NULL && *optarg != 0) {
                WrapAround = enblend::wraparoundOfString(optarg);
                if (WrapAround == UnknownWrapAround) {
                    cerr << command
                         << ": unrecognized wrap-around mode \"" << optarg << "\"\n" << endl;
                    failed = true;
                }
            } else {
                WrapAround = HorizontalStrip;
            }
            optionSet.insert(WrapAroundOption);
            break;

        case Deprecated_MinCurvatureId:
            PENALIZE_DEPRECATED_OPTION("--MinCurvature", "--contrast-min-curvature");
            // FALLTHROUGH
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
                    failed = true;
                }
            } else {
                cerr << command << ": illegal numeric format \""
                     << optarg << "\" for minimum gradient: "
                     << enblend::errorMessage(errno) << endl;
                failed = true;
            }
            optionSet.insert(MinCurvatureOption);
            break;
        }

        case Deprecated_EdgeScaleId:
            PENALIZE_DEPRECATED_OPTION("--EdgeScale", "--contrast-edge-scale");
            // FALLTHROUGH
        case EdgeScaleId: {
            char* s = new char[strlen(optarg) + 1];
            strcpy(s, optarg);
            char* save_ptr = NULL;
            char* token = enblend::strtoken_r(s, NUMERIC_OPTION_DELIMITERS, &save_ptr);
            char* tail;

            if (token == NULL || *token == 0) {
                cerr << command << ": no scale given to --contrast-edge-scale.  "
                     << "scale is required." << endl;
                failed = true;
            }
            errno = 0;
            FilterConfig.edgeScale = strtod(token, &tail);
            if (errno == 0) {
                if (*tail != 0) {
                    cerr << command << ": could not decode \"" << tail
                         << "\" in edge scale specification \""
                         << token << "\" for edge scale." << endl;
                    failed = true;
                }
            } else {
                cerr << command << ": illegal numeric format \""
                     << token << "\" for edge scale: "
                     << enblend::errorMessage(errno) << endl;
                failed = true;
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
                        failed = true;
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << token << "\" for LCE-Scale: "
                         << enblend::errorMessage(errno) << endl;
                    failed = true;
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
                        failed = true;
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << token << "\" for LCE-factor: "
                         << enblend::errorMessage(errno) << endl;
                    failed = true;
                }
            }

            if (save_ptr != NULL && *save_ptr != 0) {
                cerr << command << ": warning: ignoring trailing garbage \""
                     << save_ptr << "\" in argument to --contrast-edge-scale" << endl;
            }

            delete [] s;
            optionSet.insert(EdgeScaleOption);
            break;
        }

        case Deprecated_EntropyCutoffId:
            PENALIZE_DEPRECATED_OPTION("--EntropyCutoff", "--entropy-cutoff");
            // FALLTHROUGH
        case EntropyCutoffId: {
            char* s = new char[strlen(optarg) + 1];
            strcpy(s, optarg);
            char* save_ptr = NULL;
            char* token = enblend::strtoken_r(s, NUMERIC_OPTION_DELIMITERS, &save_ptr);
            char* tail;

            if (token == NULL || *token == 0) {
                cerr << command << ": no scale given to --entropy-cutoff.  "
                     << "lower cutoff is required." << endl;
                failed = true;
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
                    failed = true;
                }
            } else {
                cerr << command << ": illegal numeric format \""
                     << token << "\" of entropy's lower cutoff: "
                     << enblend::errorMessage(errno) << endl;
                failed = true;
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
                        failed = true;
                    }
                } else {
                    cerr << command << ": illegal numeric format \""
                         << token << "\" of entropy's upper cutoff: "
                         << enblend::errorMessage(errno) << endl;
                    failed = true;
                }
            }

            if (save_ptr != NULL && *save_ptr != 0) {
                cerr << command << ": warning: ignoring trailing garbage \""
                     << save_ptr << "\" in argument to --entropy-cutoff" << endl;
            }

            delete [] s;
            optionSet.insert(EntropyCutoffOption);
            break;
        }

        case CompressionId:
            if (optarg != NULL && *optarg != 0) {
                const std::string upper_opt(enblend::toUppercase(optarg));
                if (upper_opt == "NONE") {
                    ;           // stick with default
                } else if (upper_opt == "DEFLATE" || upper_opt == "LZW" || upper_opt == "PACKBITS" ||
                           upper_opt.find_first_not_of("0123456789") == std::string::npos) {
                    OutputCompression = upper_opt;
                } else {
                    cerr << command << ": unrecognized compression \"" << optarg << "\"" << endl;
                    failed = true;
                }
            } else {
                cerr << command << ": option \"--compression\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(CompressionOption);
            break;

        case Deprecated_GrayProjectorId:
            PENALIZE_DEPRECATED_OPTION("--GrayProjector", "--gray-projector");
            // FALLTHROUGH
        case GrayProjectorId:
            if (optarg != NULL && *optarg != 0) {
                GrayscaleProjector = optarg;
            } else {
                cerr << command << ": option \"--gray-projector\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(GrayProjectorOption);
            break;

        case 'd': // FALLTHROUGH
        case DepthId:
            if (optarg != NULL && *optarg != 0) {
                OutputPixelType = enblend::outputPixelTypeOfString(optarg);
            } else {
                cerr << command << ": options \"-d\" or \"--depth\" require arguments" << endl;
                failed = true;
            }
            optionSet.insert(DepthOption);
            break;

        case 'o': // FALLTHROUGH
        case OutputId:
            if (contains(optionSet, OutputOption)) {
                cerr << command
                     << ": warning: more than one output file specified"
                     << endl;
            }
            if (optarg != NULL && *optarg != 0) {
                OutputFileName = optarg;
            } else {
                cerr << command << ": options \"-o\" or \"--output\" require arguments" << endl;
                failed = true;
            }
            optionSet.insert(OutputOption);
            break;

        case Deprecated_SaveMasksId:
            PENALIZE_DEPRECATED_OPTION("--SaveMasks", "--save-masks");
            // FALLTHROUGH
        case SaveMasksId:
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
                         << ": warning: ignoring trailing garbage in --save-masks" << endl;
                }
                delete [] s;
            }
            SaveMasks = true;
            optionSet.insert(SaveMasksOption);
            break;

        case Deprecated_WeightExposureId:
            PENALIZE_DEPRECATED_OPTION("--wExposure", "--exposure-weight");
            // FALLTHROUGH
        case WeightExposureId:
            if (optarg != NULL && *optarg != 0) {
                WExposure = enblend::numberOfString(optarg,
                                                    _1 >= 0.0, //< src::minimum-weight-exposure 0
                                                    "exposure weight less than 0; will use 0",
                                                    0.0,
                                                    _1 <= 1.0, //< src::maximum-weight-exposure 1
                                                    "exposure weight greater than 1; will use 1",
                                                    1.0);
            } else {
                cerr << command << ": option \"--exposure-weight\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(ExposureWeightOption);
            break;

        case Deprecated_WeightContrastId:
            PENALIZE_DEPRECATED_OPTION("--wContrast", "--contrast-weight");
            // FALLTHROUGH
        case WeightContrastId:
            if (optarg != NULL && *optarg != 0) {
                WContrast = enblend::numberOfString(optarg,
                                                    _1 >= 0.0, //< src::minimum-weight-contrast 0
                                                    "contrast weight less than 0; will use 0",
                                                    0.0,
                                                    _1 <= 1.0, //< src::maximum-weight-contrast 1
                                                    "contrast weight greater than 1; will use 1",
                                                    0.0);
            } else {
                cerr << command << ": option \"--contrast-weight\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(ContrastWeightOption);
            break;

        case Deprecated_WeightSaturationId:
            PENALIZE_DEPRECATED_OPTION("--wSaturation", "--saturation-weight");
            // FALLTHROUGH
        case WeightSaturationId:
            if (optarg != NULL && *optarg != 0) {
                WSaturation = enblend::numberOfString(optarg,
                                                      _1 >= 0.0, //< src::minimum-weight-saturation 0
                                                      "saturation weight less than 0; will use 0",
                                                      0.0,
                                                      _1 <= 1.0, //< src::maximum-weight-saturation 1
                                                      "saturation weight greater than 1; will use 1",
                                                      1.0);
            } else {
                cerr << command << ": option \"--saturation-weight\" requires an argument" << endl;
                failed = true;
            }
            WSaturationIsDefault = false;
            optionSet.insert(SaturationWeightOption);
            break;

        case Deprecated_WeightMuId:
            PENALIZE_DEPRECATED_OPTION("--wMu", "--exposure-mu");
            // FALLTHROUGH
        case WeightMuId:
            if (optarg != NULL && *optarg != 0) {
                WMu = enblend::numberOfString(optarg,
                                              _1 >= 0.0, //< src::minimum-exposure-mu 0
                                              "exposure center value (mu) less than 0; will use 0",
                                              0.0,
                                              _1 <= 1.0, //< src::maximum-exposure-mu 1
                                              "exposure center value (mu) geater than 1; will use 1",
                                              1.0);
            } else {
                cerr << command << ": option \"--exposure-mu\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(ExposureMuOption);
            break;

        case Deprecated_WeightSigmaId:
            PENALIZE_DEPRECATED_OPTION("--wSigma", "--exposure-sigma");
            // FALLTHROUGH
        case WeightSigmaId:
            if (optarg != NULL && *optarg != 0) {
                WSigma = enblend::numberOfString(optarg,
                                                 _1 >= 0.0, //< src::minimum-exposure-sigma 0
                                                 "exposure standard deviation (sigma) less than 0; will use 1/1024",
                                                 1.0 / 1024.0);
            } else {
                cerr << command << ": option \"--exposure-sigma\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(ExposureSigmaOption);
            break;

        case Deprecated_WeightEntropyId:
            PENALIZE_DEPRECATED_OPTION("--wEntropy", "--entropy-weight");
            // FALLTHROUGH
        case WeightEntropyId:
            if (optarg != NULL && *optarg != 0) {
                WEntropy = enblend::numberOfString(optarg,
                                                   _1 >= 0.0, //< src::minimum-weight-entropy 0
                                                   "entropy weight less than 0; will use 0",
                                                   0.0,
                                                   _1 <= 1.0, //< src::maximum-weight-entropy 1
                                                   "entropy weight greater than 1; will use 1",
                                                   1.0);
            } else {
                cerr << command << ": option \"--entropy-weight\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(EntropyWeightOption);
            break;

        case 'v': // FALLTHROUGH
        case VerboseId:
            if (optarg != NULL && *optarg != 0) {
                Verbose = enblend::numberOfString(optarg,
                                                  _1 >= 0, //< src::minimum-verbosity-level 0
                                                  "verbosity level less than 0; will use 0",
                                                  0);
            } else {
                Verbose++;
            }
            optionSet.insert(VerboseOption);
            break;

        case Deprecated_ContrastWindowSizeId:
            PENALIZE_DEPRECATED_OPTION("--ContrastWindowSize", "--contrast-window-size");
            // FALLTHROUGH
        case ContrastWindowSizeId:
            if (optarg != NULL && *optarg != 0) {
                ContrastWindowSize =
                    enblend::numberOfString(optarg,
                                            _1 >= 3, //< src::minimum-contrast-window-size 3
                                            "contrast window size to small; will use size = 3",
                                            3);
                if (ContrastWindowSize % 2 != 1) {
                    cerr << command << ": warning: contrast window size \""
                         << ContrastWindowSize << "\" is even; increasing size to next odd number"
                         << endl;
                    ContrastWindowSize++;
                }
            } else {
                cerr << command << ": option \"--contrast-window-size\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(ContrastWindowSizeOption);
            break;

        case Deprecated_EntropyWindowSizeId:
            PENALIZE_DEPRECATED_OPTION("--EntropyWindowSize", "--entropy-window-size");
            // FALLTHROUGH
        case EntropyWindowSizeId:
            if (optarg != NULL && *optarg != 0) {
                EntropyWindowSize =
                    enblend::numberOfString(optarg,
                                            _1 >= 3, //< src::minimum-entropy-window-size 3
                                            "entropy window size to small; will use size = 3",
                                            3);
                if (EntropyWindowSize % 2 != 1) {
                    cerr << command << ": warning: entropy window size \""
                         << EntropyWindowSize << "\" is even; increasing size to next odd number"
                         << endl;
                    EntropyWindowSize++;
                }
            } else {
                cerr << command << ": option \"--entropy-window-size\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(EntropyWindowSizeOption);
            break;

        case 'b':
            if (optarg != NULL && *optarg != 0) {
                const int cache_block_size =
                    enblend::numberOfString(optarg,
                                            _1 >= 1, //< src::minimum-cache-block-size 1@dmn{KB}
                                            "cache block size must be 1 KB or more; will use 1 KB",
                                            1);
                CachedFileImageDirector::v().setBlockSize(static_cast<long long>(cache_block_size) << 10);
            } else {
                cerr << command << ": option \"-b\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(BlockSizeOption);
            break;

        case 'c':
            UseCIECAM = true;
            optionSet.insert(CIECAM02Option);
            break;

        case 'f':
            if (optarg != NULL && *optarg != 0) {
                const int n = sscanf(optarg,
                                     "%dx%d+%d+%d",
                                     &OutputWidthCmdLine, &OutputHeightCmdLine,
                                     &OutputOffsetXCmdLine, &OutputOffsetYCmdLine);
                if (n == 4) {
                    ; // ok: full geometry string
                } else if (n == 2) {
                    OutputOffsetXCmdLine = 0;
                    OutputOffsetYCmdLine = 0;
                } else {
                    cerr << command << ": option \"-f\" requires 2 or 4 arguments" << endl;
                    failed = true;
                }
            } else {
                cerr << command << ": option \"-f\" requires 2 or 4 arguments" << endl;
                failed = true;
            }
            OutputSizeGiven = true;
            optionSet.insert(SizeAndPositionOption);
            break;

        case 'g':
            GimpAssociatedAlphaHack = true;
            optionSet.insert(AssociatedAlphaOption);
            break;

        case 'l':
            // We shall take care of "too many levels" in "bounds.h".
            if (optarg != NULL && *optarg != 0) {
                ExactLevels =
                    enblend::numberOfString(optarg,
                                            _1 >= 1U, //< src::minimum-pyramid-levels 1
                                            "too few levels; will use one level",
                                            1U);
            } else {
                cerr << command << ": option \"-l\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(LevelsOption);
            break;

        case 'm':
            if (optarg != NULL && *optarg != 0) {
                const int cache_size =
                    enblend::numberOfString(optarg,
                                            _1 >= 1, //< src::minimum-cache-size 1@dmn{MB}
                                            "cache memory limit less than 1 MB; will use 1 MB",
                                            1);
                CachedFileImageDirector::v().setAllocation(static_cast<long long>(cache_size) << 20);
            } else {
                cerr << command << ": option \"-m\" requires an argument" << endl;
                failed = true;
            }
            optionSet.insert(CacheSizeOption);
            break;

        case '?':
            switch (optopt) {
            case 0: // unknown long option
                cerr << command << ": unknown option \"" << argv[optind - 1] << "\"\n";
                break;
            case 'b':           // FALLTHROUGH
            case 'd':           // FALLTHROUGH
            case 'f':           // FALLTHROUGH
            case 'l':           // FALLTHROUGH
            case 'm':           // FALLTHROUGH
            case 'o':
                cerr << command
                     << ": option \"-" << static_cast<char>(optopt) << "\" requires an argument"
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

    if (failed) {
        exit(1);
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

    if (atexit(cleanup_output) != 0) {
        cerr << command << ": warning: could not install cleanup routine\n";
    }

    sig.initialize();

    // Make sure libtiff is compiled with TIF_PLATFORM_CONSOLE
    // to avoid interactive warning dialogs.
    //TIFFSetWarningHandler(NULL);
    //TIFFSetErrorHandler(NULL);

    // List of input files.
    list<char*> inputFileNameList;
    list<char*>::iterator inputFileNameIterator;

    if (!getopt_long_works_ok())
    {
        cerr << command << ": cannot reliably parse command line; giving up\n";
        exit(1);
    }

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

#ifdef DEBUG_DUMP_GLOBAL_VARIABLES
    DUMP_GLOBAL_VARIABLES();
#endif

    sig.check();

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
        {
            std::string filename(*inputFileNameIterator);
            enblend::try_opening_file(filename);
            ImageImportInfo info(filename.c_str());
            if (layers == 0) { // OPTIMIZATION: call only once per file
                layers = info.numLayers();
            }
            if (layers >= 2) {
                filename = vigra::join_filename_layer(*inputFileNameIterator, layer);
            }
            ++layer;
            inputInfo = new ImageImportInfo(filename.c_str());
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
    OutputIsValid = false;
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
