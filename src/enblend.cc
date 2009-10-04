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
#include "signature.h"

typedef struct {
    unsigned int kmax;          // maximum number of moves for a line segment
    double tau;                 // temperature reduction factor, "cooling factor"; 0 < tau < 1
    double deltaEMax;           // maximum cost change possible by any single annealing move
    double deltaEMin;           // minimum cost change possible by any single annealing move
} anneal_para_t;

// Globals
const std::string command("enblend");
const int minimumVectorizeDistance = 4; //< src::minimum-vectorize-distance 4
const int coarseMaskVectorizeDistance = 4; //< src::coarse-mask-vectorize-distance 4
const int fineMaskVectorizeDistance = 20; //< src::fine-mask-vectorize-distance 20

// Random number generator for dithering
boost::mt19937 Twister;

// Global values from command line parameters.
int Verbose = 1;                //< src::default-verbosity-level 1
std::string OutputFileName(DEFAULT_OUTPUT_FILENAME);
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
bool Checkpoint = false;
bool UseGPU = false;
bool OptimizeMask = true;
bool CoarseMask = true;
unsigned CoarsenessFactor = 8U; //< src::default-coarseness-factor 8
double DifferenceBlurRadius = 0.0;
bool SaveMasks = false;
std::string SaveMaskTemplate("mask-%n.tif"); //< src::default-mask-template mask-%n.tif
bool LoadMasks = false;
std::string LoadMaskTemplate(SaveMaskTemplate);
std::string VisualizeTemplate("vis-%n.tif"); //< src::default-visualize-template vis-%n.tif
bool VisualizeSeam = false;
std::pair<double, double> OptimizerWeights =
    std::make_pair(8.0,      //< src::default-optimizer-weight-distance 8.0
                   1.0);     //< src::default-optimizer-weight-mismatch 1.0
anneal_para_t AnnealPara = {
    32,                         //< src::default-anneal-kmax 32
    0.75,                       //< src::default-anneal-tau 0.75
    7000.0,                     //< src::default-anneal-deltae-max 7000.0
    5.0                         //< src::default-anneal-deltae-min 5.0
};
unsigned int DijkstraRadius = 25U; //< src::default-dijkstra-radius 25
struct AlternativePercentage MaskVectorizeDistance = {0.0, false};
std::string OutputCompression;
std::string OutputPixelType;
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

Signature sig;

#include "common.h"
#include "enblend.h"
#ifdef HAVE_LIBGLEW
#include "gpu.h"
#endif

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

using enblend::enblendMain;

#ifdef _WIN32
#define strdup _strdup
#endif


void inspectGPU(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA);

    const int handle = glutCreateWindow("Enblend");

    if (handle >= 1 && glutGet(GLUT_DISPLAY_MODE_POSSIBLE)) {
        cout <<
            "  - " << GLGETSTRING(GL_VENDOR) << "\n" <<
            "  - " << GLGETSTRING(GL_RENDERER) << "\n" <<
            "  - version " << GLGETSTRING(GL_VERSION) << "\n"
            "  - extensions\n";

        const char* const extensions = GLGETSTRING(GL_EXTENSIONS);
        const char* const extensions_end = extensions + strlen(extensions);
        const unsigned extensions_per_line = 3U;
        unsigned count = 1U;

        cout << "    ";
        for (const char* c = extensions; c != extensions_end; ++c) {
            if (*c == ' ') {
                if (count % extensions_per_line == 0U) {
                    cout << "\n    ";
                } else {
                    cout << "  ";
                }
                ++count;
            } else {
                cout << *c;
            }
        }
        cout << "\n\n";
    } else {
        cout << "    <no reliable OpenGL information available>\n";
    }

    glutDestroyWindow(handle);
}


/** Print information on the current version and some configuration
 * details. */
void printVersionAndExit(int argc, char** argv) {
    cout << "enblend " << VERSION << "\n\n";

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
#ifdef CACHE_IMAGES
            "yes" <<
#else
            "no" <<
#endif
            "\n";

#ifdef HAVE_LIBGLEW
        cout << "Extra feature: GPU acceleration: yes\n";
        inspectGPU(argc, argv);
#else
        cout << "Extra feature: GPU acceleration: no\n";
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
        cout << sig.message() << "\n\n";
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
        "Usage: enblend [options] [--output=IMAGE] INPUT...\n" <<
        "Blend INPUT images into a single IMAGE.\n" <<
        "\n" <<
        "Common options:\n" <<
        "  -V, --version          output version information and exit\n" <<
        "  -a                     pre-assemble non-overlapping images\n" <<
        "  -h, --help             print this help message and exit\n" <<
        "  -l LEVELS              number of blending LEVELS to use (1 to 29)\n" <<
        "  -o, --output=FILE      write output to FILE; default: \"" << OutputFileName << "\"\n" <<
        "  -v, --verbose[=LEVEL]  verbosely report progress; repeat to\n" <<
        "                         increase verbosity or directly set to LEVEL\n" <<
        "  -w, --wrap[=MODE]      wrap around image boundary, where MODE is\n" <<
        "                         NONE, HORIZONTAL, VERTICAL, or BOTH; default: " <<
        enblend::stringOfWraparound(WrapAround) << ";\n" <<
        "                         without argument the option selects horizontal wrapping\n" <<
        "  -x                     checkpoint partial results\n" <<
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
        "  --gpu                  use graphics card to accelerate seam-line optimization\n" <<
        "  -f WIDTHxHEIGHT[+xXOFFSET+yYOFFSET]\n" <<
        "                         manually set the size and position of the output\n" <<
        "                         image; useful for cropped and shifted input\n" <<
        "                         TIFF images, such as those produced by Nona\n" <<
        "  -m CACHESIZE           set image CACHESIZE in megabytes; default: " <<
        (CachedFileImageDirector::v().getAllocation() / 1048576LL) << "MB\n" <<
        "  --visualize[=TEMPLATE] save results of optimizer in TEMPLATE; same template\n" <<
        "                         characters as \"--save-mask\"; default: \"" <<
        VisualizeTemplate << "\"\n" <<
        "\n" <<
        "Mask generation options:\n" <<
        "  --coarse-mask[=FACTOR] shrink overlap regions by FACTOR to speedup mask\n" <<
        "                         generation; this is the default; if omitted FACTOR\n" <<
        "                         defaults to " <<
        CoarsenessFactor << "\n" <<
        "  --fine-mask            generate mask at full image resolution; use e.g.\n" <<
        "                         if overlap regions are very narrow\n" <<
        "  --smooth-difference=RADIUS\n" <<
        "                         smooth the difference image prior to seam-line\n" <<
        "                         optimization with a Gaussian blur of RADIUS;\n" <<
        "                         default: " << DifferenceBlurRadius << " pixels\n" <<
        "  --optimize             turn on mask optimization; this is the default\n" <<
        "  --no-optimize          turn off mask optimization\n" <<
        "  --optimizer-weights=DISTANCEWEIGHT[:MISMATCHWEIGHT]\n" <<
        "                         set the optimizer's weigths for distance and mismatch;\n" <<
        "                         default: " << OptimizerWeights.first << ':' <<
        OptimizerWeights.second << "\n" <<
        "  --mask-vectorize=LENGTH\n" <<
        "                         set LENGTH of single seam segment; append \"%\" for\n" <<
        "                         relative value; defaults: " <<
        coarseMaskVectorizeDistance << " for coarse masks and\n" <<
        "                         " <<
        fineMaskVectorizeDistance << " for fine masks\n" <<
        "  --anneal=TAU[:DELTAEMAX[:DELTAEMIN[:KMAX]]]]\n" <<
        "                         set annealing parameters of strategy 1; defaults:\n" <<
        "                         " << AnnealPara.tau << ':' <<
        AnnealPara.deltaEMax << ':' << AnnealPara.deltaEMin << ':' << AnnealPara.kmax << "\n" <<
        "  --dijkstra=RADIUS      set search RADIUS of strategy 2; default: " <<
        DijkstraRadius << " pixels\n" <<
        "  --save-mask[=TEMPLATE] save generated masks in TEMPLATE; default: \"" <<
        SaveMaskTemplate << "\"\n" <<
        "                         conversion chars: %i: mask index, mask %n: number,\n" <<
        "                         %p: full path, %d: dirname, %b: basename,\n" <<
        "                         %f: filename, %e: extension; lowercase characters\n" <<
        "                         refer to input images uppercase to output image\n" <<
        "  --load-mask[=TEMPLATE] use existing masks in TEMPLATE instead of generating\n" <<
        "                         them; same template characters as \"--save-mask\";\n" <<
        "                         default: \"" << LoadMaskTemplate << "\"\n" <<
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


enum AllPossibleOptions {
    VersionOption, PreAssembleOption /* -a */, HelpOption, LevelsOption,
    OutputOption, VerboseOption, WrapAroundOption /* -w */,
    CheckpointOption /* -x */, CompressionOption, LZWCompressionOption,
    BlockSizeOption, CIECAM02Option /* -c */,
    DepthOption, AssociatedAlphaOption /* -g */, GPUOption,
    SizeAndPositionOption /* -f */, CacheSizeOption,
    VisualizeOption, CoarseMaskOption, FineMaskOption,
    OptimizeOption, NoOptimizeOption,
    SaveMaskOption, LoadMaskOption,
    AnnealOption, DijkstraRadiusOption, MaskVectorizeDistanceOption,
    SmoothDifferenceOption, OptimizerWeightsOption,
    // currently below the radar...
    SequentialBlendingOption
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

    if (!OptimizeMask) {
        if (contains(optionSet, VisualizeOption)) {
            cerr << command <<
                ": warning: option \"--visualize\" without mask optimization\n" <<
                command <<
                ": warning:     has no effect" <<
                endl;
        }

        if (contains(optionSet, AnnealOption)) {
            cerr << command <<
                ": warning: option \"--anneal\" without mask optimization has\n" <<
                command <<
                ": warning:     no effect" <<
                endl;
        }

        if (contains(optionSet, DijkstraRadiusOption)) {
            cerr << command <<
                ": warning: option \"--dijkstra\" without mask optimization\n" <<
                command <<
                ": warning:     has no effect" <<
                endl;
        }

        if (contains(optionSet, SmoothDifferenceOption)) {
            cerr << command <<
                ": warning: option \"--smooth-difference\" without mask optimization\n" <<
                command <<
                ": warning:     has no effect" <<
                endl;
        }

        if (contains(optionSet, OptimizerWeightsOption)) {
            cerr << command <<
                ": warning: option \"--optimizer-weights\" without mask optimization\n" <<
                command <<
                ": warning:     has no effect" <<
                endl;
        }
    }

    if (!(OptimizeMask || CoarseMask) && contains(optionSet, MaskVectorizeDistanceOption)){
        cerr << command <<
            ": warning: option \"--mask-vectorize\" without mask optimization\n" <<
            command <<
            ": warning:     or coarse mask has no effect" <<
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
        AnnealId,                   //  8
        DijkstraRadiusId,           //  9
        MaskVectorizeDistanceId,    // 10
        CompressionId,              // 11
        VerboseId,                  // 12
        HelpId,                     // 13
        VersionId,                  // 14
        DepthId,                    // 15
        OutputId,                   // 16
        WrapAroundId,               // 17
        SmoothDifferenceId,         // 18
        OptimizerWeightsId          // 19
    };

    // NOTE: See note attached to "enum OptionId" above.
    static struct option long_options[] = {
        {"gpu", no_argument, 0, NoArgument},                                   //  0
        {"coarse-mask", optional_argument, 0, IntegerArgument},                //  1
        {"fine-mask", no_argument, 0, NoArgument},                             //  2
        {"optimize", no_argument, 0, NoArgument},                              //  3
        {"no-optimize", no_argument, 0, NoArgument},                           //  4
        {"save-mask", optional_argument, 0, StringArgument},                   //  5
        {"load-mask", optional_argument, 0, StringArgument},                   //  6
        {"visualize", optional_argument, 0, StringArgument},                   //  7
        {"anneal", required_argument, 0, StringArgument},                      //  8
        {"dijkstra", required_argument, 0, IntegerArgument},                   //  9
        {"mask-vectorize", required_argument, 0, StringArgument},              // 10
        {"compression", required_argument, 0, StringArgument},                 // 11
        {"verbose", optional_argument, 0, IntegerArgument},                    // 12
        {"help", no_argument, 0, NoArgument},                                  // 13
        {"version", no_argument, 0, NoArgument},                               // 14
        {"depth", required_argument, 0, StringArgument},                       // 15
        {"output", required_argument, 0, StringArgument},                      // 16
        {"wrap", optional_argument, 0, StringArgument},                        // 17
        {"smooth-difference", required_argument, 0, StringArgument},           // 18
        {"optimizer-weights", required_argument, 0, StringArgument},           // 19
        {0, 0, 0, 0}
    };

    bool justPrintVersion = false;
    bool justPrintUsage = false;
    OptionSetType optionSet;

    // Parse command line.
    int option_index = 0;
    int c;
    opterr = 0;       // we have our own "unrecognized option" message
    while ((c = getopt_long(argc, argv, "Vab:cd:f:ghl:m:o:sv::w::xz",
                            long_options, &option_index)) != -1) {
        switch (c) {
        case NoArgument: {
            if (long_options[option_index].flag != 0) {
                break;
            }
            switch (option_index) {
            case UseGpuId:
                UseGPU = true;
                optionSet.insert(GPUOption);
                break;
            case FineMaskId:
                CoarseMask = false;
                optionSet.insert(FineMaskOption);
                break;
            case OptimizeMaskId:
                OptimizeMask = true;
                optionSet.insert(OptimizeOption);
                break;
            case NoOptimizeMaskId:
                OptimizeMask = false;
                optionSet.insert(NoOptimizeOption);
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
            case SaveMaskId:
                if (optarg != NULL && *optarg != 0) {
                    SaveMaskTemplate = optarg;
                }
                SaveMasks = true;
                optionSet.insert(SaveMaskOption);
                break;
            case LoadMaskId:
                if (optarg != NULL && *optarg != 0) {
                    LoadMaskTemplate = optarg;
                }
                LoadMasks = true;
                optionSet.insert(LoadMaskOption);
                break;
            case VisualizeId:
                if (optarg != NULL && *optarg != 0) {
                    VisualizeTemplate = optarg;
                }
                VisualizeSeam = true;
                optionSet.insert(VisualizeOption);
                break;
            case CompressionId:
                OutputCompression = optarg;
                optionSet.insert(CompressionOption);
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
            case AnnealId: {
                boost::scoped_ptr<char> s(new char[strlen(optarg) + 1]);
                strcpy(s.get(), optarg);
                char* save_ptr = NULL;
                char* token = enblend::strtoken_r(s.get(), NUMERIC_OPTION_DELIMITERS, &save_ptr);
                char* tail;

                if (token != NULL && *token != 0) {
                    errno = 0;
                    double tau = strtod(token, &tail);
                    if (errno != 0) {
                        cerr << command
                             << ": option \"--anneal\": illegal numeric format \""
                             << token << "\" of tau: " << enblend::errorMessage(errno)
                             << endl;
                        exit(1);
                    }
                    if (*tail != 0) {
                        if (*tail == '%') {
                            tau /= 100.0;
                        } else {
                            cerr << command
                                 << ": --anneal: trailing garbage \""
                                 << tail << "\" in tau: \"" << token << "\""
                                 << endl;
                            exit(1);
                        }
                    }
                    //< src::minimum-anneal-tau 0
                    if (tau <= 0.0) {
                        cerr << command
                             << ": option \"--anneal\": tau must be larger than zero"
                             << endl;
                        exit(1);
                    }
                    //< src::maximum-anneal-tau 1
                    if (tau >= 1.0) {
                        cerr << command
                             << ": option \"--anneal\": tau must be less than one"
                             << endl;
                        exit(1);
                    }
                    AnnealPara.tau = tau;
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    AnnealPara.deltaEMax = strtod(token, &tail);
                    if (errno != 0) {
                        cerr << command << ": option \"--anneal\": illegal numeric format \""
                             << token << "\" of deltaE_max: " << enblend::errorMessage(errno)
                             << endl;
                        exit(1);
                    }
                    if (*tail != 0) {
                        cerr << command
                             << ": option \"--anneal\": trailing garbage \""
                             << tail << "\" in deltaE_max: \""
                             << token << "\"" << endl;
                        exit(1);
                    }
                    //< src::minimum-anneal-deltae-max 0
                    if (AnnealPara.deltaEMax <= 0.0) {
                        cerr << command
                             << ": option \"--anneal\": deltaE_max must be larger than zero"
                             << endl;
                        exit(1);
                    }
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    AnnealPara.deltaEMin = strtod(token, &tail);
                    if (errno != 0) {
                        cerr << command
                             << ": option \"--anneal\": illegal numeric format \""
                             << token << "\" of deltaE_min: " << enblend::errorMessage(errno)
                             << endl;
                        exit(1);
                    }
                    if (*tail != 0) {
                        cerr << command
                             << ": option \"--anneal\": trailing garbage \""
                             << tail << "\" in deltaE_min: \""
                             << token << "\"" << endl;
                        exit(1);
                    }
                    //< src::minimum-anneal-deltae-min 0
                    if (AnnealPara.deltaEMin <= 0.0) {
                        cerr << command
                             << ": option \"--anneal\": deltaE_min must be larger than zero"
                             << endl;
                        exit(1);
                    }
                }
                if (AnnealPara.deltaEMin >= AnnealPara.deltaEMax) {
                    cerr << command
                         << ": option \"--anneal\": deltaE_min must be less than deltaE_max"
                         << endl;
                    exit(1);
                }

                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    errno = 0;
                    const long int kmax = strtol(token, &tail, 10);
                    if (errno != 0) {
                        cerr << command
                             << ": option \"--anneal\": illegal numeric format \""
                             << token << "\" of k_max: " << enblend::errorMessage(errno)
                             << endl;
                        exit(1);
                    }
                    if (*tail != 0) {
                        cerr << command
                             << ": option \"--anneal\": trailing garbage \""
                             << tail << "\" in k_max: \""
                             << token << "\"" << endl;
                        exit(1);
                    }
                    //< src::minimum-anneal-kmax 3
                    if (kmax < 3L) {
                        cerr << command
                             << ": option \"--anneal\": k_max must larger or equal to 3"
                             << endl;
                        exit(1);
                    }
                    AnnealPara.kmax = static_cast<unsigned int>(kmax);
                }

                optionSet.insert(AnnealOption);
                break;
            }

            case MaskVectorizeDistanceId: {
                char* tail;
                MaskVectorizeDistance.isPercentage = false;
                errno = 0;
                MaskVectorizeDistance.value = strtod(optarg, &tail);
                if (errno != 0) {
                    cerr << command
                         << ": option \"--mask-vectorize\": illegal numeric format \""
                         << optarg << "\": " << enblend::errorMessage(errno)
                         << endl;
                    exit(1);
                }
                if (*tail != 0) {
                    if (*tail == '%') {
                        MaskVectorizeDistance.isPercentage = true;
                    } else {
                        cerr << command
                             << ": option \"--mask-vectorize\": trailing garbage \""
                             << tail << "\" in \"" << optarg << "\"" << endl;
                        exit(1);
                    }
                }
                if (MaskVectorizeDistance.value <= 0.0) {
                    cerr << command
                         << ": option \"--mask-vectorize\": distance must be positive"
                         << endl;
                    exit(1);
                }

                optionSet.insert(MaskVectorizeDistanceOption);
                break;
            }

            case SmoothDifferenceId: {
                char* tail;
                errno = 0;
                const double radius = strtod(optarg, &tail);
                if (errno != 0) {
                    cerr << command
                         << ": option \"--smooth-difference\": illegal numeric format \""
                         << optarg << "\": " << enblend::errorMessage(errno)
                         << endl;
                    exit(1);
                }
                if (*tail != 0) {
                    cerr << command
                         << ": option \"--smooth-difference\": trailing garbage \""
                         << tail << "\" in \"" << optarg << "\"" << endl;
                    exit(1);
                }
                //< src::minimum-smooth-difference 0.0
                if (radius < 0.0) {
                    cerr << command
                         << ": option \"--smooth-difference\": negative radius; will not blur"
                         << endl;
                    DifferenceBlurRadius = 0.0;
                } else {
                    DifferenceBlurRadius = radius;
                }

                optionSet.insert(SmoothDifferenceOption);
                break;
            }

            case OptimizerWeightsId: {
                boost::scoped_ptr<char> s(new char[strlen(optarg) + 1]);
                strcpy(s.get(), optarg);
                char* save_ptr = NULL;
                char* token = enblend::strtoken_r(s.get(), NUMERIC_OPTION_DELIMITERS, &save_ptr);
                OptimizerWeights.first =
                    enblend::numberOfString(token,
                                            _1 >= 0.0,
                                            "negative optimizer weight; will use 0.0",
                                            0.0);
                token = enblend::strtoken_r(NULL, NUMERIC_OPTION_DELIMITERS, &save_ptr);
                if (token != NULL && *token != 0) {
                    OptimizerWeights.second =
                        enblend::numberOfString(token,
                                                _1 >= 0.0,
                                                "negative optimizer weight; will use 0.0",
                                                0.0);
                }
                if (OptimizerWeights.first == 0.0 && OptimizerWeights.second == 0.0) {
                    cerr << command
                         << ": optimizer weights cannot be both zero"
                         << endl;
                }
                break;
            }

            default:
                cerr << command
                     << ": internal error: unhandled \"StringArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case StringArgument"

        // case FloatArgument: {
        // ...
        // } // end of "case FloatArgument"

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
            case CoarseMaskId:
                CoarseMask = true;
                if (optarg != NULL && *optarg != 0) {
                    CoarsenessFactor =
                        enblend::numberOfString(optarg,
                                                _1 >= 1U,
                                                "coarseness factor less or equal to 0; will use 1",
                                                1U);
                }
                optionSet.insert(CoarseMaskOption);
                break;
            case DijkstraRadiusId:
                //< src::minimum-dijkstra-radius 1
                DijkstraRadius =
                    enblend::numberOfString(optarg,
                                            _1 >= 1U,
                                            "Dijkstra radius is 0; will use 1",
                                            1U);
                optionSet.insert(DijkstraRadiusOption);
                break;
            default:
                cerr << command
                     << ": internal error: unhandled \"IntegerArgument\" option"
                     << endl;
                exit(1);
            }
            break;
        } // end of "case IntegerArgument"

        case 'V':
            justPrintVersion = true;
            optionSet.insert(VersionOption);
            break;
        case 'a':
            OneAtATime = false;
            optionSet.insert(PreAssembleOption);
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
                     << "Try \"enblend --help\" for more information." << endl;
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
        case 'l':
            // We take care of "too many levels" in "bounds.h".
            //< src::minimum-pyramid-levels 1
            ExactLevels =
                enblend::numberOfString(optarg,
                                        _1 >= 1U,
                                        "too few levels; will use one level",
                                        1U);
            optionSet.insert(LevelsOption);
            break;
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
        case 's':
            // Deprecated sequential blending flag.
            OneAtATime = true;
            cerr << command << ": warning: flag \"-s\" is deprecated." << endl;
            optionSet.insert(SequentialBlendingOption);
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
        case 'x':
            Checkpoint = true;
            optionSet.insert(CheckpointOption);
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
            cerr << "Try \"enblend --help\" for more information." << endl;
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
        printVersionAndExit(argc, argv);
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

    sig.initialize();

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

    if (UseGPU) {
#ifdef HAVE_LIBGLEW
        initGPU(&argc, argv);
#else
        cerr << command
             << ": warning: no GPU support compiled in; option \"--gpu\" has no effect"
             << endl;
#endif
    }

    sig.check();

    //if (CachedFileImageDirector::v()->getManagedBlocks() < 4) {
    //    // Max simultaneous image access is in:
    //    // 4 in any of many calls to combineThreeImages
    //    // 4 gaussian pyramid init (src image layer, src alpha layer, dest pyramid image layer 0, dest pyramid alpha layer 0)
    //    // 4 in reduce (src image layer N, src alpha layer N, dest image layer N+1, dest alpha layer N+1)
    //    // FIXME complain or automatically adjust blocksize to get ManagedBlocks above 4?
    //}

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
    int minDim = INT_MAX;
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
                 << ": input image \""
                 << *inputFileNameIterator
                 << "\" does not have an alpha channel"
                 << endl;
            exit(1);
        }

        // Get input image's position and size.
        Rect2D imageROI(Point2D(inputInfo->getPosition()),
                        Size2D(inputInfo->width(), inputInfo->height()));

        if (inputFileNameIterator == inputFileNameList.begin()) {
            // First input image
            minDim = std::min(inputInfo->width(), inputInfo->height());
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
            if (inputInfo->width() < minDim) {
                minDim = inputInfo->width();
            }
            if (inputInfo->height() < minDim) {
                minDim = inputInfo->height();
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
             << ": warning: Enblend needs two or more overlapping input images in order to do\n"
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

    // Switch to fine mask, if the smallest coarse mask would be less
    // than 64 pixels wide or high.
    if (minDim / 8 < 64 && CoarseMask) {
        cerr << command
             << ": warning: input images to small for coarse mask; switching to fine mask"
             << endl;
        CoarseMask = false;
    }

    if (MaskVectorizeDistance.value == 0) {
        MaskVectorizeDistance.isPercentage = false;
        MaskVectorizeDistance.value =
            CoarseMask ? coarseMaskVectorizeDistance : fineMaskVectorizeDistance;
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
            if      (pixelType == "UINT8")  enblendMain<RGBValue<UInt8 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "INT8")   enblendMain<RGBValue<Int8  > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT16") enblendMain<RGBValue<UInt16> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enblendMain<RGBValue<Int16 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enblendMain<RGBValue<UInt32> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enblendMain<RGBValue<Int32 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enblendMain<RGBValue<float > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enblendMain<RGBValue<double> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                cerr << command << ": RGB images with pixel type \""
                     << pixelType
                     << "\" are not supported"
                     << endl;
                exit(1);
            }
        } else {
            if      (pixelType == "UINT8")  enblendMain<UInt8 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "INT8")   enblendMain<Int8  >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT16") enblendMain<UInt16>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enblendMain<Int16 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enblendMain<UInt32>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enblendMain<Int32 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enblendMain<float >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enblendMain<double>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
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

#ifdef HAVE_LIBGLEW
    if (UseGPU) {
        wrapupGPU();
    }
#endif

    // Success.
    return 0;
}
