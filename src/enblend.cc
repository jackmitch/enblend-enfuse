/*
 * Copyright (C) 2004-2010 Andrew Mihal
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
#define isnan _isnan
#endif // _MSC_VER

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
#include "tiff_message.h"

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
std::string OutputFileName(DEFAULT_OUTPUT_FILENAME);
int Verbose = 1;                //< src::default-verbosity-level 1
int ExactLevels = 0;
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
#include "filespec.h"
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

using enblend::FileNameList;
using enblend::TraceableFileNameList;
using enblend::enblendMain;

#ifdef _WIN32
#define strdup _strdup
#endif


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
        "+ Checkpoint = " << enblend::stringOfBool(Checkpoint) << ", option \"-x\"\n" <<
        "+ UseGPU = " << enblend::stringOfBool(UseGPU) << ", option \"--gpu\"\n" <<
        "+ OptimizeMask = " << enblend::stringOfBool(OptimizeMask) <<
        ", options \"--optimize\" and \"--no-optimize\"\n" <<
        "+ CoarseMask = " << enblend::stringOfBool(CoarseMask) <<
        ", options \"--coarse-mask\" and \"--fine-mask\"\n" <<
        "+     CoarsenessFactor = " << CoarsenessFactor << ", argument to option \"--coarse-mask\"\n" <<
        "+ DifferenceBlurRadius = " << DifferenceBlurRadius << ", option \"--smooth-difference\"\n" <<
        "+ SaveMasks = " << enblend::stringOfBool(SaveMasks) << ", option \"--save-masks\"\n" <<
        "+     SaveMaskTemplate = <" << SaveMaskTemplate << ">, argument to option \"--save-masks\"\n" <<
        "+ LoadMasks = " << enblend::stringOfBool(LoadMasks) << ", option \"--load-masks\"\n" <<
        "+     LoadMaskTemplate = <" << LoadMaskTemplate << ">, argument to option \"--load-masks\"\n" <<
        "+ VisualizeSeam = " << enblend::stringOfBool(VisualizeSeam) << ", option \"--visualize\"\n" <<
        "+     VisualizeTemplate = <" << VisualizeTemplate << ">, argument to option \"--visualize\"\n" <<
        "+ OptimizerWeights = {\n" <<
        "+     distance = " << OptimizerWeights.first << ",\n" <<
        "+     mismatch = " << OptimizerWeights.second << "\n" <<
        "+ }, arguments to option \"--visualize\"\n"
        "+ AnnealPara = {\n" <<
        "+     kmax = " << AnnealPara.kmax << ",\n" <<
        "+     tau = " << AnnealPara.tau << ",\n" <<
        "+     deltaEMax = " << AnnealPara.deltaEMax << ",\n" <<
        "+     deltaEMin = " << AnnealPara.deltaEMin << "\n" <<
        "+ }, arguments to option \"--anneal\"\n" <<
        "+ DijkstraRadius = " << DijkstraRadius << ", option \"--dijkstra\"\n" <<
        "+ MaskVectorizeDistance = {\n" <<
        "+     value = " << MaskVectorizeDistance.value << ",\n" <<
        "+     isPercentage = " << enblend::stringOfBool(MaskVectorizeDistance.isPercentage) << "\n" <<
        "+ }, arguments to option \"--mask-vectorize\"\n" <<
        "+ OutputCompression = <" << OutputCompression << ">, option \"--compression\"\n" <<
        "+ OutputPixelType = <" << OutputPixelType << ">, option \"--depth\"\n" <<
        "+ end of global variable dump\n";
}


#ifdef HAVE_LIBGLEW
void inspectGPU(int argc, char** argv)
{
#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
    cgl_init();
#else
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA);

    const int handle = glutCreateWindow("Enblend");

    if (!(handle >= 1 && glutGet(GLUT_DISPLAY_MODE_POSSIBLE))) {
        cout << "    <no reliable OpenGL information available>\n";
        return;
    }
#endif

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

#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
    CGLDestroyContext(ctx);
#else
    glutDestroyWindow(handle);
#endif
}
#endif // HAVE_LIBGLEW


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

        cout << "Supported following globbing algorithms:\n";
        const enblend::algorithm_list algos = enblend::known_globbing_algorithms();
        for (enblend::algorithm_list::const_iterator i = algos.begin(); i != algos.end(); ++i) {
            cout <<
                "  " << i->first << "\n" <<
                "    " << i->second << "\n";
        }
        cout << "\n";
    }

    if (Verbose >= VERBOSE_SIGNATURE_REPORTING) {
        cout.flush();
        std::wcout << sig.message() << L"\n\n";
        std::wcout.flush();
    }

    cout <<
        "Copyright (C) 2004-2010 Andrew Mihal.\n" <<
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
        "Usage: enblend [options] [--output=IMAGE] INPUT...\n" <<
        "Blend INPUT images into a single IMAGE.\n" <<
        "\n" <<
        "INPUT... are image filenames or response filenames.  Response\n" <<
        "filenames start with an \"" << RESPONSE_FILE_PREFIX_CHAR << "\" character.\n"
        "\n" <<
        "Common options:\n" <<
        "  -V, --version          output version information and exit\n" <<
        "  -a                     pre-assemble non-overlapping images\n" <<
        "  -h, --help             print this help message and exit\n" <<
        "  -l, --levels=LEVELS    number of blending LEVELS to use (1 to " << MAX_PYRAMID_LEVELS << ");\n" <<
        "                         negative number of LEVELS decreases maximum\n" <<
        "  -o, --output=FILE      write output to FILE; default: \"" << OutputFileName << "\"\n" <<
        "  -v, --verbose[=LEVEL]  verbosely report progress; repeat to\n" <<
        "                         increase verbosity or directly set to LEVEL\n" <<
        "  -w, --wrap[=MODE]      wrap around image boundary, where MODE is \"none\",\n" <<
        "                         \"horizontal\", \"vertical\", or \"both\"; default: " <<
        enblend::stringOfWraparound(WrapAround) << ";\n" <<
        "                         without argument the option selects horizontal wrapping\n" <<
        "  -x                     checkpoint partial results\n" <<
        "  --compression=COMPRESSION\n" <<
        "                         set compression of output image to COMPRESSION,\n" <<
        "                         where COMPRESSION is:\n" <<
        "                         \"none\", \"packbits\", \"lzw\", \"deflate\" for TIFF files and\n" <<
        "                         0 to 100 for JPEG files\n" <<
        "\n" <<
        "Extended options:\n" <<
        "  -b BLOCKSIZE           image cache BLOCKSIZE in kilobytes; default: " <<
        (CachedFileImageDirector::v().getBlockSize() / 1024LL) << "KB\n" <<
        "  -c                     use CIECAM02 to blend colors\n" <<
        "  -d, --depth=DEPTH      set the number of bits per channel of the output\n" <<
        "                         image, where DEPTH is \"8\", \"16\", \"32\", \"r32\", or \"r64\"\n" <<
        "  -g                     associated-alpha hack for Gimp (before version 2)\n" <<
        "                         and Cinepaint\n" <<
        "  --gpu                  use graphics card to accelerate seam-line optimization\n" <<
        "  -f WIDTHxHEIGHT[+xXOFFSET+yYOFFSET]\n" <<
        "                         manually set the size and position of the output\n" <<
        "                         image; useful for cropped and shifted input\n" <<
        "                         TIFF images, such as those produced by Nona\n" <<
        "  -m CACHESIZE           set image CACHESIZE in megabytes; default: " <<
        (CachedFileImageDirector::v().getAllocation() / 1048576LL) << "MB\n" <<
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
        "  --anneal=TAU[:DELTAEMAX[:DELTAEMIN[:KMAX]]]\n" <<
        "                         set annealing parameters of optimizer strategy 1;\n" <<
        "                         defaults: " << AnnealPara.tau << ':' <<
        AnnealPara.deltaEMax << ':' << AnnealPara.deltaEMin << ':' << AnnealPara.kmax << "\n" <<
        "  --dijkstra=RADIUS      set search RADIUS of optimizer strategy 2; default:\n" <<
        "                         " << DijkstraRadius << " pixels\n" <<
        "  --save-masks[=TEMPLATE]\n" <<
        "                         save generated masks in TEMPLATE; default: \"" <<
        SaveMaskTemplate << "\";\n" <<
        "                         conversion chars: %i: mask index, %n: mask number,\n" <<
        "                         %p: full path, %d: dirname, %b: basename,\n" <<
        "                         %f: filename, %e: extension; lowercase characters\n" <<
        "                         refer to input images uppercase to the output image\n" <<
        "  --load-masks[=TEMPLATE]\n" <<
        "                         use existing masks in TEMPLATE instead of generating\n" <<
        "                         them; same template characters as \"--save-masks\";\n" <<
        "                         default: \"" << LoadMaskTemplate << "\"\n" <<
        "  --visualize[=TEMPLATE] save results of optimizer in TEMPLATE; same template\n" <<
        "                         characters as \"--save-masks\"; default: \"" <<
        VisualizeTemplate << "\"\n" <<
        "\n" <<
        "Enblend accepts arguments to any option in uppercase as\n" <<
        "well as in lowercase letters.\n" <<
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

#ifdef HAVE_LIBGLEW
    if (UseGPU) {
        // FIXME what if this occurs in a GL atomic section?
        wrapupGPU();
    }
#endif

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


int process_options(int argc, char** argv)
{
    enum OptionId {
        OPTION_ID_OFFSET = 1023,    // Ids start at 1024
        UseGpuId,
        CoarseMaskId,
        FineMaskId,
        OptimizeMaskId,
        NoOptimizeMaskId,
        SaveMaskId,
        LoadMaskId,
        VisualizeId,
        AnnealId,
        DijkstraRadiusId,
        MaskVectorizeDistanceId,
        CompressionId,
        VerboseId,
        HelpId,
        VersionId,
        DepthId,
        OutputId,
        WrapAroundId,
        SmoothDifferenceId,
        OptimizerWeightsId,
        LevelsId
    };

    static struct option long_options[] = {
        {"gpu", no_argument, 0, UseGpuId},
        {"coarse-mask", optional_argument, 0, CoarseMaskId},
        {"fine-mask", no_argument, 0, FineMaskId},
        {"optimize", no_argument, 0, OptimizeMaskId},
        {"no-optimize", no_argument, 0, NoOptimizeMaskId},
        {"save-mask", optional_argument, 0, SaveMaskId}, // singular form: not documented, not deprecated
        {"save-masks", optional_argument, 0, SaveMaskId},
        {"load-mask", optional_argument, 0, LoadMaskId}, // singular form: not documented, not deprecated
        {"load-masks", optional_argument, 0, LoadMaskId},
        {"visualize", optional_argument, 0, VisualizeId},
        {"anneal", required_argument, 0, AnnealId},
        {"dijkstra", required_argument, 0, DijkstraRadiusId},
        {"mask-vectorize", required_argument, 0, MaskVectorizeDistanceId},
        {"compression", required_argument, 0, CompressionId},
        {"verbose", optional_argument, 0, VerboseId},
        {"help", no_argument, 0, HelpId},
        {"version", no_argument, 0, VersionId},
        {"depth", required_argument, 0, DepthId},
        {"output", required_argument, 0, OutputId},
        {"wrap", optional_argument, 0, WrapAroundId},
        {"smooth-difference", required_argument, 0, SmoothDifferenceId},
        {"optimizer-weights", required_argument, 0, OptimizerWeightsId},
        {"levels", required_argument, 0, LevelsId},
        {0, 0, 0, 0}
    };

    bool failed = false;
    bool justPrintVersion = false;
    bool justPrintUsage = false;
    OptionSetType optionSet;

    opterr = 0;       // we have our own "unrecognized option" message
    while (true) {
        int option_index;
        const int code = getopt_long(argc, argv, "Vab:cd:f:ghl:m:o:sv::w::x",
                                     long_options, &option_index);

        if (code == -1) {
            break;
        }

        switch (code) {
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
                    failed = true;
                }
                if (*tail != 0) {
                    if (*tail == '%') {
                        tau /= 100.0;
                    } else {
                        cerr << command
                             << ": --anneal: trailing garbage \""
                             << tail << "\" in tau: \"" << token << "\""
                             << endl;
                        failed = true;
                    }
                }
                //< src::minimum-anneal-tau 0
                if (tau <= 0.0) {
                    cerr << command
                         << ": option \"--anneal\": tau must be larger than zero"
                         << endl;
                    failed = true;
                }
                //< src::maximum-anneal-tau 1
                if (tau >= 1.0) {
                    cerr << command
                         << ": option \"--anneal\": tau must be less than one"
                         << endl;
                    failed = true;
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
                    failed = true;
                }
                if (*tail != 0) {
                    cerr << command
                         << ": option \"--anneal\": trailing garbage \""
                         << tail << "\" in deltaE_max: \""
                         << token << "\"" << endl;
                    failed = true;
                }
                //< src::minimum-anneal-deltae-max 0
                if (AnnealPara.deltaEMax <= 0.0) {
                    cerr << command
                         << ": option \"--anneal\": deltaE_max must be larger than zero"
                         << endl;
                    failed = true;
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
                    failed = true;
                }
                if (*tail != 0) {
                    cerr << command
                         << ": option \"--anneal\": trailing garbage \""
                         << tail << "\" in deltaE_min: \""
                         << token << "\"" << endl;
                    failed = true;
                }
                //< src::minimum-anneal-deltae-min 0
                if (AnnealPara.deltaEMin <= 0.0) {
                    cerr << command
                         << ": option \"--anneal\": deltaE_min must be larger than zero"
                         << endl;
                    failed = true;
                }
            }
            if (AnnealPara.deltaEMin >= AnnealPara.deltaEMax) {
                cerr << command
                     << ": option \"--anneal\": deltaE_min must be less than deltaE_max"
                     << endl;
                failed = true;
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
                    failed = true;
                }
                if (*tail != 0) {
                    cerr << command
                         << ": option \"--anneal\": trailing garbage \""
                         << tail << "\" in k_max: \""
                         << token << "\"" << endl;
                    failed = true;
                }
                //< src::minimum-anneal-kmax 3
                if (kmax < 3L) {
                    cerr << command
                         << ": option \"--anneal\": k_max must larger or equal to 3"
                         << endl;
                    failed = true;
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
                failed = true;
            }
            if (*tail != 0) {
                if (*tail == '%') {
                    MaskVectorizeDistance.isPercentage = true;
                } else {
                    cerr << command
                         << ": option \"--mask-vectorize\": trailing garbage \""
                         << tail << "\" in \"" << optarg << "\"" << endl;
                    failed = true;
                }
            }
            if (MaskVectorizeDistance.value <= 0.0) {
                cerr << command
                     << ": option \"--mask-vectorize\": distance must be positive"
                     << endl;
                failed = true;
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
                failed = true;
            }
            if (*tail != 0) {
                cerr << command
                     << ": option \"--smooth-difference\": trailing garbage \""
                     << tail << "\" in \"" << optarg << "\"" << endl;
                failed = true;
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
            optionSet.insert(OptimizerWeightsOption);
            break;
        }

        case 'v': // FALLTHROUGH
        case VerboseId:
            if (optarg != NULL && *optarg != 0) {
                Verbose =
                    enblend::numberOfString(optarg,
                                            _1 >= 0, //< src::minimum-verbosity-level 0
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
            DijkstraRadius =
                enblend::numberOfString(optarg,
                                        _1 >= 1U, //< src::minimum-dijkstra-radius 1
                                        "Dijkstra radius is 0; will use 1",
                                        1U);
            optionSet.insert(DijkstraRadiusOption);
            break;

        case 'a':
            OneAtATime = false;
            optionSet.insert(PreAssembleOption);
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

        case 'l': // FALLTHROUGH
        case LevelsId:
            if (optarg != NULL && *optarg != 0) {
                std::ostringstream oss;
                oss <<
                    "cannot use more than " << MAX_PYRAMID_LEVELS <<
                    " pyramid levels; will use at most " << MAX_PYRAMID_LEVELS << " levels";
                ExactLevels =
                    enblend::numberOfString(optarg,
                                            _1 != 0,
                                            "cannot blend with zero levels; will use at least one level",
                                            1,
                                            _1 <= MAX_PYRAMID_LEVELS,
                                            oss.str(),
                                            MAX_PYRAMID_LEVELS);
            } else {
                cerr << command << ": options \"-l\" or \"--levels\" require an argument" << endl;
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

        case 's':
            // Deprecated sequential blending flag.
            OneAtATime = true;
            cerr << command << ": warning: flag \"-s\" is deprecated." << endl;
            optionSet.insert(SequentialBlendingOption);
            break;

        case 'x':
            Checkpoint = true;
            optionSet.insert(CheckpointOption);
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
            cerr << "Try \"enblend --help\" for more information." << endl;
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

    if (atexit(cleanup_output) != 0) {
        cerr << command << ": warning: could not install cleanup routine\n";
    }

    sig.initialize();

    TIFFSetWarningHandler(tiff_warning);
    TIFFSetErrorHandler(tiff_error);

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

    TraceableFileNameList inputTraceableFileNameList;

    // Remaining parameters are input files.
    while (optind < argc) {
        TraceableFileNameList files;
        enblend::unfold_filename(files, std::string(argv[optind]));
        inputTraceableFileNameList.insert(inputTraceableFileNameList.end(),
                                          files.begin(), files.end());
        optind++;
    }

    if (inputTraceableFileNameList.empty()) {
        cerr << command << ": no input files specified\n";
        exit(1);
    }

#ifdef DEBUG_DUMP_GLOBAL_VARIABLES
    DUMP_GLOBAL_VARIABLES();
#endif

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
    int minDim = INT_MAX;
    unsigned layer = 0;
    unsigned layers = 0;
    FileNameList inputFileNameList;
    TraceableFileNameList::iterator inputFileNameIterator = inputTraceableFileNameList.begin();
    while (inputFileNameIterator != inputTraceableFileNameList.end()) {
        ImageImportInfo* inputInfo = NULL;
        std::string filename(inputFileNameIterator->first);
        if (!enblend::can_open_file(filename)) {
            enblend::unroll_trace(inputFileNameIterator->second);
            exit(1);
        }
        try {
            ImageImportInfo info(filename.c_str());
            if (layers == 0) { // OPTIMIZATION: call only once per file
                layers = info.numLayers();
            }
            if (layers >= 2) {
                filename = vigra::join_filename_layer(filename, layer);
                inputInfo = new ImageImportInfo(filename.c_str());
            } else {
                inputInfo = new ImageImportInfo(info);
            }
            ++layer;
        } catch (vigra::ContractViolation& exception) {
            cerr <<
                command << ": cannot load image \"" << filename << "\"\n" <<
                command << ":     " << exception.message() << "\n";
            if (enblend::maybe_response_file(filename)) {
                cerr <<
                    command << ": info: Maybe you meant a response file and forgot the initial '" <<
                    RESPONSE_FILE_PREFIX_CHAR << "'?\n";
            }
            exit(1);
        }

        // Save this image info in the list.
        imageInfoList.push_back(inputInfo);
        inputFileNameList.push_back(filename);

        if (Verbose >= VERBOSE_INPUT_IMAGE_INFO_MESSAGES) {
            cerr << command
                 << ": info: input image \""
                 << inputFileNameIterator->first
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
                 << ": input image \"" << inputFileNameIterator->first << "\""
                 << enblend::optional_layer_name(layer, layers)
                 << " does not have an alpha channel"
                 << endl;
            enblend::unroll_trace(inputFileNameIterator->second);
            exit(1);
        }

        // Get input image's position and size.
        Rect2D imageROI(Point2D(inputInfo->getPosition()),
                        Size2D(inputInfo->width(), inputInfo->height()));

        if (inputFileNameIterator == inputTraceableFileNameList.begin()) {
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
                         << inputFileNameIterator->first
                         << "\"" << enblend::optional_layer_name(layer, layers) << endl;
                    enblend::unroll_trace(inputFileNameIterator->second);
                    exit(1);
                }
            }
        } else {
            // Second and later images
            inputUnion |= imageROI;

            if (isColor != inputInfo->isColor()) {
                cerr << command << ": input image \""
                     << inputFileNameIterator->first << "\""
                     << enblend::optional_layer_name(layer, layers) << " is "
                     << (inputInfo->isColor() ? "color" : "grayscale") << "\n"
                     << command << ":   but previous images are "
                     << (isColor ? "color" : "grayscale")
                     << endl;
                enblend::unroll_trace(inputFileNameIterator->second);
                exit(1);
            }
            if (pixelType != inputInfo->getPixelType()) {
                cerr << command << ": input image \""
                     << inputFileNameIterator->first << "\""
                     << enblend::optional_layer_name(layer, layers) << " has pixel type "
                     << inputInfo->getPixelType() << ",\n"
                     << command << ":   but previous images have pixel type "
                     << pixelType
                     << endl;
                enblend::unroll_trace(inputFileNameIterator->second);
                exit(1);
            }
            if (resolution !=
                TiffResolution(inputInfo->getXResolution(), inputInfo->getYResolution())) {
                cerr << command << ": info: input image \""
                     << inputFileNameIterator->first << "\""
                     << enblend::optional_layer_name(layer, layers) << " has resolution "
                     << inputInfo->getXResolution() << " dpi x "
                     << inputInfo->getYResolution() << " dpi,\n"
                     << command << ": info:   but first image has resolution "
                     << resolution.x << " dpi x " << resolution.y << " dpi"
                     << endl;
                enblend::unroll_trace(inputFileNameIterator->second);
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
                             << inputFileNameIterator->first
                             << "\"" << enblend::optional_layer_name(layer, layers) << endl;
                        enblend::unroll_trace(inputFileNameIterator->second);
                        exit(1);
                    }
                }

                cerr << endl << command << ": input image \""
                     << inputFileNameIterator->first
                     << "\"" << enblend::optional_layer_name(layer, layers) << "\n";
                enblend::unroll_trace(inputFileNameIterator->second);
                cerr << command << ":     has ";
                if (newProfile) {
                    cerr << "ICC profile \""
                         << cmsTakeProductName(newProfile)
                         << " "
                         << cmsTakeProductDesc(newProfile)
                         << "\"";
                } else {
                    cerr << "no ICC profile";
                }
                cerr << ", but previous images have ";
                if (InputProfile) {
                    cerr << "ICC profile \""
                         << cmsTakeProductName(InputProfile)
                         << " "
                         << cmsTakeProductDesc(InputProfile)
                         << "\"" << endl;
                } else {
                    cerr << "no ICC profile" << endl;
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
            // inputTraceableFileNameList still lacks the filename.
            inputTraceableFileNameList.insert(inputFileNameIterator, *inputFileNameIterator);
        }
    }

    vigra_postcondition(imageInfoList.size() == inputTraceableFileNameList.size(),
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
             << inputTraceableFileNameList.begin()->first << "\";\n"
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
