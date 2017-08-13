/*
 * Copyright (C) 2004-2009 Andrew Mihal
 * Copyright (C) 2009-2017 Christoph Spiel
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

#ifndef ENBLEND_SOURCE
#error trying to compile enblend.cc without having defined ENBLEND_SOURCE
#endif
#ifdef ENFUSE_SOURCE
#error must not define ENFUSE_SOURCE when compiling enblend.cc
#endif

#ifdef _WIN32
// Make sure we bring in windows.h the right way
#define _STLP_VERBOSE_AUTO_LINK
#define _USE_MATH_DEFINES
#define NOMINMAX
#define VC_EXTRALEAN
#include <windows.h>
#undef DIFFERENCE
#endif  // _WIN32

#ifdef __GW32C__
#undef malloc
#define BOOST_NO_STDC_NAMESPACE 1
#endif

#include <algorithm>
#include <iostream>
#include <list>
#include <memory>               // std::unique_ptr
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

#include <regex>

#ifdef HAVE_EXIV2
#include <exiv2/image.hpp>
#endif

#include <lcms2.h>
#if !defined(LCMS_VERSION) || LCMS_VERSION < 2050
#error "Little CMS version 2.5 or later is required"
#endif

#include "alternativepercentage.h"
#include "global.h"
#include "layer_selection.h"
#include "optional_transitional.hpp"
#include "parameter.h"
#include "selector.h"
#include "self_test.h"
#include "signature.h"
#include "tiff_message.h"
#ifdef _MSC_VER
#include "win32helpers/delayHelper.h"
#endif

typedef enum {
    UnknownDifference,
    HueLuminanceMaxDifference,  // maximum of hue difference and luminance difference
    DeltaEDifference            // L*a*b*-based Delta E
} difference_functor_t;


typedef struct {
    unsigned int kmax;          // maximum number of moves for a line segment
    double tau;                 // temperature reduction factor, "cooling factor"; 0 < tau < 1
    double deltaEMax;           // maximum cost change possible by any single annealing move
    double deltaEMin;           // minimum cost change possible by any single annealing move
} anneal_para_t;


// Globals
const std::string command("enblend");
const int minimumVectorizeDistance = 4; //< minimum-vectorize-distance 4
const int coarseMaskVectorizeDistance = 4; //< coarse-mask-vectorize-distance 4
const int fineMaskVectorizeDistance = 20; //< fine-mask-vectorize-distance 20

// Global values from command line parameters.
std::string OutputFileName(DEFAULT_OUTPUT_FILENAME);
std::optional<std::string> OutputMaskFileName;
int Verbose = 0;                //< default-verbosity-level 0
int ExactLevels = 0;            // 0 means: automatically calculate maximum
bool OneAtATime = true;
boundary_t WrapAround = OpenBoundaries;
bool GimpAssociatedAlphaHack = false;
blend_colorspace_t BlendColorspace = UndeterminedColorspace;
bool OutputSizeGiven = false;
int OutputWidthCmdLine = 0;
int OutputHeightCmdLine = 0;
int OutputOffsetXCmdLine = 0;
int OutputOffsetYCmdLine = 0;
MainAlgo MainAlgorithm = GraphCut;
bool Checkpoint = false;
bool OptimizeMask = true;
bool CoarseMask = true;
unsigned CoarsenessFactor = 8U; //< default-coarseness-factor 8
difference_functor_t PixelDifferenceFunctor = DeltaEDifference; //< default-difference-functor Delta-E
double LuminanceDifferenceWeight = 1.0; //< default-luminance-difference-weight 1.0
double ChrominanceDifferenceWeight = 1.0; //< default-chrominance-difference-weight 1.0
bool SaveMasks = false;
bool StopAfterMaskGeneration = false;
std::string SaveMaskTemplate("mask-%n.tif"); //< default-mask-template mask-%n.tif
bool LoadMasks = false;
std::string LoadMaskTemplate(SaveMaskTemplate);
std::string VisualizeTemplate("vis-%n.tif"); //< default-visualize-template vis-%n.tif
bool VisualizeSeam = false;
std::pair<double, double> OptimizerWeights =
    std::make_pair(12.0,        //< default-optimizer-weight-distance 12.0
                   1.0);        //< default-optimizer-weight-mismatch 1.0
anneal_para_t AnnealPara = {
    32,                         //< default-anneal-kmax 32
    0.75,                       //< default-anneal-tau 0.75
    7000.0,                     //< default-anneal-deltae-max 7000.0
    5.0                         //< default-anneal-deltae-min 5.0
};
unsigned int DijkstraRadius = 25U; //< default-dijkstra-radius 25
AlternativePercentage MaskVectorizeDistance(0.0, false);
std::string OutputCompression;
std::string OutputPixelType;

TiffResolution ImageResolution;
bool OutputIsValid = true;

bool UseGPU = false;
namespace cl {class Context;}
cl::Context* GPUContext = nullptr;
namespace ocl {class BatchBuilder;}
ocl::BatchBuilder* BatchCompiler = nullptr;

// Globals related to catching SIGINT
#ifndef _WIN32
sigset_t SigintMask;
#endif

// Objects for ICC profiles
cmsHPROFILE InputProfile = nullptr;
cmsHPROFILE XYZProfile = nullptr;
cmsHPROFILE LabProfile = nullptr;
cmsHTRANSFORM InputToXYZTransform = nullptr;
cmsHTRANSFORM XYZToInputTransform = nullptr;
cmsHTRANSFORM InputToLabTransform = nullptr;
cmsHTRANSFORM LabToInputTransform = nullptr;
cmsViewingConditions ViewingConditions;
cmsHANDLE CIECAMTransform = nullptr;
cmsHPROFILE FallbackProfile = nullptr;

Signature sig;
LayerSelectionHost LayerSelection;

#include <vigra/imageinfo.hxx>
#include <vigra/impex.hxx>
#include <vigra/sized_int.hxx>

#include <tiffio.h>

#include "common.h"
#include "filespec.h"
#include "introspection.h"
#include "enblend.h"

#ifdef DMALLOC
#include "dmalloc.h"            // must be last #include
#endif

#ifdef _WIN32
#define strdup _strdup
#endif


#ifdef OPENCL
namespace GPU {
    std::unique_ptr<ocl::CalculateStateProbabilities> StateProbabilities = nullptr;
    std::unique_ptr<vigra::ocl::DistanceTransformFH> DistanceTransform = nullptr;
}
#endif


difference_functor_t
differenceFunctorOfString(const std::string& aDifferenceFunctorName)
{
    std::string name(aDifferenceFunctorName);

    name.erase(std::remove_if(name.begin(), name.end(), [](char c) {return c == '-';}), name.end());
    enblend::to_lower(name);

    if (name == "maximumhueluminance" || name == "maximumhuelum" ||
        name == "maxhueluminance" || name == "maxhuelum" || name == "max") {
        return HueLuminanceMaxDifference;
    } else if (name == "deltae" || name == "de") {
        return DeltaEDifference;
    } else {
        return UnknownDifference;
    }
}


#define DUMP_GLOBAL_VARIABLES(...) dump_global_variables(__FILE__, __LINE__, ##__VA_ARGS__)
void dump_global_variables(const char* file, unsigned line,
                           std::ostream& out = std::cout)
{
    out <<
        "+ " << file << ":" << line << ": state of global variables\n" <<
        "+ Verbose = " << Verbose << ", option \"--verbose\"\n" <<
        "+ OutputFileName = <" << OutputFileName << ">\n" <<
        "+ OutputMaskFileName = <" << OutputMaskFileName.value_or("<not defined>") << ">\n" <<
        "+ ExactLevels = " << ExactLevels << "\n" <<
        "+ UseGPU = " << UseGPU << "\n" <<
        "+ OneAtATime = " << enblend::stringOfBool(OneAtATime) << ", option \"-a\"\n" <<
        "+ WrapAround = " << enblend::stringOfWraparound(WrapAround) << ", option \"--wrap\"\n" <<
        "+ GimpAssociatedAlphaHack = " << enblend::stringOfBool(GimpAssociatedAlphaHack) <<
        ", option \"-g\"\n" <<
        "+ BlendColorspace = " << BlendColorspace << ", option \"--blend-colorspace\"\n" <<
        "+ FallbackProfile = " << (FallbackProfile ? enblend::profileDescription(FallbackProfile) : "[none]") <<
        ", option \"--fallback-profile\"\n" <<
        "+ OutputSizeGiven = " << enblend::stringOfBool(OutputSizeGiven) << ", option \"-f\"\n" <<
        "+     OutputWidthCmdLine = " << OutputWidthCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputHeightCmdLine = " << OutputHeightCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputOffsetXCmdLine = " << OutputOffsetXCmdLine << ", argument to option \"-f\"\n" <<
        "+     OutputOffsetYCmdLine = " << OutputOffsetYCmdLine << ", argument to option \"-f\"\n" <<
        "+ Checkpoint = " << enblend::stringOfBool(Checkpoint) << ", option \"-x\"\n" <<
        "+ OptimizeMask = " << enblend::stringOfBool(OptimizeMask) <<
        ", options \"--optimize\" and \"--no-optimize\"\n" <<
        "+ CoarseMask = " << enblend::stringOfBool(CoarseMask) <<
        ", options \"--coarse-mask\" and \"--fine-mask\"\n" <<
        "+     CoarsenessFactor = " << CoarsenessFactor << ", argument to option \"--coarse-mask\"\n" <<
        "+ PixelDifferenceFunctor = " << stringOfPixelDifferenceFunctor(PixelDifferenceFunctor) <<
        "+     LuminanceDifferenceWeight = " << LuminanceDifferenceWeight << "\n" <<
        "+     ChrominanceDifferenceWeight = " << ChrominanceDifferenceWeight <<
        ", option \"--image-difference\"\n" <<
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
        "+     value = " << MaskVectorizeDistance.value() << ",\n" <<
        "+     is_percentage = " << enblend::stringOfBool(MaskVectorizeDistance.is_percentage()) << "\n" <<
        "+ }, arguments to option \"--mask-vectorize\"\n" <<
        "+ OutputCompression = <" << OutputCompression << ">, option \"--compression\"\n" <<
        "+ OutputPixelType = <" << OutputPixelType << ">, option \"--depth\"\n" <<
        "+ end of global variable dump\n";
}


void
printUsage(const bool error = true)
{
    std::cout <<
        "Usage: enblend [options] [--output=IMAGE] INPUT...\n" <<
        "Blend INPUT images into a single IMAGE.\n" <<
        "\n" <<
        "INPUT... are image filenames or response filenames.  Response\n" <<
        "filenames start with an \"" << RESPONSE_FILE_PREFIX_CHAR << "\" character.\n"
        "\n" <<
        "Options:\n" <<
        "Common options:\n" <<
        "  -l, --levels=LEVELS    limit number of blending LEVELS to use (1 to " << MAX_PYRAMID_LEVELS << ");\n" <<
        "                         negative number of LEVELS decreases maximum;\n" <<
        "                         \"auto\" restores the default automatic maximization\n" <<
        "  -o, --output=FILE      write output to FILE; default: \"" << OutputFileName << "\"\n" <<
        "  -v, --verbose[=LEVEL]  verbosely report progress; repeat to\n" <<
        "                         increase verbosity or directly set to LEVEL\n" <<
        "  --compression=COMPRESSION\n" <<
        "                         set compression of output image to COMPRESSION,\n" <<
        "                         where COMPRESSION is:\n" <<
        "                         \"deflate\", \"jpeg\", \"lzw\", \"none\", \"packbits\", for TIFF files and\n" <<
        "                         0 to 100, or \"jpeg\", \"jpeg-arith\" for JPEG files,\n" <<
        "                         where \"jpeg\" and \"jpeg-arith\" accept a compression level\n" <<
#ifdef OPENCL
        "  --gpu                  employ GPU in addition to CPU for selected computations; negate\n" <<
        "                         with \"--no-gpu\"\n" <<
#endif
        "\n" <<
        "Advanced options:\n" <<
        "  --blend-colorspace=COLORSPACE\n" <<
        "                         force COLORSPACE for blending operations; Enblend uses\n" <<
        "                         \"CIELUV\" for images with ICC-profile and \"IDENTITY\" for\n" <<
        "                         those without and also for all floating-point images;\n" <<
        "                         other available blend color spaces are \"CIELAB\" and\n" <<
        "                         \"CIECAM\"\n" <<
        "  -d, --depth=DEPTH      set the number of bits per channel of the output\n" <<
        "                         image, where DEPTH is \"8\", \"16\", \"32\", \"r32\", or \"r64\"\n" <<
        "  -f WIDTHxHEIGHT[+xXOFFSET+yYOFFSET]\n" <<
        "                         manually set the size and position of the output\n" <<
        "                         image; useful for cropped and shifted input\n" <<
        "                         TIFF images, such as those produced by Nona\n" <<
        "  -g                     associated-alpha hack for Gimp (before version 2)\n" <<
        "                         and Cinepaint\n" <<
        "  --output-mask[=FILE]   write output mask to FILE if the output format does not\n" <<
        "                         support alpha-channels; default: \"" << DEFAULT_OUTPUT_MASK_FILENAME << "\"\n" <<
        "  -w, --wrap[=MODE]      wrap around image boundary, where MODE is \"none\",\n" <<
        "                         \"horizontal\", \"vertical\", or \"both\"; default: " <<
        enblend::stringOfWraparound(WrapAround) << ";\n" <<
        "                         without argument the option selects horizontal wrapping\n" <<
        "\n" <<
        "Mask generation options:\n" <<
        "  --coarse-mask[=FACTOR] shrink overlap regions by FACTOR to speedup mask\n" <<
        "                         generation; this is the default; if omitted FACTOR\n" <<
        "                         defaults to " << CoarsenessFactor << "\n" <<
        "  --fine-mask            generate mask at full image resolution; use e.g.\n" <<
        "                         if overlap regions are very narrow\n" <<
        "  --optimize             turn on mask optimization; this is the default;\n" <<
        "                         disable with \"--no-optimize\"\n" <<
        "  --save-masks[=TEMPLATE]\n" <<
        "                         save generated masks in TEMPLATE; default: \"" << SaveMaskTemplate << "\";\n" <<
        "                         conversion chars: \"%i\": mask index, \"%n\": mask number,\n" <<
        "                         \"%p\": full path, \"%d\": dirname, \"%b\": basename,\n" <<
        "                         \"%f\": filename, \"%e\": extension; lowercase characters\n" <<
        "                         refer to input images uppercase to the output image\n" <<
        "  --load-masks[=TEMPLATE]\n" <<
        "                         use existing masks in TEMPLATE instead of generating\n" <<
        "                         them; same template characters as \"--save-masks\";\n" <<
        "                         default: \"" << LoadMaskTemplate << "\"\n" <<
        "  --visualize[=TEMPLATE] save results of optimizer in TEMPLATE; same template\n" <<
        "                         characters as \"--save-masks\"; default: \"" << VisualizeTemplate << "\"\n" <<
        "\n" <<
        "Expert options:\n" <<
        "  -a, --pre-assemble     pre-assemble non-overlapping images; negate with \"--no-pre-assemble\"\n" <<
        "  -x                     checkpoint partial results\n" <<
        "  --fallback-profile=PROFILE-FILE\n" <<
        "                         use the ICC profile from PROFILE-FILE instead of sRGB\n" <<
        "  --layer-selector=ALGORITHM\n" <<
        "                         set the layer selector ALGORITHM;\n" <<
        "                         default: \"" << LayerSelection.name() << "\"; available algorithms are:\n";
    for (selector::algorithm_list::const_iterator i = selector::algorithms.begin();
         i != selector::algorithms.end();
         ++i) {
        std::cout << "                         \"" << (*i)->name() << "\": " << (*i)->description() << "\n";
    }
    std::cout <<
#ifdef OPENCL
        "  --prefer-gpu=DEVICE    select DEVICE on auto-detected platform as GPU\n" <<
        "  --prefer-gpu=PLATFORM:DEVICE\n" <<
        "                         select DEVICE on PLATFORM as GPU\n" <<
#endif
        "  --parameter=KEY1[=VALUE1][:KEY2[=VALUE2][:...]]\n" <<
        "                         set one or more KEY-VALUE pairs\n" <<
        "\n" <<
        "Expert mask generation options:\n" <<
        "  --primary-seam-generator=ALGORITHM\n" <<
        "                         use main seam finder ALGORITHM, where ALGORITHM is\n"<<
        "                         \"nearest-feature-transform\" or \"graph-cut\";\n" <<
        "                         default: \"graph-cut\"\n" <<
        "  --image-difference=ALGORITHM[:LUMINANCE-WEIGHT[:CHROMINANCE-WEIGHT]]\n" <<
        "                         use ALGORITHM for calculation of the difference image,\n" <<
        "                         where ALGORITHM is \"max-hue-luminance\" or \"delta-e\";\n" <<
        "                         LUMINANCE-WEIGHT and CHROMINANCE-WEIGHT define the\n" <<
        "                         weights of lightness and color; default: " <<
        stringOfPixelDifferenceFunctor(PixelDifferenceFunctor) << ":" << LuminanceDifferenceWeight <<
        ":" << ChrominanceDifferenceWeight << "\n" <<
        "  --optimizer-weights=DISTANCE-WEIGHT[:MISMATCH-WEIGHT]\n" <<
        "                         set the optimizer's weigths for distance and mismatch;\n" <<
        "                         default: " << OptimizerWeights.first << ':' <<
        OptimizerWeights.second << "\n" <<
        "  --mask-vectorize=LENGTH\n" <<
        "                         set LENGTH of single seam segment; append \"%\" for\n" <<
        "                         relative value; defaults: " <<
        coarseMaskVectorizeDistance << " for coarse masks and\n" <<
        "                         " <<
        fineMaskVectorizeDistance << " for fine masks\n" <<
        "  --anneal=TAU[:DELTAE-MAX[:DELTAE-MIN[:K-MAX]]]\n" <<
        "                         set annealing parameters of optimizer strategy 1;\n" <<
        "                         defaults: " << AnnealPara.tau << ':' <<
        AnnealPara.deltaEMax << ':' << AnnealPara.deltaEMin << ':' << AnnealPara.kmax << "\n" <<
        "  --dijkstra=RADIUS      set search RADIUS of optimizer strategy 2; default:\n" <<
        "                         " << DijkstraRadius << " pixels\n" <<
        "\n" <<
        "Information options:\n" <<
        "  -h, --help             print this help message and exit\n" <<
        "  -V, --version          output version information and exit\n" <<
        "  --show-globbing-algorithms\n" <<
        "                         show all globbing algorithms\n" <<
#ifdef OPENCL
        "  --show-gpu-info        list all available GPUs according to their platform and device;\n" <<
        "                         inform on current preferences\n" <<
#endif
        "  --show-image-formats   show all recognized image formats and their filename\n" <<
        "                         extensions\n" <<
        "  --show-signature       show who compiled the binary when and on which machine\n" <<
        "  --show-software-components\n" <<
        "                         show the software components with which Enblend was compiled\n" <<
        "\n" <<
        "Enblend accepts arguments to any option in uppercase as\n" <<
        "well as in lowercase letters.\n" <<
#if defined(OPENMP) || \
    (defined(OPENCL) && defined(PREFER_SEPARATE_OPENCL_SOURCE))
        "\n" <<
        "Environment:\n" <<
#endif
#ifdef OPENMP
        "  OMP_NUM_THREADS        The OMP_NUM_THREADS environment variable sets the number\n" <<
        "                         of threads to use in OpenMP parallel regions.  If unset\n" <<
        "                         Enblend uses as many threads as there are CPUs.\n" <<
        "  OMP_DYNAMIC            The OMP_DYNAMIC environment variable controls dynamic\n" <<
        "                         adjustment of the number of threads to use in executing\n" <<
        "                         OpenMP parallel regions.\n" <<
#endif
#if defined(OPENCL) && defined(PREFER_SEPARATE_OPENCL_SOURCE)
        "  ENBLEND_OPENCL_PATH    The ENBLEND_OPENCL_PATH environment variable sets the search\n" <<
        "                         path for OpenCL source files.\n" <<
#endif
        "\n" <<
        "Report bugs at <" PACKAGE_BUGREPORT ">." <<
        std::endl;

    exit(error ? 1 : 0);
}


void cleanup_output(void)
{
    if (!OutputIsValid) {
        std::cerr << command << ": info: remove invalid output image \"" << OutputFileName << "\"\n";
        errno = 0;
        if (unlink(OutputFileName.c_str()) != 0) {
            std::cerr << command <<
                ": warning: could not remove invalid output image \"" << OutputFileName << "\": " <<
                enblend::errorMessage(errno) << "\n";
        }
    }
}


void sigint_handler(int sig)
{
    std::cerr << std::endl << command << ": interrupted" << std::endl;

    cleanup_output();

#if !defined(__GW32C__) && !defined(_WIN32)
    struct sigaction action;
    action.sa_handler = SIG_DFL;
    action.sa_flags   = 0;
    sigemptyset(&(action.sa_mask));
    sigaction(sig, &action, nullptr);
#else
    signal(sig, SIG_DFL);
#endif
    raise(sig);
}


enum AllPossibleOptions {
    VersionOption, PreAssembleOption /* -a */, NoPreAssembleOption, HelpOption, LevelsOption,
    OutputOption, OutputMaskOption, VerboseOption, WrapAroundOption /* -w */,
    CheckpointOption /* -x */, CompressionOption, LZWCompressionOption,
    BlendColorspaceOption, FallbackProfileOption,
    DepthOption, AssociatedAlphaOption /* -g */,
    GPUOption, NoGPUOption, PreferredGPUOption,
    SizeAndPositionOption /* -f */,
    VisualizeOption, CoarseMaskOption, FineMaskOption,
    OptimizeOption, NoOptimizeOption,
    SaveMasksOption, LoadMasksOption,
    ImageDifferenceOption, AnnealOption, DijkstraRadiusOption, MaskVectorizeDistanceOption,
    OptimizerWeightsOption,
    LayerSelectorOption, NearestFeatureTransformOption, GraphCutOption,
    ShowImageFormatsOption, ShowSignatureOption, ShowGlobbingAlgoInfoOption, ShowSoftwareComponentsInfoOption,
    ShowGPUInfoOption,
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
    if (contains(optionSet, LoadMasksOption)) {
        if (contains(optionSet, VisualizeOption)) {
            std::cerr << command <<
                ": warning: option \"--visualize\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, CoarseMaskOption)) {
            std::cerr << command <<
                ": warning: option \"--coarse-mask\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, FineMaskOption)) {
            std::cerr << command <<
                ": warning: option \"--fine-mask\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, OptimizeOption)) {
            std::cerr << command <<
                ": warning: option \"--optimize\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, NoOptimizeOption)) {
            std::cerr << command <<
                ": warning: option \"--no-optimize\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, AnnealOption)) {
            std::cerr << command <<
                ": warning: option \"--anneal\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, DijkstraRadiusOption)) {
            std::cerr << command <<
                ": warning: option \"--dijkstra\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, MaskVectorizeDistanceOption)) {
            std::cerr << command <<
                ": warning: option \"--mask-vectorize\" has no effect with \"--load-masks\"" << std::endl;
        }
        if (contains(optionSet, OptimizerWeightsOption)) {
            std::cerr << command <<
                ": warning: option \"--optimizer-weights\" has no effect with \"--load-masks\"" << std::endl;
        }
    }

    if (contains(optionSet, SaveMasksOption) && !contains(optionSet, OutputOption)) {
        if (contains(optionSet, LevelsOption)) {
            std::cerr << command <<
                ": warning: option \"--levels\" has no effect with \"--save-masks\" and no \"--output\"" <<
                std::endl;
        }
        if (contains(optionSet, WrapAroundOption)) {
            std::cerr << command <<
                ": warning: option \"--wrap\" has no effect with \"--save-masks\" and no \"--output\"" <<
                std::endl;
        }
        if (contains(optionSet, CompressionOption)) {
            std::cerr << command <<
                ": warning: option \"--compression\" has no effect with \"--save-masks\" and no \"--output\"" <<
                std::endl;
        }
        if (contains(optionSet, DepthOption)) {
            std::cerr << command <<
                ": warning: option \"--depth\" has no effect with \"--save-masks\" and no \"--output\"" <<
                std::endl;
        }
        if (contains(optionSet, SizeAndPositionOption)) {
            std::cerr << command <<
                ": warning: option \"-f\" has no effect with \"--save-masks\" and no \"--output\"" << std::endl;
        }
    }

    if (contains(optionSet, CompressionOption) &&
        !(enblend::getFileType(OutputFileName) == "TIFF" ||
          enblend::getFileType(OutputFileName) == "JPEG")) {
        std::cerr << command <<
            ": warning: compression is not supported with output\n" <<
            command <<
            ": warning: file type \"" <<
            enblend::getFileType(OutputFileName) << "\"" <<
            std::endl;
    }

    if (contains(optionSet, AssociatedAlphaOption) &&
        enblend::getFileType(OutputFileName) != "TIFF") {
        std::cerr << command <<
            ": warning: option \"-g\" has no effect with output\n" <<
            command <<
            ": warning: file type \"" <<
            enblend::getFileType(OutputFileName) << "\"" <<
            std::endl;
    }

    if (!OptimizeMask) {

        if (contains(optionSet, AnnealOption)) {
            std::cerr << command <<
                ": warning: option \"--anneal\" without mask optimization has\n" <<
                command <<
                ": warning: no effect" <<
                std::endl;
        }

        if (contains(optionSet, DijkstraRadiusOption)) {
            std::cerr << command <<
                ": warning: option \"--dijkstra\" without mask optimization\n" <<
                command <<
                ": warning: has no effect" <<
                std::endl;
        }

        if (contains(optionSet, OptimizerWeightsOption)) {
            std::cerr << command <<
                ": warning: option \"--optimizer-weights\" without mask optimization\n" <<
                command <<
                ": warning: has no effect" <<
                std::endl;
        }
    }

    if (!(OptimizeMask || CoarseMask) && contains(optionSet, MaskVectorizeDistanceOption)) {
        std::cerr << command <<
            ": warning: option \"--mask-vectorize\" without mask optimization\n" <<
            command <<
            ": warning: or coarse mask has no effect" <<
            std::endl;
    }

    if (!CoarseMask && MainAlgorithm == GraphCut && contains(optionSet, FineMaskOption)) {
        std::cerr << command <<
            ": warning: option \"--fine-mask\" combined with option \"--primary-seam-generator=graphcut\"\n" <<
            command <<
            ": warning: incompatible with mask optimization,\n" <<
            command <<
            ": note: defaulting to no optimization" <<
            std::endl;
        OptimizeMask = false;
    }

#ifdef OPENCL
    if (!UseGPU && contains(optionSet, PreferredGPUOption)) {
        std::cerr << command << ": warning: option \"--preferred-gpu\" has no effect without enabled GPU" << std::endl;
    }
#else
    if (contains(optionSet, GPUOption)) {
        std::cerr << command << ": warning: option \"--gpu\" has no effect in this " << command << " binary,\n" <<
            command << ": warning: because " << command << " was compiled without support for OpenCL" << std::endl;
    }
    if (contains(optionSet, NoGPUOption)) {
        std::cerr << command << ": warning: option \"--no-gpu\" has no effect in this " << command << " binary,\n" <<
            command << ": warning: because " << command << " was compiled without support for OpenCL" << std::endl;
    }
    if (contains(optionSet, PreferredGPUOption)) {
        std::cerr << command << ": warning: option \"--prefer-gpu\" has no effect in this " << command << " binary,\n" <<
            command << ": warning: because " << command << " was compiled without support for OpenCL" << std::endl;
    }
#endif // OPENCL
}


const struct option*
lookup_long_option_by_id(struct option* an_option_array, int a_long_option_id)
{
    assert(an_option_array != nullptr);

    struct option* opt = an_option_array;

    while (opt->name != 0) {
        if (opt->val == a_long_option_id) {
            return opt;
        }
        ++opt;
    }

    return nullptr;
}


bool
short_option_requires_argument(const char* a_short_option_spec, char a_short_option)
{
    if (a_short_option != 0) {
        const char* c = strchr(a_short_option_spec, a_short_option);

        if (c != nullptr) {
            if (*(c + 1) == ':' && *(c + 2) != ':') {
                return true;
            }
        }
    }

    return false;
}


#ifdef OPENCL
void
initialize_gpu_subsystem(size_t a_preferred_gpu_platform, size_t a_preferred_gpu_device)
{
    try {
        cl::Platform platform = ocl::find_platform(a_preferred_gpu_platform);

        ocl::device_list_t devices;
        ocl::prefer_device(platform, a_preferred_gpu_platform, a_preferred_gpu_device, devices);

        GPUContext = ocl::create_context(platform, devices);

        if (Verbose >= VERBOSE_OPENCL_MESSAGES) {
            std::cerr << command << ": info: chose OpenCL platform #" << a_preferred_gpu_platform << ", ";

            std::string info;
            platform.getInfo(CL_PLATFORM_VENDOR, &info);
            std::cerr << info << ", ";
            platform.getInfo(CL_PLATFORM_NAME, &info);
            std::cerr << info << ", device #" << a_preferred_gpu_device << std::endl;
        }

#ifdef OPENMP
        BatchCompiler = new ocl::ThreadedBatchBuilder;
#else
        BatchCompiler = new ocl::SerialBatchBuilder;
#endif
    } catch (ocl::runtime_error& an_exception) {
        std::cerr <<
            command << ": warning: " << an_exception.what() << ";\n" <<
            command << ": warning: " << "    cannot enable GPU" << std::endl;
    }

    if (!gpu_is_ok(GPUContext)) {
        std::cerr <<
            command << ": warning: at least one of Enfuse's GPU tests has failed\n" <<
            command << ": note: refusing to activate GPU support" << std::endl;

        delete GPUContext;
        delete BatchCompiler;
        GPUContext = nullptr;
        BatchCompiler = nullptr;
    }
}
#endif // OPENCL


int
process_options(int argc, char** argv)
{
    enum OptionId {
        OPTION_ID_OFFSET = 1023,    // Ids start at 1024
        UseGpuId,
        NoUseGpuId,
        PreferGpuId,
        PreAssembleId,
        NoPreAssembleId,
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
        OutputMaskId,
        WrapAroundId,
        OptimizerWeightsId,
        LevelsId,
        BlendColorspaceId,
        FallbackProfileId,
        LayerSelectorId,
        MainAlgoId,
        ImageDifferenceId,
        ParameterId,
        NoParameterId,
        ImageFormatsInfoId,
        SignatureInfoId,
        GlobbingAlgoInfoId,
        SoftwareComponentsInfoId,
        GPUInfoId
    };

    static struct option long_options[] = {
        {"gpu", no_argument, 0, UseGpuId},
        {"no-gpu", no_argument, 0, NoUseGpuId},
        {"prefer-gpu", required_argument, 0, PreferGpuId},
        {"preferred-gpu", required_argument, 0, PreferGpuId}, // gramatically close alternative form
        {"pre-assemble", no_argument, 0, PreAssembleId},
        {"preassemble", no_argument, 0, PreAssembleId}, // dash-less form: not documented, not deprecated
        {"no-pre-assemble", no_argument, 0, NoPreAssembleId},
        {"no-preassemble", no_argument, 0, NoPreAssembleId}, // dash-less form: not documented, not deprecated
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
        {"output-mask", optional_argument, 0, OutputMaskId},
        {"wrap", optional_argument, 0, WrapAroundId},
        {"optimizer-weights", required_argument, 0, OptimizerWeightsId},
        {"levels", required_argument, 0, LevelsId},
        {"blend-colorspace", required_argument, 0, BlendColorspaceId},
        {"blend-color-space", required_argument, 0, BlendColorspaceId}, // dash form: not documented, not deprecated
        {"fallback-profile", required_argument, 0, FallbackProfileId},
        {"layer-selector", required_argument, 0, LayerSelectorId},
        {"primary-seam-generator", required_argument, 0, MainAlgoId},
        {"image-difference", required_argument, 0, ImageDifferenceId},
        {"parameter", required_argument, 0, ParameterId},
        {"no-parameter", required_argument, 0, NoParameterId},
        {"show-image-formats", no_argument, 0, ImageFormatsInfoId},
        {"show-signature", no_argument, 0, SignatureInfoId},
        {"show-globbing-algorithms", no_argument, 0, GlobbingAlgoInfoId},
        {"show-software-components", no_argument, 0, SoftwareComponentsInfoId},
        {"show-gpu-info", no_argument, 0, GPUInfoId},
        {0, 0, 0, 0}
    };

    bool failed = false;

    typedef enum {
        PROCESS_COMPLETELY,
        VERSION_ONLY, USAGE_ONLY,
        IMAGE_FORMATS_ONLY, SIGNATURE_ONLY, GLOBBING_ALGOS_ONLY, SOFTWARE_COMPONENTS_ONLY,
        GPU_INFO_COMPONENTS_ONLY
    } print_only_task_id_t;
    print_only_task_id_t print_only_task = PROCESS_COMPLETELY;

    OptionSetType optionSet;
#ifdef OPENCL
    size_t preferredGPUPlatform = 0U; // We start enumerating platforms at 1 for user convenience.
                                      // Zero means we choose the first platform found, i.e. auto-detect.
    size_t preferredGPUDevice = 1U;   // We start enumerating platforms at 1 for user convenience.
#endif

    static const char short_options[] = "Vab:cd:f:ghl:m:o:sv::w::x";
    opterr = 0;       // we have our own "unrecognized option" message

    while (true) {
        int option_index = -1;
        const int code = getopt_long(argc, argv, short_options, long_options, &option_index);

        if (code == -1) {
            break;
        }

        switch (code) {
        case UseGpuId:
            UseGPU = true;
            optionSet.insert(GPUOption);
            break;

        case NoUseGpuId:
            UseGPU = false;
            optionSet.insert(NoGPUOption);
            break;

        case PreferGpuId:
#ifdef OPENCL
        {
            char* delimiter = strpbrk(optarg, NUMERIC_OPTION_DELIMITERS);
            if (delimiter == nullptr) {
                preferredGPUDevice =
                    enblend::numberOfString(optarg, [](unsigned x) {return x >= 1U;},
                                            "preferred GPU device out of range", 1U);
            } else {
                *delimiter = 0;
                ++delimiter;
                preferredGPUPlatform =
                    enblend::numberOfString(optarg, [](unsigned x) {return x >= 1U;},
                                            "preferred GPU platform out of range", 1U);
                preferredGPUDevice =
                    enblend::numberOfString(delimiter, [](unsigned x) {return x >= 1U;},
                                            "preferred GPU device out of range", 1U);
            }
        }
#endif
            optionSet.insert(PreferredGPUOption);
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

        case 'h':
            [[fallthrough]];
        case HelpId:
            print_only_task = USAGE_ONLY;
            optionSet.insert(HelpOption);
            break;

        case 'V':
            [[fallthrough]];
        case VersionId:
            print_only_task = VERSION_ONLY;
            optionSet.insert(VersionOption);
            break;

        case ImageFormatsInfoId:
            print_only_task = IMAGE_FORMATS_ONLY;
            optionSet.insert(ShowImageFormatsOption);
            break;

        case SignatureInfoId:
            print_only_task = SIGNATURE_ONLY;
            optionSet.insert(ShowSignatureOption);
            break;

        case GlobbingAlgoInfoId:
            print_only_task = GLOBBING_ALGOS_ONLY;
            optionSet.insert(ShowGlobbingAlgoInfoOption);
            break;

        case SoftwareComponentsInfoId:
            print_only_task = SOFTWARE_COMPONENTS_ONLY;
            optionSet.insert(ShowSoftwareComponentsInfoOption);
            break;

        case GPUInfoId:
            print_only_task = GPU_INFO_COMPONENTS_ONLY;
            optionSet.insert(ShowGPUInfoOption);
            break;

        case 'w':
            [[fallthrough]];
        case WrapAroundId:
            if (optarg != nullptr && *optarg != 0) {
                WrapAround = enblend::wraparoundOfString(optarg);
                if (WrapAround == UnknownWrapAround) {
                    std::cerr << command
                              << ": unrecognized wrap-around mode \"" << optarg << "\"\n" << std::endl;
                    failed = true;
                }
            } else {
                WrapAround = HorizontalStrip;
            }
            optionSet.insert(WrapAroundOption);
            break;

        case SaveMaskId:
            if (optarg != nullptr && *optarg != 0) {
                SaveMaskTemplate = optarg;
            }
            SaveMasks = true;
            optionSet.insert(SaveMasksOption);
            break;

        case LoadMaskId:
            if (optarg != nullptr && *optarg != 0) {
                LoadMaskTemplate = optarg;
            }
            LoadMasks = true;
            optionSet.insert(LoadMasksOption);
            break;

        case VisualizeId:
            if (optarg != nullptr && *optarg != 0) {
                VisualizeTemplate = optarg;
            }
            VisualizeSeam = true;
            optionSet.insert(VisualizeOption);
            break;

        case CompressionId:
            if (optarg != nullptr && *optarg != 0) {
                std::string upper_opt(optarg);
                enblend::to_upper(upper_opt);
                if (upper_opt == "NONE" || upper_opt == "DEFLATE" || upper_opt == "LZW" || upper_opt == "PACKBITS") {
                    OutputCompression = upper_opt;
                } else if (upper_opt.find_first_not_of("0123456789") == std::string::npos) {
                    OutputCompression = "JPEG QUALITY=" + upper_opt;
                } else if (enblend::starts_with(upper_opt, "JPEG") || enblend::starts_with(upper_opt, "JPEG-ARITH")) {
                    const std::string::size_type delimiter_position = upper_opt.find_first_of(NUMERIC_OPTION_DELIMITERS);
                    if (delimiter_position == std::string::npos) {
                        if (upper_opt == "JPEG" || upper_opt == "JPEG-ARITH") {
                            OutputCompression = upper_opt;
                        } else {
                            std::cerr << command << ": trailing garbage in JPEG compression \"" << optarg << "\"" << std::endl;
                            failed = true;
                        }
                    } else {
                        const std::string algorithm(upper_opt.substr(0, delimiter_position));
                        if (algorithm == "JPEG" || algorithm == "JPEG-ARITH") {
                            const std::string level(upper_opt.substr(delimiter_position + 1U));
                            if (level.length() >= 1U && level.find_first_not_of("0123456789") == std::string::npos) {
                                upper_opt.replace(delimiter_position, 1U, " QUALITY=");
                                OutputCompression = upper_opt;
                            } else {
                                std::cerr << command << ": invalid JPEG compression level \"" << level << "\"" << std::endl;
                                failed = true;
                            }
                        } else {
                            std::cerr << command << ": unrecognized JPEG compression \"" << optarg << "\"" << std::endl;
                            failed = true;
                        }
                    }
                } else {
                    std::cerr << command << ": unrecognized compression \"" << optarg << "\"" << std::endl;
                    failed = true;
                }
            } else {
                std::cerr << command << ": option \"--compression\" requires an argument" << std::endl;
                failed = true;
            }
            optionSet.insert(CompressionOption);
            break;

        case 'd':
            [[fallthrough]];
        case DepthId:
            if (optarg != nullptr && *optarg != 0) {
                OutputPixelType = enblend::outputPixelTypeOfString(optarg);
            } else {
                std::cerr << command << ": options \"-d\" or \"--depth\" require arguments" << std::endl;
                failed = true;
            }
            optionSet.insert(DepthOption);
            break;

        case 'o':
            [[fallthrough]];
        case OutputId:
            if (contains(optionSet, OutputOption)) {
                std::cerr << command
                          << ": warning: more than one output file specified"
                          << std::endl;
            }
            if (optarg != nullptr && *optarg != 0) {
                OutputFileName = optarg;
            } else {
                std::cerr << command << ": options \"-o\" or \"--output\" require arguments" << std::endl;
                failed = true;
            }
            optionSet.insert(OutputOption);
            break;

        case OutputMaskId:
            OutputMaskFileName = std::string(optarg ? optarg : DEFAULT_OUTPUT_MASK_FILENAME);
            optionSet.insert(OutputMaskOption);
            break;

        case ImageDifferenceId: {
            std::string::size_type tail;
            const std::regex delimiterRegex(NUMERIC_OPTION_DELIMITERS_REGEX);
            const std::string arg(optarg);
            const std::sregex_token_iterator endOfSequence;
            std::sregex_token_iterator token(arg.begin(), arg.end(), delimiterRegex, -1);

            if (token == endOfSequence) {
                std::cerr << command << ": option \"--image-difference\" requires an argument" << std::endl;
                failed = true;
            } else {
                PixelDifferenceFunctor = differenceFunctorOfString(token->str());
                if (PixelDifferenceFunctor == UnknownDifference) {
                    std::cerr << command
                              << ": unknown image difference algorithm \"" << token->str() << "\"" << std::endl;
                    failed = true;
                }
                ++token;
            }

            if (token != endOfSequence) {
                try {
                    LuminanceDifferenceWeight = std::stod(token->str(), &tail);
                    if (tail != token->str().length()) {
                        std::cerr << command << ": unrecognized luminance weight \""
                                  << token->str().substr(tail) << "\" in \"" << token->str() << "\"" << std::endl;
                        failed = true;
                    }
                    if (LuminanceDifferenceWeight < 0.0) {
                        std::cerr << command << ": luminance weight must be non-negative" << std::endl;
                        failed = true;
                    }
                } catch (std::invalid_argument&) {
                    std::cerr << command << ": illegal numeric format \""
                              << token->str() << "\" of luminance weight" << std::endl;
                    failed = true;
                }
                ++token;
            }

            if (token != endOfSequence) {
                try {
                    ChrominanceDifferenceWeight = std::stod(token->str(), &tail);
                    if (tail != token->str().length()) {
                        std::cerr << command << ": unrecognized chrominance weight \""
                                  << token->str().substr(tail) << "\" in \"" << token->str() << "\"" << std::endl;
                        failed = true;
                    }
                    if (ChrominanceDifferenceWeight < 0.0) {
                        std::cerr << command << ": chrominance weight must be non-negative" << std::endl;
                        failed = true;
                    }
                } catch (std::invalid_argument&) {
                    std::cerr << command << ": illegal numeric format \""
                              << token->str() << "\" of chrominance weight" << std::endl;
                    failed = true;
                }
                ++token;
            }

            if (token != endOfSequence) {
                std::cerr << command << ": warning: ignoring trailing garbage \""
                          << token->str() << "\" in argument to \"--image-difference\"" << std::endl;
            }

            if (LuminanceDifferenceWeight + ChrominanceDifferenceWeight == 0.0) {
                std::cerr << command << ": luminance weight and chrominance weight cannot both be zero" << std::endl;
                failed = true;
            }

            optionSet.insert(ImageDifferenceOption);
            break;
        }

        case AnnealId: {
            std::string::size_type tail;
            const std::regex delimiterRegex(NUMERIC_OPTION_DELIMITERS_REGEX);
            const std::string arg(optarg);
            const std::sregex_token_iterator endOfSequence;
            std::sregex_token_iterator token(arg.begin(), arg.end(), delimiterRegex, -1);

            if (token != endOfSequence) {
                double tau = 0;
                try {
                    tau = std::stod(token->str(), &tail);
                    if (tail < token->str().length()) {
                        if (token->str().substr(tail, 1) == "%") {
                            tau /= 100.0;
                        } else {
                            std::cerr << command
                                << ": --anneal: trailing garbage \""
                                << token->str().substr(tail) << "\" in tau: \"" << token->str() << "\""
                                << std::endl;
                            failed = true;
                        }
                    }
                } catch (std::invalid_argument&) {
                    std::cerr << command
                        << ": option \"--anneal\": illegal numeric format \""
                        << token->str() << "\" of tau" << std::endl;
                    failed = true;
                }
                //< minimum-anneal-tau 0
                if (tau <= 0.0) {
                    std::cerr << command
                              << ": option \"--anneal\": tau must be larger than zero"
                              << std::endl;
                    failed = true;
                }
                //< maximum-anneal-tau 1
                if (tau >= 1.0) {
                    std::cerr << command
                              << ": option \"--anneal\": tau must be less than one"
                              << std::endl;
                    failed = true;
                }
                AnnealPara.tau = tau;
                ++token;
            }

            if (token != endOfSequence) {
                try {
                    AnnealPara.deltaEMax = std::stod(token->str(), &tail);
                    if (tail != token->str().length()) {
                        std::cerr << command
                            << ": option \"--anneal\": trailing garbage \""
                            << token->str().substr(tail) << "\" in deltaE_max: \""
                            << token->str() << "\"" << std::endl;
                        failed = true;
                    }
                } catch (std::invalid_argument&) {
                    std::cerr << command << ": option \"--anneal\": illegal numeric format \""
                              << token->str() << "\" of deltaE_max" << std::endl;
                    failed = true;
                }
                //< minimum-anneal-deltae-max 0
                if (AnnealPara.deltaEMax <= 0.0) {
                    std::cerr << command
                              << ": option \"--anneal\": deltaE_max must be larger than zero"
                              << std::endl;
                    failed = true;
                }
                ++token;
            }

            if (token != endOfSequence) {
                try {
                    AnnealPara.deltaEMin = std::stod(token->str(), &tail);
                    if (tail != token->str().length()) {
                        std::cerr << command
                            << ": option \"--anneal\": trailing garbage \""
                            << token->str().substr(tail) << "\" in deltaE_min: \""
                            << token->str() << "\"" << std::endl;
                        failed = true;
                    }
                } catch (std::invalid_argument&) {
                    std::cerr << command
                              << ": option \"--anneal\": illegal numeric format \""
                              << token->str() << "\" of deltaE_min" << std::endl;
                    failed = true;
                }
                //< minimum-anneal-deltae-min 0
                if (AnnealPara.deltaEMin <= 0.0) {
                    std::cerr << command
                              << ": option \"--anneal\": deltaE_min must be larger than zero"
                              << std::endl;
                    failed = true;
                }
                ++token;
            }
            if (AnnealPara.deltaEMin >= AnnealPara.deltaEMax) {
                std::cerr << command
                          << ": option \"--anneal\": deltaE_min must be less than deltaE_max"
                          << std::endl;
                failed = true;
            }

            if (token != endOfSequence) {
                long int kmax = 0;
                try {
                    kmax = std::stol(token->str(), &tail, 10);
                    if (tail != token->str().length()) {
                        std::cerr << command
                            << ": option \"--anneal\": trailing garbage \""
                            << token->str().substr(tail) << "\" in k_max: \""
                            << token->str() << "\"" << std::endl;
                        failed = true;
                    }

                } catch (std::invalid_argument&) {
                    std::cerr << command
                              << ": option \"--anneal\": illegal numeric format \""
                              << token->str() << "\" of k_max" << std::endl;
                    failed = true;
                }
                //< minimum-anneal-kmax 3
                if (kmax < 3L) {
                    std::cerr << command
                              << ": option \"--anneal\": k_max must larger or equal to 3"
                              << std::endl;
                    failed = true;
                }
                AnnealPara.kmax = static_cast<unsigned int>(kmax);
            }

            optionSet.insert(AnnealOption);
            break;
        }

        case MaskVectorizeDistanceId: {
            char* tail;
            MaskVectorizeDistance.set_percentage(false);
            errno = 0;
            MaskVectorizeDistance.set_value(strtod(optarg, &tail));
            if (errno != 0) {
                std::cerr << command
                          << ": option \"--mask-vectorize\": illegal numeric format \""
                          << optarg << "\": " << enblend::errorMessage(errno)
                          << std::endl;
                failed = true;
            }
            if (*tail != 0) {
                if (*tail == '%') {
                    MaskVectorizeDistance.set_percentage(true);
                } else {
                    std::cerr << command
                              << ": option \"--mask-vectorize\": trailing garbage \""
                              << tail << "\" in \"" << optarg << "\"" << std::endl;
                    failed = true;
                }
            }
            if (MaskVectorizeDistance.value() <= 0.0) {
                std::cerr << command
                          << ": option \"--mask-vectorize\": distance must be positive"
                          << std::endl;
                failed = true;
            }

            optionSet.insert(MaskVectorizeDistanceOption);
            break;
        }

        case OptimizerWeightsId: {
            const std::regex delimiterRegex(NUMERIC_OPTION_DELIMITERS_REGEX);
            const std::string arg(optarg);
            const std::sregex_token_iterator endOfSequence;
            std::sregex_token_iterator token(arg.begin(), arg.end(), delimiterRegex, -1);

            OptimizerWeights.first =
                enblend::numberOfString(token->str(),
                                        [](double x) {return x >= 0.0;},
                                        "negative optimizer weight; will use 0.0",
                                        0.0);
            ++token;
            if (token != endOfSequence) {
                OptimizerWeights.second =
                    enblend::numberOfString(token->str(),
                                            [](double x) {return x >= 0.0;},
                                            "negative optimizer weight; will use 0.0",
                                            0.0);
            }
            if (OptimizerWeights.first == 0.0 && OptimizerWeights.second == 0.0) {
                std::cerr << command
                          << ": optimizer weights cannot be both zero"
                          << std::endl;
            }
            optionSet.insert(OptimizerWeightsOption);
            break;
        }

        case 'v':
            [[fallthrough]];
        case VerboseId:
            if (optarg != nullptr && *optarg != 0) {
                Verbose = enblend::numberOfString(optarg,
                                                  [](int x) {return x >= 0;}, //< minimum-verbosity-level 0
                                                  "verbosity level less than 0; will use 0",
                                                  0);
            } else {
                Verbose++;
            }
            optionSet.insert(VerboseOption);
            break;

        case CoarseMaskId:
            CoarseMask = true;
            if (optarg != nullptr && *optarg != 0) {
                CoarsenessFactor =
                    enblend::numberOfString(optarg,
                                            [](unsigned x) {return x >= 1U;},
                                            "coarseness factor less or equal to 0; will use 1",
                                            1U);
            }
            optionSet.insert(CoarseMaskOption);
            break;

        case DijkstraRadiusId:
            DijkstraRadius =
                enblend::numberOfString(optarg,
                                        [](unsigned x) {return x >= 1U;}, //< minimum-dijkstra-radius 1
                                        "Dijkstra radius is 0; will use 1",
                                        1U);
            optionSet.insert(DijkstraRadiusOption);
            break;

        case 'a':
            [[fallthrough]];
        case PreAssembleId:
            OneAtATime = false;
            optionSet.insert(PreAssembleOption);
            break;

        case NoPreAssembleId:
            OneAtATime = true;
            optionSet.insert(NoPreAssembleOption);
            break;

        case BlendColorspaceId:
            if (optarg != nullptr && *optarg != 0) {
                std::string name(optarg);
                enblend::to_upper(name);
                if (name == "IDENTITY" || name == "ID" || name == "UNIT") {
                    BlendColorspace = IdentitySpace;
                } else if (name == "LAB" || name == "CIELAB" || name == "LSTAR" || name == "L-STAR") {
                    BlendColorspace = CIELAB;
                } else if (name == "LUV" || name == "CIELUV") {
                    BlendColorspace = CIELUV;
                } else if (name == "CIECAM" || name == "CIECAM02" || name == "JCH") {
                    BlendColorspace = CIECAM;
                } else {
                    std::cerr << command <<
                        ": unrecognized argument \"" << optarg << "\" of option \"--blend-colorspace\"" <<
                        std::endl;
                    failed = true;
                }
            } else {
                std::cerr << command << ": option \"--blend-colorspace\" requires an argument" << std::endl;
                failed = true;
            }
            optionSet.insert(BlendColorspaceOption);
            break;

        case FallbackProfileId:
            if (enblend::can_open_file(optarg)) {
                FallbackProfile = cmsOpenProfileFromFile(optarg, "r");
                if (FallbackProfile == nullptr) {
                    std::cerr << command << ": failed to open fallback ICC profile file \"" << optarg << "\"\n";
                    exit(1);
                }
            } else {
                exit(1);
            }
            optionSet.insert(FallbackProfileOption);
            break;

        case MainAlgoId:
            if (optarg != nullptr && *optarg != 0) {
                std::string algo_name(optarg);
                enblend::to_upper(algo_name);
                if (algo_name == "GRAPH-CUT" ||
                    algo_name == "GRAPHCUT" ||
                    algo_name == "GC") {
                    MainAlgorithm = GraphCut;
                    optionSet.insert(GraphCutOption);
                } else if (algo_name == "NEAREST-FEATURE-TRANSFORM" ||
                           algo_name == "NEARESTFEATURETRANSFORM" ||
                           algo_name == "NFT") {
                    MainAlgorithm = NFT;
                    optionSet.insert(NearestFeatureTransformOption);
                } else {
                    std::cerr << command <<
                        "unrecognized argument \"" << optarg << "\" of option \"--primary-seam-generator\"" <<
                        std::endl;
                    failed = true;
                }
            } else {
                std::cerr << command << ": option \"--primary-seam-generator\" requires an argument" <<
                    std::endl;
                failed = true;
            }
            break;

        case 'f':
            if (optarg != nullptr && *optarg != 0) {
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
                    std::cerr << command << ": option \"-f\" requires 2 or 4 arguments" << std::endl;
                    failed = true;
                }
            } else {
                std::cerr << command << ": option \"-f\" requires 2 or 4 arguments" << std::endl;
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
            [[fallthrough]];
        case LevelsId:
            if (optarg != nullptr && *optarg != 0) {
                std::string levels(optarg);
                enblend::to_upper(levels);
                if (levels == "AUTO" || levels == "AUTOMATIC") {
                    ExactLevels = 0;
                } else if (levels.find_first_not_of("+-0123456789") != std::string::npos) {
                    std::cerr << command <<
                        ": options \"-l\" or \"--levels\" require an integer argument or \"auto\"" << std::endl;
                    failed = true;
                } else {
                    std::ostringstream oss;
                    oss <<
                        "cannot use more than " << MAX_PYRAMID_LEVELS <<
                        " pyramid levels; will use at most " << MAX_PYRAMID_LEVELS << " levels";
                    ExactLevels =
                        enblend::numberOfString(optarg,
                                                 [](int x) {return x != 0;},
                                                "cannot blend with zero levels; will use one level",
                                                1,
                                                [](int x) {return x <= MAX_PYRAMID_LEVELS;},
                                                oss.str(),
                                                MAX_PYRAMID_LEVELS);
                }
            } else {
                std::cerr << command << ": options \"-l\" or \"--levels\" require an argument" << std::endl;
                failed = true;
            }
            optionSet.insert(LevelsOption);
            break;

        case 'x':
            Checkpoint = true;
            optionSet.insert(CheckpointOption);
            break;

        case LayerSelectorId: {
            selector::algorithm_list::const_iterator selector = selector::find_by_name(optarg);
            if (selector != selector::algorithms.end()) {
                LayerSelection.set_selector(selector->get());
            } else {
                std::cerr << command << ": unknown selector algorithm \"" << optarg << "\"";
                exit(1);
            }
            optionSet.insert(LayerSelectorOption);
            break;
        }

        case ParameterId: {
            const std::regex delimiterRegex(NUMERIC_OPTION_DELIMITERS_REGEX);
            const std::string arg(optarg);
            const std::sregex_token_iterator endOfSequence;
            std::sregex_token_iterator token(arg.begin(), arg.end(), delimiterRegex, -1);

            for (; token != endOfSequence; ++token) {
                std::string key;
                std::string value;
                size_t delimiter = token->str().find_first_of(ASSIGNMENT_CHARACTERS);

                if (delimiter == std::string::npos) {
                    key = token->str();
                    value.assign("1");
                } else {
                    key = token->str().substr(0, delimiter);
                    value = token->str().substr(delimiter + 1);
                    if (value.empty()) {
                        std::cerr <<
                            command << ": parameter key \"" << key << "\" lacks a value;\n" <<
                            command << ": note: dangling assignment operator\n";
                        exit(1);
                    }
                }
                enblend::trim(key);
                enblend::trim(value);

                if (parameter::is_valid_identifier(key)) {
                    parameter::insert(key, value);
                } else {
                    std::cerr << command << ": parameter key \"" << key << "\" is not a valid identifier\n";
                    exit(1);
                }
            }

            break;
        }

        case NoParameterId: {
            const std::regex delimiterRegex(NUMERIC_OPTION_DELIMITERS_REGEX);
            const std::string arg(optarg);
            const std::sregex_token_iterator endOfSequence;
            std::sregex_token_iterator token(arg.begin(), arg.end(), delimiterRegex, -1);

            for (; token != endOfSequence; ++token) {
                std::string key(token->str());
                enblend::trim(key);

                if (key == "*") {
                    parameter::erase_all();
                } else if (parameter::is_valid_identifier(key)) {
                    parameter::erase(key);
                } else {
                    std::cerr << command << ": warning: key \"" << key <<
                        "\" is not a valid identifier; ignoring\n";
                }
            }

            break;
        }

        case '?':
        {
            if (optopt == 0 || optopt == -1) {
                const int failing_index = optind - 1;

                std::cerr << command << ": unknown long option";
                if (failing_index >= 0 && failing_index < argc) {
                    std::cerr << " \"" << argv[failing_index] << "\"";
                }
                std::cerr << std::endl;
            } else {
                const struct option* failing_long_option = lookup_long_option_by_id(long_options, optopt);

                if (failing_long_option == nullptr) {
                    if (short_option_requires_argument(short_options, optopt)) {
                        std::cerr << command
                                  << ": option \"-" << static_cast<char>(optopt) << "\" requires an argument"
                                  << std::endl;
                    } else {
                        std::cerr << command << ": unknown option ";
                        if (isprint(optopt)) {
                            std::cerr << "\"-" << static_cast<char>(optopt) << "\"";
                        } else {
                            std::cerr << "character 0x" << std::hex << optopt;
                        }
                        std::cerr << std::endl;
                    }
                } else {
                    assert(failing_long_option->has_arg == required_argument);
                    std::cerr << command
                              << ": option \"--" << failing_long_option->name << "\" requires an argument"
                              << std::endl;
                }
            }

            std::cerr << "Try \"enblend --help\" for more information." << std::endl;
            exit(1);
        }

        default:
            std::cerr << command
                      << ": internal error: unhandled command line option"
                      << std::endl;
            exit(1);
        }
    }

    if (contains(optionSet, SaveMasksOption) && contains(optionSet, LoadMasksOption))
    {
        std::cerr << command
                  << ": options \"--load-masks\" and \"--save-masks\" are mutually exclusive" << std::endl;
        failed = true;
    }

    if (failed) {
        exit(1);
    }

    switch (print_only_task)
    {
    case VERSION_ONLY:
#ifdef _MSC_VER
        // deinstall failure hook to catch errors in printVersion directly
        // so printVersion can run to end
        __pfnDliFailureHook2 = nullptr;
#endif
        introspection::printVersion(argc, argv);
        break;                  // never reached
    case USAGE_ONLY:
        printUsage(false);
        break;                  // never reached
    case IMAGE_FORMATS_ONLY:
        introspection::printImageFormats();
        break;                  // never reached
    case SIGNATURE_ONLY:
        introspection::printSignature();
        break;                  // never reached
    case GLOBBING_ALGOS_ONLY:
        introspection::printGlobbingAlgos();
        break;                  // never reached
    case SOFTWARE_COMPONENTS_ONLY:
        introspection::printSoftwareComponents();
        break;                  // never reached
    case GPU_INFO_COMPONENTS_ONLY:
#ifdef OPENCL
        std::cout << "Available, OpenCL-compatible platform(s) and their device(s)\n";
        ocl::print_opencl_information();
        ocl::print_gpu_preference(preferredGPUPlatform, preferredGPUDevice);
        exit(0);
#else
        std::cerr <<
            command << ": option \"--show-gpu-info\" is not implemented in this binary,\n" <<
            command << ": because it was compiled without support for OpenCL" << std::endl;
        exit(1);
#endif
        break;                  // never reached

    case PROCESS_COMPLETELY:
        break;

    default:
        NEVER_REACHED("switch control expression \"print_only_task\" out of range");
    }

    StopAfterMaskGeneration = contains(optionSet, SaveMasksOption) && !contains(optionSet, OutputOption);

    warn_of_ineffective_options(optionSet);

#ifdef OPENCL
    if (UseGPU) {
        initialize_gpu_subsystem(preferredGPUPlatform, preferredGPUDevice);
    }
#endif // OPENCL

    return optind;
}


int main(int argc, char** argv)
{
#ifdef _MSC_VER
    // Make sure the FPU is set to rounding mode so that the lrint
    // functions in float_cast.h will work properly.
    // See changes in vigra numerictraits.hxx
    _controlfp(_RC_NEAR, _MCW_RC);
    // install failure hook for delayed loading of dll
    __pfnDliFailureHook2 = delayHookFailureFunc;
#else
    fesetround(FE_TONEAREST);
#endif

#ifndef _WIN32
    sigemptyset(&SigintMask);
    sigaddset(&SigintMask, SIGINT);

    struct sigaction action;
    action.sa_handler = sigint_handler;
    action.sa_flags   = 0;
    sigemptyset(&(action.sa_mask));
    sigaction(SIGINT, &action, nullptr);
#else
    signal(SIGINT, sigint_handler);
#endif

    if (atexit(cleanup_output) != 0) {
        std::cerr << command << ": warning: could not install cleanup routine\n";
    }

    sig.initialize();

    gsl_set_error_handler_off();

    TIFFSetWarningHandler(tiff_warning);
    TIFFSetErrorHandler(tiff_error);

    //< layer-selector all-layers
    LayerSelection.set_selector(selector::find_by_id(selector::id_t::AllLayersId)->get());

    if (!getopt_long_works_ok())
    {
        std::cerr << command << ": cannot reliably parse command line; giving up\n";
        exit(1);
    }

    int optind;
    try {
        optind = process_options(argc, argv);
    } catch (vigra::StdException& e) {
        std::cerr << command << ": error while processing command line options\n"
                  << command << ": " << e.what()
                  << std::endl;
        exit(1);
    }

#ifdef OPENCL
    if (GPUContext && UseGPU) {
        GPU::StateProbabilities = ocl::create_function<ocl::CalculateStateProbabilities>(GPUContext);
        GPU::DistanceTransform = ocl::create_function<vigra::ocl::DistanceTransformFH>(GPUContext);
        if (BatchCompiler) {
            BatchCompiler->submit(GPU::StateProbabilities.get());
            BatchCompiler->submit(GPU::DistanceTransform.get());
        }
    }
#endif

    enblend::TraceableFileNameList inputTraceableFileNameList;

    // Remaining parameters are input files.
    while (optind < argc) {
        enblend::TraceableFileNameList files;
        enblend::unfold_filename(files, std::string(argv[optind]));
        inputTraceableFileNameList.insert(inputTraceableFileNameList.end(),
                                          files.begin(), files.end());
        optind++;
    }

    if (inputTraceableFileNameList.empty()) {
        std::cerr << command << ": no input files specified\n";
        exit(1);
    }

    if (parameter::as_boolean("dump-global-variables", false)) {
        DUMP_GLOBAL_VARIABLES();
    }

    sig.check();

    for (enblend::TraceableFileNameList::iterator i = inputTraceableFileNameList.begin();
         i != inputTraceableFileNameList.end();
         ++i) {
        if (!enblend::can_open_file((*i)->filename())) {
            (*i)->unroll_trace();
            exit(1);
        }

        if (!vigra::isImage((*i)->filename().c_str())) {
            std::cerr <<
                command << ": cannot process \"" << (*i)->filename() << "\"; not recognized as an image\n" <<
                command << ": info: possible causes:\n" <<
                command << ": info: - An underlying image-processing library does not understand the\n" <<
                command << ": info:   particular compression, sub-format, format extension, ...\n" <<
                command << ": info: - Enblend was not compiled with support for this format, which\n" <<
                command << ": info:   can be checked with \"" << command << " --show-image-formats\".\n" <<
                command << ": info: - The image is corrupted or it is incomplete/truncated.\n" <<
                command << ": info: - It really is not an image.  Honesty, huh?\n";
            (*i)->unroll_trace();
            exit(1);
        }
    }

    LayerSelection.retrieve_image_information(inputTraceableFileNameList.begin(),
                                              inputTraceableFileNameList.end());

    // List of info structures for each input image.
    std::list<vigra::ImageImportInfo*> imageInfoList;
    std::list<vigra::ImageImportInfo*>::iterator imageInfoIterator;

    bool isColor = false;
    std::string pixelType;
    TiffResolution resolution;
    vigra::ImageImportInfo::ICCProfile iccProfile;
    vigra::Rect2D inputUnion;

    int minDim = std::numeric_limits<int>::max();
    selector::layer_ordered_list_t viable_layers;
    selector::layer_ordered_list_t::const_iterator layer;
    unsigned layers = 0;       // total number of layers in image file
    enblend::FileNameList inputFileNameList;
    enblend::TraceableFileNameList::iterator inputFileNameIterator = inputTraceableFileNameList.begin();
    while (inputFileNameIterator != inputTraceableFileNameList.end()) {
        const std::string filename((*inputFileNameIterator)->filename());
        vigra::ImageImportInfo* inputInfo = nullptr;
        try {
            vigra::ImageImportInfo info(filename.c_str());
            if (layers == 0) { // OPTIMIZATION: call only once per file
                layers = info.numImages();
                LayerSelection.set_selector((*inputFileNameIterator)->selector());
                viable_layers.clear();
                viable_layers = LayerSelection.viable_layers(filename);
                layer = viable_layers.begin();
#ifdef DEBUG_FILESPEC
                std::cout << "+ viable_layers(" << filename << ") are [ ";
                std::copy(viable_layers.begin(), viable_layers.end(),
                          std::ostream_iterator<unsigned>(std::cout, " "));
                std::cout << "]\n";
#endif
            }
            inputInfo = new vigra::ImageImportInfo(info);
        } catch (vigra::ContractViolation& exception) {
            std::cerr <<
                command << ": cannot load image \"" << filename << "\"\n" <<
                command << ": " << exception.what() << "\n";
            if (enblend::maybe_response_file(filename)) {
                std::cerr <<
                    command << ": note: maybe you meant a response file and forgot the initial '" <<
                    RESPONSE_FILE_PREFIX_CHAR << "'?\n";
            }
            exit(1);
        }

        assert(layer != viable_layers.end());
        inputInfo->setImageIndex(*layer - 1);

        if (Verbose >= VERBOSE_LAYER_SELECTION) {
            std::cerr << command << ": info: layer selector \"" << LayerSelection.name() << "\" accepts\n"
                      << command << ": info: layer " << *layer << " of " << layers << " in image \""
                      << filename << "\"\n";
        }

        // Save this image info in the list.
        imageInfoList.push_back(inputInfo);
        inputFileNameList.push_back(filename);

        if (Verbose >= VERBOSE_INPUT_IMAGE_INFO_MESSAGES) {
            std::cerr << command
                      << ": info: input image \""
                      << (*inputFileNameIterator)->filename()
                      << "\" "
                      << *layer << '/' << layers << ' ';

            if (inputInfo->isColor()) {
                std::cerr << "RGB ";
            }

            if (!inputInfo->getICCProfile().empty()) {
                std::cerr << "ICC ";
            }

            std::cerr << inputInfo->getPixelType()
                      << " position="
                      << inputInfo->getPosition().x
                      << "x"
                      << inputInfo->getPosition().y
                      << " "
                      << "size="
                      << inputInfo->width()
                      << "x"
                      << inputInfo->height()
                      << std::endl;
        }

        if (inputInfo->numExtraBands() < 1) {
            // Complain about lack of alpha channel.
            std::cerr << command
                      << ": input image \"" << (*inputFileNameIterator)->filename() << "\""
                      << enblend::optional_layer_name(*layer, layers)
                      << " does not have an alpha channel\n";
            (*inputFileNameIterator)->unroll_trace();
            exit(1);
        }

        // Get input image's position and size.
        vigra::Rect2D imageROI(vigra::Point2D(inputInfo->getPosition()),
                               vigra::Size2D(inputInfo->width(), inputInfo->height()));

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
                if (InputProfile == nullptr) {
                    std::cerr << std::endl
                              << command << ": error parsing ICC profile data from file \""
                              << (*inputFileNameIterator)->filename()
                              << "\"" << enblend::optional_layer_name(*layer, layers) << std::endl;
                    (*inputFileNameIterator)->unroll_trace();
                    exit(1);
                }
            }
        } else {
            // Second and later images
            inputUnion |= imageROI;

            if (isColor != inputInfo->isColor()) {
                std::cerr << command << ": input image \""
                          << (*inputFileNameIterator)->filename() << "\""
                          << enblend::optional_layer_name(*layer, layers) << " is "
                          << (inputInfo->isColor() ? "color" : "grayscale") << "\n"
                          << command << ": but previous images are "
                          << (isColor ? "color" : "grayscale")
                          << std::endl;
                (*inputFileNameIterator)->unroll_trace();
                exit(1);
            }
            if (pixelType != inputInfo->getPixelType()) {
                std::cerr << command << ": input image \""
                          << (*inputFileNameIterator)->filename() << "\""
                          << enblend::optional_layer_name(*layer, layers) << " has pixel type "
                          << inputInfo->getPixelType() << ",\n"
                          << command << ": but previous images have pixel type "
                          << pixelType
                          << std::endl;
                (*inputFileNameIterator)->unroll_trace();
                exit(1);
            }
            if (resolution !=
                TiffResolution(inputInfo->getXResolution(), inputInfo->getYResolution())) {
                std::cerr << command << ": info: input image \""
                          << (*inputFileNameIterator)->filename() << "\""
                          << enblend::optional_layer_name(*layer, layers) << " has resolution "
                          << inputInfo->getXResolution() << " dpi x "
                          << inputInfo->getYResolution() << " dpi,\n"
                          << command << ": info: but first image has resolution "
                          << resolution.x << " dpi x " << resolution.y << " dpi"
                          << std::endl;
                (*inputFileNameIterator)->unroll_trace();
            }
            if (iccProfile != inputInfo->getICCProfile()) {
                vigra::ImageImportInfo::ICCProfile mismatchProfile(inputInfo->getICCProfile());
                cmsHPROFILE newProfile = nullptr;
                if (!mismatchProfile.empty()) {
                    newProfile = cmsOpenProfileFromMem(mismatchProfile.data(), mismatchProfile.size());
                    if (newProfile == nullptr) {
                        std::cerr << std::endl
                                  << command << ": error parsing ICC profile data from file \""
                                  << (*inputFileNameIterator)->filename()
                                  << "\"" << enblend::optional_layer_name(*layer, layers) << std::endl;
                        (*inputFileNameIterator)->unroll_trace();
                        exit(1);
                    }
                }

                if (InputProfile == nullptr || newProfile == nullptr ||
                    enblend::profileDescription(InputProfile) != enblend::profileDescription(newProfile)) {
                    const std::string category(BlendColorspace <= IdentitySpace ? "warning" : "info");
                    std::cerr << std::endl << command << ": " << category << ": input image \""
                              << (*inputFileNameIterator)->filename()
                              << "\"" << enblend::optional_layer_name(*layer, layers) << "\n";
                    (*inputFileNameIterator)->unroll_trace();
                    std::cerr << command << ": " << category << ": has ";
                    if (newProfile) {
                        std::cerr << "ICC profile \"" << enblend::profileDescription(newProfile) << "\",\n";
                    } else {
                        std::cerr << "no ICC profile,\n";
                    }
                    std::cerr << command << ": " << category << ": but first image has ";
                    if (InputProfile) {
                        std::cerr << "ICC profile \"" << enblend::profileDescription(InputProfile) << "\";\n";
                    } else {
                        std::cerr << "no ICC profile;\n";
                    }
                    if (BlendColorspace <= IdentitySpace) {
                        std::cerr << command << ": " << category << ": blending images with different color spaces\n"
                                  << command << ": " << category << ": may have unexpected results\n";
                    }
                }
            }

            minDim = std::min(minDim, std::min(inputInfo->width(), inputInfo->height()));
        }

        ++layer;
        if (layer == viable_layers.end()) {
#ifdef DEBUG_FILESPEC
            std::cout << "+ next image\n";
#endif
            layers = 0;
            ++inputFileNameIterator;
        } else {
#ifdef DEBUG_FILESPEC
            std::cout << "+ next layer\n";
#endif
            // We are about to process the next layer in the _same_
            // image.  The imageInfoList already has been updated, but
            // inputTraceableFileNameList still lacks the filename.
            inputTraceableFileNameList.insert(inputFileNameIterator, (*inputFileNameIterator)->clone());
        }
    }

    // Check that more than one input file was given.
    if (imageInfoList.size() <= 1) {
        const size_t n = inputTraceableFileNameList.size();
        const size_t m = imageInfoList.size();

        if (n > m) {
            std::cerr << command << ": warning: selector has rejected " << n - m << " out of " << n << " images\n";
        }

        switch (m) {
        case 0:
            std::cerr << command << ": no input images given\n";
            exit(1);
            break;
        case 1:
            std::cerr << command << ": warning: only one input image given;\n"
                      << command << ": note: Enblend needs two or more overlapping input images in order to do\n"
                      << command << ": note: blending calculations.  The output will be the same as the input.\n";
            break;
        }
    }

    if (resolution == TiffResolution()) {
        std::cerr << command << ": warning: no usable resolution found in first image \""
                  << (*inputTraceableFileNameList.begin())->filename() << "\";\n"
                  << command << ": note: Enblend will assume " << DEFAULT_TIFF_RESOLUTION << " dpi\n";
        ImageResolution = TiffResolution(DEFAULT_TIFF_RESOLUTION, DEFAULT_TIFF_RESOLUTION);
    } else {
        ImageResolution = resolution;
    }

    // Switch to fine mask, if the smallest coarse mask would be less
    // than 64 pixels wide or high.
    if (minDim / CoarsenessFactor < parameter::as_unsigned("smallest-coarse-mask-size", 64) && CoarseMask) {
        std::cerr << command
                  << ": warning: input images too small for coarse mask\n"
                  << command << ": note: switching to fine mask"
                  << std::endl;
        CoarseMask = false;
        if (MainAlgorithm == GraphCut && OptimizeMask) {
            std::cerr << command
                      << ": warning: fine mask combined with graphcut incompatible with mask optimization\n"
                      << command << ": note: defaulting to no optimization"
                      << std::endl;

            OptimizeMask = false;
        }
    }

    if (MaskVectorizeDistance.value() == 0) {
        MaskVectorizeDistance.set_percentage(false);
        MaskVectorizeDistance.set_value(CoarseMask ? coarseMaskVectorizeDistance : fineMaskVectorizeDistance);
    }

    // Create the Info for the output file.
    vigra::ImageExportInfo outputImageInfo(OutputFileName.c_str());

    if (!enblend::has_known_image_extension(OutputFileName)) {
        std::string fallback_output_file_type {parameter::as_string("fallback-output-file-type",
                                                                    DEFAULT_FALLBACK_OUTPUT_FILE_TYPE)};
        if (OutputFileName == "-") {
            outputImageInfo.setFileName("/dev/stdout");
        } else {
            std::cerr <<
                command << ": warning: unknown filetype of output file \"" << OutputFileName << "\"\n" <<
                command << ": note: will fallback to type \"" << fallback_output_file_type << "\"\n";
        }
        enblend::to_upper(fallback_output_file_type);
        outputImageInfo.setFileType(fallback_output_file_type.c_str());
    }

    if (!StopAfterMaskGeneration) {
        OutputIsValid = false;

        // Make sure that inputUnion is at least as big as given by the -f paramater.
        if (OutputSizeGiven) {
            inputUnion |= vigra::Rect2D(OutputOffsetXCmdLine,
                                        OutputOffsetYCmdLine,
                                        OutputOffsetXCmdLine + OutputWidthCmdLine,
                                        OutputOffsetYCmdLine + OutputHeightCmdLine);
        }

        if (!OutputCompression.empty()) {
            outputImageInfo.setCompression(OutputCompression.c_str());
        }

        // If not overridden by the command line, the pixel type of the
        // output image is the same as the input images'.  If the pixel
        // type is not supported by the output format, replace it with the
        // best match.
        {
            const std::string outputFileType = enblend::getFileType(OutputFileName);
            const std::string neededPixelType = OutputPixelType.empty() ? std::string(pixelType) : OutputPixelType;
            const std::string bestPixelType = enblend::bestPixelType(outputFileType, neededPixelType);

            if (neededPixelType != bestPixelType) {
                std::cerr << command
                          << ": warning: "
                          << (OutputPixelType.empty() ? "deduced" : "requested")
                          << " output pixel type is \""
                          << neededPixelType
                          << "\", but image type \""
                          << outputFileType
                          << "\"\n"
                          << command << ": warning: supports \""
                          << bestPixelType
                          << "\" at best;  will use \""
                          << bestPixelType
                          << "\""
                          << std::endl;
            }
            outputImageInfo.setPixelType(bestPixelType.c_str());
            pixelType = enblend::maxPixelType(pixelType, bestPixelType);
        }

        // Set the output image ICC profile
        outputImageInfo.setICCProfile(iccProfile);

        if (BlendColorspace == UndeterminedColorspace &&
            !(iccProfile.empty() || enblend::isFloatingPoint(pixelType))) {
            BlendColorspace = CIELUV;
        }

        if (BlendColorspace == CIECAM || BlendColorspace == CIELAB || BlendColorspace == CIELUV) {
            if (enblend::isFloatingPoint(pixelType)) {
                std::cerr << command <<
                    ": warning: blend color space for floating-point images is not \"identity\"" << std::endl;
            }

            if (InputProfile == nullptr) {
                std::cerr << command << ": warning: input images do not have ICC profiles;\n";
                if (FallbackProfile == nullptr) {
                    std::cerr << command << ": warning: assuming sRGB profile" << std::endl;
                    InputProfile = cmsCreate_sRGBProfile();
                } else {
                    std::cerr << command << ": warning: using fallback profile \""
                              << enblend::profileDescription(FallbackProfile) << "\"" << std::endl;
                    InputProfile = FallbackProfile;
                    FallbackProfile = nullptr; // avoid double freeing
                }
            }
            XYZProfile = cmsCreateXYZProfile();

            const unsigned input_profile_type =
                enblend::profileChannels(InputProfile) > 1 ? TYPE_RGB_DBL : TYPE_GRAY_DBL;

            InputToXYZTransform = cmsCreateTransform(InputProfile, input_profile_type,
                                                     XYZProfile, TYPE_XYZ_DBL,
                                                     RENDERING_INTENT_FOR_BLENDING,
                                                     TRANSFORMATION_FLAGS_FOR_BLENDING);
            if (InputToXYZTransform == nullptr) {
                std::cerr << command << ": error building color transform from \""
                          << enblend::profileName(InputProfile)
                          << " "
                          << enblend::profileDescription(InputProfile)
                          << "\" to XYZ space" << std::endl;
                exit(1);
            }

            XYZToInputTransform = cmsCreateTransform(XYZProfile, TYPE_XYZ_DBL,
                                                     InputProfile, input_profile_type,
                                                     RENDERING_INTENT_FOR_BLENDING,
                                                     TRANSFORMATION_FLAGS_FOR_BLENDING);
            if (XYZToInputTransform == nullptr) {
                std::cerr << command
                          << ": error building color transform from XYZ space to \""
                          << enblend::profileName(InputProfile)
                          << " "
                          << enblend::profileDescription(InputProfile)
                          << "\"" << std::endl;
                exit(1);
            }

            // P2 Viewing Conditions: D50, 500 lumens
            ViewingConditions.whitePoint.X = XYZ_SCALE * cmsD50_XYZ()->X;
            ViewingConditions.whitePoint.Y = XYZ_SCALE * cmsD50_XYZ()->Y;
            ViewingConditions.whitePoint.Z = XYZ_SCALE * cmsD50_XYZ()->Z;
            ViewingConditions.Yb = 20.0;
            ViewingConditions.La = 31.83;
            ViewingConditions.surround = AVG_SURROUND;
            ViewingConditions.D_value = 1.0;

            CIECAMTransform = cmsCIECAM02Init(nullptr, &ViewingConditions);
            if (!CIECAMTransform) {
                std::cerr << std::endl
                          << command
                          << ": error initializing CIECAM02 transform"
                          << std::endl;
                exit(1);
            }

            cmsCIExyY white_point;
            if (cmsIsTag(InputProfile, cmsSigMediaWhitePointTag)) {
                cmsXYZ2xyY(&white_point,
                           (const cmsCIEXYZ*) cmsReadTag(InputProfile, cmsSigMediaWhitePointTag));
                if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES) {
                    double temperature;
                    cmsTempFromWhitePoint(&temperature, &white_point);
                    std::cerr << command
                              << ": info: using white point of input profile at " << temperature << "K"
                              << std::endl;
                }
            } else {
                memcpy(&white_point, cmsD50_xyY(), sizeof(cmsCIExyY));
                if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES) {
                    double temperature;
                    cmsTempFromWhitePoint(&temperature, &white_point);
                    std::cerr << command
                              << ": info: falling back to predefined (D50) white point at " << temperature << "K"
                              << std::endl;
                }
            }
            LabProfile = cmsCreateLab2Profile(&white_point);
            InputToLabTransform = cmsCreateTransform(InputProfile, input_profile_type,
                                                     LabProfile, TYPE_Lab_DBL,
                                                     RENDERING_INTENT_FOR_BLENDING,
                                                     TRANSFORMATION_FLAGS_FOR_BLENDING);
            if (!InputToLabTransform) {
                std::cerr << command << ": error building color transform from \""
                          << enblend::profileName(InputProfile)
                          << " "
                          << enblend::profileDescription(InputProfile)
                          << "\" to Lab space" << std::endl;
                exit(1);
            }
            LabToInputTransform = cmsCreateTransform(LabProfile, TYPE_Lab_DBL,
                                                     InputProfile, input_profile_type,
                                                     RENDERING_INTENT_FOR_BLENDING,
                                                     TRANSFORMATION_FLAGS_FOR_BLENDING);
            if (!LabToInputTransform) {
                std::cerr << command
                          << ": error building color transform from Lab space to \""
                          << enblend::profileName(InputProfile)
                          << " "
                          << enblend::profileDescription(InputProfile)
                          << "\"" << std::endl;
                exit(1);
            }
        } else {
            if (FallbackProfile != nullptr) {
                std::cerr << command <<
                    ": warning: blending in identity space; option \"--fallback-profile\" has no effect" <<
                    std::endl;
            }
        }

        // The size of the output image.
        if (Verbose >= VERBOSE_INPUT_UNION_SIZE_MESSAGES) {
            std::cerr << command
                      << ": info: output image size: "
                      << inputUnion
                      << std::endl;
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
        } catch (vigra::StdException & e) {
            std::cerr << std::endl
                      << command
                      << ": error opening output file \""
                      << OutputFileName
                      << "\";\n"
                      << command
                      << ": "
                      << e.what()
                      << std::endl;
            exit(1);
        }

        if (!OutputPixelType.empty()) {
            pixelType = enblend::maxPixelType(pixelType, OutputPixelType);
        }
    }

    // Invoke templatized blender.
    try {
        if (isColor) {
            if      (pixelType == "UINT8")  enblend::enblendMain<vigra::RGBValue<vigra::UInt8 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "UINT16") enblend::enblendMain<vigra::RGBValue<vigra::UInt16> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enblend::enblendMain<vigra::RGBValue<vigra::Int16 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enblend::enblendMain<vigra::RGBValue<vigra::UInt32> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enblend::enblendMain<vigra::RGBValue<vigra::Int32 > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enblend::enblendMain<vigra::RGBValue<float > >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enblend::enblendMain<vigra::RGBValue<double> >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                std::cerr << command << ": RGB images with pixel type \""
                          << pixelType
                          << "\" are not supported"
                          << std::endl;
                exit(1);
            }
        } else {
            if      (pixelType == "UINT8")  enblend::enblendMain<vigra::UInt8 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#ifndef DEBUG_8BIT_ONLY
            else if (pixelType == "UINT16") enblend::enblendMain<vigra::UInt16>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT16")  enblend::enblendMain<vigra::Int16 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "UINT32") enblend::enblendMain<vigra::UInt32>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "INT32")  enblend::enblendMain<vigra::Int32 >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "FLOAT")  enblend::enblendMain<float >(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
            else if (pixelType == "DOUBLE") enblend::enblendMain<double>(inputFileNameList, imageInfoList, outputImageInfo, inputUnion);
#endif
            else {
                std::cerr << command
                          << ": black&white images with pixel type \""
                          << pixelType
                          << "\" are not supported"
                          << std::endl;
                exit(1);
            }
        }

        for (auto x : imageInfoList) {
            delete x;
        }
        for (auto x : inputTraceableFileNameList) {
            delete x;
        }
    } catch (std::bad_alloc& e) {
        std::cerr << std::endl
                  << command << ": out of memory\n"
                  << command << ": " << e.what()
                  << std::endl;
        exit(1);
    } catch (vigra::StdException& e) {
        std::cerr << std::endl
                  << command << ": an exception occured\n"
                  << command << ": " << e.what()
                  << std::endl;
        exit(1);
    }

#ifdef OPENCL
    delete GPUContext;
#endif // OPENCL

    if (FallbackProfile) {cmsCloseProfile(FallbackProfile);}
    if (LabProfile) {cmsCloseProfile(LabProfile);}
    if (InputToLabTransform) {cmsCIECAM02Done(InputToLabTransform);}
    if (LabToInputTransform) {cmsCIECAM02Done(LabToInputTransform);}
    if (CIECAMTransform) {cmsCIECAM02Done(CIECAMTransform);}
    if (InputToXYZTransform) {cmsDeleteTransform(InputToXYZTransform);}
    if (XYZToInputTransform) {cmsDeleteTransform(XYZToInputTransform);}
    if (XYZProfile) {cmsCloseProfile(XYZProfile);}
    if (InputProfile) {cmsCloseProfile(InputProfile);}

    // Success.
    return 0;
}
