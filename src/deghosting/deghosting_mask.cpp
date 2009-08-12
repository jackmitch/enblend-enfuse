
#include <iostream>
#include <string>
#include <cstring>

#include <hugin_config.h>
#include <hugin_version.h>
// for stripExtension
#include <hugin_utils/utils.h>
// for exportImage
#include <vigra_ext/impexalpha.hxx>

#include "deghosting.h"
#include "support.h"
#include "threshold.h"

// deghosting algorithms
#include "khan.h"

//#ifdef WIN32
 #include <getopt.h>
//#else
 #include <unistd.h>
//#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace deghosting;

// options for otherFlags
static uint16_t SAVE_GENWEIGHTS = 1;

// globals containing settings
static int iterations = 1;
static double sigma = 30;
//static uint16_t flags = ADV_ONLYP;
static uint16_t flags = ADV_ONLYP + ADV_MULTIRES;
static uint16_t otherSaveFlags = 0;
static uint16_t otherThresholdFlags = 0;
static uint16_t debugFlags = 0;
static EMoR response(0.0f);
static int verbosity = 0;
static double thresholdLim = 150;
static double contrast = 1.3;

/** function handling advanced options
 */
void parseOptions_advanced(char* optarg) {
    for(char *c = optarg; *c; c++) {
        switch(*c) {
            case 'g':
                if (flags & ADV_LOGARITHM) {
                    cerr << "g cannot be combined with l" << endl;
                    exit(1);
                }
                flags += ADV_GAMMA;
                break;
            case 'l':
                if (flags & ADV_GAMMA) {
                    cerr << "l cannot be combined with g" << endl;
                    exit(1);
                }
                flags += ADV_LOGARITHM;
                break;
            case 'm':
                flags -= ADV_MULTIRES;
                break;
            case 'n':
                flags += ADV_NOINITWEIGHTS;
                break;
            case 't':
                otherThresholdFlags += THRESHOLD_DONTCARE;
                break;
            case 'w':
                flags -= ADV_ONLYP;
                break;
            default:
                cerr<< "Error: unknown option" << endl;
                exit(1);
        }
    }
}

/** function handling save options
 */
void parseOptions_save(char* optarg) {
    for(char *c = optarg; *c; c++) {
        switch(*c) {
            case 'i':
                debugFlags += SAVE_INITWEIGHTS;
                break;
            case 'w':
                otherSaveFlags += SAVE_GENWEIGHTS;
                break;
            default:
                cerr<< "Error: unknown option" << endl;
                exit(1);
        }
    }
}

static void usage()
{
    cerr << "deghosting_mask: creates mask for removing ghosting in images" << endl
         << "deghosting_mask version " << PACKAGE_VERSION << endl
         << endl
         << "Usage: deghosting_mask [options] inputfile(s) " << endl
         << "   option are: " << endl
         << "     -o, --output=PREFIX       prefix for output masks" << endl
         << "     -i, --iterations=ITER     number of iterations, default is (ITER > 0)" << endl
         << "                               default: 1" << endl
         << "     -s, --sigma=SIGMA         standard deviation of Gaussian weighting" << endl
         << "                               function (SIGMA > 0); default: 30" << endl
         << "     -r, --response=E:M:o:R    use camera response specified in EMoR format" << endl
         << "     -t, --threshold=THRESH    threshold; default: 150" << endl
         << "     -c, --contrast=CONTR      change constrast before applying threshold;" << endl
         << "                               default: 1.3" << endl
         << "     -a, --advanced=SET        advanced settings. Possible options are:" << endl
         << "                               l   use logarithm on input if input images are HDR," << endl
         << "                                   cannot be combined with g" << endl
         << "                               g   use gamma 2.2 correction on input if input images are HDR," << endl
         << "                                   cannot be combined with l" << endl
         << "                               m   do not scale image, NOTE: slows down process" << endl
         << "                               n   do not compute initial weigths" << endl
         << "                               t   use simple threshold, may result in holes in images" << endl
         << "                               w   compute \"complete\" weights, not only probabilities" << endl
         << "     -w, --save=SET            advanced save settings" << endl
         << "                               i   save initial weigths" << endl
         << "                               w   save generated weigths" << endl
         << "     -b BLOCKSIZE              image cache BLOCKSIZE in kilobytes; default: " <<
            (CachedFileImageDirector::v().getBlockSize() / 1024LL) << "KB" << endl
         << "     -m CACHESIZE              set image CACHESIZE in megabytes; default: " << 
            (CachedFileImageDirector::v().getAllocation() / 1048576LL) << "MB" << endl
         << "     -h, --help                display this help" << endl
         << "     -v, --verbose             verbose, repeat for more verbose output" << endl;
}

int main(int argc, char *argv[]) {
    const char * optstring = "o:i:s:r:t:c:a:w:b:m:hv";
    opterr = 0;
    int c;
    
    string outputPrefix = "weight";
    
    enum optionArgumentKind {
        NoArgument,
        StringArgument,
        DoubleArgument,
        IntegerArgument,
        ArrayArgument
    };
    
    enum optionId {
        outputID,
        iterationsID,
        sigmaID,
        responseID,
        thresholdID,
        constrastID,
        advancedID,
        saveID,
        helpID,
        verboseID
    };
      
    static struct option longOptions[] = {
        {"output", 1, 0, StringArgument},
        {"iterations", 1, 0, IntegerArgument},
        {"sigma", 1, 0, DoubleArgument},
        {"response", 1, 0, ArrayArgument},
        {"threshold", 1, 0, DoubleArgument},
        {"contrast", 1, 0, DoubleArgument},
        {"advanced", 1, 0, StringArgument},
        {"save", 1, 0, StringArgument},
        {"help", 0, 0, NoArgument},
        {"verbose", 0, 0, NoArgument},
        {0, 0, 0, 0}
    };
    
    // TEST
    // response for testing
    response.resize(5);
    response[0] = -3.59f;
    response[1] = -0.93f;
    response[2] =  0.11f;
    response[3] = -0.22f;
    response[4] =  0.34f;
       
    int optionIndex = 0;
    
    while ((c = getopt_long(argc, argv, optstring, longOptions, &optionIndex)) != -1) {
        switch (c) {
            case NoArgument: {
                if (longOptions[optionIndex].flag != 0) break;
                switch (optionIndex) {
                    case helpID:
                        usage();
                        return 0;
                    case verboseID:
                        verbosity++;
                        break;
                    default:
                        cerr << "There's a problem with parsing options" << endl;
                        return 1;
                }
                break;
            }
            
            case StringArgument: {
                if (longOptions[optionIndex].flag != 0) break;
                switch (optionIndex) {
                    case outputID:
                        outputPrefix = optarg;
                        break;
                    case advancedID:
                        parseOptions_advanced(optarg);
                        break;
                    case saveID:
                        parseOptions_save(optarg);
                        break;
                    default:
                        cerr << "There's a problem with parsing options" << endl;
                        return 1;
                }
                break;
            }
            
            case IntegerArgument: {
                if (longOptions[optionIndex].flag != 0) break;
                switch (optionIndex) {
                    case iterationsID:
                        iterations = atoi(optarg);
                        break;
                    default:
                        cerr << "There's a problem with parsing options" << endl;
                        return 1;
                }
                break;
            }
            
            case DoubleArgument: {
                if (longOptions[optionIndex].flag != 0) break;
                switch (optionIndex) {
                    case sigmaID:
                        sigma = atof(optarg);
                        break;
                    case thresholdID:
                        thresholdLim = atof(optarg);
                        break;
                    case constrastID:
                        contrast = atof(optarg);
                        break;
                }
                break;
            }
            
            case ArrayArgument: {
                if (longOptions[optionIndex].flag != 0) break;
                switch (optionIndex) {
                    case responseID:
                        // TODO
                        break;
                }
                break;
            }
            
            case 'o':
                outputPrefix = optarg;
                break;
            case 'i':
                iterations = atoi(optarg);
                break;
            case 's':
                sigma = atof(optarg);
                break;
            case 'r':
                // TODO
                break;
            case 't':
                thresholdLim = atof(optarg);
                break;
            case 'c':
                contrast = atof(optarg);
                break;
            case 'a':
                parseOptions_advanced(optarg);
                break;
            case 'w':
                parseOptions_save(optarg);
                break;
            case 'b': {
                const int kilobytes = atoi(optarg);
                if (kilobytes < 1) {
                    cerr << "cache block size must be 1 or more." << endl;
                    return 1;
                }
                CachedFileImageDirector::v().setBlockSize(static_cast<long long>(kilobytes) << 10);
                break;
            }
            case 'm': {
                const int megabytes = atoi(optarg);
                if (megabytes < 1) {
                    cerr << "memory limit must be 1 or more." << endl;
                    return 1;
                }
                CachedFileImageDirector::v().setAllocation(static_cast<long long>(megabytes) << 20);
                break;
            }
            case 'h':
                usage();
                return 0;
            case 'v':
                verbosity++;
                break;
        }
    }
    
    cout << endl;
    
    unsigned nFiles = argc - optind;
    if (nFiles == 0) {
        cerr << std::endl << "Error: at least three input images needed" << std::endl <<std::endl;
        usage();
        return 1;
    } else if (nFiles < 3) {
        std::cout << std::endl << "Error: You have to specify at least three images." << std::endl;
        return 1;
    }
    
    // load all images
    vector<string> inputFiles;
    for (size_t i=optind; i < (size_t)argc; i++)
    {
        inputFiles.push_back(argv[i]);
    }
    
    Deghosting* deghoster = NULL;
    
    Khan khanDeghoster(inputFiles, flags, debugFlags, iterations, verbosity);
    khanDeghoster.setSigma(sigma);
    
    deghoster = &khanDeghoster;
    deghoster->setCameraResponse(response);
    
    vector<FImagePtr> weights = deghoster->createWeightMasks();
    
    // save weights
    if (otherSaveFlags & SAVE_GENWEIGHTS) {
        for (unsigned int i=0; i<weights.size(); ++i) {
            char tmpfn[100];
            snprintf(tmpfn, 99, "%s_%d.tif", outputPrefix.c_str(), i);
            ImageExportInfo exWeights(tmpfn);
            exportImage(srcImageRange(*weights[i]), exWeights.setPixelType("UINT8"));
        }
    }
    
    vector<BImagePtr> thresholded = threshold(weights, thresholdLim, BrightnessContrastFunctor<FImage::PixelType>(1, contrast, 0, 255), otherThresholdFlags);
    
    // save masks with treshold applied
    for (unsigned int i=0; i<weights.size(); ++i) {
        char tmpfn[100];
        string fileName = hugin_utils::stripExtension(inputFiles[i]);
        fileName = hugin_utils::stripPath(fileName);
        snprintf(tmpfn, 99, "%s_mask.tif", fileName.c_str());
        ImageExportInfo exWeights(tmpfn);
        exportImage(srcImageRange(*thresholded[i]), exWeights.setPixelType("UINT8"));
    }
}
