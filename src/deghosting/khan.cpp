
#include <signal.h>

#include "khan.h"
#include <vigra/impex.hxx>
// for importImageAlpha, it's not really necessary
#include <vigra_ext/impexalpha.hxx>

#ifdef ATAN_KH
// neded for atan and abs
#include <cmath>
#endif

// for GammaFunctor
#include <vigra/transformimage.hxx>

// for resampleImage
#include <vigra/resizeimage.hxx>

// for snprintf
#include <cstdio>
using std::snprintf;

using std::endl;
using std::cout;

namespace deghosting {
    Khan::Khan(std::vector<std::string>& setInputFiles, const uint16_t setFlags, const uint16_t setDebugFlags,
                int setIterations, int setVerbosity) {
        inputFiles = setInputFiles;
        flags = setFlags;
        debugFlags = setDebugFlags;
        iterations = setIterations;
        for (unsigned int i=0; i<5; i++)
            response.push_back(0);
        verbosity = setVerbosity;
        sigma = 30;
        PIPOW = sigma*std::sqrt(2*PI);
        denom = 1/PIPOW;
    }

    void Khan::loadImages(std::vector<std::string>& newInputFiles) {
        inputFiles = newInputFiles;
    }

    void Khan::setFlags(const uint16_t newFlags) {
        flags = newFlags;
    }

    void Khan::setDebugFlags(const uint16_t newFlags) {
        debugFlags = newFlags;
    }

    void Khan::setIterationNum(const int newIterations) {
        iterations = newIterations;
    }
    
    void Khan::setCameraResponse(EMoR newResponse) {
        response = newResponse;
    }
    
    void Khan::setVerbosity(int newVerbosity) {
        verbosity = newVerbosity;
    }
    
    void Khan::setSigma(double newSigma) {
        sigma = newSigma;
    }

    float Khan::hat(RGBValue<float> pixel) {
        double value = (pixel[0]+pixel[1]+pixel[2])/3;
        double t = (value/127.5 -1);
        t *= t; // ^2
        t *= t; // ^4
        t *= t; // ^8
        t *= t; // ^16
        return 1.0 - t; 
    }
    
    double Khan::Kh(deghosting::AlgTinyVector< float, 3 > x) {
        #ifdef ATAN_KH
            // good choice for sigma for this function is around 600
            return std::atan(-(x*x)+sigma)/PI + 0.5;
        #else
            // good choice for sigma for this function is around 30
            return (std::exp(-(x*x)/(2*sigma*sigma)) * denom);
        #endif
    }
    
    void Khan::preprocessImage(unsigned int i, FImagePtr &weight, FLabImagePtr &LabImage) {
        cout << "Loading image number " << i << endl;
        ImageImportInfo imgInfo(inputFiles[i].c_str());
        FRGBImage * pInputImg =  new FRGBImage(imgInfo.size());
        BImage imgAlpha(imgInfo.size());
        weight = FImagePtr(new FImage(imgInfo.size()));
        
        // import alpha
        importImageAlpha(imgInfo, destImage(*pInputImg),destImage(imgAlpha));
        
        // take logarithm or gamma correction if the input images are HDR
        // I'm not sure if it's the right way how to
        // find out if they are HDR
        if ( (!strcmp(imgInfo.getFileType(),"TIFF") && strcmp(imgInfo.getPixelType(),"UINT8")) ||
                !strcmp(imgInfo.getFileType(),"EXR") ||
                !strcmp(imgInfo.getPixelType(),"FLOAT"))
        {
            // take logarithm
            if (flags & ADV_LOGARITHM) {
                transformImage(srcImageRange(*pInputImg),destImage(*pInputImg),LogarithmFunctor<RGBValue<float> >(1.0));
            }
            // use gamma 2.2
            else if (flags & ADV_GAMMA) {
                // GammaFunctor is only in vigra 1.6 GRRR
                // I have to use BrightnessContrastFunctor
                // TODO: change to the GammaFunctor in the future
                vigra::FindMinMax<float> minmax;
                vigra::inspectImage(srcImageRange(*pInputImg), minmax);
                transformImage(srcImageRange(*pInputImg),destImage(*pInputImg),BrightnessContrastFunctor<RGBValue<float> >(0.45,1.0,minmax.min, minmax.max));
            }
        }
        
        // generate initial weights using hat function
        if (flags & ADV_NOINITWEIGHTS) {
            // convert image to grayscale (and take logarithm)
            RGBToGrayAccessor<FRGBImage::PixelType> color2gray;
            transformImage(srcImageRange(*pInputImg, color2gray), destImage(*weight), log(Arg1()+Param(1.0f)));
        }
        else {
            transformImage(srcImageRange(*pInputImg),destImage(*weight),hat);
        }
        
        // convert from linear RGB to L*a*b //
        RGB2LabFunctor<float> RGB2Lab;
        LabImage = FLabImagePtr(new FLabImage(imgInfo.size()));
        transformImage(srcImageRange(*pInputImg), destImage(*LabImage), RGB2Lab);
        
        delete pInputImg;
        pInputImg = 0;
    }
    
    std::vector<FImagePtr> Khan::createWeightMasks() {
        for (unsigned int i = 0; i < inputFiles.size(); i++) {
            FImagePtr weight;
            FLabImagePtr LabImage;
            preprocessImage(i, weight, LabImage);
            LabImages.push_back(LabImage);
            weights.push_back(weight);
            
            // save init weights
            if (debugFlags & SAVE_INITWEIGHTS) {
                char tmpfn[100];
                snprintf(tmpfn, 99, "init_weights_%d.tiff", i);
                ImageExportInfo exWeights(tmpfn);
                exportImage(srcImageRange(*weight), exWeights.setPixelType("UINT8"));
            }
        }
        
        float maxWeight = 0;
        // image size
        const int origWidth = weights[0]->width();
        const int origHeight = weights[0]->height();
        
        // if we doing scaling, we have to backup L*a*b images of original size
        std::vector<FLabImagePtr> backupLab;
        if (flags & ADV_MULTIRES) {
            for (unsigned int i = 0; i < LabImages.size(); i++) {
                // backup original size L*a*b
                backupLab.push_back(LabImages[i]);
            }
        }
        
        cout << endl << "Running khan algorithm" << endl;
        // and we can run khan algorithm
        // khan iteration
        for (int it = 0; it < iterations; it++) {
            if (verbosity > 0)
                cout << "iteration " << it+1 << endl;
            // copy weights from previous iteration
            if (verbosity > 1)
                cout << "copying weights from previous iteration" << endl;
            
            std::vector<FImagePtr> prevWeights;
            for (unsigned int i = 0; i < weights.size(); i++) {
                // scale weights to the requied size
                if (flags & ADV_MULTIRES) {
                    // it would be better to use resampleImage, but it seems to be present only in VIGRA 1.6
                    // so let's use some of the resizeImageINTERPOLATION() functions
                    
                    // compute width
                    int resized_width = origWidth / ( iterations/(it+1) );
                    //compute height
                    int resized_height = origHeight / ( iterations/(it+1) );
                    // destination images
                    FImage resizedWeight;
                    FLabImage resizedLab;
                    // it's not worthy to scale to less than 100px per side
                    if (resized_width > 100 && resized_height > 100) {
                        // create destination image of desired size
                        resizedWeight = FImage(Size2D(resized_width,resized_height));
                        resizedLab = FLabImage(Size2D(resized_width,resized_height));
                    } else if (origWidth >= 100 && origHeight >= 100) {
                        // resize it to the smallest value (ie 100px for the shorter side)
                        if (origWidth >= origHeight) {
                            resizedWeight = FImage(Size2D(100*origWidth/origHeight, 100));
                            resizedLab = FLabImage(Size2D(100*origWidth/origHeight, 100));
                        } else {
                            resizedWeight = FImage(Size2D(100, 100*origHeight/origWidth));
                            resizedLab = FLabImage(Size2D(100, 100*origHeight/origWidth));
                        }
                    } else {
                        // don't scale at all
                        // just copy weights as if no scaling seting was applied
                        goto DONTSCALE;
                    }
                    
                    // No interpolation – only for testing
                    resizeImageNoInterpolation(srcImageRange(*weights[i]), destImageRange(resizedWeight));
                    resizeImageNoInterpolation(srcImageRange(*backupLab[i]), destImageRange(resizedLab));
                    
                    FImagePtr tmp(new FImage(resizedWeight));
                    prevWeights.push_back(tmp);
                    LabImages[i] = FLabImagePtr(new FLabImage(resizedLab));
                    weights[i] = FImagePtr(new FImage(resizedWeight));
                } else {
                    DONTSCALE:
                    FImagePtr tmp(new FImage(*weights[i]));
                    prevWeights.push_back(tmp);
                }
            }
            
            // loop through all images
            for (unsigned int i = 0; i < LabImages.size(); i++) {
                if (verbosity > 1)
                    cout << "processing image " << i+1 << endl;
                
                // vector storing pixel data
                AlgTinyVector<float, 3> X;
                // sums for eq. 6
                double wpqssum = 0;
                double wpqsKhsum = 0;
                // image size
                const int width = LabImages[i]->width();
                const int height = LabImages[i]->height();

                // iterator to the upper left corner
                FLabImage::traverser sy = LabImages[i]->upperLeft();
                // iterator to the lower right corner
                FLabImage::traverser send = LabImages[i]->lowerRight();
                // iterator to the weight image left corner
                FImage::traverser wy = weights[i]->upperLeft();
                // loop through the row
                for (int y=0; sy.y != send.y; ++sy.y, ++wy.y, ++y) {
                    // iterator to the source (L*a*b image)
                    FLabImage::traverser sx = sy;
                    // iterator to the weight
                    FImage::traverser wx = wy;
                    // loop over the pixels
                    for (int x=0; sx.x != send.x; ++sx.x, ++wx.x, ++x) {
                        if (verbosity > 2)
                            cout << "processing pixel (" << x+1 << "," << y+1 << ")" << endl;
                        // set pixel vector
                        X = *sx;
                        
                        // loop through all layers
                        for (unsigned int j = 0; j < LabImages.size(); j++) {
                            if (verbosity > 2)
                                cout << "processing layer " << j << endl;
                            
                            // iterator to the neighbour
                            FLabImage::traverser neighby = LabImages[j]->upperLeft();
                            // iterator to the weight
                            FImage::traverser weighty = prevWeights[j]->upperLeft();
                            // pixel offset
                            int ndy = -NEIGHB_DIST;
                            // move them to the upper bound
                            // find the best upper bound
                            if (y-NEIGHB_DIST < 0) {
                                ndy = -y;
                            }
                            else {
                                neighby.y += y-NEIGHB_DIST;
                                weighty.y += y-NEIGHB_DIST;
                            }
                            
                            // iterate through neighbourhoods y axis
                            int maxDisty = (height - y) > NEIGHB_DIST ? NEIGHB_DIST : (height - y-1);
                            for (; ndy <= maxDisty; ++neighby.y, ++weighty.y, ++ndy) {
                                FLabImage::traverser neighbx = neighby;
                                FImage::traverser weightx = weighty;
                                // pixel offset
                                int ndx = -NEIGHB_DIST;
                                // move them to the upper bound
                                // find the best upper bound
                                if (x-NEIGHB_DIST < 0) {
                                    ndx = -x;
                                }
                                else {
                                    neighbx.x += x-NEIGHB_DIST;
                                    weightx.x += x-NEIGHB_DIST;
                                }
                                // iterate through neighbourhoods x axis
                                int maxDistx = (width - x) > NEIGHB_DIST ? NEIGHB_DIST : (width - x-1);
                                for (; ndx <= maxDistx; ++neighbx.x, ++weightx.x, ++ndx) {
                                    if (verbosity > 3)
                                        cout << "(" << ndx << "," << ndy << ")";
                                    // now we can construct pixel vector
                                    // should omit the middle pixel, ie use only neighbours
                                    if (ndx != ndy != 0) {
                                        wpqsKhsum += (*weightx * Kh(X-(*neighbx)));
                                        wpqssum += *weightx;
                                    }
                                    
                                    maxDistx = (width - x) > NEIGHB_DIST ? NEIGHB_DIST : (width - x-1);
                                }
                                if (verbosity > 3)
                                    cout << endl;
                                
                                maxDisty = (height - y) > NEIGHB_DIST ? NEIGHB_DIST : (height - y-1);
                            }
                        }
                        
                        if (verbosity > 2)
                            cout << "computing new weight" << endl;
                        // compute probability and set weight
                        //cout << "P=" << (float) wpqsKhsum/wpqssum << endl;
                        if (flags & ADV_ONLYP)
                            *wx = (float) wpqsKhsum/wpqssum;
                        else
                            *wx *= (float) wpqsKhsum/wpqssum;
                        if (maxWeight < *wx)
                            maxWeight = *wx;
                        wpqsKhsum = wpqssum = 0;
                        
                    }
                }
            }
        }
        
        if (verbosity > 1)
                cout << "normalizing weights" << endl;
        double factor = 255.0f/maxWeight;
        for (unsigned int i=0; i<weights.size(); ++i) {
            transformImage(srcImageRange(*(weights[i])), destImage(*(weights[i])), NormalizeFunctor<float>(factor));
        }
        return weights;
    }
}
