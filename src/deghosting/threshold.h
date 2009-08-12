/**
 * thresshold
 */

#include <vector>
#include <stdint.h>

#include <boost/shared_ptr.hpp>
#include <vigra/stdimage.hxx>
#include <vigra/transformimage.hxx>

#include "deghosting.h"

using std::vector;
using namespace vigra;
using namespace deghosting;

const uint16_t THRESHOLD_DONTCARE = 1;

/** Threshold function
 * used for creating alpha masks for images
 * @param const vector<FImagePtr> vector of images
 * @param const int threshold all pixels above this thresshold are set to 255, others to 0
 * @param const uint16_t flags flags for setting the behavior
 *          possible values are:
 *              0 – applies only threshold
 *              ONE_UNMASKED – if pixel should be black in all images after applying threshold
 *                             leave it in one image (where the pixel value is highest) white
 *              ONE_MASKED – make the pixel black only in one image where it's value is lowest
 */
template <class Functor>
vector<BImagePtr> threshold(const vector<FImagePtr> &inputImages, const double threshold, const Functor &f, uint16_t flags) {
    vector<BImagePtr> retVal;
    const uint8_t minValue = 0;
    const uint8_t maxValue = 255;
    
    // don't care about masking
    if (flags & THRESHOLD_DONTCARE) {
        for (unsigned int i=0; i < inputImages.size(); ++i) {
            BImagePtr tmpImg(new BImage(inputImages[i]->size()));
            // preprocess image using functor f
            transformImage(srcImageRange(*(inputImages[i])), destImage(*tmpImg), f);
            transformImage(srcImageRange(*tmpImg), destImage(*tmpImg),
                            Threshold<FImage::PixelType, BImage::PixelType>(threshold, 255, 0, 255));
            retVal.push_back(tmpImg);
        }
        return retVal;
    }
    
    // arrays with iterators
    FImage::traverser siterators[inputImages.size()];
    BImage::traverser diterators[inputImages.size()];
    // iterator to the end
    FImage::traverser send = inputImages[0]->lowerRight();
    // fill iterators and retVal
    for (unsigned int i=0; i < inputImages.size(); ++i) {
        // fill retVal
        BImagePtr tmpImg(new BImage(inputImages[i]->size()));
        retVal.push_back(tmpImg);
        // fill iterators
        siterators[i] = inputImages[i]->upperLeft();
        diterators[i] = retVal[i]->upperLeft();
    }
    
    // leave pixels not masked in at least one image
    // this is DEFAULT
    // loop over row
    while (siterators[0].y != send.y) {
        // array with column iterators
        FImage::traverser siteratorsX[inputImages.size()];
        BImage::traverser diteratorsX[inputImages.size()];
        for (unsigned int i=0; i < inputImages.size(); ++i) {
            siteratorsX[i] = siterators[i];
            diteratorsX[i] = diterators[i];
        }
        // now we can loop over the same pixels in all images at once
        while (siteratorsX[0].x != send.x) {
            // image with highest weight for the pixel
            unsigned int highestI = 0;
            for (unsigned int i=0; i<inputImages.size(); ++i) {
                // apply threshold
                if (*(siteratorsX[i]) < threshold)
                    *(diteratorsX[i]) = minValue;
                else
                    *(diteratorsX[i]) = maxValue;
                // highest search
                highestI = (*(siteratorsX[highestI]) > *(siteratorsX[i])) ? highestI : i;
            }
            // set the highest one to maxValue;
            *(diteratorsX[highestI]) = maxValue;
            
            // move all iterators by one pixel to the right
            for (unsigned int i=0; i<inputImages.size(); ++i) {
                ++siteratorsX[i].x;
                ++diteratorsX[i].x;
            }
        }
        // move all iterators to the next row
        for (unsigned int i=0; i<inputImages.size(); ++i) {
            ++siterators[i].y;
            ++diterators[i].y;
        }
    }
}
