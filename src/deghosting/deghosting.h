/*
 * deghosting.h
 *
 *  Created on: May 23, 2009
 *      Author: Lukas "stativ" Jirkovsky
 */

#ifndef DEGHOSTING_H_
#define DEGHOSTING_H_

#include <vector>
#include <string>
#include <stdint.h>

#include <boost/shared_ptr.hpp>
#include <vigra/stdimage.hxx>

// define if you want to use image cache
#define DEGHOSTING_CACHE_IMAGES

namespace deghosting {
    typedef boost::shared_ptr<vigra::BImage> BImagePtr;
    typedef boost::shared_ptr<vigra::FImage> FImagePtr;
    // type for camera response
    typedef std::vector<float> EMoR;

    // constants for advanced modes
    const uint16_t ADV_LOGARITHM     = 1;
    const uint16_t ADV_GAMMA         = 2;
    const uint16_t ADV_NOINITWEIGHTS = 4;
    const uint16_t ADV_ONLYP         = 8;
    const uint16_t ADV_MULTIRES      = 16;

    // constants for debug modes
    const uint16_t SAVE_INITWEIGHTS   = 1;

    class Deghosting
    {
    public:
        Deghosting() {}
        
        /** create weight masks
         * create weight masks for masking out ghosting regions
         */
        virtual std::vector<FImagePtr> createWeightMasks() {}

        /** load images for processing
         * @param inputFiles images to be processed
         */
        virtual void loadImages(std::vector<std::string>& inputFiles) {}
        
        /** set advanced flags
         * Allows to change behavior of used algorithm
         * @param flags one of the constants describing advanced mode
         */
        virtual void setFlags(const uint16_t flags) {}
        
        /** set flags for debugging purposes
         * @param debugFlags one of the constants describing action which should be done
         */
        virtual void setDebugFlags(const uint16_t debugFlags) {}
        
        /** set number of iterations
         */
        virtual void setIterationNum(const int iterations) {}
        
        /** set camera response function
         * set camera response function in EMoR format
         * @param response array of five floats representing response
         */
        virtual void setCameraResponse(EMoR response) {}
        
        /** set verbosity level
         * @param verbosity the higher the number is, the more verbose algorithm will be
         */
        virtual void setVerbosity(int verbosity) {}
        virtual ~Deghosting() {}
    };

}

#endif /* DEGHOSTING_H_ */
