/*
 * Copyright (C) 2004 Andrew Mihal
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
#ifndef __ENBLEND_H__
#define __ENBLEND_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>
#include <stdio.h>

#include <boost/static_assert.hpp>

#include "common.h"
#include "assemble.h"
#include "blend.h"
#include "bounds.h"
#include "mask.h"
#include "pyramid.h"

#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::BasicImage;
using vigra::BCFImage;
using vigra::BImage;
using vigra::CachedFileImage;
using vigra::CachedFileImageDirector;
using vigra::FindMinMax;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::initImage;
using vigra::initImageIf;
using vigra::inspectImage;
using vigra::NumericTraits;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

namespace enblend {

/** Enblend's main blending loop. Templatized to handle different image types.
 */
template <typename ImageType, typename MaskPyramidType, typename ImagePyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    #ifdef ENBLEND_CACHE_IMAGES
    // FIXME mask caching turned off
    typedef BImage MaskType;
    typedef BCFImage AlphaType;
    #else
    typedef BImage MaskType;
    typedef BImage AlphaType;
    #endif

    typedef typename ImageType::value_type ImageValueType;
    typedef typename AlphaType::value_type AlphaValueType;
    typedef typename MaskPyramidType::value_type MaskPyramidValueType;
    typedef typename ImagePyramidType::value_type ImagePyramidValueType;

    // Create the initial black image.
    EnblendROI blackBB;
    pair<ImageType*, AlphaType*> blackPair =
            assemble<ImageType, AlphaType>(imageInfoList, inputUnion, blackBB);
    exportImageAlpha(srcImageRange(*(blackPair.first)),
                     srcImage(*(blackPair.second)),
                     outputImageInfo);
    // mem usage before = 0
    // mem xsection = up to 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
    // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType

    #ifdef ENBLEND_CACHE_IMAGES
    if (Verbose > VERBOSE_CFI_MESSAGES) {
        CachedFileImageDirector &v = CachedFileImageDirector::v();
        cout << "Image cache statistics after loading black image:" << endl;
        v.printStats("blackImage", blackPair.first);
        v.printStats("blackAlpha", blackPair.second);
        v.printStats();
        v.resetCacheMisses();
        cout << "--------------------------------------------------------------------------------" << endl;
    }
    #endif

    // Main blending loop.
    while (!imageInfoList.empty()) {

        // Create the white image.
        EnblendROI whiteBB;
        pair<ImageType*, AlphaType*> whitePair =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, whiteBB);
        // mem usage before = inputUnion*ImageValueType + inputUnion*AlphaValueType
        // mem xsection = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after loading white image:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats("whiteImage", whitePair.first);
            v.printStats("whiteAlpha", whitePair.second);
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        //ImageExportInfo whiteInfo("enblend_white.tif");
        //exportImageAlpha(srcImageRange(*(whitePair.first)),
        //                 srcImage(*(whitePair.second)),
        //                 whiteInfo);

        // Union bounding box of whiteImage and blackImage.
        EnblendROI uBB;
        whiteBB.unite(blackBB, uBB);

        if (Verbose > VERBOSE_UBB_MESSAGES) {
            cout << "image union bounding box: ("
                 << uBB.getUL().x
                 << ", "
                 << uBB.getUL().y
                 << ") -> ("
                 << uBB.getLR().x
                 << ", "
                 << uBB.getLR().y
                 << ")" << endl;
        }

        // Intersection bounding box of whiteImage and blackImage.
        EnblendROI iBB;
        bool iBBValid = whiteBB.intersect(blackBB, iBB);

        if (Verbose > VERBOSE_IBB_MESSAGES) {
            if (iBBValid) {
                cout << "image intersection bounding box: ("
                     << iBB.getUL().x
                     << ", "
                     << iBB.getUL().y
                     << ") -> ("
                     << iBB.getLR().x
                     << ", "
                     << iBB.getLR().y
                     << ")" << endl;
            } else {
                cout << "image intersection bounding box: (no intersection)"
                     << endl;
            }
        }

        // Determine what kind of overlap we have.
        Overlap overlap = inspectOverlap(uBB.apply(srcImageRange(*(blackPair.second))),
                uBB.apply(srcImage(*(whitePair.second))));

        // If white image is redundant, skip it and go to next images.
        if (overlap == CompleteOverlap) {
            // White image is redundant.
            delete whitePair.first;
            delete whitePair.second;
            cerr << "enblend: some images are redundant and will not be blended."
                 << endl;
            continue;
        }
        else if (overlap == NoOverlap && ExactLevels == 0) {
            // Images do not actually overlap.
            cerr << "enblend: images do not overlap - they will be combined without blending."
                 << endl;
            cerr << "enblend: use the -l flag to force blending with a certain number of levels."
                 << endl;

            // Copy white image into black image verbatim.
            copyImageIf(srcImageRange(*(whitePair.first)),
                        maskImage(*(whitePair.second)),
                        destImage(*(blackPair.first)));
            copyImageIf(srcImageRange(*(whitePair.second)),
                        maskImage(*(whitePair.second)),
                        destImage(*(blackPair.second)));

            delete whitePair.first;
            delete whitePair.second;

            if (Verbose > VERBOSE_CHECKPOINTING_MESSAGES) {
                cout << "Checkpointing..." << endl;
            }
            // Checkpoint results.
            exportImageAlpha(srcImageRange(*(blackPair.first)),
                             srcImage(*(blackPair.second)),
                             outputImageInfo);

            blackBB = uBB;
            continue;
        }

        // Estimate memory requirements for mask generation
        if (Verbose > VERBOSE_MEMORY_ESTIMATION_MESSAGES) {
            // mem usage = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
            //           + 2*BImage*ubb + 2*UIImage*ubb
            long long inputUnionPixels = inputUnion.size().x * inputUnion.size().y;
            long long uBBPixels = uBB.size().x * uBB.size().y;
            long long bytes =
                    inputUnionPixels
                        * (2*sizeof(ImageValueType) + 2*sizeof(AlphaValueType))
                    + uBBPixels
                        * (2*sizeof(typename BImage::value_type) + 2*sizeof(typename UIImage::value_type));
            int mbytes = (int)ceil(bytes / 1000000.0);
            cout << "Estimated space required for mask generation: "
                 << mbytes
                 << "MB" << endl;
        }

        // Create the blend mask.
        bool wraparoundForMask = Wraparound && 
                (uBB.size().x == inputUnion.size().x);
        EnblendROI mBB;
        MaskType *mask = createMask<ImageType, AlphaType, MaskType>(whitePair, blackPair,
                iBB, uBB, wraparoundForMask, mBB);
        // mem usage before = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem xsection = 2*BImage*ubb + 2*UIImage*ubb
        // mem usage after = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after mask generation:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats("whiteImage", whitePair.first);
            v.printStats("whiteAlpha", whitePair.second);
            //v.printStats("mask", mask);
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        ImageExportInfo maskInfo("enblend_mask.tif");
        maskInfo.setPosition(uBB.getUL());
        exportImage(srcImageRange(*mask), maskInfo);
return;
        // Calculate ROI bounds and number of levels from mBB.
        // ROI bounds must be at least mBB but not to extend uBB.
        EnblendROI roiBB;
        unsigned int numLevels = roiBounds<MaskPyramidValueType>(inputUnion, iBB, mBB, uBB, roiBB, wraparoundForMask);
        bool wraparoundForBlend = Wraparound && (roiBB.size().x == inputUnion.size().x);
        //cout << "Wraparound = " << wraparoundForBlend << endl;

        // Estimate memory requirements for this blend iteration
        if (Verbose > VERBOSE_MEMORY_ESTIMATION_MESSAGES) {
            // Maximum utilization is when all three pyramids have been built
            // inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
            // + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType
            long long inputUnionPixels = inputUnion.size().x * inputUnion.size().y;
            long long roiBBPixels = roiBB.size().x * roiBB.size().y;
            long long bytes =
                    inputUnionPixels * (sizeof(ImageValueType) + 2*sizeof(AlphaValueType))
                    + (4/3) * roiBBPixels *
                            (sizeof(MaskPyramidValueType) + 2*sizeof(ImagePyramidValueType));
            int mbytes = (int)ceil(bytes / 1000000.0);
            cout << "Estimated space required for this blend step: "
                 << mbytes
                 << "MB" << endl;
        }

        // Create a version of roiBB relative to uBB upperleft corner.
        // This is to access roi within images of size uBB.
        // For example, the mask.
        EnblendROI roiBB_uBB;
        roiBB_uBB.setCorners(roiBB.getUL() - uBB.getUL(), roiBB.getLR() - uBB.getUL());

        // Build Gaussian pyramid from mask.
        vector<MaskPyramidType*> *maskGP = gaussianPyramid<MaskType, MaskPyramidType>(
                numLevels, wraparoundForBlend, roiBB_uBB.apply(srcImageRange(*mask)));
        // mem usage before = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem usage after = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after calculating mask pyramid:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats("whiteImage", whitePair.first);
            v.printStats("whiteAlpha", whitePair.second);
            //v.printStats("mask", mask);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats("maskGP", i, (*maskGP)[i]);
            }
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        // Now it is safe to make changes to mask image.
        // Black out the ROI in the mask.
        // Make an roiBounds relative to uBB origin.
        initImage(roiBB_uBB.apply(destImageRange(*mask)),
                NumericTraits<MaskPyramidValueType>::zero());

        // Copy pixels inside whiteBB and inside white part of mask into black image.
        // These are pixels where the white image contributes outside of the ROI.
        // We cannot modify black image inside the ROI yet because we haven't built the
        // black pyramid.
        copyImageIf(uBB.apply(srcImageRange(*(whitePair.first))),
                    maskImage(*mask),
                    uBB.apply(destImage(*(blackPair.first))));

        // We no longer need the mask.
        delete mask;
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType

        // Build Laplacian pyramid from white image.
        vector<ImagePyramidType*> *whiteLP =
                laplacianPyramid<ImageType, AlphaType, ImagePyramidType>(
                        numLevels, wraparoundForBlend,
                        roiBB.apply(srcImageRange(*(whitePair.first))),
                        roiBB.apply(maskImage(*(whitePair.second))));
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after calculating white pyramid:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats("whiteImage", whitePair.first);
            v.printStats("whiteAlpha", whitePair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats("maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats("whiteLP", i, (*whiteLP)[i]);
            }
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType

        // We no longer need the white rgb data.
        delete whitePair.first;
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType

        // Build Laplacian pyramid from black image.
        vector<ImagePyramidType*> *blackLP =
                laplacianPyramid<ImageType, AlphaType, ImagePyramidType>(
                        numLevels, wraparoundForBlend,
                        roiBB.apply(srcImageRange(*(blackPair.first))),
                        roiBB.apply(maskImage(*(blackPair.second))));
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after calculating black pyramid:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats("whiteAlpha", whitePair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats("maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats("whiteLP", i, (*whiteLP)[i]);
            }
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats("blackLP", i, (*blackLP)[i]);
            }
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif
        // Peak memory xsection is here!
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType
        //exportPyramid(blackLP, "enblend_black_lp");

        // Make the black image alpha equal to the union of the
        // white and black alpha channels.
        initImageIf(whiteBB.apply(destImageRange(*(blackPair.second))),
                whiteBB.apply(maskImage(*(whitePair.second))),
                NumericTraits<AlphaValueType>::max());

        // We no longer need the white alpha data.
        delete whitePair.second;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType

        // Blend pyramids
        blend<MaskType, MaskPyramidType, ImagePyramidType>(maskGP, whiteLP, blackLP);
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after blending pyramids:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats("maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats("whiteLP", i, (*whiteLP)[i]);
            }
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats("blackLP", i, (*blackLP)[i]);
            }
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        // delete mask pyramid
        //exportPyramid(maskGP, "enblend_mask_gp");
        for (unsigned int i = 0; i < maskGP->size(); i++) {
            delete (*maskGP)[i];
        }
        delete maskGP;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + 2*(4/3)*roiBB*ImagePyramidType

        // delete white pyramid
        //exportPyramid(whiteLP, "enblend_white_lp");
        for (unsigned int i = 0; i < whiteLP->size(); i++) {
            delete (*whiteLP)[i];
        }
        delete whiteLP;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + (4/3)*roiBB*ImagePyramidType

        //exportPyramid(blackLP, "enblend_blend_lp");

        // collapse black pyramid
        collapsePyramid(wraparoundForBlend, blackLP);
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after collapsing black pyramid:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats("blackLP", i, (*blackLP)[i]);
            }
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        // copy collapsed black pyramid into black image ROI, using black alpha mask.
        copyFromPyramidImageIf<ImageType, ImagePyramidType, AlphaType>(
                srcImageRange(*((*blackLP)[0])),
                roiBB.apply(maskImage(*(blackPair.second))),
                roiBB.apply(destImage(*(blackPair.first))));

        // delete black pyramid
        for (unsigned int i = 0; i < blackLP->size(); i++) {
            delete (*blackLP)[i];
        }
        delete blackLP;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType

        if (Verbose > VERBOSE_CHECKPOINTING_MESSAGES) {
            cout << "Checkpointing..." << endl;
        }

        // Checkpoint results.
        exportImageAlpha(srcImageRange(*(blackPair.first)),
                         srcImage(*(blackPair.second)),
                         outputImageInfo);
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after checkpointing:" << endl;
            v.printStats("blackImage", blackPair.first);
            v.printStats("blackAlpha", blackPair.second);
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif

        // Now set blackBB to uBB.
        blackBB = uBB;
    }

    delete blackPair.first;
    delete blackPair.second;

};

template <typename ImageType, typename ImagePyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion,
        VigraTrueType) {

    // ImagePyramidType::value_type is a scalar.
    enblendMain<ImageType, ImagePyramidType, ImagePyramidType>(
            imageInfoList, outputImageInfo, inputUnion);
}

template <typename ImageType, typename ImagePyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion,
        VigraFalseType) {

    // ImagePyramidType::value_type is a vector.
    typedef typename ImagePyramidType::value_type VectorType;
    #ifdef ENBLEND_CACHE_IMAGES
    typedef CachedFileImage<typename VectorType::value_type> MaskPyramidType;
    #else
    typedef BasicImage<typename VectorType::value_type> MaskPyramidType;
    #endif

    enblendMain<ImageType, MaskPyramidType, ImagePyramidType>(
            imageInfoList, outputImageInfo, inputUnion);
}

template <typename ImageType, typename ImagePyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    // This indicates if ImagePyramidType is RGB or grayscale.
    typedef typename NumericTraits<typename ImagePyramidType::value_type>::isScalar pyramid_is_scalar;
    typedef typename NumericTraits<typename ImageType::value_type>::isScalar image_is_scalar;

    // If the image is RGB, the pyramid must also be RGB.
    // If the image is scalar, the pyramid must also be scalar.
    BOOST_STATIC_ASSERT(pyramid_is_scalar::asBool == image_is_scalar::asBool);

    enblendMain<ImageType, ImagePyramidType>(
            imageInfoList, outputImageInfo, inputUnion, pyramid_is_scalar());
}

} // namespace enblend

#endif /* __ENBLEND_H__ */
