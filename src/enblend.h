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

//#include <boost/static_assert.hpp>
#include <iostream>
#include <list>
#include <stdio.h>

#include "assemble.h"
#include "blend.h"
#include "bounds.h"
#include "mask.h"
#include "pyramid.h"

#include "common.h"
#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::BasicImage;
using vigra::BImage;
using vigra::FindMinMax;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::initImageIf;
using vigra::initImage;
using vigra::inspectImage;
using vigra::NumericTraits;
using vigra::VigraTrueType;
using vigra::VigraFalseType;

namespace enblend {

template <typename ImageType, typename MaskPyramidType, typename ImagePyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    typedef typename ImageType::value_type ImageValueType;
    typedef BImage MaskType;
    typedef BImage AlphaType;
    typedef typename AlphaType::value_type AlphaValueType;
    typedef typename MaskPyramidType::value_type MaskPyramidValueType;
    typedef typename ImagePyramidType::value_type ImagePyramidValueType;

    //cout << "sizeof(ImageValueType) = " << sizeof(ImageValueType) << endl;
    //cout << "sizeof(AlphaValueType) = " << sizeof(AlphaValueType) << endl;
    //cout << "sizeof(MaskPyramidValueType) = " << sizeof(MaskPyramidValueType) << endl;
    //cout << "sizeof(ImagePyramidValueType) = " << sizeof(ImagePyramidValueType) << endl;

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

    // Main blending loop.
    while (!imageInfoList.empty()) {

        // Create the white image.
        EnblendROI whiteBB;
        pair<ImageType*, AlphaType*> whitePair =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, whiteBB);
        // mem usage before = inputUnion*ImageValueType + inputUnion*AlphaValueType
        // mem xsection = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType

        //ImageExportInfo whiteInfo("enblend_white.tif");
        //exportImageAlpha(srcImageRange(*(whitePair.first)),
        //                 srcImage(*(whitePair.second)),
        //                 whiteInfo);

        // Union bounding box of whiteImage and blackImage.
        EnblendROI uBB;
        whiteBB.unite(blackBB, uBB);

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

        // Intersection bounding box of whiteImage and blackImage.
        EnblendROI iBB;
        whiteBB.intersect(blackBB, iBB);

        // If images are disjoint, make the iBB the full size of the uBB.
        if (overlap == NoOverlap) {
            iBB = uBB;
        }

        // Calculate ROI bounds and number of levels from iBB.
        // ROI bounds not to extend uBB.
        EnblendROI roiBB;
        unsigned int numLevels = roiBounds<MaskPyramidValueType>(inputUnion, iBB, uBB, roiBB);
        bool wraparoundThisIteration = Wraparound && (roiBB.size().x == inputUnion.size().x);

        // Estimate memory requirements for this blend iteration
        if (Verbose > 0) {
            // Maximum utilization is when all three pyramids have been built
            // inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
            // + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType
            long long inputUnionPixels = inputUnion.size().x * inputUnion.size().y;
            long long roiBBPixels = roiBB.size().x * roiBB.size().y;
            long long bytes =
                    inputUnionPixels * (sizeof(ImageValueType) + 2*sizeof(AlphaValueType))
                    + (4/3) * roiBBPixels *
                            (sizeof(MaskPyramidValueType) + 2*sizeof(ImagePyramidValueType));
            long long mbytes = bytes / 1000000;
            cout << "Estimated space required for this blend step: "
                 << mbytes
                 << "MB" << endl;
        }

        // Create a version of roiBB relative to uBB upperleft corner.
        // This is to access roi within images of size uBB.
        // For example, the mask.
        EnblendROI roiBB_uBB;
        roiBB_uBB.setCorners(roiBB.getUL() - uBB.getUL(), roiBB.getLR() - uBB.getUL());

        // Create the blend mask.
        MaskType *mask = createMask<AlphaType, MaskType>(whitePair.second, blackPair.second,
                uBB, wraparoundThisIteration);
        // mem usage before = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem xsection = 2*BImage*ubb + 2*UIImage*ubb
        // mem usage after = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType

        //ImageExportInfo maskInfo("enblend_mask.tif");
        //maskInfo.setPosition(uBB.getUL());
        //exportImage(srcImageRange(*mask), maskInfo);

        // Build Gaussian pyramid from mask.
        vector<MaskPyramidType*> *maskGP = gaussianPyramid<MaskType, MaskPyramidType>(
                numLevels, wraparoundThisIteration, roiBB_uBB.apply(srcImageRange(*mask)));
        // mem usage before = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem usage after = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType

        //for (unsigned int i = 0; i < (numLevels - 1); i++) {
        //    // Clear all levels except last.
        //    //ImagePyramidValueType zero = NumericTraits<ImagePyramidValueType>::zero();
        //    initImage(destImageRange(*((*maskGP)[i])), NumericTraits<MaskPyramidValueType>::zero());
        //}
        //collapsePyramid(wraparoundThisIteration, maskGP);
        //for (unsigned int i = 0; i < numLevels; i++) {
        //    char filenameBuf[512];
        //    snprintf(filenameBuf, 512, "enblend_mask_gp%04u.tif", i);
        //    cout << filenameBuf << endl;
        //    ImageExportInfo gpInfo(filenameBuf);
        //    vigra::USImage usLabPyramid((*maskGP)[i]->width(), (*maskGP)[i]->height());
        //    vigra::transformImage(srcImageRange(*((*maskGP)[i])), destImage(usLabPyramid),
        //            vigra::linearRangeMapping(NumericTraits<MaskPyramidValueType>::min(),
        //                                 NumericTraits<MaskPyramidValueType>::max(),
        //                                 NumericTraits<unsigned short>::min(),
        //                                 NumericTraits<unsigned short>::max()));
        //    exportImage(srcImageRange(usLabPyramid), gpInfo);
        //}

        //for (unsigned int i = 0; i < numLevels; i++) {
        //    char filenameBuf[512];
        //    snprintf(filenameBuf, 512, "enblend_mask_gp%04u.tif", i);
        //    cout << filenameBuf << endl;
        //    ImageExportInfo mgpInfo(filenameBuf);
        //    exportImage(srcImageRange(*((*maskGP)[i])), mgpInfo);
        //}

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
                        numLevels, wraparoundThisIteration,
                        roiBB.apply(srcImageRange(*(whitePair.first))),
                        roiBB.apply(maskImage(*(whitePair.second))));
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType

        //for (unsigned int i = 0; i < numLevels; i++) {
        //    char filenameBuf[512];
        //    snprintf(filenameBuf, 512, "enblend_white_lp%04u.tif", i);
        //    cout << filenameBuf << endl;
        //    ImageExportInfo wlpInfo(filenameBuf);
        //    exportImage(srcImageRange(*((*whiteLP)[i])), wlpInfo);
        //}

        // We no longer need the white rgb data.
        delete whitePair.first;
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType

        // Build Laplacian pyramid from black image.
        vector<ImagePyramidType*> *blackLP =
                laplacianPyramid<ImageType, AlphaType, ImagePyramidType>(
                        numLevels, wraparoundThisIteration,
                        roiBB.apply(srcImageRange(*(blackPair.first))),
                        roiBB.apply(maskImage(*(blackPair.second))));
        // Peak memory xsection is here!
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType

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

        // delete mask pyramid
        for (unsigned int i = 0; i < maskGP->size(); i++) {
            delete (*maskGP)[i];
        }
        delete maskGP;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + 2*(4/3)*roiBB*ImagePyramidType

        // delete white pyramid
        for (unsigned int i = 0; i < whiteLP->size(); i++) {
            delete (*whiteLP)[i];
        }
        delete whiteLP;
        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + (4/3)*roiBB*ImagePyramidType

        //for (unsigned int i = 0; i < (numLevels - 1); i++) {
        //    // Clear all levels except last.
        //    //ImagePyramidValueType zero = NumericTraits<ImagePyramidValueType>::zero();
        //    initImage(destImageRange(*((*blackLP)[i])), NumericTraits<ImagePyramidValueType>::zero());
        //}
        //for (unsigned int i = 0; i < numLevels; i++) {
        //    char filenameBuf[512];
        //    snprintf(filenameBuf, 512, "enblend_blend_lp%04u.tif", i);
        //    cout << filenameBuf << endl;
        //    ImageExportInfo lpInfo(filenameBuf);
        //    vigra::USRGBImage usLabPyramid((*blackLP)[i]->width(), (*blackLP)[i]->height());
        //    vigra::transformImage(srcImageRange(*((*blackLP)[i])), destImage(usLabPyramid),
        //            vigra::linearRangeMapping(ImagePyramidValueType(NumericTraits<MaskPyramidValueType>::min()),
        //                                 ImagePyramidValueType(NumericTraits<MaskPyramidValueType>::max()),
        //                                 typename vigra::USRGBImage::value_type(NumericTraits<unsigned short>::min()),
        //                                 typename vigra::USRGBImage::value_type(NumericTraits<unsigned short>::max())));
        //    exportImage(srcImageRange(usLabPyramid), lpInfo);
        //}

        // collapse black pyramid
        collapsePyramid(wraparoundThisIteration, blackLP);

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

        if (Verbose > 0) {
            cout << "Checkpointing..." << endl;
        }

        // Checkpoint results.
        exportImageAlpha(srcImageRange(*(blackPair.first)),
                         srcImage(*(blackPair.second)),
                         outputImageInfo);

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
    //typedef BasicImage<typename ImagePyramidType::value_type> MaskPyramidType;

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
    typedef BasicImage<typename VectorType::value_type> MaskPyramidType;

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
