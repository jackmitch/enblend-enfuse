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

#include "assemble.h"
#include "bounds.h"
#include "mask.h"
#include "pyramid.h"

#include "common.h"
#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::initImageIf;
using vigra::NumericTraits;

namespace enblend {

template <typename ImageType, typename AlphaType, typename PyramidType, typename PyramidRGBType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    typedef typename ImageType::value_type ImageValueType;
    typedef typename AlphaType::value_type AlphaValueType;
    typedef typename PyramidType::value_type PyramidValueType;
    typedef typename PyramidRGBType::value_type PyramidRGBValueType;

    cout << "sizeof(ImageValueType) = " << sizeof(ImageValueType) << endl;
    cout << "sizeof(AlphaValueType) = " << sizeof(AlphaValueType) << endl;
    cout << "sizeof(PyramidValueType) = " << sizeof(PyramidValueType) << endl;
    cout << "sizeof(PyramidRGBValueType) = " << sizeof(PyramidRGBValueType) << endl;

    // Create the initial black image.
    EnblendROI blackBB;
    pair<ImageType*, AlphaType*> blackPair =
            assemble<ImageType, AlphaType>(imageInfoList, inputUnion, blackBB);
    exportImageAlpha(srcImageRange(*(blackPair.first)),
                     srcImage(*(blackPair.second)),
                     outputImageInfo);
    // mem xsection = up to 2*inputUnion*ImageValueType

    // Main blending loop.
    while (!imageInfoList.empty()) {

        // Create the white image.
        EnblendROI whiteBB;
        pair<ImageType*, AlphaType*> whitePair =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, whiteBB);
        ImageExportInfo whiteInfo("enblend_white.tif");
        exportImageAlpha(srcImageRange(*(whitePair.first)),
                         srcImage(*(whitePair.second)),
                         whiteInfo);
        // mem xsection = up to 2*inputUnion*ImageValueType

        // Union bounding box of whiteImage and blackImage.
        EnblendROI uBB;
        whiteBB.unite(blackBB, uBB);

        // Intersection bounding box of whiteImage and blackImage.
        EnblendROI iBB;
        bool overlap = whiteBB.intersect(blackBB, iBB);

        // Calculate ROI bounds and number of levels from iBB.
        // ROI bounds not to extend uBB.
        // FIXME consider case where overlap==false
        EnblendROI roiBB;
        unsigned int numLevels = roiBounds<PyramidValueType>(inputUnion, iBB, uBB, roiBB);

        // Create the blend mask.
        PyramidType *mask =
                createMask<AlphaType, PyramidType>(whitePair.second, blackPair.second, uBB);
        ImageExportInfo maskInfo("enblend_mask.tif");
        maskInfo.setPosition(uBB.getUL());
        exportImage(srcImageRange(*mask), maskInfo);

        //// Get the first level of the gaussian pyramid for black image.
        //PyramidType *blackGP0 = new PyramidType(roiBB.size());
        //AlphaType *blackGP0a = new AlphaType(roiBB.size());
        //copyImage(roiBB.apply(srcImageRange(*(blackPair.first))),
        //        destImage(*blackGP0));
        //copyImage(roiBB.apply(srcImageRange(*(blackPair.second))),
        //        destImage(*blackGP0a));

        // Build Gaussian pyramid from mask.
        vector<int*> *maskGP = gaussianPyramid(numLevels,
                mask->upperLeft() + (roiBB.getUL() - uBB.getUL()),
                mask->upperLeft() + (roiBB.getLR() - uBB.getUL()),
                mask->accessor());

        // Now it is safe to make changes to mask image.
        // Black out the ROI in the mask.
        // Make an roiBounds relative to uBB origin.
        initImage(mask->upperLeft() + (roiBB.getUL() - uBB.getUL()),
                mask->upperLeft() + (roiBB.getLR() - uBB.getUL()),
                mask->accessor(),
                NumericTraits<PyramidValueType>::zero());

        // Copy pixels inside whiteBB and inside white part of mask into black image.
        // These are pixels where the white image contributes outside of the ROI.
        // We cannot modify black image inside the ROI yet because we haven't built the
        // black pyramid.
        copyImageIf(uBB.apply(srcImageRange(*(whitePair.first))),
                    maskImage(*mask),
                    uBB.apply(destImage(*(blackPair.first))));

        // We no longer need the mask.
        delete mask;

        // Build Laplacian pyramid from white image.
        //vector<PyramidType*> *whiteLP = whitePair.first, whitePair.second

        // We no longer need the white rgb data.
        delete whitePair.first;

        // Build Laplacian pyramid from black image.

        // Make the black image alpha equal to the union of the
        // white and black alpha channels.
        initImageIf(whiteBB.apply(destImageRange(*(blackPair.second))),
                whiteBB.apply(maskImage(*(whitePair.second))),
                NumericTraits<AlphaValueType>::max());

        // We no longer need the white alpha data.
        delete whitePair.second;

        // Blend pyramids
        // delete mask pyramid
        // delete white pyramid
        // collapse black pyramid

        // copy collapsed black pyramid into black image ROI, using black alpha mask.

        // delete black pyramid

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

} // namespace enblend

#endif /* __ENBLEND_H__ */
