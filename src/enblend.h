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

#include "common.h"
#include "vigra/impex.hxx"

using std::cout;
using std::endl;
using std::list;

using vigra::ImageExportInfo;
using vigra::ImageImportInfo;

namespace enblend {

template <typename ImageType, typename AlphaType, typename MaskType, typename PyramidType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    typedef typename ImageType::value_type ImageValueType;
    typedef typename AlphaType::value_type AlphaValueType;
    typedef typename MaskType::value_type MaskValueType;
    typedef typename PyramidType::value_type PyramidValueType;

    cout << "sizeof(ImageValueType) = " << sizeof(ImageValueType) << endl;
    cout << "sizeof(AlphaValueType) = " << sizeof(AlphaValueType) << endl;
    cout << "sizeof(MaskValueType) = " << sizeof(MaskValueType) << endl;
    cout << "sizeof(PyramidValueType) = " << sizeof(PyramidValueType) << endl;

    // Create the initial white image. This should go into outputImageInfo.
    // This way the output file always contains the results of the last completed
    // blend iteration.
    EnblendROI whiteBB;
    ImageImportInfo *whiteImageInfo =
            assemble<ImageType, AlphaType>(imageInfoList, inputUnion, whiteBB,
                    &outputImageInfo);

    // Main blending loop.
    while (!imageInfoList.empty()) {

        // Create the black image.
        EnblendROI blackBB;
        ImageImportInfo *blackImageInfo =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, blackBB);

        // Union bounding box of whiteImage and blackImage.
        EnblendROI uBB;
        whiteBB.unite(blackBB, uBB);

        // Intersection bounding box of whiteImage and blackImage.
        EnblendROI iBB;
        bool overlap = whiteBB.intersect(blackBB, iBB);

        // Calculate ROI bounds and number of levels from iBB.
        // ROI bounds not to extend uBB.
        EnblendROI roiBB;
        unsigned int numLevels =
                roiBounds<PyramidValueType, MaskValueType>(inputUnion, iBB, uBB, roiBB);

        // Create the blend mask and the union mask.
        ImageImportInfo *maskImageInfo =
                mask<ImageType, AlphaType, MaskType>(whiteImageInfo, blackImageInfo,
                        inputUnion, uBB, iBB, overlap);

        // Copy pixels inside blackBB and outside ROI into white image.

        // Build Laplacian pyramid from blackImage.
        //ImageImportInfo *blackLPInfo = 

        // Now we no longer need the blackImage temp file, delete it.
        unlink(blackImageInfo->getFileName());
        delete blackImageInfo;

        // Build Gaussian pyramid from mask.
        //ImageImportInfo *maskGPInfo =

        // We no longer need the maskImage temp file, delete it.
        unlink(maskImageInfo->getFileName());
        delete maskImageInfo;

        // Build Laplacian pyramid from whiteImage
        //vector<PyramidType*> *whiteLP =

        // Blend pyramids

        // We no longer need the blackLPInfo file or the maskGPInfo files.
        // Delete them.
        //unlink(blackLPInfo->getFileName());
        //delete blackLPInfo;
        //unlink(maskGPInfo->getFileName());
        //delete maskGPInfo;

        // Collapse result back into whiteImageInfo.
        //collapsePyramid

        // Copy result into whiteImageFile using unionMaskFile as a template.

        // Done with unionMaskFile.

        // done with whiteLP.
        //for (unsigned int i = 0; i < whiteLP->size(); i++) {
        //    delete whiteLP[i];
        //}
        //delete whiteLP;

        // Now set whiteBB to uBB.
        whiteBB = uBB;
    }

    // We no longer need whiteImageInfo.
    // Final results are already in
    delete whiteImageInfo;
};

} // namespace enblend

#endif /* __ENBLEND_H__ */
