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
#ifndef __MASK_H__
#define __MASK_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "common.h"
#include "nearest.h"

#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/impexalpha.hxx"
#include "vigra/initimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"

using vigra::BImage;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::exportImage;
using vigra::initImageIf;
using vigra::linearIntensityTransform;
using vigra::NumericTraits;
using vigra::transformImage;
using vigra::transformImageIf;

using vigra::functor::Arg1;
using vigra::functor::ifThenElse;
using vigra::functor::Param;

namespace enblend {

/** Calculate a blending mask between whiteImage and blackImage.
 */
template <typename AlphaType, typename MaskType>
MaskType *createMask(AlphaType *whiteAlpha,
        AlphaType *blackAlpha,
        EnblendROI &uBB,
        bool wraparound) {

    typedef typename MaskType::PixelType MaskPixelType;

    // Mask initializer pixel values:
    // 0 = outside both black and white image, or inside both images.
    // 1 = inside white image only.
    // 255 = inside black image only.
    BImage *maskInit = new BImage(uBB.size());

    // Set maskInit = 1 at all pixels where whiteImage contributes.
    initImageIf(destImageRange(*maskInit),
            maskIter(whiteAlpha->upperLeft() + uBB.getUL()),
            1);

    // maskInit = maskInit + 255 at all pixels where blackImage contributes.
    // if whiteImage also contributes, this wraps around to zero.
    transformImageIf(srcImageRange(*maskInit),
            maskIter(blackAlpha->upperLeft() + uBB.getUL()),
            destImage(*maskInit),
            linearIntensityTransform(1, 255));

    // Do mask transform here.
    // replaces 0 areas with either 1 or 255.
    BImage *maskTransform = new BImage(uBB.size());
    nearestFeatureTransform(wraparound,
            srcImageRange(*maskInit),
            destImage(*maskTransform));

    delete maskInit;

    MaskType *mask = new MaskType(uBB.size());

    // Dump maskTransform into mask
    // maskTransform = 1, then mask = max value (white image)
    // maskTransform != 1, then mask = zero - (black image)
    transformImage(srcImageRange(*maskTransform),
            destImage(*mask),
            ifThenElse(Arg1() == Param(1),
                    Param(NumericTraits<MaskPixelType>::max()),
                    Param(NumericTraits<MaskPixelType>::zero())));

    delete maskTransform;

    //char tmpFilename[] = "enblend_mask.tif";
    //ImageImportInfo maskImageInfo(tmpFilename);
    //importImage(maskImageInfo, destImage(*mask));

    return mask;

    ////char tmpFilename[] = ".enblend_mask_XXXXXX";
    //char tmpFilename[] = "enblend_mask.tif";
    ////int tmpFD = mkstemp(tmpFilename);
    //ImageExportInfo maskImageInfo(tmpFilename);
    //maskImageInfo.setPosition(uBB.getUL());
    //maskImageInfo.setFileType("TIFF");
    //exportImage(srcImageRange(mask), maskImageInfo);
    ////close(tmpFD);

    //return new ImageImportInfo(tmpFilename);
}

} // namespace enblend

#endif /* __MASK_H__ */
    // mem xsection = BImage*uBB + ImageType*inputUnion + AlphaType*inputUnion
    // mem xsection = 2*BImage*uBB + 2*UIImage*uBB
