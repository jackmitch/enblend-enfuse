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
#include "vigra/stdcachedfileimage.hxx"

using vigra::BCFImage;
using vigra::BImage;
using vigra::exportImage;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
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

    //// Read mask from a file instead of calculating it.
    //MaskType *fileMask = new MaskType(uBB.size());
    //ImageImportInfo fileMaskInfo("enblend_mask.tif");
    //importImage(fileMaskInfo, destImage(*fileMask));
    //return fileMask;

    // Mask initializer pixel values:
    // 0 = outside both black and white image, or inside both images.
    // 1 = inside white image only.
    // 255 = inside black image only.
    #ifdef ENBLEND_CACHE_IMAGES
    BCFImage *maskInit = new BCFImage(uBB.size());
    #else
    BImage *maskInit = new BImage(uBB.size());
    #endif
    // mem xsection = BImage*ubb

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

    // Mask transform replaces 0 areas with either 1 or 255.
    #ifdef ENBLEND_CACHE_IMAGES
    BCFImage *maskTransform = new BCFImage(uBB.size());
    #else
    BImage *maskTransform = new BImage(uBB.size());
    #endif
    // mem xsection = 2*BImage*ubb
    nearestFeatureTransform(wraparound,
            srcImageRange(*maskInit),
            destImage(*maskTransform));
    // mem xsection = 2*BImage*ubb + 2*UIImage*ubb

    delete maskInit;
    // mem xsection = BImage*ubb

    MaskType *mask = new MaskType(uBB.size());
    // mem xsection = BImage*ubb + MaskType*ubb

    // Dump maskTransform into mask
    // maskTransform = 1, then mask = max value (white image)
    // maskTransform != 1, then mask = zero - (black image)
    transformImage(srcImageRange(*maskTransform),
            destImage(*mask),
            ifThenElse(Arg1() == Param(1),
                    Param(NumericTraits<MaskPixelType>::max()),
                    Param(NumericTraits<MaskPixelType>::zero())));

    delete maskTransform;
    // mem xsection = MaskType*ubb

    return mask;
}

} // namespace enblend

#endif /* __MASK_H__ */
