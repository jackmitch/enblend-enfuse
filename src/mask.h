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

#include "vigra/error.hxx"
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
        bool wraparound,
        EnblendROI &mBB) {

    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

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

    // Find the bounding box of the mask transition line and put it in mBB.
    MaskIteratorType firstMulticolorColumn = mask->lowerRight();
    MaskIteratorType lastMulticolorColumn = mask->upperLeft();
    MaskIteratorType firstMulticolorRow = mask->lowerRight();
    MaskIteratorType lastMulticolorRow = mask->upperLeft();

    MaskIteratorType myPrev = mask->upperLeft();
    MaskIteratorType my = mask->upperLeft() + Diff2D(0,1);
    MaskIteratorType mend = mask->lowerRight();
    for (; my.y != mend.y; ++my.y, ++myPrev.y) {
        MaskIteratorType mxLeft = my;
        MaskIteratorType mx = my + Diff2D(1,0);
        MaskIteratorType mxUpLeft = myPrev;
        MaskIteratorType mxUp = myPrev + Diff2D(1,0);

        if (*mxUpLeft != *mxLeft) {
            // Transition line is between mxUpLeft and mxLeft.
            if (firstMulticolorRow.y > mxUpLeft.y) firstMulticolorRow = mxUpLeft;
            if (lastMulticolorRow.y < mxLeft.y) lastMulticolorRow = mxLeft;
        }

        for (; mx.x != mend.x; ++mx.x, ++mxLeft.x, ++mxUp.x) {
            if (*mxLeft != *mx || *mxUp != *mx) {
                // Transition line is between mxLeft and mx and between mx and mxUp
                if (firstMulticolorColumn.x > mxLeft.x) firstMulticolorColumn = mxLeft;
                if (lastMulticolorColumn.x < mx.x) lastMulticolorColumn = mx;
                if (firstMulticolorRow.y > mxUp.y) firstMulticolorRow = mxUp;
                if (lastMulticolorRow.y < mx.y) lastMulticolorRow = mx;
            }
        }
    }

    // Check that mBB is well-defined.
    if ((firstMulticolorColumn.x >= lastMulticolorColumn.x)
            || (firstMulticolorRow.y >= lastMulticolorRow.y)) {
        // No transition pixels were found in the mask at all.
        // This means that one image has no contribution.
        vigra_fail("Mask transition line bounding box undefined.");
    }

    // Move mBB lower right corner out one pixel, per VIGRA convention.
    ++lastMulticolorColumn.x;
    ++lastMulticolorRow.y;

    // mBB is defined relative to the inputUnion origin.
    mBB.setCorners(uBB.getUL() + Diff2D(firstMulticolorColumn.x, firstMulticolorRow.y),
            uBB.getUL() + Diff2D(lastMulticolorColumn.x, lastMulticolorRow.y));

    if (Verbose > VERBOSE_ROIBB_SIZE_MESSAGES) {
        cout << "Mask transition line bounding box: ("
             << mBB.getUL().x
             << ", "
             << mBB.getUL().y
             << ") -> ("
             << mBB.getLR().x
             << ", "
             << mBB.getLR().y
             << ")" << endl;
    }

    return mask;
}

} // namespace enblend

#endif /* __MASK_H__ */
