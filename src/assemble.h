/*
 * Copyright (C) 2004-2005 Andrew Mihal
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
#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>
#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#include "common.h"
#include "vigra/copyimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"
#include "vigra_ext/FunctorAccessor.h"
#include "vigra_ext/impexalpha.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::copyImageIf;
using vigra::Diff2D;
using vigra::exportImageAlpha;
using vigra::FindBoundingRectangle;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::inspectImageIf;
using vigra::NumericTraits;
using vigra::Threshold;
using vigra::transformImage;

using vigra_ext::WriteFunctorAccessor;

namespace enblend {

/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of ImageImportInfos.
 *  Returns an ImageImportInfo for the temporary file.
 *  memory xsection = 2 * (ImageType*inputUnion + AlphaType*inputUnion)
 */
template <typename ImageType, typename AlphaType>
pair<ImageType*, AlphaType*> assemble(list<ImageImportInfo*> &imageInfoList,
        EnblendROI &inputUnion,
        EnblendROI &bb) {

    typedef typename AlphaType::PixelType AlphaPixelType;
    typedef typename AlphaType::traverser AlphaIteratorType;
    typedef typename AlphaType::Accessor AlphaAccessor;

    // No more images to assemble?
    if (imageInfoList.empty()) return pair<ImageType*, AlphaType*>(NULL, NULL);

    // Create an image to assemble input images into.
    ImageType *image = new ImageType(inputUnion.size());
    AlphaType *imageA = new AlphaType(inputUnion.size());

    if (Verbose > VERBOSE_ASSEMBLE_MESSAGES) {
        if (OneAtATime) {
            cout << "Loading next image: "
                 << imageInfoList.front()->getFileName()
                 << endl;
        } else {
            cout << "Combining non-overlapping images: "
                 << imageInfoList.front()->getFileName();
            cout.flush();
        }
    }

    // Load the first image into the destination.
    // Use a thresholding accessor to write to the alpha image.
    typedef WriteFunctorAccessor<
            Threshold<AlphaPixelType, AlphaPixelType>, AlphaAccessor>
            ThresholdingAccessor;
    // Threshold the alpha mask so that all pixels are either contributing
    // or not contributing.
    ThresholdingAccessor imageATA(
            Threshold<AlphaPixelType, AlphaPixelType>(
                    NumericTraits<AlphaPixelType>::max() / 2,
                    NumericTraits<AlphaPixelType>::max(),
                    NumericTraits<AlphaPixelType>::zero(),
                    NumericTraits<AlphaPixelType>::max()
            ),
            imageA->accessor());

    Diff2D imagePos = imageInfoList.front()->getPosition();
    importImageAlpha(*imageInfoList.front(),
            destIter(image->upperLeft() + imagePos - inputUnion.getUL()),
            destIter(imageA->upperLeft() + imagePos - inputUnion.getUL(), imageATA));
    imageInfoList.erase(imageInfoList.begin());

    if (!OneAtATime) {
        // Attempt to assemble additional non-overlapping images.

        // List of ImageImportInfos we decide to assemble.
        list<list<ImageImportInfo*>::iterator> toBeRemoved;

        list<ImageImportInfo*>::iterator i;
        for (i = imageInfoList.begin(); i != imageInfoList.end(); i++) {
            ImageImportInfo *info = *i;

            // Load the next image.
            ImageType *src = new ImageType(info->size());
            AlphaType *srcA = new AlphaType(info->size());

            // Use a thresholding accessor to write to the alpha image.
            ThresholdingAccessor srcATA(
                    Threshold<AlphaPixelType, AlphaPixelType>(
                            NumericTraits<AlphaPixelType>::max() / 2,
                            NumericTraits<AlphaPixelType>::max(),
                            NumericTraits<AlphaPixelType>::zero(),
                            NumericTraits<AlphaPixelType>::max()
                    ),
                    srcA->accessor());
            importImageAlpha(*info, destImage(*src), destImage(*srcA, srcATA));

            // Check for overlap.
            bool overlapFound = false;
            AlphaIteratorType dy =
                    imageA->upperLeft() - inputUnion.getUL() + info->getPosition();
            AlphaAccessor da = imageA->accessor();
            AlphaIteratorType sy = srcA->upperLeft();
            AlphaIteratorType send = srcA->lowerRight();
            AlphaAccessor sa = srcA->accessor();
            for(; sy.y != send.y; ++sy.y, ++dy.y) {
                AlphaIteratorType sx = sy;
                AlphaIteratorType dx = dy;
                for(; sx.x != send.x; ++sx.x, ++dx.x) {
                    if (sa(sx) && da(dx)) {
                        overlapFound = true;
                        break;
                    }
                }
                if (overlapFound) break;
            }

            if (!overlapFound) {
                // Copy src and srcA into image and imageA.

                if (Verbose > VERBOSE_ASSEMBLE_MESSAGES) {
                    cout << " " << info->getFileName();
                    cout.flush();
                }

                Diff2D srcPos = info->getPosition();
                copyImageIf(srcImageRange(*src),
                        maskImage(*srcA),
                        destIter(image->upperLeft() - inputUnion.getUL() + srcPos));
                copyImageIf(srcImageRange(*srcA),
                        maskImage(*srcA),
                        destIter(imageA->upperLeft() - inputUnion.getUL() + srcPos));

                // Remove info from list later.
                toBeRemoved.push_back(i);
            }

            delete src;
            delete srcA;
        }
        
        // Erase the ImageImportInfos we used.
        list<list<ImageImportInfo*>::iterator>::iterator r;
        for (r = toBeRemoved.begin(); r != toBeRemoved.end(); r++) {
            imageInfoList.erase(*r);
        }
    }

    if (Verbose > VERBOSE_ASSEMBLE_MESSAGES && !OneAtATime) cout << endl;

    // Calculate bounding box of image.
    FindBoundingRectangle unionRect;
    inspectImageIf(srcIterRange(Diff2D(), Diff2D() + image->size()),
            srcImage(*imageA), unionRect);
    bb.setCorners(unionRect.upperLeft, unionRect.lowerRight);

    if (Verbose > VERBOSE_ABB_MESSAGES) {
        cout << "assembled images bounding box: ("
             << unionRect.upperLeft.x
             << ", "
             << unionRect.upperLeft.y
             << ") -> ("
             << unionRect.lowerRight.x
             << ", "
             << unionRect.lowerRight.y
             << ")" << endl;
    }

    return pair<ImageType*, AlphaType*>(image, imageA);

};

} // namespace enblend

#endif /* __ASSEMBLE_H__ */
