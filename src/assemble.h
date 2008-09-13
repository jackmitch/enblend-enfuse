/*
 * Copyright (C) 2004-2007 Andrew Mihal
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
#include "fixmath.h"
#include "vigra/copyimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"
#include "vigra_ext/FunctorAccessor.h"
#include "vigra_ext/impexalpha.hxx"

using std::cerr;
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
using vigra::Rect2D;
using vigra::Threshold;
using vigra::transformImage;

using vigra_ext::ReadFunctorAccessor;
using vigra_ext::WriteFunctorAccessor;

namespace enblend {

/** Write output images.
 */
template <typename ImageType, typename AlphaType>
void checkpoint(pair<ImageType*, AlphaType*> &p, ImageExportInfo &outputImageInfo) {

    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType ImagePixelComponentType;
    typedef typename AlphaType::Accessor AlphaAccessor;
    typedef typename AlphaType::PixelType AlphaPixelType;

    typedef ReadFunctorAccessor<
            Threshold<AlphaPixelType, ImagePixelComponentType>, AlphaAccessor>
            ThresholdingAccessor;
 
    ThresholdingAccessor ata(
            Threshold<AlphaPixelType, ImagePixelComponentType>(
                    NumericTraits<AlphaPixelType>::zero(),
                    NumericTraits<AlphaPixelType>::zero(),
                    AlphaTraits<ImagePixelComponentType>::max(),
                    AlphaTraits<ImagePixelComponentType>::zero()
            ),
            (p.second)->accessor());

    try {
        exportImageAlpha(srcImageRange(*(p.first)),
                         srcIter((p.second)->upperLeft(), ata),
                         outputImageInfo);
    } catch (std::exception & ) {
        // try to export without alpha channel
        exportImage(srcImageRange(*(p.first)),
                    outputImageInfo);
    }
};

template <typename DestIterator, typename DestAccessor,
          typename AlphaIterator, typename AlphaAccessor>
void import(const ImageImportInfo &info,
            const pair<DestIterator, DestAccessor> &image,
            const pair<AlphaIterator, AlphaAccessor> &alpha) {

    typedef typename DestIterator::PixelType ImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType ImagePixelComponentType;
    typedef typename AlphaIterator::PixelType AlphaPixelType;

    // Use a thresholding accessor to write to the alpha image.
    typedef WriteFunctorAccessor<
            Threshold<ImagePixelComponentType, AlphaPixelType>, AlphaAccessor>
            ThresholdingAccessor;

    // Threshold the alpha mask so that all pixels are either contributing
    // or not contributing.
    ThresholdingAccessor ata(
            Threshold<ImagePixelComponentType, AlphaPixelType>(
                    AlphaTraits<ImagePixelComponentType>::max() / 2,
                    AlphaTraits<ImagePixelComponentType>::max(),
                    AlphaTraits<AlphaPixelType>::zero(),
                    AlphaTraits<AlphaPixelType>::max()
            ),
            alpha.second);

    if (info.numExtraBands() > 0) {
        importImageAlpha(info, image, destIter(alpha.first, ata));
    }
    else {
        // Import image without alpha for enfuse. Init the alpha image to 100%.
        importImage(info, image.first, image.second);
        initImage(alpha.first,
                  (alpha.first + Diff2D(info.width(), info.height())),
                  alpha.second,
                  AlphaTraits<AlphaPixelType>::max());
    }

};

/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of ImageImportInfos.
 *  Returns an ImageImportInfo for the temporary file.
 *  memory xsection = 2 * (ImageType*inputUnion + AlphaType*inputUnion)
 */
template <typename ImageType, typename AlphaType>
pair<ImageType*, AlphaType*> assemble(list<ImageImportInfo*> &imageInfoList,
        Rect2D &inputUnion,
        Rect2D &bb) {

    typedef typename AlphaType::traverser AlphaIteratorType;
    typedef typename AlphaType::Accessor AlphaAccessor;

    // No more images to assemble?
    if (imageInfoList.empty())
        return pair<ImageType*, AlphaType*>(static_cast<ImageType*>(NULL),
                                            static_cast<AlphaType*>(NULL));

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

    Diff2D imagePos = imageInfoList.front()->getPosition();
    import(*imageInfoList.front(),
            destIter(image->upperLeft() + imagePos - inputUnion.upperLeft()),
            destIter(imageA->upperLeft() + imagePos - inputUnion.upperLeft()));
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

            import(*info, destImage(*src), destImage(*srcA));

            // Check for overlap.
            bool overlapFound = false;
            AlphaIteratorType dy =
                    imageA->upperLeft() - inputUnion.upperLeft() + info->getPosition();
            AlphaAccessor da = imageA->accessor();
            AlphaIteratorType sy = srcA->upperLeft();
            AlphaIteratorType send = srcA->lowerRight();
            AlphaAccessor sa = srcA->accessor();
            for(; sy.y < send.y; ++sy.y, ++dy.y) {
                AlphaIteratorType sx = sy;
                AlphaIteratorType dx = dy;
                for(; sx.x < send.x; ++sx.x, ++dx.x) {
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
                        destIter(image->upperLeft() - inputUnion.upperLeft() + srcPos));
                copyImageIf(srcImageRange(*srcA),
                        maskImage(*srcA),
                        destIter(imageA->upperLeft() - inputUnion.upperLeft() + srcPos));

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
    bb = unionRect();

    if (Verbose > VERBOSE_ABB_MESSAGES) {
        cout << "assembled images bounding box: " << unionRect() << endl;
    }

    return pair<ImageType*, AlphaType*>(image, imageA);

};

} // namespace enblend

#endif /* __ASSEMBLE_H__ */
