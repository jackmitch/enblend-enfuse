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
#include "vigra/imageinfo.hxx"
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
using vigra::linearRangeMapping;

using vigra_ext::ReadFunctorAccessor;
using vigra_ext::WriteFunctorAccessor;

namespace enblend {

template <typename T> struct IntegralSelect {typedef T Result;};
template <> struct IntegralSelect<float> {typedef vigra::UInt32 Result;};
template <> struct IntegralSelect<double> {typedef vigra::UInt32 Result;};
template <> struct IntegralSelect<vigra::RGBValue<float> > {
    typedef vigra::RGBValue<vigra::UInt32> Result;
};
template <> struct IntegralSelect<vigra::RGBValue<double> > {
    typedef vigra::RGBValue<vigra::UInt32> Result;
};


template <typename ImageType, typename AlphaType, typename Accessor>
void
exportImagePreferablyWithAlpha(const ImageType* image,
                               const AlphaType* mask,
                               const Accessor& accessor,
                               const ImageExportInfo& outputImageInfo)
{
    try {
        exportImageAlpha(srcImageRange(*image),
                         srcIter(mask->upperLeft(), accessor),
                         outputImageInfo);
    } catch (std::exception&) {
        // Oh well, there is no alpha-channel.  So we export without it.
        exportImage(srcImageRange(*image), outputImageInfo);
    }
    OutputIsValid = true;
}


/** Write output images.
 */
template <typename ImageType, typename AlphaType>
void
checkpoint(const pair<ImageType*, AlphaType*>& p,
           const ImageExportInfo& outputImageInfo)
{
    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType
        ImagePixelComponentType;
    typedef typename AlphaType::Accessor AlphaAccessor;
    typedef typename AlphaType::PixelType AlphaPixelType;

    ImageType* image = p.first;
    AlphaType* mask = p.second;

    ReadFunctorAccessor<Threshold<AlphaPixelType, ImagePixelComponentType>, AlphaAccessor>
        ata(Threshold<AlphaPixelType, ImagePixelComponentType>
            (AlphaTraits<AlphaPixelType>::zero(),
             AlphaTraits<AlphaPixelType>::zero(),
             AlphaTraits<ImagePixelComponentType>::max(),
             AlphaTraits<ImagePixelComponentType>::zero()),
            mask->accessor());

    const pair<double, double> outputRange =
        enblend::rangeOfPixelType(outputImageInfo.getPixelType());
    const ImagePixelComponentType inputMin =
        NumericTraits<ImagePixelComponentType>::isIntegral::asBool ?
        NumericTraits<ImagePixelComponentType>::min() :
        0.0;
    const ImagePixelComponentType inputMax =
        NumericTraits<ImagePixelComponentType>::isIntegral::asBool ?
        NumericTraits<ImagePixelComponentType>::max() :
        1.0;
#ifdef DEBUG
    cerr << "+ checkpoint: input range:  ("
         << static_cast<double>(inputMin) << ", "
         << static_cast<double>(inputMax) << ")\n"
         << "+ checkpoint: output range: ("
         << outputRange.first << ", " << outputRange.second << ")" << endl;
#endif

    if (inputMin <= outputRange.first && inputMax >= outputRange.second) {
        if (inputMin == outputRange.first && inputMax == outputRange.second) {
            // No rescaling is necessary here: We skip the redundant
            // transformation of the input to the output range and
            // leave the channel width alone.
            ;
#ifdef DEBUG
            cerr << "+ checkpoint: leaving channel width alone" << endl;
#endif
            exportImagePreferablyWithAlpha(image, mask, ata, outputImageInfo);
        } else {
            cerr << command
                 << ": info: narrowing channel width for output as \""
                 << toLowercase(outputImageInfo.getPixelType()) << "\"" << endl;

            ImageType lowDepthImage(image->width(), image->height());
            transformImageMP(srcImageRange(*image),
                             destImage(lowDepthImage),
                             linearRangeMapping(ImagePixelType(inputMin),
                                                ImagePixelType(inputMax),
                                                ImagePixelType(outputRange.first),
                                                ImagePixelType(outputRange.second)));
            exportImagePreferablyWithAlpha(&lowDepthImage, mask, ata, outputImageInfo);
        }
    } else {
        cerr << command
             << ": info: rescaling floating-point data for output as \""
             << toLowercase(outputImageInfo.getPixelType()) << "\"" << endl;

        typedef typename IntegralSelect<ImagePixelType>::Result IntegralPixelType;
        typedef typename EnblendNumericTraits<IntegralPixelType>::ImagePixelComponentType
            IntegralPixelComponentType;

        IMAGETYPE<IntegralPixelType> integralImage(image->width(), image->height());

        ReadFunctorAccessor<Threshold<AlphaPixelType, IntegralPixelComponentType>, AlphaAccessor>
            ata(Threshold<AlphaPixelType, IntegralPixelComponentType>
                (AlphaTraits<AlphaPixelType>::zero(),
                 AlphaTraits<AlphaPixelType>::zero(),
                 AlphaTraits<IntegralPixelComponentType>::max(),
                 AlphaTraits<IntegralPixelComponentType>::zero()),
                mask->accessor());

        transformImageMP(srcImageRange(*image),
                         destImage(integralImage),
                         linearRangeMapping(ImagePixelType(inputMin),
                                            ImagePixelType(inputMax),
                                            IntegralPixelType(outputRange.first),
                                            IntegralPixelType(outputRange.second)));
        exportImagePreferablyWithAlpha(&integralImage, mask, ata, outputImageInfo);
    }
}


template <typename DestIterator, typename DestAccessor,
          typename AlphaIterator, typename AlphaAccessor>
void
import(const ImageImportInfo& info,
       const pair<DestIterator, DestAccessor>& image,
       const pair<AlphaIterator, AlphaAccessor>& alpha)
{
    typedef typename DestIterator::PixelType ImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType
        ImagePixelComponentType;
    typedef typename AlphaIterator::PixelType AlphaPixelType;

    const Diff2D extent = Diff2D(info.width(), info.height());
    const std::string pixelType = info.getPixelType();
    const range_t inputRange = enblend::rangeOfPixelType(pixelType);

    if (info.numExtraBands() > 0) {
        // Threshold the alpha mask so that all pixels are either
        // contributing or not contributing.
        WriteFunctorAccessor<Threshold<ImagePixelComponentType, AlphaPixelType>, AlphaAccessor>
            ata(Threshold<ImagePixelComponentType, AlphaPixelType>
                (inputRange.second / 2,
                 inputRange.second,
                 AlphaTraits<AlphaPixelType>::zero(),
                 AlphaTraits<AlphaPixelType>::max()),
                alpha.second);

        importImageAlpha(info, image, destIter(alpha.first, ata));
    } else {
        // Import image without alpha.  Initialize the alpha image to 100%.
        importImage(info, image.first, image.second);
        initImage(srcIterRange(alpha.first, alpha.first + extent, alpha.second),
                  AlphaTraits<AlphaPixelType>::max());
    }

    // Performance Optimization: Transform only if ranges do not
    // match.
    const double min =
        NumericTraits<ImagePixelComponentType>::isIntegral::asBool ?
        static_cast<double>(NumericTraits<ImagePixelComponentType>::min()) :
        0.0;
    const double max =
        NumericTraits<ImagePixelComponentType>::isIntegral::asBool ?
        static_cast<double>(NumericTraits<ImagePixelComponentType>::max()) :
        1.0;
    if (inputRange.first != min || inputRange.second != max) {
        transformImageMP(srcIterRange(image.first, image.first + extent, image.second),
                         destIter(image.first, image.second),
                         linearRangeMapping(ImagePixelType(inputRange.first),
                                            ImagePixelType(inputRange.second),
                                            ImagePixelType(min),
                                            ImagePixelType(max)));
    }
}


/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of ImageImportInfos.
 *  Returns an ImageImportInfo for the temporary file.
 *  memory xsection = 2 * (ImageType*inputUnion + AlphaType*inputUnion)
 */
template <typename ImageType, typename AlphaType>
pair<ImageType*, AlphaType*>
assemble(list<ImageImportInfo*>& imageInfoList, Rect2D& inputUnion, Rect2D& bb)
{
    typedef typename AlphaType::traverser AlphaIteratorType;
    typedef typename AlphaType::Accessor AlphaAccessor;

    // No more images to assemble?
    if (imageInfoList.empty()) {
        return pair<ImageType*, AlphaType*>(static_cast<ImageType*>(NULL),
                                            static_cast<AlphaType*>(NULL));
    }

    // Create an image to assemble input images into.
    ImageType* image = new ImageType(inputUnion.size());
    AlphaType* imageA = new AlphaType(inputUnion.size());

    if (Verbose >= VERBOSE_ASSEMBLE_MESSAGES) {
        const vigra::FilenameLayerPair file_layer =
            vigra::split_filename(imageInfoList.front()->getFileName());
        const int layers = imageInfoList.front()->numLayers();
        if (OneAtATime) {
            cerr << command
                 << ": info: loading next image: "
                 << file_layer.first << " "
                 << file_layer.second + 1 << '/' << layers
                 << endl;
        } else {
            cerr << command
                 << ": info: combining non-overlapping images: "
                 << file_layer.first << " "
                 << file_layer.second + 1 << '/' << layers;
            cerr.flush();
        }
    }

    const Diff2D imagePos = imageInfoList.front()->getPosition();
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
            ImageImportInfo* info = *i;

            // Load the next image.
            ImageType* src = new ImageType(info->size());
            AlphaType* srcA = new AlphaType(info->size());

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

                if (Verbose >= VERBOSE_ASSEMBLE_MESSAGES) {
                    cerr << " " << info->getFileName();
                    cerr.flush();
                }

                const Diff2D srcPos = info->getPosition();
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
        for (list<list<ImageImportInfo*>::iterator>::iterator r = toBeRemoved.begin();
             r != toBeRemoved.end();
             ++r) {
            imageInfoList.erase(*r);
        }
    }

    if (Verbose >= VERBOSE_ASSEMBLE_MESSAGES && !OneAtATime) {
        cerr << endl;
    }

    // Calculate bounding box of image.
    FindBoundingRectangle unionRect;
    inspectImageIf(srcIterRange(Diff2D(), Diff2D() + image->size()),
                   srcImage(*imageA), unionRect);
    bb = unionRect();

    if (Verbose >= VERBOSE_ABB_MESSAGES) {
        cerr << command
             << ": info: assembled images bounding box: "
             << unionRect()
             << endl;
    }

    return pair<ImageType*, AlphaType*>(image, imageA);
}

} // namespace enblend

#endif /* __ASSEMBLE_H__ */

// Local Variables:
// mode: c++
// End:
