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
#ifndef __ENFUSE_H__
#define __ENFUSE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>
#include <stdio.h>

#include <boost/static_assert.hpp>

#include "common.h"
#include "numerictraits.h"
#include "fixmath.h"
#include "assemble.h"
#include "blend.h"
#include "bounds.h"
#include "pyramid.h"

#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/transformimage.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::functor::Arg1;
using vigra::functor::Arg2;
using vigra::functor::Param;
using vigra::BasicImage;
using vigra::CachedFileImage;
using vigra::CachedFileImageDirector;
using vigra::FImage;
using vigra::FindMinMax;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::initImage;
using vigra::initImageIf;
using vigra::inspectImage;
using vigra::NumericTraits;
using vigra::Size2D;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

namespace enblend {

template <typename MaskPixelType>
class ImageMaskMultiplyFunctor {
public:
    ImageMaskMultiplyFunctor(MaskPixelType d) : divisor(NumericTraits<MaskPixelType>::toRealPromote(d)) {}

    template <typename ImagePixelType>
    ImagePixelType operator()(const ImagePixelType &iP, const MaskPixelType &maskP) const {

        typedef typename NumericTraits<ImagePixelType>::RealPromote RealImagePixelType;

        // Convert mask pixel to blend coefficient in range [0.0, 1.0].
        double maskCoeff = NumericTraits<MaskPixelType>::toRealPromote(maskP) / divisor;

        RealImagePixelType riP = NumericTraits<ImagePixelType>::toRealPromote(iP);

        RealImagePixelType blendP = riP * maskCoeff;

        return NumericTraits<ImagePixelType>::fromRealPromote(blendP);
    }

protected:
    double divisor;
};

template <typename InputType, typename ResultType>
class ExposureFunctor {
public:
    ExposureFunctor(double w) : weight(w) {}

    inline ResultType operator()(const InputType a) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(a, srcIsScalar());
    }

protected:

    template <typename T>
    inline ResultType f(const T a, VigraTrueType) const {
        const T b = NumericTraits<T>::max() / 2;
        const T c = NumericTraits<T>::max() / 5;
        typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(weight * exp(-1 * (ra-b) * (ra-b) / (2 * c * c)));
    }

    template <typename T>
    inline ResultType f(const T a, VigraFalseType) const {
        return f(a.luminance(), VigraTrueType());
    }

    double weight;
};

template <typename ImageType, typename AlphaType, typename MaskType>
void enfuseMask(triple<typename ImageType::const_traverser, typename ImageType::const_traverser, typename ImageType::ConstAccessor> src,
                pair<typename AlphaType::const_traverser, typename AlphaType::ConstAccessor> mask,
                pair<typename MaskType::traverser, typename MaskType::Accessor> result) {

    // Exposure
    transformImageIf(src, mask, result, ExposureFunctor<typename ImageType::value_type, typename MaskType::value_type>(WExposure));

    // TODO: contrast

    // TODO: saturation

};

/** Enfuse's main blending loop. Templatized to handle different image types.
 */
template <typename ImagePixelType>
void enfuseMain(list<ImageImportInfo*> &imageInfoList,
        ImageExportInfo &outputImageInfo,
        Rect2D &inputUnion) {

    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType ImagePixelComponentType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImageType ImageType;
    typedef typename EnblendNumericTraits<ImagePixelType>::AlphaPixelType AlphaPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::AlphaType AlphaType;
    typedef typename FImage::value_type MaskPixelType;
    typedef FImage MaskType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePyramidPixelType ImagePyramidPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePyramidType ImagePyramidType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskPyramidPixelType MaskPyramidPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskPyramidType MaskPyramidType;

    enum {ImagePyramidIntegerBits = EnblendNumericTraits<ImagePixelType>::ImagePyramidIntegerBits};
    enum {ImagePyramidFractionBits = EnblendNumericTraits<ImagePixelType>::ImagePyramidFractionBits};
    enum {MaskPyramidIntegerBits = EnblendNumericTraits<ImagePixelType>::MaskPyramidIntegerBits};
    enum {MaskPyramidFractionBits = EnblendNumericTraits<ImagePixelType>::MaskPyramidFractionBits};
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMImagePixelType SKIPSMImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMAlphaPixelType SKIPSMAlphaPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMMaskPixelType SKIPSMMaskPixelType;

    // List of input image / input alpha / mask triples
    list <triple<ImageType*, AlphaType*, MaskType*> > imageList;

    // Sum of all masks
    MaskType *normImage = new MaskType(inputUnion.size());

    while (!imageInfoList.empty()) {

        Rect2D imageBB;
        pair<ImageType*, AlphaType*> imagePair =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, imageBB);

        MaskType *mask = new MaskType(inputUnion.size());

        enfuseMask<ImageType, AlphaType, MaskType>(srcImageRange(*(imagePair.first)),
                                                   srcImage(*(imagePair.second)),
                                                   destImage(*mask));

        // Add the mask to the norm image.
        combineTwoImages(srcImageRange(*mask), srcImage(*normImage), destImage(*normImage), Arg1() + Arg2());

        imageList.push_back(make_triple(imagePair.first, imagePair.second, mask));

    }

    // Result image. Alpha will be union of all input alphas.
    pair<ImageType*, AlphaType*> outputPair(NULL, new AlphaType(inputUnion.size()));

    Rect2D junkBB;
    unsigned int numLevels = roiBounds<ImagePixelComponentType>(inputUnion, inputUnion, inputUnion, inputUnion, junkBB, Wraparound);

    vector<ImagePyramidType*> *resultLP = NULL;

    while (!imageList.empty()) {
        triple<ImageType*, AlphaType*, MaskType*> imageTriple = imageList.front();
        imageList.erase(imageList.begin());

        // Make output alpha the union of all input alphas.
        copyImageIf(srcImageRange(*(imageTriple.second)),
                    maskImage(*(imageTriple.second)),
                    destImage(*(outputPair.second)));

        vector<ImagePyramidType*> *imageLP =
                laplacianPyramid<ImageType, AlphaType, ImagePyramidType,
                                 ImagePyramidIntegerBits, ImagePyramidFractionBits,
                                 SKIPSMImagePixelType, SKIPSMAlphaPixelType>(
                        "imageGP",
                        numLevels, Wraparound,
                        srcImageRange(*(imageTriple.first)),
                        maskImage(*(imageTriple.second)));

        delete imageTriple.first;
        delete imageTriple.second;

        typename EnblendNumericTraits<ImagePixelType>::MaskPixelType maxMaskPixelType =
            NumericTraits<typename EnblendNumericTraits<ImagePixelType>::MaskPixelType>::max();

        // Normalize the mask coefficients.
        // Scale to the range expected by the MaskPyramidPixelType.
        combineTwoImagesIf(srcImageRange(*(imageTriple.third)),
                           srcImage(*normImage),
                           maskImage(*normImage),
                           destImage(*(imageTriple.third)),
                           Param(maxMaskPixelType) * Arg1() / Arg2());

        vector<MaskPyramidType*> *maskGP =
                gaussianPyramid<MaskType, MaskPyramidType,
                                MaskPyramidIntegerBits, MaskPyramidFractionBits,
                                SKIPSMMaskPixelType>(
                        numLevels, Wraparound, srcImageRange(*(imageTriple.third)));

        delete imageTriple.third;

        ConvertScalarToPyramidFunctor<typename EnblendNumericTraits<ImagePixelType>::MaskPixelType,
                                      MaskPyramidPixelType,
                                      MaskPyramidIntegerBits,
                                      MaskPyramidFractionBits> maskConvertFunctor;
        MaskPyramidPixelType maxMaskPyramidPixelValue = maskConvertFunctor(maxMaskPixelType);

        for (unsigned int i = 0; i < maskGP->size(); ++i) {
            // Multiply image lp with the mask gp.
            combineTwoImages(srcImageRange(*((*imageLP)[i])),
                             srcImage(*((*maskGP)[i])),
                             destImage(*((*imageLP)[i])),
                             ImageMaskMultiplyFunctor<MaskPyramidPixelType>(maxMaskPyramidPixelValue));

            // Done with maskGP.
            delete (*maskGP)[i];
        }
        delete maskGP;

        if (resultLP != NULL) {
            // Add imageLP to resultLP.
            for (unsigned int i = 0; i < imageLP->size(); ++i) {
                combineTwoImages(srcImageRange(*((*imageLP)[i])),
                                 srcImage(*((*resultLP)[i])),
                                 destImage(*((*resultLP)[i])),
                                 Arg1() + Arg2());
                delete (*imageLP)[i];
            }
            delete imageLP;
        }
        else {
            resultLP = imageLP;
        }

    }

    collapsePyramid<SKIPSMImagePixelType>(Wraparound, resultLP);

    outputPair.first = new ImageType(inputUnion.size());

    copyFromPyramidImageIf<ImagePyramidType, AlphaType, ImageType,
                           ImagePyramidIntegerBits, ImagePyramidFractionBits>(
            srcImageRange(*((*resultLP)[0])),
            maskImage(*(outputPair.second)),
            destImage(*(outputPair.first)));

    // Delete result pyramid.
    for (unsigned int i = 0; i < resultLP->size(); ++i) {
        delete (*resultLP)[i];
    }
    delete resultLP;

    checkpoint(outputPair, outputImageInfo);

    delete outputPair.first;
    delete outputPair.second;

};

} // namespace enblend

#endif /* __ENFUSE_H__ */
