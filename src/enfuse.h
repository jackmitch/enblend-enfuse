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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

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
using vigra::Kernel1D;
using vigra::VectorNormFunctor;

using boost::lambda::_1;
using boost::lambda::_2;
using boost::lambda::bind;
using boost::lambda::const_parameters;

namespace enblend {


// compute the local variance inside a window
// TODO: respect alpha mask and properly calculate borders
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void localVariance(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                   DestIterator dest_ul, DestAccessor dest_acc, 
                   Size2D size)
{
    vigra_precondition(size.x > 1 && size.y > 1,
                       "localVariance(): window for local variance must be at least 2x2");

    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SrcSumType;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    // calculate width and height of the image
    int w = src_lr.x - src_ul.x;
    int h = src_lr.y - src_ul.y;

    vigra_precondition(w >= size.x && h >= size.y,
                       "localVariance(): kernel larger than image.");

    int x,y;
    x = 0;
    y = 0;

    // create iterators for the interior part of the image (where the kernel always fits into the image)
    Diff2D border(size.x/2, size.y/2);
    DestIterator yd = dest_ul + border;
    SrcIterator ys = src_ul + border;
    SrcIterator send = src_lr - border;
    int n = size.x*size.y;
     
    // iterate over the interior part
    for(; ys.y < send.y; ++ys.y, ++yd.y)
    {
        // create x iterators
        DestIterator xd(yd);
        SrcIterator xs(ys);

        for(; xs.x < send.x; ++x, ++xs.x, ++xd.x)
        {
            // init the sum
            SrcSumType sum = NumericTraits<SrcSumType>::zero();
            SrcSumType sum_sqr = NumericTraits<SrcSumType>::zero();

            // calculate the window, required for border case
            SrcIterator yys = xs - border;
            SrcIterator yyend = xs + border;

            for(; yys.y <= yyend.y; ++yys.y)
            {
                typename SrcIterator::row_iterator xxs = yys.rowIterator();
                typename SrcIterator::row_iterator xxe = xxs + size.x;

                for(; xxs < xxe; ++xxs)
                {
                    sum += src_acc(xxs);
                    sum_sqr += src_acc(xxs) * src_acc(xxs);
                }
            }

            // store convolution result in destination pixel
            dest_acc.set(DestTraits::fromRealPromote( (sum_sqr - sum*sum/n) / (n-1)), xd);
        }
    }
}


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
    typedef ResultType result_type;

    ExposureFunctor(double w) : weight(w) {}

    inline ResultType operator()(const InputType &a) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(a, srcIsScalar());
    }

protected:

    template <typename T>
    inline ResultType f(const T &a, VigraTrueType) const {
        const T b = NumericTraits<T>::max() / 2;
        const T c = NumericTraits<T>::max() / 5;
        typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(weight * exp(-1 * (ra-b) * (ra-b) / (2 * c * c)));
    }

    template <typename T>
    inline ResultType f(const T &a, VigraFalseType) const {
        return f(a.luminance(), VigraTrueType());
    }

    double weight;
};

template <typename InputType, typename ResultType>
class SaturationFunctor {
public:
    typedef ResultType result_type;

    SaturationFunctor(double w) : weight(w) {}

    inline ResultType operator()(const InputType &a) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(a, srcIsScalar());
    }

protected:

    template <typename T>
    inline ResultType f(const T &a, VigraTrueType) const {
        return NumericTraits<ResultType>::zero();
    }

    template <typename T>
    inline ResultType f(const T &a, VigraFalseType) const {
        typedef typename T::value_type TComponentType;
        typename NumericTraits<TComponentType>::RealPromote rsa = NumericTraits<TComponentType>::toRealPromote(a.saturation());
        return NumericTraits<ResultType>::fromRealPromote(weight * rsa / NumericTraits<TComponentType>::max());
    }

    double weight;
};


template <typename InputType, typename ScaleType, typename ResultType>
class ContrastFunctor {
public:
    typedef ResultType result_type;

    ContrastFunctor(double w) : weight(w) {}

    inline ResultType operator()(const InputType &a) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(a, srcIsScalar());
    }

protected:

    template <typename T>
    inline ResultType f(const T &a, VigraTrueType) const {
        typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(sqrt(weight * ra) / NumericTraits<ScaleType>::max());
    }

    template <typename T>
    inline ResultType f(const T &a, VigraFalseType) const {
        typedef typename T::value_type TComponentType;
        typedef typename ScaleType::value_type ScaleComponentType;
        typename NumericTraits<TComponentType>::RealPromote norm = sqrt( (double) a[0] + a[1] + a[2]);
        return NumericTraits<ResultType>::fromRealPromote(weight * norm / NumericTraits<ScaleComponentType>::max());
    }

    double weight;
};


template <typename ImageType, typename AlphaType, typename MaskType>
void enfuseMask(triple<typename ImageType::const_traverser, typename ImageType::const_traverser, typename ImageType::ConstAccessor> src,
                pair<typename AlphaType::const_traverser, typename AlphaType::ConstAccessor> mask,
                pair<typename MaskType::traverser, typename MaskType::Accessor> result) {

    // Exposure
    if (WExposure) {
        transformImageIf(src, mask, result, ExposureFunctor<typename ImageType::value_type, typename MaskType::value_type>(WExposure));
    }

    // contrast criteria
    if (WContrast) {
        typedef typename ImageType::value_type ImageValueType;
        typedef typename NumericTraits<ImageValueType>::Promote GradientType;

        #ifdef ENBLEND_CACHE_IMAGES
            typedef CachedFileImage<GradientType > GradImage;
        #else
            typedef BasicImage<GradientType> GradImage;
        #endif

        GradImage grad(src.second - src.first);

        localVariance(src.first, src.second, src.third, grad.upperLeft(), grad.accessor(), Size2D(ContrastWindowSize ,ContrastWindowSize));

        // use the gray value standart deviation / norm of the color standart deviation as a contrast measure
        // The standart deviation is scaled by the max pixel value, and tends to be in the
        // range between 0.1 .. 0.3
        // use a heuristic multiplier of 5 to bring it into a range similar to the other
        // criteria
        ContrastFunctor<GradientType, ImageValueType, typename MaskType::value_type> cf(5.0*WContrast);
        combineTwoImagesIf(srcImageRange(grad), result, mask, result, const_parameters(bind(cf, _1) + _2));
    }

    // Saturation
    if (WSaturation) {
        combineTwoImagesIf(src, result, mask, result, const_parameters(bind(SaturationFunctor<typename ImageType::value_type, typename MaskType::value_type>(WSaturation), _1) + _2));
    }

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
    #ifdef ENBLEND_CACHE_IMAGES
        typedef CachedFileImage<float> MaskType;
    #else
        typedef FImage MaskType;
    #endif
    typedef typename MaskType::value_type MaskPixelType;
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

    // Result image. Alpha will be union of all input alphas.
    pair<ImageType*, AlphaType*> outputPair(NULL, new AlphaType(inputUnion.size()));

    int m = 0;
    while (!imageInfoList.empty()) {

        Rect2D imageBB;
        pair<ImageType*, AlphaType*> imagePair =
                assemble<ImageType, AlphaType>(imageInfoList, inputUnion, imageBB);

        MaskType *mask = new MaskType(inputUnion.size());

        enfuseMask<ImageType, AlphaType, MaskType>(srcImageRange(*(imagePair.first)),
                                                   srcImage(*(imagePair.second)),
                                                   destImage(*mask));

        if (Debug) {
            std::ostringstream oss;
            oss << "mask" << m << ".tif";
            ImageExportInfo maskInfo(oss.str().c_str());
            exportImage(srcImageRange(*mask), maskInfo);
        }

        // Make output alpha the union of all input alphas.
        copyImageIf(srcImageRange(*(imagePair.second)),
                    maskImage(*(imagePair.second)),
                    destImage(*(outputPair.second)));

        // Add the mask to the norm image.
        combineTwoImages(srcImageRange(*mask), srcImage(*normImage), destImage(*normImage), Arg1() + Arg2());

        imageList.push_back(make_triple(imagePair.first, imagePair.second, mask));

        ++m;
    }

    typename EnblendNumericTraits<ImagePixelType>::MaskPixelType maxMaskPixelType =
        NumericTraits<typename EnblendNumericTraits<ImagePixelType>::MaskPixelType>::max();

    if (HardMask) {
        Size2D sz = normImage->size();
        typename list< triple <ImageType*, AlphaType* , MaskType* > >::iterator imageIter;
        for (int x=0; x<sz.x; ++x) {
            for (int y=0; y < sz.y; ++y) {
                float max = 0;
                double maxi = 0;
                int i=0;
                for(imageIter=imageList.begin(); imageIter != imageList.end(); ++imageIter) {
                    float w = (*(*imageIter).third)(x,y);
                    if (w > max) {
                        max = w;
                        maxi = i;
                    }
		    i++;
                }
                i=0;
                for(imageIter=imageList.begin(); imageIter != imageList.end(); ++imageIter) {
                    if (i == maxi) {
                        (*(*imageIter).third)(x,y) = maxMaskPixelType;
                    } else {
                        (*(*imageIter).third)(x,y) = 0.0f;
                    }
                    i++;
                }
            }
        }
        int i = 0;
        if (Debug) {
            for(imageIter=imageList.begin(); imageIter != imageList.end(); ++imageIter) {
	        std::ostringstream oss;
        	oss << "mask" << i << "_wta.tif";
	        ImageExportInfo maskInfo(oss.str().c_str());
        	exportImage(srcImageRange(*(imageIter->third)), maskInfo);
                i++;
            }
	}
    }

    Rect2D junkBB;
    unsigned int numLevels = roiBounds<ImagePixelComponentType>(inputUnion, inputUnion, inputUnion, inputUnion, junkBB, Wraparound);

    vector<ImagePyramidType*> *resultLP = NULL;

    m = 0;
    while (!imageList.empty()) {
        triple<ImageType*, AlphaType*, MaskType*> imageTriple = imageList.front();
        imageList.erase(imageList.begin());

        // imageLP is constructed using the image's own alpha channel
        // as the boundary for extrapolation.
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

        if (!HardMask) {
            // Normalize the mask coefficients.
            // Scale to the range expected by the MaskPyramidPixelType.
            combineTwoImagesIf(srcImageRange(*(imageTriple.third)),
                               srcImage(*normImage),
                               maskImage(*normImage),
                               destImage(*(imageTriple.third)),
                               Param(maxMaskPixelType) * Arg1() / Arg2());
        }

        // maskGP is constructed using the union of the input alpha channels
        // as the boundary for extrapolation.
        vector<MaskPyramidType*> *maskGP =
                gaussianPyramid<MaskType, AlphaType, MaskPyramidType,
                                MaskPyramidIntegerBits, MaskPyramidFractionBits,
                                SKIPSMMaskPixelType, SKIPSMAlphaPixelType>(
                        numLevels, Wraparound, srcImageRange(*(imageTriple.third)), maskImage(*(outputPair.second)));

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

        ++m;
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
