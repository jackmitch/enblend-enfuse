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
#include <iomanip>
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

#include "vigra/colorconversions.hxx"
#include "vigra/flatmorphology.hxx"
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


// compute the local standard deviation inside a window
// TODO: respect alpha mask and properly calculate borders
template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
void localStdDevIf(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                   MaskIterator mask_ul, MaskAccessor mask_acc,
                   DestIterator dest_ul, DestAccessor dest_acc,
                   Size2D size)
{
    vigra_precondition(size.x > 1 && size.y > 1,
                       "localStdDevIf(): window for local variance must be at least 2x2");

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote SrcSumType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;

    // calculate width and height of the image
    const int w = src_lr.x - src_ul.x;
    const int h = src_lr.y - src_ul.y;

    vigra_precondition(w >= size.x && h >= size.y,
                       "localStdDevIf(): kernel larger than image.");

    // create iterators for the interior part of the image (where the kernel always fits into the image)
    Diff2D border(size.x / 2, size.y / 2);
    DestIterator yd = dest_ul + border;
    SrcIterator ys = src_ul + border;
    MaskIterator ym = mask_ul + border;
    SrcIterator send = src_lr - border;

    // iterate over the interior part
    for(; ys.y < send.y; ++ys.y, ++yd.y, ++ym.y)
    {
        // create x iterators
        DestIterator xd(yd);
        SrcIterator xs(ys);
        MaskIterator xm(ym);

        for(; xs.x < send.x; ++xs.x, ++xd.x, ++xm.x)
        {
            // init the sum
            SrcSumType sum = NumericTraits<SrcSumType>::zero();
            SrcSumType sum_sqr = NumericTraits<SrcSumType>::zero();
            int n = 0;

            // calculate the window, required for border case
            // TODO: move border cases into an own loop.
            SrcIterator yys = xs - border;
            MaskIterator yym = xm - border;

/*
            if (yys.x - ys.x < 0) {
                // left border
                yys.x = ys.x;
                yym.x = ym.x;
            }
            if (yys.y - ys.y < 0) {
                // top border
                yys.y = ys.y;
                yym.y = ym.y;
            }
*/
            SrcIterator yyend = xs + border + Diff2D(1, 1);
/*
            if (send.x - yyend.x < 0) {
                // right border
                yyend.x = send.x;
            }
            if (send.y - yyend.y < 0) {
                // bottom border
                yyend.y = send.y;
            }
*/
            for(; yys.y < yyend.y; ++yys.y, ++yym.y)
            {
                typename SrcIterator::row_iterator xxs = yys.rowIterator();
                typename SrcIterator::row_iterator xxe = yyend.rowIterator();
                typename MaskIterator::row_iterator xxm = yym.rowIterator();

                for(; xxs < xxe; ++xxs, ++xxm)
                {
                    if (mask_acc(xxm)) {
                        sum += src_acc(xxs);
                        sum_sqr += src_acc(xxs) * src_acc(xxs);
                        n++;
                    }
                }
            }

            // store convolution result in destination pixel
            // s^2 = (1 / (n-1)) * sum_sqr - (n / (n-1)) * (sum/n)^2
            //     = (1 / (n-1)) * (sum_sqr - sum^2 / n)
            SrcSumType ss = (n > 1) ? sqrt((sum_sqr - sum*sum/n) / (n-1))
                                    : NumericTraits<SrcSumType>::zero();
            dest_acc.set(DestTraits::fromRealPromote(ss), xd);
        }
    }
}


template <typename SrcIterator, typename SrcAccessor,
          typename MaskIterator, typename MaskAccessor,
          typename DestIterator, typename DestAccessor>
inline void
localStdDevIf(triple<SrcIterator, SrcIterator, SrcAccessor> src,
              pair<MaskIterator, MaskAccessor> mask,
              pair<DestIterator, DestAccessor> dest,
              Size2D size)
{
    localStdDevIf(src.first, src.second, src.third,
                  mask.first, mask.second,
                  dest.first, dest.second,
                  size);
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

    ExposureFunctor(double w,double m, double s) : weight(w), mu(m), sigma(s){}

    inline ResultType operator()(const InputType &a) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(a, srcIsScalar());
    }

protected:
    template <typename T>
    inline ResultType f(const T &a, VigraTrueType) const {
        const double b = NumericTraits<T>::max() * mu;
        const double c = NumericTraits<T>::max() * sigma;
        typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(weight * exp(-1 * (ra-b) * (ra-b) / (2 * c * c)));
    }

    template <typename T>
    inline ResultType f(const T &a, VigraFalseType) const {
        return f(a.luminance(), VigraTrueType());
    }

    double weight,mu,sigma;
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
        typedef typename NumericTraits<ScaleType>::isIntegral scaleIsIntegral;
        return f(a, srcIsScalar(), scaleIsIntegral());
    }

protected:
    // grayscale, integral
    template <typename T>
    inline ResultType f(const T &a, VigraTrueType, VigraTrueType) const {
        const typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(weight * ra / NumericTraits<ScaleType>::max());
    }

    // grayscale, floating-point
    template <typename T>
    inline ResultType f(const T &a, VigraTrueType, VigraFalseType) const {
        const typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
        return NumericTraits<ResultType>::fromRealPromote(weight * ra);
    }

    // RGB, integral
    template <typename T>
    inline ResultType f(const T &a, VigraFalseType, VigraTrueType) const {
        typedef typename T::value_type TComponentType;
        typedef typename NumericTraits<TComponentType>::RealPromote RealTComponentType;
        typedef typename ScaleType::value_type ScaleComponentType;
        const RealTComponentType ra = static_cast<RealTComponentType>(a.lightness());
        return NumericTraits<ResultType>::fromRealPromote(weight * ra / NumericTraits<ScaleComponentType>::max());
    }

    // RGB, floating-point
    template <typename T>
    inline ResultType f(const T &a, VigraFalseType, VigraFalseType) const {
        typedef typename T::value_type TComponentType;
        typedef typename NumericTraits<TComponentType>::RealPromote RealTComponentType;
        const RealTComponentType ra = static_cast<RealTComponentType>(a.lightness());
        return NumericTraits<ResultType>::fromRealPromote(weight * ra);
    }

    double weight;
};


template <typename ValueType>
struct MagnitudeAccessor
{
    typedef ValueType value_type;

    template <class Iterator>
    ValueType operator()(const Iterator& i) const {return std::abs(*i);}

    ValueType operator()(const ValueType* i) const {return std::abs(*i);}

    template <class Iterator, class Difference>
    ValueType operator()(const Iterator& i, Difference d) const {return std::abs(i[d]);}

    template <class Value, class Iterator>
    void set(const Value& v, const Iterator& i) const {
        *i = vigra::detail::RequiresExplicitCast<ValueType>::cast(std::abs(v));
    }

    template <class Value, class Iterator>
    void set(const Value& v, Iterator& i) const {
        *i = vigra::detail::RequiresExplicitCast<ValueType>::cast(std::abs(v));
    }

    template <class Value, class Iterator, class Difference>
    void set(const Value& v, const Iterator& i, const Difference& d) const {
        i[d] = vigra::detail::RequiresExplicitCast<ValueType>::cast(std::abs(v));
    }
};


template <typename InputType, typename ResultType>
class ClampingFunctor
{
public:
    ClampingFunctor(InputType lower, ResultType lowerValue,
                    InputType upper, ResultType upperValue):
        lo(lower), up(upper), loval(lowerValue), upval(upperValue)
        {}

    ResultType operator()(const InputType& x) const {
        if (x <= lo) return loval;
        else if (x >= up) return upval;
        else return x;
    }

private:
    InputType lo, up;
    ResultType loval, upval;
};


// If the first argument is lower than THRESHOLD return the second
// argument, i.e. the "fill-in" value multiplied with SCALE2,
// otherwise return the first argument multiplied with SCALE1.
template <typename InputType, typename ResultType>
class FillInFunctor
{
public:
    FillInFunctor(InputType thr, double s1, double s2):
        threshold(thr), scale1(s1), scale2(s2)
        {}

    ResultType operator()(const InputType& x, const InputType& y) const {
        if (x >= threshold) return NumericTraits<ResultType>::fromRealPromote(scale1 * x);
        else return NumericTraits<ResultType>::fromRealPromote(scale2 * y);
    }

private:
    InputType threshold;
    double scale1, scale2;
};


template <typename InputType, typename ResultType>
class MultiGrayscaleAccessor
{
public:
    typedef ResultType value_type;

    MultiGrayscaleAccessor(const char* accessorName) {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        initializeTypeSpecific(srcIsScalar());
        initialize(accessorName);
    }

    template <class Iterator>
    ResultType operator()(const Iterator& i) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(i, srcIsScalar());
    }

    template <class Iterator, class Difference>
    ResultType operator()(const Iterator& i, Difference d) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(i, d, srcIsScalar());
    }

private:
#define CHANNEL_MIXER "channel-mixer"
    void initialize(const char* accessorName) {
        if (accessorName == NULL) kind = AVERAGE;
        else
        {
            char dummy;
            double red, green, blue;
            std::string name(accessorName);

            std::transform(name.begin(), name.end(), name.begin(), tolower);
            if (name == "average") kind = AVERAGE;
            else if (name == "l-star") kind = LSTAR;
            else if (name == "lightness") kind = LIGHTNESS;
            else if (name == "value") kind = VALUE;
            else if (name == "luminance") kind = LUMINANCE;
            else if (name == CHANNEL_MIXER)
            {
                cerr <<
                    "enfuse: \"" CHANNEL_MIXER "\" is a grayscale projector requiring\n" <<
                    "enfuse:      arguments like e.g. \"channel-mixer:0.30:0.59:0.11\".";
                exit(1);
            }
            else if (sscanf(name.c_str(),
                            CHANNEL_MIXER "%["
                            OPTION_DELIMITERS "]%lf%["
                            OPTION_DELIMITERS "]%lf%["
                            OPTION_DELIMITERS "]%lf",
                            &dummy, &red, &dummy, &green, &dummy, &blue) == 6)
            {
                const double sum = red + green + blue;

                redWeight = red / sum;
                greenWeight = green / sum;
                blueWeight = blue / sum;
                kind = MIXER;
            }
            else
            {
                cerr << "enfuse: unknown grayscale projector \"" << name << "\".\n";
                exit(1);
            }
        }
    }

    void initializeTypeSpecific(VigraTrueType) {}

    void initializeTypeSpecific(VigraFalseType) {
        typedef typename InputType::value_type ValueType;
        labfun = vigra::RGB2LabFunctor<double>(NumericTraits<ValueType>::max());
    }

    template <class Iterator>
    ResultType f(const Iterator& i, VigraFalseType) const {
        typedef typename InputType::value_type ValueType;
        switch (kind)
        {
        case AVERAGE:
            return NumericTraits<ResultType>::fromRealPromote
                ((NumericTraits<ValueType>::toRealPromote((*i).red()) +
                  NumericTraits<ValueType>::toRealPromote((*i).green()) +
                  NumericTraits<ValueType>::toRealPromote((*i).blue())) /
                 3.0);
        case LSTAR:
        {
            typedef typename vigra::RGB2LabFunctor<double>::result_type LABResultType;
            const LABResultType y = labfun.operator()(*i) / 100.0;
            return NumericTraits<ResultType>::fromRealPromote(NumericTraits<ValueType>::max() * y[0]);
        }
        case LIGHTNESS:
            return NumericTraits<ResultType>::fromRealPromote
                ((std::min((*i).red(), std::min((*i).green(), (*i).blue())) +
                  std::max((*i).red(), std::max((*i).green(), (*i).blue()))) /
                 2.0);
        case VALUE:
            return std::max((*i).red(), std::max((*i).green(), (*i).blue()));
        case LUMINANCE:
            return NumericTraits<ResultType>::fromRealPromote((*i).luminance());
        case MIXER:
            return NumericTraits<ResultType>::fromRealPromote
                (redWeight * NumericTraits<ValueType>::toRealPromote((*i).red()) +
                 greenWeight * NumericTraits<ValueType>::toRealPromote((*i).green()) +
                 blueWeight * NumericTraits<ValueType>::toRealPromote((*i).blue()));
        }

        // never reached
        return ResultType();
    }

    template <class Iterator, class Difference>
    ResultType f(const Iterator& i, Difference d, VigraFalseType) const {
        typedef typename InputType::value_type ValueType;
        switch (kind)
        {
        case AVERAGE:
            return NumericTraits<ResultType>::fromRealPromote
                ((NumericTraits<ValueType>::toRealPromote(i[d].red()) +
                  NumericTraits<ValueType>::toRealPromote(i[d].green()) +
                  NumericTraits<ValueType>::toRealPromote(i[d].blue())) /
                 3.0);
        case LSTAR:
        {
            typedef typename vigra::RGB2LabFunctor<double>::result_type LABResultType;
            const LABResultType y = labfun.operator()(i[d]) / 100.0;
            return NumericTraits<ResultType>::fromRealPromote(NumericTraits<ValueType>::max() * y[0]);
        }
        case LIGHTNESS:
            return NumericTraits<ResultType>::fromRealPromote
                ((std::min(i[d].red(), std::min(i[d].green(), i[d].blue())) +
                  std::max(i[d].red(), std::max(i[d].green(), i[d].blue()))) /
                 2.0);
        case VALUE:
            return std::max(i[d].red(), std::max(i[d].green(), i[d].blue()));
        case LUMINANCE:
            return NumericTraits<ResultType>::fromRealPromote(i[d].luminance());
        case MIXER:
            return NumericTraits<ResultType>::fromRealPromote
                (redWeight * NumericTraits<ValueType>::toRealPromote(i[d].red()) +
                 greenWeight * NumericTraits<ValueType>::toRealPromote(i[d].green()) +
                 blueWeight * NumericTraits<ValueType>::toRealPromote(i[d].blue()));
        }

        // never reached
        return ResultType();
    }

    template <class Iterator>
    ResultType f(const Iterator& i, VigraTrueType) const {return *i;}

    template <class Iterator, class Difference>
    ResultType f(const Iterator& i, Difference d, VigraTrueType) const {return i[d];}

    enum AccKind {AVERAGE, LSTAR, LIGHTNESS, VALUE, LUMINANCE, MIXER} kind;
    double redWeight, greenWeight, blueWeight;
    vigra::RGB2LabFunctor<double> labfun;
};


template <typename ImageType, typename AlphaType, typename MaskType>
void enfuseMask(triple<typename ImageType::const_traverser, typename ImageType::const_traverser, typename ImageType::ConstAccessor> src,
                pair<typename AlphaType::const_traverser, typename AlphaType::ConstAccessor> mask,
                pair<typename MaskType::traverser, typename MaskType::Accessor> result) {

    // Exposure
    if (WExposure > 0.0) {
        transformImageIf(src, mask, result, ExposureFunctor<typename ImageType::value_type,
                typename MaskType::value_type>(WExposure, WMu, WSigma));
    }

    // contrast criteria
    if (WContrast > 0.0) {
        typedef typename ImageType::PixelType PixelType;
        typedef typename NumericTraits<PixelType>::ValueType ScalarType;
        typedef typename NumericTraits<ScalarType>::Promote LongScalarType;
#ifdef ENBLEND_CACHE_IMAGES
        typedef CachedFileImage<LongScalarType> GradImage;
#else
        typedef BasicImage<LongScalarType> GradImage;
#endif
        typedef typename GradImage::iterator GradIterator;

        const typename ImageType::difference_type imageSize = src.second - src.first;
        GradImage grad(imageSize);
        MultiGrayscaleAccessor<PixelType, LongScalarType> ga(GrayscaleProjector);

        if (FilterConfig.edgeScale > 0.0)
        {
#ifdef DEBUG_LOG
            cout << "+ Laplacian Edge Detection, scale = " << FilterConfig.edgeScale << " pixels" << endl;
#endif
            GradImage laplacian(imageSize);

            if (FilterConfig.lceScale > 0.0)
            {
#ifdef DEBUG_LOG
                cout <<
                    "+ Local Contrast Enhancement, (scale, amount) = " << FilterConfig.lceScale <<
                    " pixels, " << 100.0 * FilterConfig.lceFactor << "%" << endl;
#endif
                GradImage lce(imageSize);
                gaussianSharpening(src.first, src.second, ga,
                                   lce.upperLeft(), lce.accessor(),
                                   FilterConfig.lceFactor, FilterConfig.lceScale);
                laplacianOfGaussian(lce.upperLeft(), lce.lowerRight(), lce.accessor(),
                                    laplacian.upperLeft(), MagnitudeAccessor<LongScalarType>(),
                                    FilterConfig.edgeScale);
            }
            else
            {
                laplacianOfGaussian(src.first, src.second, ga,
                                    laplacian.upperLeft(), MagnitudeAccessor<LongScalarType>(),
                                    FilterConfig.edgeScale);
            }

#ifdef DEBUG_LOG
            {
                vigra::FindMinMax<LongScalarType> minmax;
                inspectImage(srcImageRange(laplacian), minmax);
                cout << "+ after Laplacian and Magnitude: min = " <<
                    minmax.min << ", max = " << minmax.max << endl;
            }
#endif

            const double minCurve =
                MinCurvature.isPercentage ?
                static_cast<double>(NumericTraits<ScalarType>::max()) * MinCurvature.value / 100.0 :
                MinCurvature.value;
            if (minCurve <= 0.0)
            {
#ifdef DEBUG_LOG
                cout << "+ truncate values below " << -minCurve << endl;;
#endif
                transformImageIf(laplacian.upperLeft(), laplacian.lowerRight(), laplacian.accessor(),
                                 mask.first, mask.second,
                                 grad.upperLeft(), grad.accessor(),
                                 ClampingFunctor<LongScalarType, LongScalarType>
                                 (static_cast<LongScalarType>(-minCurve), LongScalarType(),
                                  NumericTraits<LongScalarType>::max(), NumericTraits<LongScalarType>::max()));
            }
            else
            {
#ifdef DEBUG_LOG
                cout << "+ merge local contrast and edges - switch at " << minCurve << endl;
#endif
                GradImage localContrast(imageSize);
                // TODO: use localStdDev
                localStdDevIf(src.first, src.second, ga,
                              mask.first, mask.second,
                              localContrast.upperLeft(), localContrast.accessor(),
                              Size2D(ContrastWindowSize, ContrastWindowSize));

                combineTwoImagesIf(laplacian.upperLeft(), laplacian.lowerRight(), laplacian.accessor(),
                                   localContrast.upperLeft(), localContrast.accessor(),
                                   mask.first, mask.second,
                                   grad.upperLeft(), grad.accessor(),
                                   FillInFunctor<LongScalarType, LongScalarType>
                                   (static_cast<LongScalarType>(minCurve), // threshold
                                    1.0, // scale factor for "laplacian"
                                    minCurve / NumericTraits<ScalarType>::max())); // scale factor for "localContrast"
            }
        }
        else
        {
#ifdef DEBUG_LOG
            cout << "+ Variance of Local Contrast" << endl;
#endif
            localStdDevIf(src.first, src.second, ga,
                          mask.first, mask.second,
                          grad.upperLeft(), grad.accessor(),
                          Size2D(ContrastWindowSize, ContrastWindowSize));
        }

#ifdef DEBUG_LOG
        {
            vigra::FindMinMax<LongScalarType> minmax;
            inspectImage(srcImageRange(grad), minmax);
            cout << "+ final grad: min = " << minmax.min << ", max = " << minmax.max << endl;
        }
#endif
        ContrastFunctor<LongScalarType, ScalarType, typename MaskType::value_type> cf(WContrast);
        combineTwoImagesIf(srcImageRange(grad), result, mask, result, const_parameters(bind(cf, _1) + _2));
    }

    // Saturation
    if (WSaturation > 0.0) {
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
            oss << "mask" << std::setw(4) << std::setfill('0') << m << ".tif";
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

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after loading image " << m << " :" << endl;
            v.printStats("image", imagePair.first);
            v.printStats("alpha", imagePair.second);
            v.printStats("weight", mask);
            v.printStats("normImage", normImage);
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif


        ++m;
    }

    const int totalImages = imageList.size();

    typename EnblendNumericTraits<ImagePixelType>::MaskPixelType maxMaskPixelType =
        NumericTraits<typename EnblendNumericTraits<ImagePixelType>::MaskPixelType>::max();

    if (HardMask) {
        if (Verbose) {
            cout << "Creating hard blend mask" << std::endl;
        }
        Size2D sz = normImage->size();
        typename list< triple <ImageType*, AlphaType* , MaskType* > >::iterator imageIter;
        for (int y = 0; y < sz.y; ++y) {
            for (int x = 0; x < sz.x; ++x) {
                float max = 0.0f;
                int maxi = 0;
                int i = 0;
                for(imageIter = imageList.begin(); imageIter != imageList.end(); ++imageIter) {
                    const float w = static_cast<float>((*(*imageIter).third)(x, y));
                    if (w > max) {
                        max = w;
                        maxi = i;
                    }
		    i++;
                }
                i = 0;
                for(imageIter = imageList.begin(); imageIter != imageList.end(); ++imageIter) {
                    if (max == 0.0f) {
                        (*(*imageIter).third)(x, y) = static_cast<MaskPixelType>(maxMaskPixelType) / totalImages;
                    } else if (i == maxi) {
                        (*(*imageIter).third)(x, y) = maxMaskPixelType;
                    } else {
                        (*(*imageIter).third)(x, y) = 0.0f;
                    }
                    i++;
                }
            }
        }
        int i = 0;
        if (Debug) {
            for(imageIter = imageList.begin(); imageIter != imageList.end(); ++imageIter) {
	        std::ostringstream oss;
        	oss << "mask" << std::setw(4) << std::setfill('0') << i << "_wta.tif";
	        ImageExportInfo maskInfo(oss.str().c_str());
        	exportImage(srcImageRange(*(imageIter->third)), maskInfo);
                i++;
            }
	}
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector &v = CachedFileImageDirector::v();
            cout << "Image cache statistics after creating hard mask:" << endl;
            v.printStats();
            v.resetCacheMisses();
            cout << "--------------------------------------------------------------------------------" << endl;
        }
        #endif
    }

    Rect2D junkBB;
    unsigned int numLevels = roiBounds<ImagePixelComponentType>(inputUnion, inputUnion, inputUnion, inputUnion, junkBB, Wraparound);

    vector<ImagePyramidType*> *resultLP = NULL;

    m = 0;
    while (!imageList.empty()) {
        triple<ImageType*, AlphaType*, MaskType*> imageTriple = imageList.front();
        imageList.erase(imageList.begin());

        std::ostringstream oss0;
        oss0 << "imageGP" << m << "_";

        // imageLP is constructed using the image's own alpha channel
        // as the boundary for extrapolation.
        vector<ImagePyramidType*> *imageLP =
                laplacianPyramid<ImageType, AlphaType, ImagePyramidType,
                                 ImagePyramidIntegerBits, ImagePyramidFractionBits,
                                 SKIPSMImagePixelType, SKIPSMAlphaPixelType>(
                        oss0.str().c_str(),
                        numLevels, Wraparound,
                        srcImageRange(*(imageTriple.first)),
                        maskImage(*(imageTriple.second)));

        delete imageTriple.first;
        delete imageTriple.second;

        //std::ostringstream oss1;
        //oss1 << "imageLP" << m << "_";
        //exportPyramid<ImagePyramidType>(imageLP, oss1.str().c_str());

        if (!HardMask) {
            // Normalize the mask coefficients.
            // Scale to the range expected by the MaskPyramidPixelType.
            combineTwoImages(srcImageRange(*(imageTriple.third)),
                             srcImage(*normImage),
                             destImage(*(imageTriple.third)),
                             ifThenElse(Arg2() > Param(0.0),
                                        Param(maxMaskPixelType) * Arg1() / Arg2(),
                                        Param(maxMaskPixelType / totalImages)));
        }

        // maskGP is constructed using the union of the input alpha channels
        // as the boundary for extrapolation.
        vector<MaskPyramidType*> *maskGP =
                gaussianPyramid<MaskType, AlphaType, MaskPyramidType,
                                MaskPyramidIntegerBits, MaskPyramidFractionBits,
                                SKIPSMMaskPixelType, SKIPSMAlphaPixelType>(
                        numLevels, Wraparound, srcImageRange(*(imageTriple.third)), maskImage(*(outputPair.second)));

        delete imageTriple.third;

        //std::ostringstream oss2;
        //oss2 << "maskGP" << m << "_";
        //exportPyramid<MaskPyramidType>(maskGP, oss2.str().c_str());

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

        //std::ostringstream oss3;
        //oss3 << "multLP" << m << "_";
        //exportPyramid<ImagePyramidType>(imageLP, oss3.str().c_str());

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

        //std::ostringstream oss4;
        //oss4 << "resultLP" << m << "_";
        //exportPyramid<ImagePyramidType>(resultLP, oss4.str().c_str());

        ++m;
    }

    delete normImage;

    //exportPyramid<ImagePyramidType>(resultLP, "resultLP");

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