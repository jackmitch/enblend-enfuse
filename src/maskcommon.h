/* 
 * File:   maskcommon.h
 * Author: rosomack
 *
 * Created on June 2, 2011, 3:12 PM
 */

#ifndef MASKCOMMON_H
#define	MASKCOMMON_H

using vigra::LinearIntensityTransform;
namespace enblend{
    
    template <typename PixelType, typename ResultType>
    class PixelDifferenceFunctor
    {
        typedef typename EnblendNumericTraits<PixelType>::ImagePixelComponentType PixelComponentType;
        typedef typename EnblendNumericTraits<ResultType>::ImagePixelComponentType ResultPixelComponentType;
        typedef LinearIntensityTransform<ResultType> RangeMapper;

    public:
        PixelDifferenceFunctor() :
            rm(linearRangeMapping(NumericTraits<PixelComponentType>::min(),
                                  NumericTraits<PixelComponentType>::max(),
                                  ResultType(NumericTraits<ResultPixelComponentType>::min()),
                                  ResultType(NumericTraits<ResultPixelComponentType>::max()))) {}

        ResultType operator()(const PixelType& a, const PixelType& b) const {
            typedef typename NumericTraits<PixelType>::isScalar src_is_scalar;
            return diff(a, b, src_is_scalar());
        }

    protected:
        ResultType diff(const PixelType& a, const PixelType& b, VigraFalseType) const {
            PixelComponentType aLum = a.luminance();
            PixelComponentType bLum = b.luminance();
            PixelComponentType aHue = a.hue();
            PixelComponentType bHue = b.hue();
            PixelComponentType lumDiff = (aLum > bLum) ? (aLum - bLum) : (bLum - aLum);
            PixelComponentType hueDiff = (aHue > bHue) ? (aHue - bHue) : (bHue - aHue);
            if (hueDiff > (NumericTraits<PixelComponentType>::max() / 2)) {
                hueDiff = NumericTraits<PixelComponentType>::max() - hueDiff;
            }
            return rm(std::max(hueDiff, lumDiff));
        }

        ResultType diff(const PixelType& a, const PixelType& b, VigraTrueType) const {
            typedef typename NumericTraits<PixelType>::isSigned src_is_signed;
            return scalar_diff(a, b, src_is_signed());
        }

        ResultType scalar_diff(const PixelType& a, const PixelType& b, VigraTrueType) const {
            return rm(std::abs(a - b));
        }

        // This appears necessary because NumericTraits<unsigned int>::Promote
        // is an unsigned int instead of an int.
        ResultType scalar_diff(const PixelType& a, const PixelType& b, VigraFalseType) const {
            return rm(std::abs(static_cast<int>(a) - static_cast<int>(b)));
        }

        RangeMapper rm;
    };

} //namespace enblend

#endif	/* MASKCOMMON_H */

