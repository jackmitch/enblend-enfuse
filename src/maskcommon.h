/*
 * Copyright (C) 2004-2011 Andrew Mihal, Mikolaj Leszczynski
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
#ifndef MASKCOMMON_H
#define	MASKCOMMON_H

using vigra::LinearIntensityTransform;
namespace enblend {
    
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
    
    template <typename PixelType, typename ResultType>
    class PixelSumFunctor
    {
        typedef typename EnblendNumericTraits<PixelType>::ImagePixelComponentType PixelComponentType;
        typedef typename EnblendNumericTraits<ResultType>::ImagePixelComponentType ResultPixelComponentType;
        typedef LinearIntensityTransform<ResultType> RangeMapper;

    public:
        PixelSumFunctor() :
            rm(linearRangeMapping(NumericTraits<PixelComponentType>::min(),
                                  NumericTraits<PixelComponentType>::max(),
                                  ResultType(NumericTraits<ResultPixelComponentType>::min()),
                                  ResultType(NumericTraits<ResultPixelComponentType>::max()))) {}

        ResultType operator()(const PixelType& a, const PixelType& b) const {
            typedef typename NumericTraits<PixelType>::isScalar src_is_scalar;
            return sum(a, b, src_is_scalar());
        }

    protected:
        ResultType sum(const PixelType& a, const PixelType& b, VigraFalseType) const {
            PixelComponentType aLum = a.luminance();
            PixelComponentType bLum = b.luminance();
            PixelComponentType lumDiff = (aLum + bLum) / 2;
            return rm(lumDiff);
        }

        ResultType sum(const PixelType& a, const PixelType& b, VigraTrueType) const {
            typedef typename NumericTraits<PixelType>::isSigned src_is_signed;
            return scalar_sum(a, b, src_is_signed());
        }

        ResultType scalar_sum(const PixelType& a, const PixelType& b, VigraTrueType) const {
            return rm(a + b);
        }

        // This appears necessary because NumericTraits<unsigned int>::Promote
        // is an unsigned int instead of an int.
        ResultType scalar_sum(const PixelType& a, const PixelType& b, VigraFalseType) const {
            return rm(std::abs(static_cast<int>(a) + static_cast<int>(b)));
        }

        RangeMapper rm;
    };
    
    template <typename PixelType, typename ResultType>
    class MapFunctor
    {
        typedef typename EnblendNumericTraits<PixelType>::ImagePixelComponentType PixelComponentType;
        typedef typename EnblendNumericTraits<ResultType>::ImagePixelComponentType ResultPixelComponentType;
        typedef LinearIntensityTransform<ResultType> RangeMapper;

    public:
        MapFunctor() :
            rm(linearRangeMapping(NumericTraits<PixelComponentType>::min(),
                                  NumericTraits<PixelComponentType>::max(),
                                  ResultType(NumericTraits<ResultPixelComponentType>::min()),
                                  ResultType(NumericTraits<ResultPixelComponentType>::max()))) {}

        ResultType operator()(const PixelType& a) const {
            typedef typename NumericTraits<PixelType>::isScalar src_is_scalar;
            return map(a, src_is_scalar());
        }

    protected:
        ResultType map(const PixelType& a, VigraFalseType) const {
            PixelComponentType aLum = a.luminance();
            return rm(aLum);
        }

        ResultType map(const PixelType& a, VigraTrueType) const {
            typedef typename NumericTraits<PixelType>::isSigned src_is_signed;
            return scalar_map(a, src_is_signed());
        }

        ResultType scalar_map(const PixelType& a, VigraTrueType) const {
            return rm(a);
        }

        // This appears necessary because NumericTraits<unsigned int>::Promote
        // is an unsigned int instead of an int.
        ResultType scalar_map(const PixelType& a, VigraFalseType) const {
            return rm(std::abs(static_cast<int>(a)));
        }

        RangeMapper rm;
    };
    
} //namespace enblend

#endif	/* MASKCOMMON_H */

// Local Variables:
// mode: c++
// End:
