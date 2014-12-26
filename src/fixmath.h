/*
 * Copyright (C) 2004-2014 Andrew Mihal
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

#ifndef FIXMATH_H_INCLUDED_
#define FIXMATH_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include <time.h>

#ifdef _WIN32
#include <boost/math/special_functions.hpp>
using namespace boost::math;
#endif

#include <gsl/gsl_rng.h>

#include <vigra/basicimage.hxx>
#include <vigra/colorconversions.hxx>
#include <vigra/mathutil.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/utilities.hxx>

#include "minimizer.h"
#include "muopt.h"
#include "parameter.h"


#define XYZ_SCALE 100.0


namespace enblend {

static inline double
wrap_cyclically(double x, double modulus)
{
    while (x < 0.0)
    {
        x += modulus;
    }

    return fmod(x, modulus);
}


static inline double
limit(double x, double lower_limit, double upper_limit)
{
    if (EXPECT_RESULT(std::isnan(x), false))
    {
        throw std::range_error("limit: not a number");
    }

    return std::min(std::max(lower_limit, x), upper_limit);
}


static inline void
jch_to_rgb(const cmsJCh* jch, double* rgb)
{
    cmsCIEXYZ scaled_xyz;
    cmsCIECAM02Reverse(CIECAMTransform, jch, &scaled_xyz);
    // xyz values *approximately* in range [0, 100]

    // scale xyz values to range [0, 1]
    const double xyz[] = {
        scaled_xyz.X / XYZ_SCALE,
        scaled_xyz.Y / XYZ_SCALE,
        scaled_xyz.Z / XYZ_SCALE
    };

    cmsDoTransform(XYZToInputTransform, xyz, rgb, 1U);
    // rgb values *approximately* in range [0, 1]
}


/** A functor for converting scalar pixel values to the number representation used
 *  for pyramids.  These are either fixed-point integers or floating-point numbers.
 */
template <typename SrcPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertScalarToPyramidFunctor
{
    typedef vigra::NumericTraits<SrcPixelType> SrcTraits;
    typedef typename SrcTraits::isIntegral SrcIsIntegral;

    typedef vigra::NumericTraits<PyramidPixelType> PyramidTraits;
    typedef typename PyramidTraits::isIntegral PyramidIsIntegral;

public:
    ConvertScalarToPyramidFunctor() :
        pyramid_scale(double(1U << PyramidFractionBits))
    {
#ifdef DEBUG_COLORSPACE_STATISTICS
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            "+ ConvertScalarToPyramidFunctor::ConvertScalarToPyramidFunctor: " <<
            (SrcIsIntegral().value ? "integral" : "floating-point") << " source  =>  " <<
            (PyramidIsIntegral().value? "integral" : "floating-point") << " pyramid\n";
#endif
    }

    PyramidPixelType operator()(const SrcPixelType& v) const
    {
        return doConvert(v, SrcIsIntegral(), PyramidIsIntegral());
    }

protected:
    // Convert an integral pixel type to an integral pyramid value type.
    PyramidPixelType doConvert(const SrcPixelType& v, vigra::VigraTrueType, vigra::VigraTrueType) const
    {
        return convertIntegerToFixedPoint(v);
    }

    // Convert an integral pixel type to a real pyramid value type.
    PyramidPixelType doConvert(const SrcPixelType& v, vigra::VigraTrueType, vigra::VigraFalseType) const
    {
        return SrcTraits::toRealPromote(v);
    }

    // Convert a real pixel type to an integral pyramid value type.
    PyramidPixelType doConvert(const SrcPixelType& v, vigra::VigraFalseType, vigra::VigraTrueType) const
    {
        return convertDoubleToFixedPoint(v);
    }

    // Convert a real pixel type to a real pyramid value type.
    PyramidPixelType doConvert(const SrcPixelType& v, vigra::VigraFalseType, vigra::VigraFalseType) const
    {
        // Convert real data using a log transform.  These achieves
        // two purposes:
        //  1. During blending, even completely non-negative images
        //     can result in negative pixels.  A log transform
        //     followed by the exp inverse guarantees all-positive
        //     output.
        //  2. For HDR data, the log transform put the samples closer
        //     to a perceptual space making the blending a little more
        //     pleasing.  Ideally, all blending should be done in a
        //     strictly perceptually-linear space, such as Luv or Lab.
        //
        // See ConvertPyramidToScalarFunctor::doConvert for the
        // inverse transform.
        //
        // Check for non-positive values -- they should not be in the
        // input, but if they are we need to handle them or log will
        // return a NaN.
        //
        // v >= 0.0 ? 1.0 + log(v + 1.0) : 1.0 / (1.0 - v)
        return v >= 0.0 ? 1.0 + log1p(v) : 1.0 / (1.0 - v);
    }

    PyramidPixelType convertDoubleToFixedPoint(const double& v) const
    {
        // Shift v to get the appropriate number of fraction bits into
        // the integer part, then fromRealPromote this value into the
        // fixed-point type.
        return PyramidTraits::fromRealPromote(v * pyramid_scale);
    }

    PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType& v) const
    {
        // Shift v left to move the decimal point and set the fraction
        // bits to zero.
        return static_cast<PyramidPixelType>(v) << PyramidFractionBits;
    }

private:
    const double pyramid_scale;
};


class MersenneTwister
{
public:
    typedef unsigned long result_type;

    MersenneTwister() : generator_(gsl_rng_alloc(gsl_rng_mt19937)) {assert(generator_);}
    MersenneTwister(const MersenneTwister& another_generator) :
        generator_(gsl_rng_clone(another_generator.generator_)) {assert(generator_);}
    ~MersenneTwister() {gsl_rng_free(generator_);}

    MersenneTwister& operator=(const MersenneTwister& another_generator)
    {
        if (this != &another_generator)
        {
            gsl_rng_free(generator_);
            generator_ = gsl_rng_clone(another_generator.generator_);
            assert(generator_);
        }
        return *this;
    }

    result_type min() const {return gsl_rng_min(generator_);}
    result_type max() const {return gsl_rng_max(generator_);}

    void seed() {gsl_rng_set(generator_, gsl_rng_default_seed);}
    void seed(result_type a_seed) {gsl_rng_set(generator_, a_seed);}

    result_type operator()() {return gsl_rng_get(generator_);}

private:
    gsl_rng* generator_;
};


inline static unsigned
non_deterministic_seed()
{
    unsigned seed = static_cast<unsigned>(1 + omp_get_thread_num());

    const clock_t now = clock();
    if (now != static_cast<clock_t>(-1))
    {
        seed ^= static_cast<unsigned>(now);
    }

    return seed;
}


// A functor for converting numbers stored in the pyramid number
// representation back into normal pixel values.

template <typename DestPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToScalarFunctor
{
    typedef vigra::NumericTraits<DestPixelType> DestTraits;
    typedef typename DestTraits::isIntegral DestIsIntegral;

    typedef vigra::NumericTraits<PyramidPixelType> PyramidTraits;
    typedef typename PyramidTraits::isIntegral PyramidIsIntegral;

public:
    ConvertPyramidToScalarFunctor()
    {
#ifdef DEBUG
        random_number_generator_.seed(); // apply default seed in all threads for reproducibility
#else
        random_number_generator_.seed(non_deterministic_seed());
#endif
    }

    DestPixelType operator()(const PyramidPixelType& v) const
    {
        return doConvert(v, DestIsIntegral(), PyramidIsIntegral());
    }

protected:
    // test time with floating-point dithering: 100.01 sec
    // test time with integer dithering: 94.89 sec
    // Convert an integral pyramid pixel to an integral image pixel.
    DestPixelType doConvert(const PyramidPixelType& v, vigra::VigraTrueType, vigra::VigraTrueType) const
    {
        // Integer Dithering
        constexpr PyramidPixelType half_m1 = (1U << (PyramidFractionBits - 1)) - 1;
        constexpr PyramidPixelType quarter = 1U << (PyramidFractionBits - 2);
        constexpr PyramidPixelType three_quarter = 3U << (PyramidFractionBits - 2);

        const PyramidPixelType vShifted = v >> PyramidFractionBits;
        const PyramidPixelType vFraction = v & ((1U << PyramidFractionBits) - 1U);

        if (vFraction >= quarter && vFraction < three_quarter)
        {
            const PyramidPixelType random =
                (PyramidPixelType(random_number_generator_()) & half_m1) + quarter;

            if (random <= vFraction)
            {
                return DestPixelType(DestTraits::fromPromote(vShifted + 1));
            }
            else
            {
                return DestPixelType(DestTraits::fromPromote(vShifted));
            }
        }
        else if (vFraction >= quarter)
        {
            return DestPixelType(DestTraits::fromPromote(vShifted + 1));
        }
        else
        {
            return DestPixelType(DestTraits::fromPromote(vShifted));
        }
    }

    // Convert a real pyramid pixel to an integral image pixel.
    DestPixelType doConvert(const PyramidPixelType& v, vigra::VigraTrueType, vigra::VigraFalseType) const
    {
        return DestTraits::fromRealPromote(dither(v));
    }

    // Convert an integral pyramid pixel to a real image pixel.
    DestPixelType doConvert(const PyramidPixelType& v, vigra::VigraFalseType, vigra::VigraTrueType) const
    {
        return convertFixedPointToDouble(v);
    }

    // Convert a real pyramid pixel to a real image pixel.
    DestPixelType doConvert(const PyramidPixelType& v, vigra::VigraFalseType, vigra::VigraFalseType) const
    {
        // Undo logarithmic/rational mapping that was done in building
        // the pyramid.  See ConvertScalarToPyramidFunctor::doConvert
        // for the forward transformation.
        return v >= 1.0 ? expm1(v - 1.0) : 1.0 - 1.0 / v;
    }

    // Dithering is used to fool the eye into seeing gradients that are finer
    // than the precision of the pixel type.
    // This prevents the occurence of cleanly-bordered regions in the output where
    // the pixel values suddenly change from N to N+1.
    // Such regions are especially objectionable in the green channel of 8-bit images.
    double dither(const double& v) const
    {
        const double vFraction = v - floor(v);

        // Only dither values within a certain range of the rounding cutoff point.
        if (vFraction > 0.25 && vFraction <= 0.75)
        {
            if (vFraction - 0.25 >= 0.5 * random())
            {
                return ceil(v);
            }
            else
            {
                return floor(v);
            }
        }

        return v;
    }

    double convertFixedPointToDouble(const PyramidPixelType& v) const
    {
        return PyramidTraits::toRealPromote(v) / static_cast<double>(1U << PyramidFractionBits);
    }

private:
    double random() const
    {
        return static_cast<double>(random_number_generator_()) / static_cast<double>(random_number_generator_.max());
    }

    mutable MersenneTwister random_number_generator_;
};


//
// Trivial RGB wrapper for vector pixel types
//

template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToPyramidFunctor
{
    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertScalarToPyramidFunctor<SrcComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;
public:
    ConvertVectorToPyramidFunctor() : converter() {}

    PyramidVectorType operator()(const SrcVectorType& v) const
    {
        return PyramidVectorType(converter(v.red()), converter(v.green()), converter(v.blue()));
    }

protected:
    ConvertFunctorType converter;
};


//
// Trivial RGB wrapper for vector pixel types
//

template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToVectorFunctor
{
    typedef typename DestVectorType::value_type DestComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertPyramidToScalarFunctor<DestComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertPyramidToVectorFunctor() : converter() {}

    DestVectorType operator()(const PyramidVectorType& v) const
    {
        return DestVectorType(converter(v.red()), converter(v.green()), converter(v.blue()));
    }

protected:
    ConvertFunctorType converter;
};


//
//  Configuration of CIELAB and CIELUV forward and backward transforms
//

template <int PyramidIntegerBits, int PyramidFractionBits>
class PyramidScale
{
    enum
    {
        // Maximum value of the `L'-channel for both L*a*b* and L*u*v*
        // color spaces.
        MAXIMUM_LIGHTNESS = 100,

        // PYRAMID_HEADROOM_BITS is the number of bits we must reserve
        // when we scale to the fixed-point based representation
        // inside the pyramid.  They buffer soft overflows, like
        // e.g. a lightness of 100.001 and they are required to hold
        // the result of an addition or subtraction.
        PYRAMID_HEADROOM_BITS = 2
    };

    //                  Vigra Manual             |                 Our Findings
    //    ==============================================================================
    //                                           |
    //    Lab                                    |    Lab / copy into pyramid
    //           0       <=  L*  <=  100         |           0      <=  L*  <=  100
    //         -86.1813  <=  a*  <=   98.2352    |         -79.275  <=  a*  <=   93.547
    //        -107.862   <=  b*  <=   94.4758    |        -112.033  <=  b*  <=   93.391
    //                                           |    Lab / copy from pyramid
    //                                           |           0      <=  L*  <=  108.539
    //                                           |         -79.375  <=  a*  <=   91.281
    //                                           |        -111.125  <=  b*  <=  182.562
    //                                           |
    //    Luv                                    |    Luv / copy into pyramid
    //           0       <=  L*  <=  100         |           0      <=  L*  <=  100.001
    //         -83.077   <=  u*  <=  175.015     |         -34.796  <=  u*  <=  155.726
    //        -134.101   <=  v*  <=  107.393     |        -114.366  <=  v*  <=  109.423
    //                                           |    Luv / copy from pyramid
    //                                           |           0      <=  L*  <=  136.719
    //                                           |         -37.703  <=  u*  <=  101.203
    //                                           |        -115.094  <=  v*  <=  200.422

public:
    PyramidScale() :
        color_limit(parameter::as_double("lab-color-limit", 127.0)),
        pyramid_scale(double(1U << (PyramidIntegerBits - PYRAMID_HEADROOM_BITS)))
    {
#ifdef DEBUG_COLORSPACE_STATISTICS
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            "+ PyramidScale::PyramidScale: PyramidIntegerBits = " << PyramidIntegerBits <<
            ", PyramidFractionBits = " << PyramidFractionBits << "  =>  pyramid-scale = " << pyramid_scale <<
            "\n";
#endif // DEBUG_COLORSPACE_STATISTICS
    }

    double scale_lightness_for_pyramid(double a_lightness) const
    {
        const double scaled_lightness = pyramid_scale * a_lightness / static_cast<double>(MAXIMUM_LIGHTNESS);
        const double result = std::round(scaled_lightness);

#ifdef DEBUG_COLORSPACE_STATISTICS
        // Note: For we rescale color channels with respect to
        // `color_limit' and color channels as well as luminance ends
        // up in the same pyramid, we can check for overflow with the
        // same limit.
        if (result < 0.0 || result > 1.0 + color_limit)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            std::cout <<
                "+ PyramidScale::scale_lightness_for_pyramid: out-of-range L = " <<
                a_lightness << ", scaled = " << scaled_lightness << ", rounded = " << result << "\n";
        }
#endif // DEBUG_COLORSPACE_STATISTICS

        assert(result >= 0.0);
        return result;
    }

    double scale_color_difference_for_pyramid(double a_color_difference) const
    {
        const double scaled_color_difference =
            pyramid_scale * (color_limit + a_color_difference) / (2.0 * color_limit);
        const double result = std::round(scaled_color_difference);

#ifdef DEBUG_COLORSPACE_STATISTICS
        if (std::abs(result) > color_limit)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            std::cout <<
                "+ PyramidScale::scale_color_difference_for_pyramid: out-of-range a,b|u,v = " <<
                a_color_difference << ", scaled = " << scaled_color_difference << ", rounded = " << result <<
                "\n";
        }
#endif // DEBUG_COLORSPACE_STATISTICS

        assert(result >= -color_limit);
        return result;
    }

    double scale_lightness_of_pyramid(double a_lightness) const
    {
        return a_lightness * static_cast<double>(MAXIMUM_LIGHTNESS) / pyramid_scale;
    }

    double scale_color_difference_of_pyramid(double a_color_difference) const
    {
        return a_color_difference * (2.0 * color_limit) / pyramid_scale - color_limit;
    }

private:
    const double color_limit;
    const double pyramid_scale;
};


#ifdef DEBUG_COLORSPACE_STATISTICS
class TriplePeakHold
{
    enum {n = 3};

public:
    TriplePeakHold()
    {
        for (int i = 0; i < n; ++i)
        {
            std::ostringstream s;
            s << '#' << (i + 1);
            labels_[i] = s.str();
        }

        for (auto& p : peaks_)
        {
            p.first = std::numeric_limits<double>::max();
            p.second = std::numeric_limits<double>::min();
        }
    }

    TriplePeakHold(const std::string& label1, const std::string& label2, const std::string& label3)
    {
        labels_[0] = label1;
        labels_[1] = label2;
        labels_[2] = label3;

        for (auto& p : peaks_)
        {
            p.first = std::numeric_limits<double>::max();
            p.second = std::numeric_limits<double>::min();
        }
    }

    void update(double value1, double value2, double value3)
    {
        peaks_[0].first = std::min(peaks_[0].first, value1);
        peaks_[0].second = std::max(peaks_[0].second, value1);
        peaks_[1].first = std::min(peaks_[1].first, value2);
        peaks_[1].second = std::max(peaks_[1].second, value2);
        peaks_[2].first = std::min(peaks_[2].first, value3);
        peaks_[2].second = std::max(peaks_[2].second, value3);
    }

    std::string as_string(const std::string a_prefix = std::string()) const
    {
        std::ostringstream s;

        for (int i = 0; i < n; ++i)
        {
            s <<
                a_prefix <<
                peaks_[i].first << " <= " << labels_[i] << " <= " << peaks_[i].second <<
                "\n";
        }

        return s.str();
    }

private:
    std::array<std::string, n> labels_;
    std::array<std::pair<double, double>, n> peaks_;
};
#endif // DEBUG_COLORSPACE_STATISTICS


//
//  Fixed point converter that uses ICC profile transformation and L*a*b* color space
//

template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToLabPyramidFunctor :
    private PyramidScale<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScale<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef vigra::NumericTraits<SrcComponentType> SrcTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertScalarToPyramidFunctor<SrcComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertVectorToLabPyramidFunctor() :
        converter(),
        rgb_source_scale(1.0 / SrcTraits::toRealPromote(SrcTraits::max()))
#ifdef DEBUG_COLORSPACE_STATISTICS
        , range("L", "a", "b")
#endif // DEBUG_COLORSPACE_STATISTICS
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            "+ ConvertVectorToLabPyramidFunctor::ConvertVectorToLabPyramidFunctor source range [" <<
            static_cast<double>(SrcTraits::min()) << ", " <<
            static_cast<double>(SrcTraits::max()) << "]\n" <<
            "+ ConvertVectorToLabPyramidFunctor::ConvertVectorToLabPyramidFunctor: pyramid range [" <<
            static_cast<double>(vigra::NumericTraits<PyramidComponentType>::min()) << ", " <<
            static_cast<double>(vigra::NumericTraits<PyramidComponentType>::max()) << "]\n";
    }

#ifdef DEBUG_COLORSPACE_STATISTICS
    ~ConvertVectorToLabPyramidFunctor()
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            range.as_string("+ ConvertVectorToLabPyramidFunctor::~ConvertVectorToLabPyramidFunctor: ") <<
            std::endl;
    }
#endif // DEBUG_COLORSPACE_STATISTICS

    PyramidVectorType operator()(const SrcVectorType& v) const
    {
        const double rgb[] = {
            rgb_source_scale * SrcTraits::toRealPromote(v.red()),
            rgb_source_scale * SrcTraits::toRealPromote(v.green()),
            rgb_source_scale * SrcTraits::toRealPromote(v.blue())
        };
        cmsCIELab lab;

        cmsDoTransform(InputToLabTransform, rgb, &lab, 1U);
#ifdef DEBUG_COLORSPACE_STATISTICS
        range.update(lab.L, lab.a, lab.b);
#endif // DEBUG_COLORSPACE_STATISTICS

        return PyramidVectorType(converter(Scale::scale_lightness_for_pyramid(lab.L)),
                                 converter(Scale::scale_color_difference_for_pyramid(lab.a)),
                                 converter(Scale::scale_color_difference_for_pyramid(lab.b)));
    }

protected:
    ConvertFunctorType converter;
    const double rgb_source_scale;
#ifdef DEBUG_COLORSPACE_STATISTICS
    mutable TriplePeakHold range;
#endif // DEBUG_COLORSPACE_STATISTICS
};


//
// Fixed point converter that uses ICC profile transformation and L*a*b* color space
//

template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertLabPyramidToVectorFunctor :
    private PyramidScale<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScale<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename DestVectorType::value_type DestComponentType;
    typedef vigra::NumericTraits<DestComponentType> DestTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertPyramidToScalarFunctor<DestComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertLabPyramidToVectorFunctor() :
        converter(),
        rgb_dest_scale(DestTraits::toRealPromote(DestTraits::max()))
#ifdef DEBUG_COLORSPACE_STATISTICS
        , range("L", "a", "b")
#endif // DEBUG_COLORSPACE_STATISTICS
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            "+ ConvertLabPyramidToVectorFunctor::ConvertLabPyramidToVectorFunctor: dest range [" <<
            static_cast<double>(DestTraits::min()) << ", " <<
            static_cast<double>(DestTraits::max()) << "]\n" <<
            "+ ConvertLabPyramidToVectorFunctor::ConvertLabPyramidToVectorFunctor: pyramid range [" <<
            static_cast<double>(vigra::NumericTraits<PyramidComponentType>::min()) << ", " <<
            static_cast<double>(vigra::NumericTraits<PyramidComponentType>::max()) << "]\n";
    }

#ifdef DEBUG_COLORSPACE_STATISTICS
    ~ConvertLabPyramidToVectorFunctor()
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            range.as_string("+ ConvertLabPyramidToVectorFunctor::~ConvertLabPyramidToVectorFunctor: ") <<
            std::endl;
    }
#endif // DEBUG_COLORSPACE_STATISTICS

    DestVectorType operator()(const PyramidVectorType& v) const
    {
        const cmsCIELab lab = {
            Scale::scale_lightness_of_pyramid(converter(v.red())),
            Scale::scale_color_difference_of_pyramid(converter(v.green())),
            Scale::scale_color_difference_of_pyramid(converter(v.blue()))
        };
        double rgb[3];

#ifdef DEBUG_COLORSPACE_STATISTICS
        range.update(lab.L, lab.a, lab.b);
#endif // DEBUG_COLORSPACE_STATISTICS
        cmsDoTransform(LabToInputTransform, &lab, rgb, 1U);

        return DestVectorType(DestTraits::fromRealPromote(rgb_dest_scale * rgb[0]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[1]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[2]));
    }

protected:
    ConvertFunctorType converter;
    const double rgb_dest_scale;
#ifdef DEBUG_COLORSPACE_STATISTICS
    mutable TriplePeakHold range;
#endif // DEBUG_COLORSPACE_STATISTICS
};


//
//  Fixed point converter that uses ICC profile transformation and L*u*v* color space
//

template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToLuvPyramidFunctor :
    private PyramidScale<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScale<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef vigra::NumericTraits<SrcComponentType> SrcTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertScalarToPyramidFunctor<SrcComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;
    typedef vigra::XYZ2LuvFunctor<double> XYZ2LuvFunctor;

public:
    ConvertVectorToLuvPyramidFunctor() :
        converter(),
        rgb_source_scale(1.0 / SrcTraits::toRealPromote(SrcTraits::max()))
#ifdef DEBUG_COLORSPACE_STATISTICS
        , range("L", "u", "v")
#endif // DEBUG_COLORSPACE_STATISTICS
    {}

#ifdef DEBUG_COLORSPACE_STATISTICS
    ~ConvertVectorToLuvPyramidFunctor()
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            range.as_string("+ ConvertVectorToLuvPyramidFunctor::~ConvertVectorToLuvPyramidFunctor: ") <<
            std::endl;
    }
#endif // DEBUG_COLORSPACE_STATISTICS

    PyramidVectorType operator()(const SrcVectorType& v) const
    {
        const double rgb[] = {
            rgb_source_scale * SrcTraits::toRealPromote(v.red()),
            rgb_source_scale * SrcTraits::toRealPromote(v.green()),
            rgb_source_scale * SrcTraits::toRealPromote(v.blue())
        };
        XYZ2LuvFunctor::argument_type xyz;

        cmsDoTransform(InputToXYZTransform, rgb, &xyz[0], 1U);
        const XYZ2LuvFunctor::result_type luv {xyz2luv(xyz)};
#ifdef DEBUG_COLORSPACE_STATISTICS
        range.update(luv[0], luv[1], luv[2]);
#endif // DEBUG_COLORSPACE_STATISTICS

        return PyramidVectorType(converter(Scale::scale_lightness_for_pyramid(luv[0])),
                                 converter(Scale::scale_color_difference_for_pyramid(luv[1])),
                                 converter(Scale::scale_color_difference_for_pyramid(luv[2])));
    }

protected:
    XYZ2LuvFunctor xyz2luv;
    ConvertFunctorType converter;
    const double rgb_source_scale;
#ifdef DEBUG_COLORSPACE_STATISTICS
    mutable TriplePeakHold range;
#endif // DEBUG_COLORSPACE_STATISTICS
};


//
// Fixed point converter that uses ICC profile transformation and L*u*v* color space
//

template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertLuvPyramidToVectorFunctor :
    private PyramidScale<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScale<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename DestVectorType::value_type DestComponentType;
    typedef vigra::NumericTraits<DestComponentType> DestTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertPyramidToScalarFunctor<DestComponentType, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;
    typedef vigra::Luv2XYZFunctor<double> Luv2XYZFunctor;

public:
    ConvertLuvPyramidToVectorFunctor() :
        converter(),
        rgb_dest_scale(DestTraits::toRealPromote(DestTraits::max()))
#ifdef DEBUG_COLORSPACE_STATISTICS
        , range("L", "u", "v")
#endif // DEBUG_COLORSPACE_STATISTICS
    {}

#ifdef DEBUG_COLORSPACE_STATISTICS
    ~ConvertLuvPyramidToVectorFunctor()
    {
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            range.as_string("+ ConvertLuvPyramidToVectorFunctor::~ConvertLuvPyramidToVectorFunctor: ") <<
            std::endl;
    }
#endif // DEBUG_COLORSPACE_STATISTICS

    DestVectorType operator()(const PyramidVectorType& v) const
    {
        const Luv2XYZFunctor::value_type luv {
            Scale::scale_lightness_of_pyramid(converter(v.red())),
            Scale::scale_color_difference_of_pyramid(converter(v.green())),
            Scale::scale_color_difference_of_pyramid(converter(v.blue()))
        };
#ifdef DEBUG_COLORSPACE_STATISTICS
        range.update(luv[0], luv[1], luv[2]);
#endif // DEBUG_COLORSPACE_STATISTICS
        const Luv2XYZFunctor::result_type xyz {luv2xyz(luv)};
        double rgb[3];

        cmsDoTransform(XYZToInputTransform, &xyz[0], rgb, 1U);

        return DestVectorType(DestTraits::fromRealPromote(rgb_dest_scale * rgb[0]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[1]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[2]));
    }

protected:
    Luv2XYZFunctor luv2xyz;
    ConvertFunctorType converter;
    const double rgb_dest_scale;
#ifdef DEBUG_COLORSPACE_STATISTICS
    mutable TriplePeakHold range;
#endif // DEBUG_COLORSPACE_STATISTICS
};


template <int PyramidIntegerBits, int PyramidFractionBits>
class PyramidScaleJCh
{
public:
    enum
    {
        MAXIMUM_LIGHTNESS = 100,    // J
        MAXIMUM_CHROMA = 120,       // C
        MAXIMUM_HUE = 360           // h
    };

    PyramidScaleJCh() :
        pyramid_scale(double(1U << (PyramidIntegerBits - 1 - 7)))
    {
#ifdef DEBUG_COLORSPACE_STATISTICS
#ifdef OPENMP
#pragma omp critical
#endif
        std::cout <<
            "+ PyramidScaleJCh::PyramidScaleJCh: PyramidIntegerBits = " << PyramidIntegerBits <<
            ", PyramidFractionBits = " << PyramidFractionBits << "  =>  pyramid-scale = " << pyramid_scale <<
            "\n";
#endif // DEBUG_COLORSPACE_STATISTICS
    }

    double scale_lightness_for_pyramid(double a_lightness) const
    {
        return a_lightness * pyramid_scale;
    }

    double scale_chroma_hue_x_for_pyramid(double a_chroma, double a_hue_angle) const
    {
        return a_chroma * sin(a_hue_angle) * pyramid_scale;
    }

    double scale_chroma_hue_y_for_pyramid(double a_chroma, double a_hue_angle) const
    {
        return a_chroma * cos(a_hue_angle) * pyramid_scale;
    }

    double scale_lightness_of_pyramid(double a_scaled_lightness) const
    {
        return a_scaled_lightness / pyramid_scale;
    }

    double scale_chroma_of_pyramid(double a_scaled_x, double a_scaled_y) const
    {
        return hypot(a_scaled_x / pyramid_scale, a_scaled_y / pyramid_scale);
    }

    double scale_hue_of_pyramid(double a_scaled_x, double a_scaled_y) const
    {
        return wrap_cyclically(degree_of_radian(atan2(a_scaled_x / pyramid_scale, a_scaled_y / pyramid_scale)),
                               MAXIMUM_HUE);
    }

private:
    double degree_of_radian(double x) const {return x * (180.0 / M_PI);}

    const double pyramid_scale;
};


//
// Fixed point converter that uses ICC profile transformation and JCh color space
//

template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToJCHPyramidFunctor :
    private PyramidScaleJCh<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScaleJCh<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef vigra::NumericTraits<SrcComponentType> SrcTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertScalarToPyramidFunctor<double, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertVectorToJCHPyramidFunctor() :
        converter(),
        rgb_source_scale(1.0 / SrcTraits::toRealPromote(SrcTraits::max()))
    {}

    PyramidVectorType operator()(const SrcVectorType& v) const
    {
        // rgb values must be in range [0, 1]
        const double rgb[] = {
            rgb_source_scale * SrcTraits::toRealPromote(v.red()),
            rgb_source_scale * SrcTraits::toRealPromote(v.green()),
            rgb_source_scale * SrcTraits::toRealPromote(v.blue())
        };
        cmsJCh jch;

        rgb_to_jch(rgb, &jch);

        const double theta = radian_of_degree(jch.h);

        return PyramidVectorType(converter(Scale::scale_lightness_for_pyramid(jch.J)),
                                 converter(Scale::scale_chroma_hue_x_for_pyramid(jch.C, theta)),
                                 converter(Scale::scale_chroma_hue_y_for_pyramid(jch.C, theta)));
    }

protected:
    double radian_of_degree(double x) const
    {
        return x * (M_PI / 180.0);
    }

    void rgb_to_jch(const double* rgb, cmsJCh* jch) const
    {
        double xyz[3];
        cmsDoTransform(InputToXYZTransform, rgb, xyz, 1U);

        const cmsCIEXYZ scaled_xyz = {XYZ_SCALE * xyz[0], XYZ_SCALE * xyz[1], XYZ_SCALE * xyz[2]};
        cmsJCh jch_unlimited;
        cmsCIECAM02Forward(CIECAMTransform, &scaled_xyz, &jch_unlimited);

        jch->J = EXPECT_RESULT(std::isnan(jch_unlimited.J), false) ? 0.0 : jch_unlimited.J;
        jch->C = EXPECT_RESULT(std::isnan(jch_unlimited.C), false) ? 0.0 : jch_unlimited.C;
        jch->h = wrap_cyclically(jch_unlimited.h, Scale::MAXIMUM_HUE);
    }

    ConvertFunctorType converter;
    const double rgb_source_scale;
};


namespace ciecam_detail
{
    template <typename forward_iterator>
    static inline void
    limit_sequence(forward_iterator first, forward_iterator last, double lower_limit, double upper_limit)
    {
        while (first != last)
        {
            *first = limit(*first, lower_limit, upper_limit);
            ++first;
        }
    }


    static inline double
    uniform_random(unsigned* seed)
    {
        return static_cast<double>(enblend::rand_r(seed)) / static_cast<double>(RAND_MAX);
    }


    static inline bool
    bracket_minimum(const gsl_function& cost,
                    double& x_initial, double x_lower, double x_upper,
                    unsigned maximum_tries)
    {
        const double y_minimum_bound =
            std::min(cost.function(x_lower, cost.params), cost.function(x_upper, cost.params));
        double y_initial = cost.function(x_initial, cost.params);

        if (y_initial < y_minimum_bound)
        {
            return true;
        }

        unsigned i = 0U;
        const double lower = std::max(0.001, 1.001 * x_lower);
        const double upper = 0.999 * x_upper;
        unsigned seed = 1000003U; // fixed seed for reproducibility

        while (y_initial >= y_minimum_bound && i < maximum_tries)
        {
            x_initial = uniform_random(&seed) * (upper - lower) + lower;
            y_initial = cost.function(x_initial, cost.params);
            ++i;
        }

        return i < maximum_tries;
    }


    static inline int
    alternating_power_spacing(int i, int n,
                              double a, double b, double c,
                              double p)
    {
        const bool is_even_n = n % 2 == 0;
        const double left_unit_stride = 1.0 / static_cast<double>(n - (is_even_n ? 1 : 2));
        const double right_unit_stride = 1.0 / static_cast<double>(n - (is_even_n ? 2 : 1));
        const double left_width = c - a;
        const double right_width = b - c;

        const double x = static_cast<double>(i);
        double y = c;

        if (i % 2 == 1)         // 1, 3, 5, ...
        {
            y -= left_width * std::pow(x * left_unit_stride, 1.0 / p);
        }
        else                    // 0, 2, 4, 6, ...
        {
            y += right_width * std::pow(x * right_unit_stride, p);
        }

        return y;
    }


    struct extra_minimizer_parameter
    {
        explicit extra_minimizer_parameter(const cmsJCh& out_of_box_jch) : jch(out_of_box_jch)
        {
            double rgb[3];

            jch_to_rgb(&jch, rgb);
            cmsDoTransform(InputToLabTransform, rgb, &bad_lab, 1);
        }

        cmsJCh jch;
        cmsCIELab bad_lab;
    };


    inline double
    delta_e_of_lab_and_rgb(const cmsCIELab* lab, const double* rgb)
    {
        cmsCIELab lab_of_rgb;

        cmsDoTransform(InputToLabTransform, rgb, &lab_of_rgb, 1);

        return cmsCMCdeltaE(lab, &lab_of_rgb, 2.0, 1.0);
    }


    inline double
    out_of_box_penalty(const double* rgb)
    {
        const double infinite_badness = 10000.0;
        double result = 0.0;

        for (int i = 0; i < 3; ++i)
        {
            if (rgb[i] > 1.0)
            {
                result += rgb[i] * infinite_badness;
            }
            else if (rgb[i] < 0.0)
            {
                result += (1.0 - rgb[i]) * infinite_badness;
            }
        }

        return result;
    }


    inline double
    delta_e_cost(const cmsJCh* jch, const extra_minimizer_parameter* parameter)
    {
        double rgb[3];
        jch_to_rgb(jch, rgb);

        return delta_e_of_lab_and_rgb(&parameter->bad_lab, rgb) + out_of_box_penalty(rgb);
    }


    double
    delta_e_min_cost(double luminance, void* data)
    {
        const extra_minimizer_parameter* parameter = static_cast<const extra_minimizer_parameter*>(data);
        const cmsJCh jch = {luminance, parameter->jch.C, parameter->jch.h};

        return delta_e_cost(&jch, parameter);
    }


    double
    delta_e_multimin_cost(const gsl_vector* x, void* data)
    {
        const extra_minimizer_parameter* parameter = static_cast<const extra_minimizer_parameter*>(data);
        const cmsJCh jch = {gsl_vector_get(x, 0), gsl_vector_get(x, 1), parameter->jch.h};

        return delta_e_cost(&jch, parameter);
    }


    void
    show_jch_rgb(const std::string& label, const cmsJCh* jch)
    {
        double rgb[3];

        jch_to_rgb(jch, rgb);

        std::cout <<
            label << " J = " << jch->J << ", C = " << jch->C << ", h = " << jch->h << "\n" <<
            label << " RGB = (" << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << ")" <<
            std::endl;
    }
} // namespace ciecam_detail


//
// Fixed point converter that uses ICC profile transformation and JCh color space
//

template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertJCHPyramidToVectorFunctor :
    private PyramidScaleJCh<PyramidIntegerBits, PyramidFractionBits>
{
    typedef PyramidScaleJCh<PyramidIntegerBits, PyramidFractionBits> Scale;
    typedef typename DestVectorType::value_type DestComponentType;
    typedef vigra::NumericTraits<DestComponentType> DestTraits;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertPyramidToScalarFunctor<double, PyramidComponentType,
                                          PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertJCHPyramidToVectorFunctor() :
        converter(),
        rgb_dest_scale(DestTraits::toRealPromote(DestTraits::max())),

        // Parameters for highlight optimizer only
        highlight_lightness_guess_1d_factor(limit(parameter::as_double("highlight-recovery-lightness-guess-factor", 0.975),
                                                  0.25, 4.0)),
        highlight_lightness_guess_1d_offset(parameter::as_double("highlight-recovery-lightness-guess-offset", 0.0)),

        maximum_highlight_iterations(limit(parameter::as_unsigned("highlight-recovery-maximum-iterations", 100U),
                                           10U, 1000U)),
        maximum_highlight_bracket_tries(limit(parameter::as_unsigned("highlight-recovery-bracket-maximum-tries", 1000U),
                                              10U, 1000000U)),
        highlight_simplex_lightness_step_length(limit(parameter::as_double("highlight-recovery-lightness-step-length", 12.5),
                                                      1.0 / 65536.0, 100.0)),
        highlight_simplex_chroma_step_length(limit(parameter::as_double("highlight-recovery-chroma-step-length", 6.25),
                                                   1.0 / 65536.0, 120.0)),
        highlight_iterations_per_leg(limit(parameter::as_unsigned("highlight-recovery-iterations-per-leg", 50U),
                                           5U, 500U)),
        maximum_highlight_leg(limit(parameter::as_unsigned("highlight-recovery-maximum-legs", 10U), 1U, 100U)),
        shadow_disguised_as_highlight_j(limit(parameter::as_double("shadow-disguised-as-highlight-lightness", 1.0),
                                              0.0001, 10.0)),

        // Parameters for shadow optimizer only
        shadow_lightness_lightness_guess_factor(parameter::as_double("shadow-recovery-lightness-lightness-guess-factor", 1.24)),
        shadow_lightness_chroma_guess_factor(parameter::as_double("shadow-recovery-lightness-chroma-guess-factor", -0.136)),
        shadow_lightness_guess_offset(parameter::as_double("shadow-recovery-lightness-guess-offset", 0.0)),
        shadow_chroma_lightness_guess_factor(parameter::as_double("shadow-recovery-chroma-lightness-guess-factor", -0.604)),
        shadow_chroma_chroma_guess_factor(parameter::as_double("shadow-recovery-chroma-chroma-guess-factor", 1.33)),
        shadow_chroma_guess_offset(parameter::as_double("shadow-recovery-chroma-guess-offset", 0.0)),

        shadow_simplex_lightness_step_length(limit(parameter::as_double("shadow-recovery-lightness-step-length", 0.625),
                                                   1.0 / 65536.0, 100.0)),
        shadow_simplex_chroma_step_length(limit(parameter::as_double("shadow-recovery-chroma-step-length", 1.25),
                                                1.0 / 65536.0, 120.0)),
        shadow_iterations_per_leg(limit(parameter::as_unsigned("shadow-recovery-iterations-per-leg", 40U),
                                        4U, 400U)),
        maximum_shadow_leg(limit(parameter::as_unsigned("shadow-recovery-maximum-legs", 5U), 1U, 50U)),
        maximum_multistart_tries(limit(parameter::as_unsigned("shadow-recovery-maximum_tries", 20U), 1U, 500U)),

        // Parameters for both optimizers
        // Desired error limits: LoFi: 0.5/2^8, HiFi: 0.5/2^16, Super-HiFi: 0.5/2^24
        optimizer_error(limit(parameter::as_double("ciecam-optimizer-error", 0.5 / 65536.0),
                              0.5 / 16777216.0, 1.0)),
        // Delta-E goals: LoFi: 1.0, HiFi: 0.5, Super-HiFi: 0.0
        optimizer_goal(limit(parameter::as_double("ciecam-optimizer-deltae-goal", 0.5), 0.0, 10.0))
    {}

    double highlight_lightness_guess_1d(const cmsJCh& jch) const
    {
        return std::min( // heuristic function with fitted parameter
                        highlight_lightness_guess_1d_factor * jch.J + highlight_lightness_guess_1d_offset,
                        0.995 * Scale::MAXIMUM_LIGHTNESS); // backstop such that our guess is less than the maximum
    }

    double highlight_lightness_guess_2d(const cmsJCh& jch) const
    {
        return std::min(0.99609375 * Scale::MAXIMUM_LIGHTNESS, jch.J);
    }

    double highlight_chroma_guess_2d(const cmsJCh& jch) const
    {
        return std::min(0.99609375 * Scale::MAXIMUM_CHROMA, jch.C);
    }

    double shadow_lightness_guess_2d(const cmsJCh& jch) const
    {
        return std::max(shadow_lightness_lightness_guess_factor * jch.J +
                        shadow_lightness_chroma_guess_factor * jch.C +
                        shadow_lightness_guess_offset,
                        0.0);
    }

    double shadow_chroma_guess_2d(const cmsJCh& jch) const
    {
        return std::max(shadow_chroma_lightness_guess_factor * jch.J +
                        shadow_chroma_chroma_guess_factor * jch.C +
                        shadow_chroma_guess_offset,
                        0.0);
    }

    double optimize_1d(cmsJCh initial_jch,
                              double initial_lightness,
                              unsigned maximum_iterations,
                              cmsJCh& final_jch) const
    {
        ciecam_detail::extra_minimizer_parameter extra(initial_jch);
        gsl_function cost = {ciecam_detail::delta_e_min_cost, &extra};

        GoldenSectionMinimizer1D
            optimizer(cost, initial_lightness,
                      0.0, std::max(static_cast<double>(Scale::MAXIMUM_LIGHTNESS), initial_jch.J));

        optimizer.set_absolute_error(optimizer_error)->
            set_goal(optimizer_goal)->
            set_maximum_number_of_iterations(maximum_iterations);
        optimizer.run();

        final_jch.J = optimizer.x_minimum();

#ifdef DEBUG_OPTIMIZE_1D_STATISTICS
#ifdef OPENMP
#pragma omp critical
#endif
        {
            ciecam_detail::show_jch_rgb("+ optimize_1d: initial", &initial_jch);
            std::cout <<
                "+ optimize_1d: reached after " << optimizer.number_of_iterations() << " iterations\n" <<
                "+ optimize_1d: delta-E = " << optimizer.f_minimum() << "\n";
            ciecam_detail::show_jch_rgb("+ optimize_1d: final", &final_jch);
        }
#endif

        return optimizer.f_minimum();
    }

    double optimize_2d(cmsJCh initial_jch,
                       double initial_lightness, double initial_chroma,
                       double initial_lightness_step_length, double initial_chroma_step_length,
                       unsigned maximum_leg, unsigned iterations_per_leg,
                       cmsJCh& final_jch) const
    {
        ciecam_detail::extra_minimizer_parameter extra(initial_jch);
        gsl_multimin_function cost = {ciecam_detail::delta_e_multimin_cost, 2U, &extra};
        const MinimizerMultiDimensionSimplex::array_type initial = {
            initial_lightness, initial_chroma
        };
        MinimizerMultiDimensionSimplex::array_type step = {
            initial_lightness_step_length,
            initial_chroma_step_length
        };
        MinimizerMultiDimensionSimplex2Randomized optimizer(cost, initial, step);

        optimizer.set_absolute_error(optimizer_error)->set_goal(optimizer_goal);
        for (unsigned leg = 1U; leg <= maximum_leg; ++leg)
        {
            optimizer.set_maximum_number_of_iterations(leg * iterations_per_leg);
            optimizer.run();
            if (optimizer.has_reached_goal())
            {
                break;
            }

            step[0] = optimizer.characteristic_size();
            step[1] = optimizer.characteristic_size();
            optimizer.set_step_sizes(step);
        }

        MinimizerMultiDimensionSimplex::array_type minimum_parameter(2U);
        optimizer.x_minimum(minimum_parameter.begin());

        final_jch = {minimum_parameter[0], minimum_parameter[1], initial_jch.h};

        return optimizer.f_minimum();
    }

    void multistart_optimize_2d(const cmsJCh* initial_jch, double* rgb) const
    {
        const double guessed_j = shadow_lightness_guess_2d(*initial_jch);
        const double guessed_c = shadow_chroma_guess_2d(*initial_jch);
        double best_deltae = std::numeric_limits<double>::infinity();
        cmsJCh best_jch = *initial_jch;

        unsigned n = 0U;

        while (n < maximum_multistart_tries)
        {
            cmsJCh jch;
            const double opt_deltae =
                optimize_2d(*initial_jch,
                            ciecam_detail::alternating_power_spacing(n, maximum_multistart_tries,
                                                                     0.0, Scale::MAXIMUM_LIGHTNESS, guessed_j,
                                                                     2.0),
                            guessed_c,
                            shadow_simplex_lightness_step_length, shadow_simplex_chroma_step_length,
                            maximum_shadow_leg, shadow_iterations_per_leg,
                            jch);
            if (opt_deltae < best_deltae)
            {
                best_deltae = opt_deltae;
                best_jch = jch;
            }
            if (opt_deltae <= optimizer_goal)
            {
                best_deltae = opt_deltae;
                best_jch = jch;
                break;
            }

            ++n;
        }

        jch_to_rgb(&best_jch, rgb);

#ifdef DEBUG_SHADOW_HIGHLIGHT_STATISTICS
        if (best_deltae > optimizer_goal)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            std::cout <<
                "\n" <<
                "+ multistart_optimize_2d: recovery failure: deltaE = " << best_deltae <<
                " after " << n + 1U << " iteration[s]\n";
            ciecam_detail::show_jch_rgb("+ multistart_optimize_2d: initial", initial_jch);
            ciecam_detail::show_jch_rgb("+ multistart_optimize_2d: final", &best_jch);
        }
#endif
    }

    void flexible_optimize_1d_2d(const cmsJCh* jch, double* rgb) const
    {
        double guessed_j = highlight_lightness_guess_1d(*jch);
        ciecam_detail::extra_minimizer_parameter extra(*jch);
        gsl_function cost = {ciecam_detail::delta_e_min_cost, &extra};
        cmsJCh opt_jch = *jch;
        double delta_e_1d = 0.0;
        double delta_e_2d __attribute__((unused)) = 0.0;

        if (EXPECT_RESULT(ciecam_detail::bracket_minimum(cost, guessed_j,
                                                         0.0, std::max(static_cast<double>(Scale::MAXIMUM_LIGHTNESS), jch->J),
                                                         maximum_highlight_bracket_tries),
                          true))
        {
            delta_e_1d = optimize_1d(*jch, guessed_j, maximum_highlight_iterations, opt_jch);
            if (delta_e_1d > optimizer_goal)
            {
#ifdef DEBUG_HIGHLIGHT_FALLBACK_STATISTICS
#ifdef OPENMP
#pragma omp critical
#endif
                std::cout <<
                    "+ flexible_optimize_1d_2d: falling back from J-optimizer to (J, C)-optimizer for J = " <<
                    jch->J << ", {C = " << jch->C << ", h = " << jch->h << "}\n" <<
                    "+          1d opt J = " << opt_jch.J << " and 1d delta-E = " << delta_e_1d <<
                    std::endl;
#endif
                delta_e_2d =
                    optimize_2d(opt_jch,
                                highlight_lightness_guess_2d(opt_jch), highlight_chroma_guess_2d(opt_jch),
                                highlight_simplex_lightness_step_length, highlight_simplex_chroma_step_length,
                                maximum_highlight_leg, highlight_iterations_per_leg,
                                opt_jch);
            }
        }
        else
        {
            delta_e_2d =
                optimize_2d(*jch,
                            highlight_lightness_guess_2d(*jch), highlight_chroma_guess_2d(*jch),
                            highlight_simplex_lightness_step_length, highlight_simplex_chroma_step_length,
                            maximum_highlight_leg, highlight_iterations_per_leg,
                            opt_jch);
        }

        jch_to_rgb(&opt_jch, rgb);

#ifdef DEBUG_SHADOW_HIGHLIGHT_STATISTICS
        if (delta_e_1d > optimizer_goal && delta_e_2d > optimizer_goal)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            {
                std::cout <<
                    "\n" <<
                    "+ flexible_optimize_1d_2d: failed to reach optimizer goal " <<
                    optimizer_goal << " -- only achived deltaE{1d} = " << delta_e_1d <<
                    ", deltaE{2d} = " << delta_e_2d << "\n";
                ciecam_detail::show_jch_rgb("+ flexible_optimize_1d_2d: initial", jch);
                ciecam_detail::show_jch_rgb("+ flexible_optimize_1d_2d: final", &opt_jch);
            }
        }
#endif
    }

    DestVectorType operator()(const PyramidVectorType& v) const
    {
        const double j = converter(v.red());
        const double ch_x = converter(v.green());
        const double ch_y = converter(v.blue());
        const cmsJCh jch = {
            Scale::scale_lightness_of_pyramid(j),
            Scale::scale_chroma_of_pyramid(ch_x, ch_y),
            Scale::scale_hue_of_pyramid(ch_x, ch_y)
        };

        if (jch.J <= 0.0)
        {
            // Lasciate ogne speranza, voi ch'intrate.
            return DestVectorType(0, 0, 0);
        }

        double rgb[3];
        jch_to_rgb(&jch, rgb);

        // Implementation Notes
        //
        //         New LittleCMS versions use "open color space" arithmetics, which means color
        // coordinates can end up outside their domains, e.g. JCh cylinder or RGB cube.  We just
        // let LittleCMS chug along freewheeling as long as possible.  Right here, we must take
        // care of out-of-cube RGB values, because the array `rgb[3]' ends up as pixel in the
        // user's image.  We exert great care on the renegade pixels to preserve as much
        // information as possible; for example we *always* preserve the color component (h).
        //
        // (1) Overflow - at least one component is larger than one.
        //     Function: flexible_optimize_1d_2d()
        //     First try to find a visually similar pixel with different luminance (J).  If no
        //     pixel is close enough extend the search to luminance-saturation space (J, C).
        // (2) Underflow - at least one component less than zero.
        //     Function: multistart_optimize_2d()
        //     Search luminance-saturation space (J, C) for a visually similar pixel using many
        //     different initial luminance (J) values.
        //
        // Bear in mind that cases (1) and (2) are *not* disjoint!
        //
        // (3) Low-J pixels and the special cases rgb[i] == 0.0, i = 0, 1, 2.
        //     Some, but not all, RGB pixels that overflow, i.e. case (1), underflow, i.e. case
        //     (2), or have one or more components equal to zero are in fact dim shadows and
        //     neither bright highlights nor fancy sparkling green -- and infrequently --
        //     blue-ish dots.  We call them "shadows disguised as highlights".
        //
        //     These kind of deceivers *cannot* be roped in with our optimization strategies
        //     because their RGB values diverge more and more the closer we get to reasonable
        //     luminances.  The reason is that we always assume a non-zero saturation (C), the
        //     JCh model wants to compensate the low luminance, and the RGB components go
        //     haywire.  We solve the problem by expunging the saturation before launching the
        //     shadow optimizer.
        if (rgb[0] > 1.0 || rgb[1] > 1.0 || rgb[2] > 1.0)
        {
            if (jch.J <= shadow_disguised_as_highlight_j)
            {
                const cmsJCh ich {jch.J, 0.0, jch.h};
                multistart_optimize_2d(&ich, rgb);
            }
            else
            {
                flexible_optimize_1d_2d(&jch, rgb);
            }
        }
        else if (rgb[0] <= 0.0 || rgb[1] <= 0.0 || rgb[2] <= 0.0)
        {
            if (jch.J <= shadow_disguised_as_highlight_j)
            {
                const cmsJCh ich {jch.J, 0.0, jch.h};
                multistart_optimize_2d(&ich, rgb);
            }
            else
            {
                multistart_optimize_2d(&jch, rgb);
            }
        }

#ifdef DEBUG_STUBBORN_SHADOW_HIGHLIGHT
        if (rgb[0] > 1.0 || rgb[1] > 1.0 || rgb[2] > 1.0)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            {
                std::cout << "\n";
                ciecam_detail::show_jch_rgb("+ stubborn highlight:", &jch);
            }
        }

        if (rgb[0] < 0.0 || rgb[1] < 0.0 || rgb[2] < 0.0)
        {
#ifdef OPENMP
#pragma omp critical
#endif
            {
                std::cout << "\n";
                ciecam_detail::show_jch_rgb("+ stubborn shadow:", &jch);
            }
        }
#endif

        ciecam_detail::limit_sequence(rgb, rgb + 3U, 0.0, 1.0);

        return DestVectorType(DestTraits::fromRealPromote(rgb_dest_scale * rgb[0]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[1]),
                              DestTraits::fromRealPromote(rgb_dest_scale * rgb[2]));
    }

protected:
    ConvertFunctorType converter;
    const double rgb_dest_scale;

    const double highlight_lightness_guess_1d_factor;
    const double highlight_lightness_guess_1d_offset;
    const unsigned maximum_highlight_iterations;
    const unsigned maximum_highlight_bracket_tries;
    const double highlight_simplex_lightness_step_length;
    const double highlight_simplex_chroma_step_length;
    const unsigned highlight_iterations_per_leg;
    const unsigned maximum_highlight_leg;
    const double shadow_disguised_as_highlight_j;

    const double shadow_lightness_lightness_guess_factor;
    const double shadow_lightness_chroma_guess_factor;
    const double shadow_lightness_guess_offset;
    const double shadow_chroma_lightness_guess_factor;
    const double shadow_chroma_chroma_guess_factor;
    const double shadow_chroma_guess_offset;

    const double shadow_simplex_lightness_step_length;
    const double shadow_simplex_chroma_step_length;
    const unsigned shadow_iterations_per_leg;
    const unsigned maximum_shadow_leg;
    const unsigned maximum_multistart_tries;

    const double optimizer_error;
    const double optimizer_goal;
};


// Copy a scalar image into a scalar pyramid image.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyToPyramidImage(typename SrcImageType::const_traverser src_upperleft,
                   typename SrcImageType::const_traverser src_lowerright,
                   typename SrcImageType::ConstAccessor sa,
                   typename PyramidImageType::traverser dest_upperleft,
                   typename PyramidImageType::Accessor da,
                   vigra::VigraTrueType)
{
    typedef typename SrcImageType::value_type SrcPixelType;
    typedef typename PyramidImageType::value_type PyramidPixelType;

    vigra::omp::transformImage(src_upperleft, src_lowerright, sa,
                               dest_upperleft, da,
                               ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType,
                                                             PyramidIntegerBits, PyramidFractionBits>());
}


// Copy a vector image into a vector pyramid image.  Uses an optional color space conversion.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
void
copyToPyramidImage(typename SrcImageType::const_traverser src_upperleft,
                   typename SrcImageType::const_traverser src_lowerright,
                   typename SrcImageType::ConstAccessor sa,
                   typename PyramidImageType::traverser dest_upperleft,
                   typename PyramidImageType::Accessor da,
                   vigra::VigraFalseType)
{
    typedef typename SrcImageType::value_type SrcVectorType;
    typedef typename PyramidImageType::value_type PyramidVectorType;

    switch (BlendColorspace)
    {
    case UndeterminedColorspace:
    case IdentitySpace:
        vigra::omp::transformImage(src_upperleft, src_lowerright, sa,
                                   dest_upperleft, da,
                                   ConvertVectorToPyramidFunctor<SrcVectorType, PyramidVectorType,
                                                                 PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIELAB:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIELAB color conversion";
            if (!enblend::profileName(InputProfile).empty())
            {
                std::cerr << " from/to \"" << enblend::profileName(InputProfile) << "\" profile";
            }
            std::cerr << "\n";
        }
        vigra::omp::transformImage(src_upperleft, src_lowerright, sa,
                                   dest_upperleft, da,
                                   ConvertVectorToLabPyramidFunctor<SrcVectorType, PyramidVectorType,
                                                                    PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIELUV:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIELUV color conversion";
            if (!enblend::profileName(InputProfile).empty())
            {
                std::cerr << " from/to \"" << enblend::profileName(InputProfile) << "\" profile";
            }
            std::cerr << "\n";
        }
        vigra::omp::transformImage(src_upperleft, src_lowerright, sa,
                                   dest_upperleft, da,
                                   ConvertVectorToLuvPyramidFunctor<SrcVectorType, PyramidVectorType,
                                                                    PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIECAM:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIECAM02 color conversion";
            if (!enblend::profileName(InputProfile).empty())
            {
                std::cerr << " from/to \"" << enblend::profileName(InputProfile) << "\" profile";
            }
            std::cerr << "\n";
        }
        vigra::omp::transformImage(src_upperleft, src_lowerright, sa,
                                   dest_upperleft, da,
                                   ConvertVectorToJCHPyramidFunctor<SrcVectorType, PyramidVectorType,
                                                                    PyramidIntegerBits, PyramidFractionBits>());
        break;

    default:
        NEVER_REACHED("switch control expression \"BlendColorspace\" out of range");
    }
}


// Compile-time switch based on scalar or vector image type.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyToPyramidImage(typename SrcImageType::const_traverser src_upperleft,
                   typename SrcImageType::const_traverser src_lowerright,
                   typename SrcImageType::ConstAccessor sa,
                   typename PyramidImageType::traverser dest_upperleft,
                   typename PyramidImageType::Accessor da)
{
    typedef typename vigra::NumericTraits<typename SrcImageType::value_type>::isScalar src_is_scalar;

    copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>
        (src_upperleft, src_lowerright, sa,
         dest_upperleft, da,
         src_is_scalar());
}


// Version using argument object factories.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyToPyramidImage(vigra::triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
                   vigra::pair<typename PyramidImageType::traverser, typename PyramidImageType::Accessor> dest)
{
    copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>
        (src.first, src.second, src.third,
         dest.first, dest.second);
}


// Copy a scalar pyramid image into a scalar image.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType,
          int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyFromPyramidImageIf(typename PyramidImageType::const_traverser src_upperleft,
                       typename PyramidImageType::const_traverser src_lowerright,
                       typename PyramidImageType::ConstAccessor sa,
                       typename MaskImageType::const_traverser mask_upperleft,
                       typename MaskImageType::ConstAccessor ma,
                       typename DestImageType::traverser dest_upperleft,
                       typename DestImageType::Accessor da,
                       vigra::VigraTrueType)
{
    typedef typename DestImageType::value_type DestPixelType;
    typedef typename PyramidImageType::value_type PyramidPixelType;

    vigra::omp::transformImageIf(src_upperleft, src_lowerright, sa,
                                 mask_upperleft, ma,
                                 dest_upperleft, da,
                                 ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType,
                                                               PyramidIntegerBits, PyramidFractionBits>());
}


// Copy a vector pyramid image into a vector image.  Uses an optional
// color space conversion.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType,
          int PyramidIntegerBits, int PyramidFractionBits>
void
copyFromPyramidImageIf(typename PyramidImageType::const_traverser src_upperleft,
                       typename PyramidImageType::const_traverser src_lowerright,
                       typename PyramidImageType::ConstAccessor sa,
                       typename MaskImageType::const_traverser mask_upperleft,
                       typename MaskImageType::ConstAccessor ma,
                       typename DestImageType::traverser dest_upperleft,
                       typename DestImageType::Accessor da,
                       vigra::VigraFalseType)
{
    typedef typename DestImageType::value_type DestVectorType;
    typedef typename PyramidImageType::value_type PyramidVectorType;

    switch (BlendColorspace)
    {
    case UndeterminedColorspace:
    case IdentitySpace:
        // OpenMP changes the result here!  The maximum absolute
        // difference is 1 of 255 for 8-bit images.  -- cls
        vigra::omp::transformImageIf(src_upperleft, src_lowerright, sa,
                                     mask_upperleft, ma,
                                     dest_upperleft, da,
                                     ConvertPyramidToVectorFunctor<DestVectorType, PyramidVectorType,
                                                                   PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIELAB:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIELAB color conversion" << std::endl;
        }
        vigra::omp::transformImageIf(src_upperleft, src_lowerright, sa,
                                     mask_upperleft, ma,
                                     dest_upperleft, da,
                                     ConvertLabPyramidToVectorFunctor<DestVectorType, PyramidVectorType,
                                                                      PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIELUV:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIELUV color conversion" << std::endl;
        }
        vigra::omp::transformImageIf(src_upperleft, src_lowerright, sa,
                                     mask_upperleft, ma,
                                     dest_upperleft, da,
                                     ConvertLuvPyramidToVectorFunctor<DestVectorType, PyramidVectorType,
                                                                      PyramidIntegerBits, PyramidFractionBits>());
        break;

    case CIECAM:
        if (Verbose >= VERBOSE_COLOR_CONVERSION_MESSAGES)
        {
            std::cerr << command << ": info: CIECAM02 color conversion" << std::endl;
        }
        vigra::omp::transformImageIf(src_upperleft, src_lowerright, sa,
                                     mask_upperleft, ma,
                                     dest_upperleft, da,
                                     ConvertJCHPyramidToVectorFunctor<DestVectorType, PyramidVectorType,
                                                                      PyramidIntegerBits, PyramidFractionBits>());
        break;

    default:
        NEVER_REACHED("switch control expression \"BlendColorspace\" out of range");
    }
}


// Compile-time switch based on scalar or vector image type.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType,
          int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyFromPyramidImageIf(typename PyramidImageType::const_traverser src_upperleft,
                       typename PyramidImageType::const_traverser src_lowerright,
                       typename PyramidImageType::ConstAccessor sa,
                       typename MaskImageType::const_traverser mask_upperleft,
                       typename MaskImageType::ConstAccessor ma,
                       typename DestImageType::traverser dest_upperleft,
                       typename DestImageType::Accessor da)
{
    typedef typename vigra::NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

    copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType,
                           PyramidIntegerBits, PyramidFractionBits>
        (src_upperleft, src_lowerright, sa,
         mask_upperleft, ma,
         dest_upperleft, da,
         src_is_scalar());
}


// Version using argument object factories.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType,
          int PyramidIntegerBits, int PyramidFractionBits>
inline void
copyFromPyramidImageIf(vigra::triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
                       vigra::pair<typename MaskImageType::const_traverser, typename MaskImageType::ConstAccessor> mask,
                       vigra::pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest)
{
    copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType,
                           PyramidIntegerBits, PyramidFractionBits>
        (src.first, src.second, src.third,
         mask.first, mask.second,
         dest.first, dest.second);
}

} // namespace enblend

#endif // FIXMATH_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
