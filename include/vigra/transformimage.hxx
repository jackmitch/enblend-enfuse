/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.2.0, Aug 07 2003 )                                    */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_TRANSFORMIMAGE_HXX
#define VIGRA_TRANSFORMIMAGE_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/rgbvalue.hxx"

namespace vigra {

/** \addtogroup TransformAlgo Algorithms to Transform Images
    Apply functor to calculate a pixelwise transformation of one image

    @{
*/

/********************************************************/
/*                                                      */
/*                      transformLine                   */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class Functor>
void
transformLine(SrcIterator s,
              SrcIterator send, SrcAccessor src,
              DestIterator d, DestAccessor dest,
              Functor const & f)
{
    for(; s != send; ++s, ++d)
        dest.set(f(src(s)), d);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor,
          class Functor>
void
transformLineIf(SrcIterator s,
                SrcIterator send, SrcAccessor src,
                MaskIterator m, MaskAccessor mask,
                DestIterator d, DestAccessor dest,
                Functor const & f)
{
    for(; s != send; ++s, ++d, ++m)
        if(mask(m))
            dest.set(f(src(s)), d);
}

/********************************************************/
/*                                                      */
/*                      transformImage                  */
/*                                                      */
/********************************************************/

/** \brief Apply unary point transformation to each pixel.

    The transformation given by the functor is applied to every source
    pixel and the result written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    Note that the unary functors of the STL can be used in addition to
    the functors specifically defined in \ref TransformFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor, class Functor>
        void
        transformImage(SrcImageIterator src_upperleft,
               SrcImageIterator src_lowerright, SrcAccessor sa,
               DestImageIterator dest_upperleft, DestAccessor da,
               Functor const & f)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor, class Functor>
        void
        transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
               pair<DestImageIterator, DestAccessor> dest,
               Functor const & f)
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code

    #include <math.h>         // for sqrt()

    vigra::transformImage(srcImageRange(src),
                          destImage(dest),
                          &::sqrt );

    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator      dest_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    Functor functor;

    dest_accessor.set(functor(src_accessor(sx)), dx);

    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
void
transformImage(SrcImageIterator src_upperleft,
               SrcImageIterator src_lowerright, SrcAccessor sa,
               DestImageIterator dest_upperleft, DestAccessor da,
               Functor const & f)
{
    int w = src_lowerright.x - src_upperleft.x;

    for(; src_upperleft.y < src_lowerright.y; ++src_upperleft.y, ++dest_upperleft.y)
    {
        transformLine(src_upperleft.rowIterator(),
                      src_upperleft.rowIterator() + w, sa,
                      dest_upperleft.rowIterator(), da, f);
    }
}

template <class SrcImageIterator, class SrcAccessor,
      class DestImageIterator, class DestAccessor, class Functor>
inline
void
transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
               pair<DestImageIterator, DestAccessor> dest,
               Functor const & f)
{
    transformImage(src.first, src.second, src.third,
                   dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                   transformImageIf                   */
/*                                                      */
/********************************************************/

/** \brief Apply unary point transformation to each pixel within the ROI
    (i.e., where the mask is non-zero).

    The transformation given by the functor is applied to every source
    pixel in the ROI (i.e. when the return vlaue of the mask's accessor
    is not zero)
    and the result is written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    Note that the unary functors of the STL can be used in addition to
    the functors specifically defined in \ref TransformFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        transformImageIf(SrcImageIterator src_upperleft,
            SrcImageIterator src_lowerright, SrcAccessor sa,
            MaskImageIterator mask_upperleft, MaskAccessor ma,
            DestImageIterator dest_upperleft, DestAccessor da,
            Functor const & f)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                 pair<MaskImageIterator, MaskAccessor> mask,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code
    #include <math.h>         // for sqrt()

    vigra::transformImageIf(srcImageRange(src),
                            maskImage(mask),
                            destImage(dest),
                            &::sqrt );

    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator  dest_upperleft;
    MaskImageIterator mask_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    MaskAccessor mask_accessor;
    Functor functor;

    if(mask_accessor(mx))
       dest_accessor.set(functor(src_accessor(sx)), dx);

    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
transformImageIf(SrcImageIterator src_upperleft,
                 SrcImageIterator src_lowerright, SrcAccessor sa,
                 MaskImageIterator mask_upperleft, MaskAccessor ma,
                 DestImageIterator dest_upperleft, DestAccessor da,
                 Functor const & f)
{
    int w = src_lowerright.x - src_upperleft.x;

    for(; src_upperleft.y < src_lowerright.y;
             ++src_upperleft.y, ++mask_upperleft.y, ++dest_upperleft.y)
    {
        transformLineIf(src_upperleft.rowIterator(),
                        src_upperleft.rowIterator() + w, sa,
                        mask_upperleft.rowIterator(), ma,
                        dest_upperleft.rowIterator(), da, f);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                 pair<MaskImageIterator, MaskAccessor> mask,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
{
    transformImageIf(src.first, src.second, src.third,
                     mask.first, mask.second,
                     dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*               gradientBasedTransform                 */
/*                                                      */
/********************************************************/

/** \brief Calculate a function of the image gradient.

    The gradient and the function represented by <TT>Functor f</TT>
    are calculated in one go: for each location, the symmetric
    difference in x- and y-directions (asymmetric difference at the
    image borders) are passed to the given functor, and the result is
    written the destination image. Functors to be used with this
    function include \ref MagnitudeFunctor and \ref
    RGBGradientMagnitudeFunctor.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        gradientBasedTransform(SrcImageIterator srcul, SrcImageIterator srclr, SrcAccessor sa,
                      DestImageIterator destul, DestAccessor da, Functor const & f)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        gradientBasedTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                       pair<DestImageIterator, DestAccessor> dest, Functor const & const & f)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"


    \code
    vigra::FImage src(w,h), magnitude(w,h);
    ...

    gradientBasedTransform(srcImageRange(src), destImage(magnitude),
                                vigra::MagnitudeFunctor<float>());
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator is, isend;
    DestImageIterator id;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    typename NumericTraits<typename SrcAccessor::value_type>::RealPromote
        diffx, diffy;

    diffx = src_accessor(is, Diff2D(-1,0)) - src_accessor(is, Diff2D(1,0));
    diffy = src_accessor(is, Diff2D(0,-1)) - src_accessor(is, Diff2D(0,1));

    Functor f;

    dest_accessor.set(f(diffx, diffy), id);

    \endcode

*/

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
void
gradientBasedTransform(SrcImageIterator srcul, SrcImageIterator srclr, SrcAccessor sa,
                       DestImageIterator destul, DestAccessor da, Functor const & grad)
{
    int w = srclr.x - srcul.x;
    int h = srclr.y - srcul.y;
    int x,y;

    SrcImageIterator sy = srcul;
    DestImageIterator dy = destul;

    static const Diff2D left(-1,0);
    static const Diff2D right(1,0);
    static const Diff2D top(0,-1);
    static const Diff2D bottom(0,1);

    typename NumericTraits<typename SrcAccessor::value_type>::RealPromote
             diffx, diffy;

    SrcImageIterator sx = sy;
    DestImageIterator dx = dy;

    diffx = sa(sx) - sa(sx, right);
    diffy = sa(sx) - sa(sx, bottom);
    da.set(grad(diffx, diffy), dx);

    for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
    {
        diffx = (sa(sx, left) - sa(sx, right)) / 2.0;
        diffy = sa(sx) - sa(sx, bottom);
        da.set(grad(diffx, diffy), dx);
    }

    diffx = sa(sx, left) - sa(sx);
    diffy = sa(sx) - sa(sx, bottom);
    da.set(grad(diffx, diffy), dx);

    ++sy.y;
    ++dy.y;

    for(y=2; y<h; ++y, ++sy.y, ++dy.y)
    {
        sx = sy;
        dx = dy;

        diffx = sa(sx) - sa(sx, right);
        diffy = (sa(sx, top) - sa(sx, bottom)) / 2.0;
        da.set(grad(diffx, diffy), dx);

        for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
        {
            diffx = (sa(sx, left) - sa(sx, right)) / 2.0;
            diffy = (sa(sx, top) - sa(sx, bottom)) / 2.0;
            da.set(grad(diffx, diffy), dx);
        }

        diffx = sa(sx, left) - sa(sx);
        diffy = (sa(sx, top) - sa(sx, bottom)) / 2.0;
        da.set(grad(diffx, diffy), dx);
    }

    sx = sy;
    dx = dy;

    diffx = sa(sx) - sa(sx, right);
    diffy = sa(sx, top) - sa(sx);
    da.set(grad(diffx, diffy), dx);

    for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
    {
        diffx = (sa(sx, left) - sa(sx, right)) / 2.0;
        diffy = sa(sx, top) - sa(sx);
        da.set(grad(diffx, diffy), dx);
    }

    diffx = sa(sx, left) - sa(sx);
    diffy = sa(sx, top) - sa(sx);
    da.set(grad(diffx, diffy), dx);
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
inline
void
gradientBasedTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                       pair<DestImageIterator, DestAccessor> dest, Functor const & grad)
{
    gradientBasedTransform(src.first, src.second, src.third,
                           dest.first, dest.second, grad);
}

/** @} */
/** \addtogroup TransformFunctor Functors to Transform Images

    Note that the unary functors of the STL can also be used in
    connection with \ref transformImage().
*/
//@{

template <class DestValueType, class Multiplier = double>
class LinearIntensityTransform
{
  public:
        /* the functors argument type (actually, since 
           <tt>operator()</tt> is a template, much more types are possible)
        */
    typedef DestValueType argument_type;

        /* the functors result type
        */
    typedef DestValueType result_type;

        /* \deprecated use argument_type and result_type
        */
    typedef DestValueType value_type;

        /* type of the offset (used in internal calculations to prevent
            overflows and minimize round-off errors).
        */
    typedef typename
            NumericTraits<DestValueType>::RealPromote argument_promote;

        /* type of the scale factor
        */
    typedef Multiplier scalar_multiplier_type;

        /* init scale and offset
        */
    LinearIntensityTransform(scalar_multiplier_type scale, argument_promote offset)
    : scale_(scale), offset_(offset)
    {}

        /* calculate transform
        */
    template <class SrcValueType>
    result_type operator()(SrcValueType const & s) const
    {
        return NumericTraits<result_type>::fromRealPromote(scale_ * (s + offset_));
    }

  private:

    scalar_multiplier_type scale_;
    argument_promote offset_;
};


template <class DestValueType, class Multiplier = double>
class ScalarIntensityTransform
{
  public:
        /* the functors argument type (actually, since 
           <tt>operator()</tt> is a template, much more types are possible)
        */
    typedef DestValueType argument_type;

        /* the functors result type
        */
    typedef DestValueType result_type;

        /* \deprecated use argument_type and result_type
        */
    typedef DestValueType value_type;

        /* type of the scale factor
        */
    typedef Multiplier scalar_multiplier_type;

        /* init scale
        */
    ScalarIntensityTransform(scalar_multiplier_type scale)
    : scale_(scale)
    {}

        /* calculate transform
        */
    template <class SrcValueType>
    result_type operator()(SrcValueType const & s) const
    {
        return NumericTraits<result_type>::fromRealPromote(scale_ * s);
    }

  private:
    scalar_multiplier_type scale_;
};

/********************************************************/
/*                                                      */
/*              linearIntensityTransform                */
/*                                                      */
/********************************************************/

/** \brief Apply a linear transform to the source pixel values

    Factory function for a functor that linearly transforms the
    source pixel values. The functor applies the transform
    '<TT>destvalue = scale * (srcvalue + offset)</TT>' to every pixel.
    This can, for example, be used to transform images into the visible
    range 0...255 or to invert an image.

    If you leave out the second parameter / offset, you will get an
    optimized version of the functor which only scales by the given
    factor, however you have to make the template parameter (pixel
    type) explicit.

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <class Multiplier, class DestValueType>
        LinearIntensityTransform<DestValueType, Multiplier>
        linearIntensityTransform(Multiplier scale, DestValueType offset);

        template <class DestValueType, class Multiplier>
        ScalarIntensityTransform<DestValueType, Multiplier>
        linearIntensityTransform(Multiplier scale);
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::IImage src(width, height);
    vigra::BImage dest(width, height);
    ...
    vigra::FindMinMax<IImage::PixelType> minmax;   // functor to find range

    vigra::inspectImage(srcImageRange(src), minmax); // find original range

    // transform to range 0...255
    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearIntensityTransform(
                            255.0 / (minmax.max - minmax.min), // scaling
                          - minmax.min));                    // offset
    \endcode

	The one-parameter version can be used like this:

    \code
	// scale from 0..255 to 0..1.0
	FImage dest(src.size());

    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearIntensityTransform<float>(1.0 / 255));
    \endcode

    <b> Required Interface:</b>

    The source and destination value types must be models of \ref LinearSpace in both cases.

*/
template <class Multiplier, class DestValueType>
LinearIntensityTransform<DestValueType, Multiplier>
linearIntensityTransform(Multiplier scale, DestValueType offset)
{
    return LinearIntensityTransform<DestValueType, Multiplier>(scale, offset);
}

template <class DestValueType, class Multiplier>
ScalarIntensityTransform<DestValueType, Multiplier>
linearIntensityTransform(Multiplier scale)
{
    return ScalarIntensityTransform<DestValueType, Multiplier>(scale);
}

/********************************************************/
/*                                                      */
/*                   linearRangeMapping                 */
/*                                                      */
/********************************************************/

/** \brief Map a source intensity range linearly to a destination range.

    Factory function for a functor that linearly transforms the
    source pixel values. The functor applies the transform
    '<TT>destvalue = scale * (srcvalue + offset)</TT>' to every pixel,
    where <tt>scale = (dest_max - dest_min) / (src_max - src_min)</tt>
    and <tt>offset = dest_min / scale - src_min</tt>. As a result,
    the pixel values <tt>src_max</tt>, <tt>src_min</tt> in the source image 
    are mapped onto <tt>dest_max</tt>, <tt>dest_min</tt> respectively. 
    This works for scalar as well as vector pixel types.

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <class SrcValueType, class DestValueType>
        LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
        linearRangeMapping(SrcValueType src_min, SrcValueType src_max, 
                           DestValueType dest_min, DestValueType dest_max );
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::IImage src(width, height);
    vigra::BImage dest(width, height);
    ...
    vigra::FindMinMax<IImage::PixelType> minmax;   // functor to find range

    vigra::inspectImage(srcImageRange(src), minmax); // find original range

    // transform to range 0...255
    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearRangeTransform(
                            minmax.min, minmax.max,               // src range
                            (unsigned char)0, (unsigned char)255) // dest range
                          );
    \endcode

    <b> Required Interface:</b>

    The source and destination value types must be models of \ref LinearSpace in both cases.

*/
template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(SrcValueType src_min, SrcValueType src_max, 
                   DestValueType dest_min, DestValueType dest_max )
{
    return linearRangeMapping(src_min, src_max, dest_min, dest_max,
            typename NumericTraits<DestValueType>::isScalar());
} 

template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(
    SrcValueType src_min, SrcValueType src_max, 
    DestValueType dest_min, DestValueType dest_max,
    VigraTrueType /* isScalar */ )
{
    typedef typename NumericTraits<DestValueType>::RealPromote Multiplier;
    Multiplier diff = src_max - src_min;
    Multiplier scale = diff == NumericTraits<Multiplier>::zero()
                     ? NumericTraits<Multiplier>::one() 
                     : (dest_max - dest_min) / diff;
    return LinearIntensityTransform<DestValueType, Multiplier>(
                                   scale, dest_min / scale - src_min );
} 

template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(
    SrcValueType src_min, SrcValueType src_max, 
    DestValueType dest_min, DestValueType dest_max,
    VigraFalseType /* isScalar */ )
{
    typedef typename NumericTraits<DestValueType>::RealPromote Multiplier;
    typedef typename Multiplier::value_type MComponent;
    Multiplier scale(dest_max), offset(dest_max);
    for(int i=0; i<src_min.size(); ++i)
    { 
        MComponent diff = src_max[i] - src_min[i];
        scale[i] = diff == NumericTraits<MComponent>::zero()
                     ? NumericTraits<MComponent>::one() 
                     : (dest_max[i] - dest_min[i]) / diff;
        offset[i] = dest_min[i] / scale[i] - src_min[i];
    }
    return LinearIntensityTransform<DestValueType, Multiplier>(scale, offset);
} 

/********************************************************/
/*                                                      */
/*                      Threshold                       */
/*                                                      */
/********************************************************/

/** \brief Threshold an image.

    If a source pixel is above or equal the lower and below
    or equal the higher threshold (i.e. within the closed interval
    [lower, heigher]) the destination pixel is set to 'yesresult',
    otherwise to 'noresult'.

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage src(width, height), dest(width, height);
    ...
    vigra::transformImage(src.upperLeft(), src.lowerRight(), src.accessor(),
       dest.upperLeft(), dest.accessor(),
       vigra::Threshold<
          vigra::BImage::PixelType, vigra::BImage::PixelType>(10, 100, 0, 255));

    \endcode

    <b> Required Interface:</b>

    \code

    SrcValueType   src;
    DestValueType  dest, yesresult, noresult;

    dest = ((src < lower) || (higher < src)) ? noresult : yesresult;

    \endcode

*/
template <class SrcValueType, class DestValueType>
class Threshold
{
   public:

        /** the functor's argument type
        */
    typedef SrcValueType argument_type;

        /** the functor's result type
        */
    typedef DestValueType result_type;

        /** init thresholds and return values
        */
    Threshold(argument_type lower, argument_type higher,
              result_type noresult, result_type yesresult)
    : lower_(lower), higher_(higher),
      yesresult_(yesresult), noresult_(noresult)
    {}

        /** calculate transform
        */
    result_type operator()(argument_type s) const
    {
        return ((s < lower_) || (higher_ < s)) ? noresult_ : yesresult_;
    }

  private:

    argument_type lower_, higher_;
    result_type yesresult_, noresult_;
};

/********************************************************/
/*                                                      */
/*                BrightnessContrastFunctor             */
/*                                                      */
/********************************************************/

/** \brief Adjust brightness and contrast of an image.

    This functor applies a gamma correction to each pixel in order to
    modify the brightness of the image. To the result of the gamma correction,
    another transform is applied that modifies the contrast. The brightness and
    contrast parameters must be positive. Values greater than 1 will increase image
    brightness and contrast, values smaller than 1 decrease them. A value = 1 will
    have no effect.
    For \ref RGBValue "RGBValue's", the transforms are applied component-wise. The pixel
    values are assumed to lie between the given minimum and maximum
    values. In case of RGB, this is again understood component-wise. In case
    of <TT>unsigned char</TT>, min and max default to 0 and 255 respectively.
    Precisely, the following transform is applied to each <em> PixelValue</em>:

    \f[
    \begin{array}{rcl}
    V_1 & = & \frac{PixelValue - min}{max - min} \\
    V_2 & = & V_1^\frac{1}{brightness} \\
    V_3 & = & 2 V_2 - 1 \\
    V_4 & = & \left\lbrace
        \begin{array}{l}
         V_3^\frac{1}{contrast} \mbox{\rm \quad if  } V_3 \ge 0 \\
         - (-V_3)^\frac{1}{contrast} \mbox{\rm \quad otherwise}
        \end{array} \right. \\
    Result & = & \frac{V_4 + 1}{2} (max - min) + min
    \end{array}
    \f]

    If the <TT>PixelType</TT> is <TT>unsigned char</TT>, a look-up-table is used
    for faster computation.

    <b> Usage:</b>

        <b>\#include</b> "<a href="transformimage_8hxx-source.html">vigra/transformimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage bimage(width, height);
    double brightness, contrast;
    ...
    vigra::transformImage(srcImageRange(bimage), destImage(bimage),
       vigra::BrightnessContrastFunctor<unsigned char>(brightness, contrast));



    vigra::FImage fimage(width, height);
    ...

    vigra::FindMinmax<float> minmax;
    vigra::inspectImage(srcImageRange(fimage), minmax);

    vigra::transformImage(srcImageRange(fimage), destImage(fimage),
       vigra::BrightnessContrastFunctor<float>(brightness, contrast, minmax.min, minmax.max));


    \endcode

    <b> Required Interface:</b>

    Scalar types: must be a linear algebra (+, - *, NumericTraits),
    strict weakly ordered (<), and <TT>pow()</TT> must be defined.

    RGB values: the component type must meet the above requirements.
*/
template <class PixelType>
class BrightnessContrastFunctor
{
    typedef typename
        NumericTraits<PixelType>::RealPromote promote_type;

 public:

        /** the functor's argument type
        */
    typedef PixelType argument_type;

        /** the functor's result type
        */
    typedef PixelType result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef PixelType value_type;

        /** Init functor for argument range <TT>[min, max]</TT>.
            <TT>brightness</TT> and <TT>contrast</TT> values > 1 will
            increase brightness and contrast, < 1 will decrease them, and == 1 means
            no change.
        */
    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              argument_type const & min, argument_type const & max)
    : b_(1.0/brightness),
      c_(1.0/contrast),
      min_(min),
      diff_(max - min),
      zero_(NumericTraits<promote_type>::zero()),
      one_(NumericTraits<promote_type>::one())
    {}

        /** Calculate modified gray or color value
        */
    result_type operator()(argument_type const & v) const
    {
        promote_type v1 = (v - min_) / diff_;
        promote_type brighter = pow(v1, b_);
        promote_type v2 = 2.0 * brighter - one_;
        promote_type contrasted = (v2 < zero_) ?
                                     -pow(-v2, c_) :
                                      pow(v2, c_);
        return result_type(0.5 * diff_ * (contrasted + one_) + min_);
    }

  private:
    promote_type b_, c_;
    argument_type min_;
    promote_type diff_, zero_, one_;
};

template <>
class BrightnessContrastFunctor<unsigned char>
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
     unsigned char lut[256];

 public:

    typedef unsigned char value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min = 0, value_type const & max = 255)
    {
        BrightnessContrastFunctor<promote_type> f(brightness, contrast, min, max);

        for(int i = min; i <= max; ++i)
        {
            lut[i] = static_cast<unsigned char>(f(i)+0.5);
        }
    }

    value_type operator()(value_type const & v) const
    {

        return lut[v];
    }
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class ComponentType>
class BrightnessContrastFunctor<RGBValue<ComponentType> >
{
    typedef typename
        NumericTraits<ComponentType>::RealPromote promote_type;
    BrightnessContrastFunctor<ComponentType> red, green, blue;

 public:

    typedef RGBValue<ComponentType> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class BrightnessContrastFunctor<RGBValue<int> >
{
    typedef NumericTraits<int>::RealPromote promote_type;
    BrightnessContrastFunctor<int> red, green, blue;

 public:

    typedef RGBValue<int> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

template <>
class BrightnessContrastFunctor<RGBValue<float> >
{
    typedef NumericTraits<float>::RealPromote promote_type;
    BrightnessContrastFunctor<float> red, green, blue;

 public:

    typedef RGBValue<float> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class BrightnessContrastFunctor<RGBValue<unsigned char> >
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
    BrightnessContrastFunctor<unsigned char> red, green, blue;

 public:

    typedef RGBValue<unsigned char> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
       value_type const & min = value_type(0,0,0),
       value_type const & max = value_type(255, 255, 255))
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

//@}

} // namespace vigra

#endif // VIGRA_TRANSFORMIMAGE_HXX