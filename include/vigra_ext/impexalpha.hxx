/** @file impexalpha.hxx
 *
 *  Routines to save images with alpha masks.
 *
 *  These routines handle the conversion of byte alpha
 *  channels into the final output types.
 *
 *  @author Pablo d'Angelo <pablo.dangelo@web.de>
 *
 *  $Id: impexalpha.hxx,v 1.2 2007/01/27 05:00:36 acmihal Exp $
 *
 *  This is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this software; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef IMPEXALPHA_HXX_
#define IMPEXALPHA_HXX_

#include <vigra/imageiterator.hxx>
#include <vigra/transformimage.hxx>
#include <vigra/initimage.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/impex.hxx>


namespace vigra_ext {

// define values for mask true value. max for integers, 1 for floats
template <class T1> struct GetMaskTrue;

#define VIGRA_EXT_GETMASKTRUE(T1, S)            \
    template<>                                  \
    struct GetMaskTrue<T1>                      \
    {                                           \
        static T1 get()                         \
        {                                       \
            return S;                           \
        }                                       \
    }

#define VIGRA_EXT_GETMASKMAX(T1)                        \
    template<>                                          \
    struct GetMaskTrue<T1>                              \
    {                                                   \
        static T1 get()                                 \
        {                                               \
            return vigra::NumericTraits<T1>::max();     \
        }                                               \
    }

VIGRA_EXT_GETMASKMAX(vigra::UInt8);
VIGRA_EXT_GETMASKMAX(vigra::Int16);
VIGRA_EXT_GETMASKMAX(vigra::UInt16);
VIGRA_EXT_GETMASKMAX(vigra::Int32);
VIGRA_EXT_GETMASKMAX(vigra::UInt32);

VIGRA_EXT_GETMASKTRUE(float, 1.0f);
VIGRA_EXT_GETMASKTRUE(double, 1.0);


template <class ValueIterator, class ValueAccessor, class AlphaIterator, class AlphaAccessor>
class MultiImageMaskAccessor2
{
public:
    enum {STATIC_SIZE = 2};

    typedef typename ValueAccessor::value_type component_type;
    typedef typename AlphaAccessor::value_type alpha_type;
    typedef vigra::TinyVector<component_type, STATIC_SIZE> value_type;

    //typedef vigra::VectorElementAccessor<vigra::VectorAccessor<value_type> > ElementAccessor;
    //operator vigra::VectorAccessor<value_type>() {return a_;}

    MultiImageMaskAccessor2(ValueIterator i1, ValueAccessor a1, AlphaIterator i2, AlphaAccessor a2) :
        i1_(i1), a1_(a1), i2_(i2), a2_(a2)
    {}

    template <class Difference>
    value_type operator()(const Difference& d) const
    {
        return value_type(a1_(i1_, d), a2_(i2_, d));
    }

    template <class ValueDifference, class AlphaDifference>
    value_type operator()(ValueDifference d, const AlphaDifference& d2) const
    {
        d += d2;
        return value_type(a1_(i1_, d), a2_(i2_, d));
    }

    template <class Difference>
    value_type set(const value_type& v, const Difference& d) const
    {
        a1_.set(v[0], i1_, d);
        a2_.set(v[1], i2_, d);
    }

    template <class Value, class Iterator>
    void setComponent(const Value& value, const Iterator& i, int index) const
    {
        switch (index) {
        case 0:
            a1_.set(value, i1_, *i);
            break;
        case 1:
            a2_.set(value, i2_, *i);
            break;
        default:
            vigra_fail("MultiImageMaskAccessor2::setComponent: index out of range");
        }
    }

    template <class Iterator>
    component_type getComponent(const Iterator& i, int index) const
    {
        switch (index) {
        case 0:
            return a1_(i1_, *i);
        case 1:
            return a2_(i2_, *i);
        default:
            vigra_fail("MultiImageMaskAccessor2::getComponent: index out of range");
        }
    }

    template <class Iterator>
    unsigned size(const Iterator& i) const
    {
        return STATIC_SIZE;
    }

private:
    ValueIterator i1_;
    ValueAccessor a1_;
    AlphaIterator i2_;
    AlphaAccessor a2_;
};


template <class ValueIterator, class ValueAccessor, class AlphaIterator, class AlphaAccessor>
class MultiImageVectorMaskAccessor4
{
public:
    enum {STATIC_SIZE = 4};

    typedef typename ValueAccessor::value_type VT1;
    typedef typename AlphaAccessor::value_type alpha_type;

    typedef vigra::TinyVector<typename VT1::value_type, STATIC_SIZE> value_type;
    typedef typename value_type::value_type component_type;

    MultiImageVectorMaskAccessor4(ValueIterator i1, ValueAccessor a1, AlphaIterator i2, AlphaAccessor a2) :
        i1_(i1), a1_(a1), i2_(i2), a2_(a2)
    {}

    template <class Difference>
    value_type operator()(const Difference& d) const
    {
        const VT1& v1 = a1_.get(i1_, d);
        return value_type(v1[0], v1[1], v1[2],
                          a2_(i2_, d));
    }

    template <class ValueDifference, class AlphaDifference>
    value_type operator()(ValueDifference d, const AlphaDifference& d2) const
    {
        d += d2;
        const VT1& v1 = a1_.get(i1_, d);
        return value_type(v1[0], v1[1], v1[2],
                          a2_(i2_, d));
    }

    template <class Difference>
    value_type set(const value_type& v, const Difference& d) const
    {
        for (int index = 0; index != STATIC_SIZE; ++index) {
            a1_.setComponent(v[index], i1_, index);
        }
        a2_.set(v[3], i2_, d);
    }

    template <class Value, class Iterator>
    void setComponent(const Value& value, const Iterator& i, int index) const
    {
        if (index < STATIC_SIZE - 1) {
            a1_.setComponent(value, i1_, *i, index);
        } else if (index == STATIC_SIZE - 1) {
            a2_.set(value, i2_, *i);
        } else {
            vigra_fail("MultiImageVectorMaskAccessor4::setComponent: index out of range");
        }
    }

    template <class Iterator>
    component_type getComponent(const Iterator& i, int index) const
    {
        if (index < STATIC_SIZE - 1) {
            return a1_.getComponent(i1_, *i, index);
        } else if (index == STATIC_SIZE - 1) {
            return a2_(i2_, *i);
        } else {
            vigra_fail("MultiImageVectorMaskAccessor4::getComponent: index out of range");
        }
    }

    template <class Iterator>
    unsigned size (const Iterator& i) const
    {
        return STATIC_SIZE;
    }

private:
    ValueIterator i1_;
    ValueAccessor a1_;
    AlphaIterator i2_;
    AlphaAccessor a2_;
};


// scalar image
template<class SrcIterator, class SrcAccessor,
         class AlphaIterator, class AlphaAccessor>
void
exportImageAlpha(const vigra::triple<SrcIterator, SrcIterator, SrcAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha,
                 const vigra::ImageExportInfo& info,
                 vigra::VigraTrueType)
{
    typedef MultiImageMaskAccessor2<SrcIterator, SrcAccessor, AlphaIterator, AlphaAccessor> MultiAccessor;

    exportImage(vigra::CoordinateIterator(),
                vigra::CoordinateIterator(image.second - image.first),
                MultiAccessor(image.first, image.third, alpha.first, alpha.second),
                info);
}


// vector image
template<class SrcIterator, class SrcAccessor,
         class AlphaIterator, class AlphaAccessor>
void
exportImageAlpha(const vigra::triple<SrcIterator, SrcIterator, SrcAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha,
                 const vigra::ImageExportInfo& info,
                 vigra::VigraFalseType)
{
    typedef MultiImageVectorMaskAccessor4<SrcIterator, SrcAccessor, AlphaIterator, AlphaAccessor> MultiAccessor;

    exportImage(vigra::CoordinateIterator(),
                vigra::CoordinateIterator(image.second - image.first),
                MultiAccessor(image.first, image.third, alpha.first, alpha.second),
                info);
}


/** export an image with a differently typed alpha channel.
 *
 *  This function handles the merging of the images and the
 *  scales the alpha channel to the correct values.
 *
 *  can write to all output formats that support 4 channel images.
 *  (currently only png and tiff).
 */
template<class SrcIterator, class SrcAccessor,
         class AlphaIterator, class AlphaAccessor>
void
exportImageAlpha(const vigra::triple<SrcIterator, SrcIterator, SrcAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha,
                 const vigra::ImageExportInfo& info)
{
    typedef typename vigra::NumericTraits<typename SrcAccessor::value_type>::isScalar is_scalar;

    // Select function for scalar, or vector image, depending on
    // source type.  The alpha image has to be scalar all the time;
    // stuff will break with strange compile error if it is not.
    exportImageAlpha(image, alpha, info, is_scalar());
}


// vector image
template<class DestIterator, class DestAccessor,
         class AlphaIterator, class AlphaAccessor>
void
importImageAlpha(const vigra::ImageImportInfo& info,
                 const std::pair<DestIterator, DestAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha,
                 vigra::VigraFalseType)
{
    vigra_precondition(image.second(image.first).size() == 3,
                       "only scalar and 3-channel (i.e. vector) images supported by vigra_ext/impexalpha.hxx");

    typedef MultiImageVectorMaskAccessor4<DestIterator, DestAccessor, AlphaIterator, AlphaAccessor> MultiAccessor;

    importImage(info,
                vigra::CoordinateIterator(),
                MultiAccessor(image.first, image.second, alpha.first, alpha.second));
}

// scalar image
template<class DestIterator, class DestAccessor,
         class AlphaIterator, class AlphaAccessor>
void
importImageAlpha(const vigra::ImageImportInfo& info,
                 const std::pair<DestIterator, DestAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha,
                 vigra::VigraTrueType)
{
    typedef MultiImageMaskAccessor2<DestIterator, DestAccessor, AlphaIterator, AlphaAccessor> MultiAccessor;

    importImage(info,
                vigra::CoordinateIterator(),
                MultiAccessor(image.first, image.second, alpha.first, alpha.second));
}


/** import an image with a differently typed alpha channel.
 *
 *  This function loads an image, and splits it into a
 *  color image and a separate alpha channel, the alpha channel
 *  should be a 8 bit image.
 *
 *  If the image doesn't contain any alpha channel, a completely
 *  white one is created.
 *
 *  can write to all output formats that support 4 channel images.
 *  (currently only png and tiff).
 */
template<class DestIterator, class DestAccessor,
         class AlphaIterator, class AlphaAccessor>
void
importImageAlpha(const vigra::ImageImportInfo& info,
                 const vigra::pair<DestIterator, DestAccessor>& image,
                 const std::pair<AlphaIterator, AlphaAccessor>& alpha)
{
    typedef typename vigra::NumericTraits<typename DestAccessor::value_type>::isScalar is_scalar;

    if (info.numExtraBands() == 1) {
	// import image and alpha channel
	importImageAlpha(info, image, alpha, is_scalar());
    } else if (info.numExtraBands() == 0) {
	// no alpha channel in file, import as usual.
	importImage(info, image);
	// fill alpha image
	vigra::initImage(alpha.first,
                         alpha.first + vigra::Diff2D(info.width(), info.height()),
                         alpha.second,
                         255);
    } else {
	vigra_fail("Images with two or more alpha channels are not supported");
    }
}

} // namespace vigra_ext

#endif // IMPEXALPHA_HXX
