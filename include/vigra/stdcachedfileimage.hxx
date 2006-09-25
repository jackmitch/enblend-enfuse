/*
 *  Copyright (C) 2004 Andrew Mihal
 *
 *  This software is an extension of the VIGRA computer vision library.
 *  ( Version 1.2.0, Aug 07 2003 )
 *  You may use, modify, and distribute this software according
 *  to the terms stated in the LICENSE file included in
 *  the VIGRA distribution.
 *
 *  VIGRA is Copyright 1998-2002 by Ullrich Koethe
 *  Cognitive Systems Group, University of Hamburg, Germany
 *
 *  The VIGRA Website is
 *      http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/
 *  Please direct questions, bug reports, and contributions to
 *      koethe@informatik.uni-hamburg.de
 *
 *  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR
 *  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */
#ifndef VIGRA_EXT_STDCACHEDFILEIMAGE_HXX
#define VIGRA_EXT_STDCACHEDFILEIMAGE_HXX

#include <vigra/tuple.hxx>
#include <vigra/iteratortraits.hxx>
#include <vigra/accessor.hxx>
#include <vigra/rgbvalue.hxx>
#include <vigra/sized_int.hxx>
#include <vigra/cachedfileimage.hxx>

namespace vigra {

#define CFI_DEFINE_ITERATORTRAITS(VALUETYPE, ACCESSOR, CONSTACCESSOR) \
    template<> \
    struct IteratorTraits< \
        CachedFileImageIterator<VALUETYPE > > \
    { \
        typedef CachedFileImageIterator<VALUETYPE >  Iterator; \
        typedef Iterator                             iterator; \
        typedef iterator::iterator_category          iterator_category; \
        typedef iterator::value_type                 value_type; \
        typedef iterator::reference                  reference; \
        typedef iterator::index_reference            index_reference; \
        typedef iterator::pointer                    pointer; \
        typedef iterator::difference_type            difference_type; \
        typedef iterator::row_iterator               row_iterator; \
        typedef iterator::column_iterator            column_iterator; \
        typedef ACCESSOR<VALUETYPE >                 default_accessor; \
        typedef ACCESSOR<VALUETYPE >                 DefaultAccessor; \
    }; \
    template<> \
    struct IteratorTraits< \
        ConstCachedFileImageIterator<VALUETYPE > > \
    { \
        typedef \
          ConstCachedFileImageIterator<VALUETYPE >   Iterator; \
        typedef Iterator                             iterator; \
        typedef iterator::iterator_category          iterator_category; \
        typedef iterator::value_type                 value_type; \
        typedef iterator::reference                  reference; \
        typedef iterator::index_reference            index_reference; \
        typedef iterator::pointer                    pointer; \
        typedef iterator::difference_type            difference_type; \
        typedef iterator::row_iterator               row_iterator; \
        typedef iterator::column_iterator            column_iterator; \
        typedef CONSTACCESSOR<VALUETYPE >            default_accessor; \
        typedef CONSTACCESSOR<VALUETYPE >            DefaultAccessor; \
    }; \
    template<> \
    struct IteratorTraits< \
        StridedCachedFileImageIterator<VALUETYPE > > \
    { \
        typedef StridedCachedFileImageIterator<VALUETYPE > Iterator; \
        typedef Iterator                                   iterator; \
        typedef iterator::iterator_category                iterator_category; \
        typedef iterator::value_type                       value_type; \
        typedef iterator::reference                        reference; \
        typedef iterator::index_reference                  index_reference; \
        typedef iterator::pointer                          pointer; \
        typedef iterator::difference_type                  difference_type; \
        typedef iterator::row_iterator                     row_iterator; \
        typedef iterator::column_iterator                  column_iterator; \
        typedef ACCESSOR<VALUETYPE >                       default_accessor; \
        typedef ACCESSOR<VALUETYPE >                       DefaultAccessor; \
    }; \
    template<> \
    struct IteratorTraits< \
        ConstStridedCachedFileImageIterator<VALUETYPE > > \
    { \
        typedef ConstStridedCachedFileImageIterator<VALUETYPE > Iterator; \
        typedef Iterator                                        iterator; \
        typedef iterator::iterator_category                     iterator_category; \
        typedef iterator::value_type                            value_type; \
        typedef iterator::reference                             reference; \
        typedef iterator::index_reference                       index_reference; \
        typedef iterator::pointer                               pointer; \
        typedef iterator::difference_type                       difference_type; \
        typedef iterator::row_iterator                          row_iterator; \
        typedef iterator::column_iterator                       column_iterator; \
        typedef CONSTACCESSOR<VALUETYPE >                       default_accessor; \
        typedef CONSTACCESSOR<VALUETYPE >                       DefaultAccessor; \
    };

CFI_DEFINE_ITERATORTRAITS(UInt8, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<UInt8> UInt8CFImage;

CFI_DEFINE_ITERATORTRAITS(Int8, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<Int8> Int8CFImage;

CFI_DEFINE_ITERATORTRAITS(UInt16, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<UInt16> UInt16CFImage;

CFI_DEFINE_ITERATORTRAITS(Int16, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<Int16> Int16CFImage;

CFI_DEFINE_ITERATORTRAITS(UInt32, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<UInt32> UInt32CFImage;

CFI_DEFINE_ITERATORTRAITS(Int32, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<Int32> Int32CFImage;

CFI_DEFINE_ITERATORTRAITS(float, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<float> FCFImage;

CFI_DEFINE_ITERATORTRAITS(double, StandardValueAccessor, StandardConstValueAccessor)
typedef CachedFileImage<double> DCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<UInt8>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<UInt8> > UInt8RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<Int8>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<Int8> > Int8RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<UInt16>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<UInt16> > UInt16RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<Int16>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<Int16> > Int16RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<UInt32>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<UInt32> > UInt32RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<Int32>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<Int32> > Int32RGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<float>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<float> > FRGBCFImage;

CFI_DEFINE_ITERATORTRAITS(RGBValue<double>, RGBAccessor, RGBAccessor)
typedef CachedFileImage<RGBValue<double> > DRGBCFImage;

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

// define traits for BasicImageIterator instanciations that
// were not explicitly defined above
template <class T>
struct IteratorTraits<CachedFileImageIterator<T> >
{
    typedef CachedFileImageIterator<T>           Iterator;
    typedef Iterator                             iterator;
    typedef typename iterator::iterator_category iterator_category;
    typedef typename iterator::value_type        value_type;
    typedef typename iterator::reference         reference;
    typedef typename iterator::index_reference   index_reference;
    typedef typename iterator::pointer           pointer;
    typedef typename iterator::difference_type   difference_type;
    typedef typename iterator::row_iterator      row_iterator;
    typedef typename iterator::column_iterator   column_iterator;
    typedef StandardAccessor<T>                  DefaultAccessor;
    typedef StandardAccessor<T>                  default_accessor;
};

template <class T>
struct IteratorTraits<ConstCachedFileImageIterator<T> >
{
    typedef ConstCachedFileImageIterator<T>        Iterator;
    typedef Iterator                               iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
    typedef StandardConstAccessor<T>               DefaultAccessor;
    typedef StandardConstAccessor<T>               default_accessor;
};

#endif

} // namespace vigra

#endif /* VIGRA_EXT_STDCACHEDFILEIMAGE_HXX */
