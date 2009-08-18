/*
 * Copyright (C) 2009 Christoph L. Spiel
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
#ifndef __OPENMP_H__
#define __OPENMP_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vigra/diff2d.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/convolution.hxx"
#include "vigra/distancetransform.hxx"


#if _OPENMP >= 200505 // at least OpenMP version 2.5

#include <omp.h>


#define OPENMP
#define OPENMP_YEAR (_OPENMP / 100)
#define OPENMP_MONTH (_OPENMP % 100)


// These are the image sizes (measured in pixels) where we switch from
// serial (single thread) to multi processing.  The crossover points
// can be different for scalar, i.e. black-and-white images and
// non-scalar, i.e. RGB images.

#define CROSSOVER_COMBINETWOIMAGES_SCALAR 65536
#define CROSSOVER_COMBINETWOIMAGES_NON_SCALAR 16384

#define CROSSOVER_COMBINETHREEIMAGES_SCALAR 46656
#define CROSSOVER_COMBINETHREEIMAGES_NON_SCALAR 12544

#define CROSSOVER_TRANSFORMIMAGE_SCALAR 57600
#define CROSSOVER_TRANSFORMIMAGE_NON_SCALAR 32768

#define CROSSOVER_DISTANCE_TRANSFORM_SCALAR 1
#define CROSSOVER_DISTANCE_TRANSFORM_NON_SCALAR 1


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                   DestImageIterator dest_upperleft, DestAccessor dest_acc,
                   const Functor& func)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src1_lowerright - src1_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_COMBINETWOIMAGES_SCALAR : CROSSOVER_COMBINETWOIMAGES_NON_SCALAR))
    {
#pragma omp parallel
        {
            const int n = omp_get_num_threads();
            const int i = omp_get_thread_num();
            const vigra::Diff2D begin(0, (i * size.y) / n);
            const vigra::Diff2D end(size.x, ((i + 1) * size.y) / n);

            vigra::combineTwoImages(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                    src2_upperleft + begin, src2_acc,
                                    dest_upperleft + begin, dest_acc,
                                    func);
        } // omp parallel
    }
    else
    {
        vigra::combineTwoImages(src1_upperleft, src1_lowerright, src1_acc,
                                src2_upperleft, src2_acc,
                                dest_upperleft, dest_acc,
                                func);
    }
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesIfMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                     SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                     MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                     DestImageIterator dest_upperleft, DestAccessor dest_acc,
                     const Functor& func)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src1_lowerright - src1_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_COMBINETWOIMAGES_SCALAR : CROSSOVER_COMBINETWOIMAGES_NON_SCALAR))
    {
#pragma omp parallel
        {
            const int n = omp_get_num_threads();
            const int i = omp_get_thread_num();
            const vigra::Diff2D begin(0, (i * size.y) / n);
            const vigra::Diff2D end(size.x, ((i + 1) * size.y) / n);

            vigra::combineTwoImagesIf(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                      src2_upperleft + begin, src2_acc,
                                      mask_upperleft + begin, mask_acc,
                                      dest_upperleft + begin, dest_acc,
                                      func);
        } // omp parallel
    }
    else
    {
        vigra::combineTwoImagesIf(src1_upperleft, src1_lowerright, src1_acc,
                                  src2_upperleft, src2_acc,
                                  mask_upperleft, mask_acc,
                                  dest_upperleft, dest_acc,
                                  func);
    }
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineThreeImagesMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                     SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                     SrcImageIterator3 src3_upperleft, SrcAccessor3 src3_acc,
                     DestImageIterator dest_upperleft, DestAccessor dest_acc,
                     const Functor& func)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src1_lowerright - src1_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_COMBINETHREEIMAGES_SCALAR : CROSSOVER_COMBINETHREEIMAGES_NON_SCALAR))
    {
#pragma omp parallel
        {
            const int n = omp_get_num_threads();
            const int i = omp_get_thread_num();
            const vigra::Diff2D begin(0, (i * size.y) / n);
            const vigra::Diff2D end(size.x, ((i + 1) * size.y) / n);

            vigra::combineThreeImages(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                      src2_upperleft + begin, src2_acc,
                                      src3_upperleft + begin, src3_acc,
                                      dest_upperleft + begin, dest_acc,
                                      func);
        } // omp parallel
    }
    else
    {
        vigra::combineThreeImages(src1_upperleft, src1_lowerright, src1_acc,
                                  src2_upperleft, src2_acc,
                                  src3_upperleft, src3_acc,
                                  dest_upperleft, dest_acc,
                                  func);
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                 DestImageIterator dest_upperleft, DestAccessor dest_acc,
                 const Functor& func)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src_lowerright - src_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_TRANSFORMIMAGE_SCALAR : CROSSOVER_TRANSFORMIMAGE_NON_SCALAR))
    {
#pragma omp parallel
        {
            const int n = omp_get_num_threads();
            const int i = omp_get_thread_num();
            const vigra::Diff2D begin(0, (i * size.y) / n);
            const vigra::Diff2D end(size.x, ((i + 1) * size.y) / n);

            vigra::transformImage(src_upperleft + begin, src_upperleft + end, src_acc,
                                  dest_upperleft + begin, dest_acc,
                                  func);
        } // omp parallel
    }
    else
    {
        vigra::transformImage(src_upperleft, src_lowerright, src_acc,
                              dest_upperleft, dest_acc,
                              func);
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageIfMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                   MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                   DestImageIterator dest_upperleft, DestAccessor dest_acc,
                   const Functor& func)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src_lowerright - src_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_TRANSFORMIMAGE_SCALAR : CROSSOVER_TRANSFORMIMAGE_NON_SCALAR))
    {
#pragma omp parallel
        {
            const int n = omp_get_num_threads();
            const int i = omp_get_thread_num();
            const vigra::Diff2D begin(0, (i * size.y) / n);
            const vigra::Diff2D end(size.x, ((i + 1) * size.y) / n);

            vigra::transformImageIf(src_upperleft + begin, src_upperleft + end, src_acc,
                                    mask_upperleft + begin, mask_acc,
                                    dest_upperleft + begin, dest_acc,
                                    func);
        } // omp parallel
    }
    else
    {
        vigra::transformImageIf(src_upperleft, src_lowerright, src_acc,
                                mask_upperleft, mask_acc,
                                dest_upperleft, dest_acc,
                                func);
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class ValueType>
void
distanceTransformMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor sa,
                    DestImageIterator dest_upperleft, DestAccessor da,
                    ValueType background, int norm)
{
    typedef typename DestAccessor::value_type value_type;
    typedef typename vigra::NumericTraits<value_type>::isScalar isScalar;

    const vigra::Diff2D size(src_lowerright - src_upperleft);

    if (size.x * size.y >=
        (isScalar().asBool ? CROSSOVER_DISTANCE_TRANSFORM_SCALAR : CROSSOVER_DISTANCE_TRANSFORM_NON_SCALAR))
    {
        distanceTransform(src_upperleft, src_lowerright, sa,
                          dest_upperleft, da,
                          background, norm);
    }
    else
    {
        distanceTransform(src_upperleft, src_lowerright, sa,
                          dest_upperleft, da,
                          background, norm);
    }
}


#else


#undef OPENMP
#define OPENMP_YEAR 0
#define OPENMP_MONTH 0

inline void omp_set_num_threads(int) {}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_max_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
inline int omp_get_num_procs() {return 1;}
inline void omp_set_dynamic(int) {}
inline int omp_get_dynamic() {return 0;}
inline int omp_in_parallel() {return 0;}
inline void omp_set_nested(int) {}
inline int omp_get_nested() {return 0;}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesMP(SrcImageIterator1 src1_upperleft,
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                   DestImageIterator dest_upperleft, DestAccessor dest_acc,
                   const Functor& func)
{
    vigra::combineTwoImages(src1_upperleft, src1_lowerright, src1_acc,
                            src2_upperleft, src2_acc,
                            dest_upperleft, dest_acc,
                            func);
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesIfMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                     SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                     MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                     DestImageIterator dest_upperleft, DestAccessor dest_acc,
                     const Functor& func)
{
    vigra::combineTwoImagesIf(src1_upperleft, src1_lowerright, src1_acc,
                              src2_upperleft, src2_acc,
                              mask_upperleft, mask_acc,
                              dest_upperleft, dest_acc,
                              func);
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineThreeImagesMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                     SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                     SrcImageIterator3 src3_upperleft, SrcAccessor3 src3_acc,
                     DestImageIterator dest_upperleft, DestAccessor dest_acc,
                     const Functor& func)
{
    vigra::combineThreeImages(src1_upperleft, src1_lowerright, src1_acc,
                              src2_upperleft, src2_acc,
                              src3_upperleft, src3_acc,
                              dest_upperleft, dest_acc,
                              func);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                 DestImageIterator dest_upperleft, DestAccessor dest_acc,
                 const Functor& func)
{
    vigra::transformImage(src_upperleft, src_lowerright, src_acc,
                          dest_upperleft, dest_acc,
                          func);
}


template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageIfMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                   MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                   DestImageIterator dest_upperleft, DestAccessor dest_acc,
                   const Functor& func)
{
    vigra::transformImageIf(src_upperleft, src_lowerright, src_acc,
                            mask_upperleft, mask_acc,
                            dest_upperleft, dest_acc,
                            func);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class ValueType>
void
distanceTransformMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor sa,
                    DestImageIterator dest_upperleft, DestAccessor da,
                    ValueType background, int norm)
{
    vigra::distanceTransform(src_upperleft, src_lowerright, src_acc,
                             dest_upperleft, dest_acc,
                             background, norm);
}

#endif // _OPENMP >= 200505


//
// Argument Object Factory versions
//

template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesMP(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                   vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                   vigra::pair<DestImageIterator, DestAccessor> dest,
                   const Functor& func)
{
    combineTwoImagesMP(src1.first, src1.second, src1.third,
                       src2.first, src2.second,
                       dest.first, dest.second,
                       func);
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineTwoImagesIfMP(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                     vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                     vigra::pair<MaskImageIterator, MaskAccessor> mask,
                     vigra::pair<DestImageIterator, DestAccessor> dest,
                     const Functor& func)
{
    combineTwoImagesIfMP(src1.first, src1.second, src1.third,
                         src2.first, src2.second,
                         mask.first, mask.second,
                         dest.first, dest.second,
                         func);
}


template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
combineThreeImagesMP(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                     vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                     vigra::pair<SrcImageIterator3, SrcAccessor3> src3,
                     vigra::pair<DestImageIterator, DestAccessor> dest,
                     const Functor& func)
{
    combineThreeImagesMP(src1.first, src1.second, src1.third,
                         src2.first, src2.second,
                         src3.first, src3.second,
                         dest.first, dest.second,
                         func);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageMP(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                 vigra::pair<DestImageIterator, DestAccessor> dest,
                 const Functor& func)
{
    transformImageMP(src.first, src.second, src.third,
                     dest.first, dest.second,
                     func);
}


template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageIfMP(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                   vigra::pair<MaskImageIterator, MaskAccessor> mask,
                   vigra::pair<DestImageIterator, DestAccessor> dest,
                   const Functor& func)
{
    transformImageIfMP(src.first, src.second, src.third,
                       mask.first, mask.second,
                       dest.first, dest.second,
                       func);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class ValueType>
inline void
distanceTransformMP(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    vigra::pair<DestImageIterator, DestAccessor> dest,
                    ValueType background, int norm)
{
    distanceTransformMP(src.first, src.second, src.third,
                        dest.first, dest.second,
                        background, norm);
}


#endif // __OPENMP_H__

// Local Variables:
// mode: c++
// End:
