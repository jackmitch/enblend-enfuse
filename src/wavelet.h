/*
 * Copyright (C) 2005 Andrew Mihal
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
#ifndef __WAVELET_H__
#define __WAVELET_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vigra/basicimage.hxx"

using std::max;

using vigra::BasicImage;

namespace enblend {

// Forward CDF(4,2) wavelet transform with integer lifting
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void wavelet(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;
    typedef BasicImage<DestPixelType> TempImage;
    typedef typename TempImage::traverser TempImageIterator;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in wavelet");

    TempImage tempS((max(src_w, src_h) + 1) >> 1, 1);
    TempImage tempD((max(src_w, src_h) + 1) >> 1, 1);

    // Transform rows from source to dest
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dy = dest_upperleft;
    for (; sy.y != send.y; ++sy.y, ++dy.y) {

        SrcImageIterator sx = sy;

        // Splitting from source into tempS and tempD.
        TempImageIterator tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; sx.x < (send.x - 1); ++sx.x, ++tsi.x, ++tdi.x) {
            *tsi = sa(sx);
            ++sx.x;
            *tdi = sa(sx);
        }
        // Odd sample
        if (sx.x != send.x) {
            *tsi = sa(sx);
        }

        // From tempS and tempD to tempS
        // s = s - (1/4)(dp + d)
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        if (wraparound) {
            *tsi = *tsi - ((*tdi + tdi((src_w>>1)-1, 0)) >> 2);
        } else {
            *tsi = *tsi - (*tdi >> 1);
        }
        ++tsi.x;
        ++tdi.x;
        for (int i = 1; i < (src_w >> 1); ++i, ++tsi.x, ++tdi.x) {
            *tsi = *tsi - ((tdi(-1, 0) + *tdi) >> 2);
        }
        if ((src_w & 1) == 1) {
            if (wraparound) {
                *tsi = *tsi - ((tdi(-1, 0) + *(tempD.upperLeft())) >> 2);
            } else {
                *tsi = *tsi - (tdi(-1, 0) >> 1);
            }
        }
        
        // Dual lifting, from tempS and tempD to second half of destImage
        // d = d - (s + sn)
        DestImageIterator dx = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator dend = dest_lowerright;
        if ((src_w & 1) == 0) dend += Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++tdi.x) {
            //*dx = *tdi - ((*tsi + tsi(1, 0)) >> 1);
            *dx = *tdi - (*tsi + tsi(1, 0));
        }
        if ((src_w & 1) == 0){
            if (wraparound) {
                //*dx = *tdi - ((*tsi + *(tempS.upperLeft())) >> 1);
                *dx = *tdi - (*tsi + *(tempS.upperLeft()));
            } else {
                // Mirror edge treatment
                //*dx = *tdi - *tsi;
                *dx = *tdi - (*tsi << 1);
            }
        }

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        // s = s + (3/16)(d + dp)
        dx = dy;
        dend = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator dualLiftResult = dend;
        tsi = tempS.upperLeft();
        if (wraparound) {
            //*dx = *tsi + ((*dualLiftResult + dy(src_w - 1, 0)) >> 2);
            DestPixelType tmp = *dualLiftResult + dy(src_w - 1, 0);
            tmp += (tmp << 1);
            *dx = *tsi + (tmp >> 4);
        } else {
            // Mirror edge treatment
            //*dx = *tsi + (*dualLiftResult >> 1);
            DestPixelType tmp = *dualLiftResult;
            tmp += (tmp << 1);
            *dx = *tsi + (tmp >> 3);
        }
        ++dx.x;
        ++tsi.x;
        ++dualLiftResult.x;
        if ((src_w & 1) == 1) dend += Diff2D(-1, 0);
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++dualLiftResult.x) {
            //*dx = *tsi + ((dualLiftResult(-1, 0) + *dualLiftResult) >> 2);
            DestPixelType tmp = dualLiftResult(-1, 0) + *dualLiftResult;
            tmp += (tmp << 1);
            *dx = *tsi + (tmp >> 4);
        }
        if ((src_w & 1) == 1) {
            if (wraparound) {
                //*dx = *tsi + ((dualLiftResult(-1, 0) + dy((src_w + 1) >> 1, 0)) >> 2);
                DestPixelType tmp = dualLiftResult(-1, 0) + dy((src_w + 1) >> 1, 0);
                tmp += (tmp << 1);
                *dx = *tsi + (tmp >> 4);
            } else {
                //*dx = *tsi + (dualLiftResult(-1, 0) >> 1);
                DestPixelType tmp = dualLiftResult(-1, 0);
                tmp += (tmp << 1);
                *dx = *tsi + (tmp >> 3);
            }
        }

    }

    // Transform columns from dest to dest
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    //// If lowpass only, do half the columns.
    //if (lowpassOnly) dend += Diff2D((src_w + 1) >> 1 - dend.x, 0);
    for (; dx.x != dend.x; ++dx.x) {

        dy = dx;

        // Split from dest into tempS and tempD.
        TempImageIterator tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; dy.y < (dend.y - 1); ++dy.y, ++tsi.x, ++tdi.x) {
            *tsi = da(dy);
            ++dy.y;
            *tdi = da(dy);
        }
        // Odd sample
        if (dy.y != dend.y) {
            *tsi = da(dy);
        }

        // From tempS and tempD to tempS
        // s = s - (1/4)(dp + d)
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        *tsi = *tsi - ((*tdi + tdi((src_h>>1)-1, 0)) >> 2);
        ++tsi.x;
        ++tdi.x;
        for (int i = 1; i < (src_h >> 1); ++i, ++tsi.x, ++tdi.x) {
            *tsi = *tsi - ((tdi(-1, 0) + *tdi) >> 2);
        }
        if ((src_h & 1) == 1) {
            *tsi = *tsi - (tdi(-1, 0) >> 1);
        }
        
        // Dual lifting, from tempS and tempD to second half of destImage
        dy = dx + Diff2D(0, (src_h + 1) >> 1);
        DestImageIterator dend = dest_lowerright;
        if ((src_h & 1) == 0) dend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++tdi.x) {
            //*dy = *tdi - ((*tsi + tsi(1, 0)) >> 1);
            *dy = *tdi - (*tsi + tsi(1, 0));
        }
        if ((src_h & 1) == 0) {
            // Mirror edge treatment
            //*dy = *tdi - *tsi;
            *dy = *tdi - (*tsi << 1);
        }

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        dy = dx;
        dend = dx + Diff2D(0, (src_h + 1) >> 1);
        DestImageIterator dualLiftResult = dend;
        if ((src_h & 1) == 1) dend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        // Mirror edge treatment
        //*dy = *tsi + (*dualLiftResult >> 1);
        {
            DestPixelType tmp = *dualLiftResult;
            tmp += (tmp << 1);
            *dy = *tsi + (tmp >> 3);
        }
        ++dy.y;
        ++tsi.x;
        ++dualLiftResult.y;
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++dualLiftResult.y) {
            //*dy = *tsi + ((dualLiftResult(0, -1) + *dualLiftResult) >> 2);
            DestPixelType tmp = dualLiftResult(0, -1) + *dualLiftResult;
            tmp += (tmp << 1);
            *dy = *tsi + (tmp >> 4);
        }
        if ((src_h & 1) == 1) {
            //*dy = *tsi + (dualLiftResult(0, -1) >> 1);
            DestPixelType tmp = dualLiftResult(0, -1);
            tmp += (tmp << 1);
            *dy = *tsi + (tmp >> 3);
        }

    }

};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void wavelet(unsigned int levels, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    if (levels == 0) return;

    wavelet(wraparound,
            src_upperleft, src_lowerright, sa,
            dest_upperleft, dest_lowerright, da);

    for (unsigned int i = 1; i < levels; i++) {
        int dest_w = dest_lowerright.x - dest_upperleft.x;
        int dest_h = dest_lowerright.y - dest_upperleft.y;
        dest_lowerright = dest_upperleft + Diff2D((dest_w + 1) >> 1, (dest_h + 1) >> 1);

        wavelet(wraparound,
                dest_upperleft, dest_lowerright, da,
                dest_upperleft, dest_lowerright, da);
    }
};

// Version using argument object factories
template <typename SrcImageType, typename DestImageType>
inline void wavelet(unsigned int levels, bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        triple<typename DestImageType::traverser, typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    wavelet(levels, wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

// Inverse CDF(4,2) wavelet transform with integer lifting
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void iwavelet(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;
    typedef BasicImage<DestPixelType> TempImage;
    typedef typename TempImage::traverser TempImageIterator;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in wavelet");

    TempImage tempS((max(src_w, src_h) + 1) >> 1, 1);
    TempImage tempD((max(src_w, src_h) + 1) >> 1, 1);

    // Transform columns from source to dest.
    SrcImageIterator sx = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; sx.x != send.x; ++sx.x, ++dx.x) {

        // Inverse primal lifting from source to tempS
        SrcImageIterator si = sx;
        SrcImageIterator siend = sx + Diff2D(0, (src_h + 1) >> 1);
        SrcImageIterator di = siend;
        if ((src_h & 1) == 1) siend += Diff2D(0, -1);
        TempImageIterator tsi = tempS.upperLeft();
        // Mirror edge treatment
        //*tsi = *si - (*di >> 1);
        {
            DestPixelType tmp = *di;
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 3);
        }
        ++si.y;
        ++di.y;
        ++tsi.x;
        for (; si.y != siend.y; ++si.y, ++di.y, ++tsi.x) {
            //*tsi = *si - ((di(0, -1) + *di) >> 2);
            DestPixelType tmp = di(0, -1) + *di;
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 4);
        }
        if ((src_h & 1) == 1) {
            //*tsi = *si - (di(0, -1) >> 1);
            DestPixelType tmp = di(0, -1);
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 3);
        }

        // Inverse dual lifting from source and tempS to tempD
        di = sx + Diff2D(0, (src_h + 1) >> 1);
        SrcImageIterator diend = src_lowerright;
        if ((src_h & 1) == 0) diend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.y != diend.y; ++di.y, ++tsi.x, ++tdi.x) {
            //*tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
            *tdi = *di + (*tsi + tsi(1, 0));
        }
        if ((src_h & 1) == 0) {
            // Mirror edge treatment
            //*tdi = *di + *tsi;
            *tdi = *di + (*tsi << 1);
        }

        // From tempS and tempD to tempS
        // s = s + (1/4)(dp + d)
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        *tsi = *tsi + ((*tdi + tdi((src_h>>1)-1, 0)) >> 2);
        ++tsi.x;
        ++tdi.x;
        for (int i = 1; i < (src_h >> 1); ++i, ++tsi.x, ++tdi.x) {
            *tsi = *tsi + ((tdi(-1, 0) + *tdi) >> 2);
        }
        if ((src_h & 1) == 1) {
            *tsi = *tsi + (tdi(-1, 0) >> 1);
        }

        // Merging from tempS and tempD to destImage
        DestImageIterator dy = dx;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y < (dend.y - 1); ++dy.y, ++tsi.x, ++tdi.x) {
            *dy = *tsi;
            ++dy.y;
            *dy = *tdi;
        }
        // Odd sample
        if (dy.y != dend.y) {
            *dy = *tsi;
        }

    }

    // Transform rows from dest to dest.
    DestImageIterator dy = dest_upperleft;
    for (; dy.y != dend.y; ++dy.y) {

        // Inverse primal lifting from dest to tempS
        DestImageIterator si = dy;
        DestImageIterator siend = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator di = siend;
        TempImageIterator tsi = tempS.upperLeft();
        if (wraparound) {
            //*tsi = *si - ((*di + dy(src_w - 1, 0)) >> 2);
            DestPixelType tmp = *di + dy(src_w - 1, 0);
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 4);
        } else {
            // Mirror edge treatment
            //*tsi = *si - (*di >> 1);
            DestPixelType tmp = *di;
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 3);
        }
        ++si.x;
        ++di.x;
        ++tsi.x;
        if ((src_w & 1) == 1) siend += Diff2D(-1, 0);
        for (; si.x != siend.x; ++si.x, ++di.x, ++tsi.x) {
            //*tsi = *si - ((di(-1, 0) + *di) >> 2);
            DestPixelType tmp = di(-1, 0) + *di;
            tmp += (tmp << 1);
            *tsi = *si - (tmp >> 4);
        }
        if ((src_w & 1) == 1) {
            if (wraparound) {
                //*tsi = *si - ((di(-1, 0) + dy((src_w + 1) >> 1, 0)) >> 2);
                DestPixelType tmp = di(-1, 0) + dy((src_w + 1) >> 1, 0);
                tmp += (tmp << 1);
                *tsi = *si - (tmp >> 4);
            } else {
                //*tsi = *si - (di(-1, 0) >> 1);
                DestPixelType tmp = di(-1, 0);
                tmp += (tmp << 1);
                *tsi = *si - (tmp >> 3);
            }
        }

        // Inverse dual lifting from dest and tempS to tempD
        di = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator diend = dend;
        if ((src_w & 1) == 0) diend += Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.x != diend.x; ++di.x, ++tsi.x, ++tdi.x) {
            //*tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
            *tdi = *di + (*tsi + tsi(1, 0));
        }
        if ((src_w & 1) == 0) {
            if (wraparound) {
                //*tdi = *di + ((*tsi + *(tempS.upperLeft())) >> 1);
                *tdi = *di + (*tsi + *(tempS.upperLeft()));
            } else {
                // Mirror edge treatment
                //*tdi = *di + *tsi;
                *tdi = *di + (*tsi << 1);
            }
        }

        // From tempS and tempD to tempS
        // s = s + (1/4)(dp + d)
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        if (wraparound) {
            *tsi = *tsi + ((*tdi + tdi((src_w>>1)-1, 0)) >> 2);
        } else {
            *tsi = *tsi + (*tdi >> 1);
        }
        ++tsi.x;
        ++tdi.x;
        for (int i = 1; i < (src_w >> 1); ++i, ++tsi.x, ++tdi.x) {
            *tsi = *tsi + ((tdi(-1, 0) + *tdi) >> 2);
        }
        if ((src_w & 1) == 1) {
            if (wraparound) {
                *tsi = *tsi + ((tdi(-1, 0) + *(tempD.upperLeft())) >> 2);
            } else {
                *tsi = *tsi + (tdi(-1, 0) >> 1);
            }
        }

        // Merging from tempS and tempD to destImage
        dx = dy;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x < (dend.x - 1); ++dx.x, ++tsi.x, ++tdi.x) {
            *dx = *tsi;
            ++dx.x;
            *dx = *tdi;
        }
        // Odd sample
        if (dx.x != dend.x) {
            *dx = *tsi;
        }
    }

};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void iwavelet(unsigned int levels, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    copyImage(src_upperleft, src_lowerright, sa, dest_upperleft, da);

    if (levels == 0) return;

    DestImageIterator *lowerRightArray = new DestImageIterator[levels];
    lowerRightArray[0] = dest_lowerright;
    for (unsigned int i = 1; i < levels; i++) {
        int dest_w = dest_lowerright.x - dest_upperleft.x;
        int dest_h = dest_lowerright.y - dest_upperleft.y;
        dest_lowerright = dest_upperleft + Diff2D((dest_w + 1) >> 1, (dest_h + 1) >> 1);
        lowerRightArray[i] = dest_lowerright;
    }

    for (unsigned int i = levels; i > 0; --i) {
        iwavelet(wraparound,
                dest_upperleft, lowerRightArray[i-1], da,
                dest_upperleft, lowerRightArray[i-1], da);
    }

    delete [] lowerRightArray;

}

// Version using argument object factories
template <typename SrcImageType, typename DestImageType>
inline void iwavelet(unsigned int levels, bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        triple<typename DestImageType::traverser, typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    iwavelet(levels, wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void wavelet(int levels, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        MaskIterator mask_upperleft,
        MaskAccessor ma,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da,
        DestMaskIterator dest_mask_upperleft,
        DestMaskIterator dest_mask_lowerright,
        DestMaskAccessor dma) {

};

} // namespace enblend

#endif /* __WAVELET_H__ */
