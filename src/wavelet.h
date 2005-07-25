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

template <typename T>
inline T three16(const T & v) {
    T res = v + (v << 1);
    return (res >> 4);
}

template <typename T>
inline T three8(const T & v) {
    T res = v + (v << 1);
    return (res >> 3);
}

// Forward cubic B-spline wavelet transform with integer lifting
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void _wavelet(unsigned int srcLevel, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;

    // Size of input image
    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;

    // Distance between s pixels
    int stride = 2 << srcLevel;

    // Distance between s-d pair
    int adjacent = 1 << srcLevel;

    // Number of s and d pixels
    int sCount = src_w;
    int dCount;
    for (unsigned int i = 0; i <= srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }

    // First do the rows
    {
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dy = dest_upperleft;
    for (; sy.y != send.y; ++sy.y, ++dy.y) {

        SrcImageIterator sx = sy;
        DestImageIterator dx = dy;

        // Lifting 1:
        // ds = ss - 1/4 * (sd_prev + sd)
        {
            SrcImageIterator ss = sx;
            SrcImageIterator sd = ss + Diff2D(adjacent, 0);

            DestImageIterator ds = dx;
            //DestImageIterator dd = ds + Diff2D(adjacent, 0);

            // In case sCount > dCount and wraparound:
            // last s will depend on first s
            DestPixelType firstPixel;

            // Left border
            if (wraparound) {
                //cout << "l1 left ww" << endl;
                firstPixel = *ss - ( (sx(adjacent + stride * (dCount - 1), 0) + *sd) >> 2 );
            } else {
                firstPixel = *ss - ( *sd >> 1 );
            }
            ds.x += stride;
            ss.x += stride;
            sd.x += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds = *ss - ( (sd(-stride, 0) + *sd) >> 2 );
                ds.x += stride;
                ss.x += stride;
                sd.x += stride;
            }

            // Right border
            if (i < sCount) {
                if (wraparound) {
                    //cout << "l1 right ww" << endl;
                    *ds = *ss - ( (sd(-stride, 0) + sx(adjacent, 0)) >> 2 );
                } else {
                    *ds = *ss - ( sd(-stride, 0) >> 1 );
                }
            }

            // First pixel
            *dx = firstPixel;
        }

        // Lifting 2:
        // dd = sd - (ds + ds_next)
        {
            SrcImageIterator ss = sx;
            SrcImageIterator sd = ss + Diff2D(adjacent, 0);

            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            int i;
            for (i = 0; i < (sCount - 1); i++) {
                *dd = *sd - (*ds + ds(stride, 0));
                dd.x += stride;
                sd.x += stride;
                ds.x += stride;
            }

            // Right border
            if (i < dCount) {
                if (wraparound) {
                    //cout << "l2 right ww" << endl;
                    *dd = *sd - (*ds + *dx);
                } else {
                    *dd = *sd - (*ds << 1);
                }
            }
        }

        // Lifting 3:
        // ds = ds + 3/16 * (dd_prev + dd)
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            // In case sCount > dCount and wraparound:
            // last s will depend on first s
            DestPixelType firstPixel;

            // Left border
            if (wraparound) {
                //cout << "l3 left ww" << endl;
                firstPixel = *ds + three16( dx(adjacent + stride * (dCount - 1), 0) + *dd );
            } else {
                firstPixel = *ds + three8( *dd );
            }
            ds.x += stride;
            dd.x += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds += three16( dd(-stride, 0) + *dd );
                ds.x += stride;
                dd.x += stride;
            }

            // Right border
            if (i < sCount) {
                if (wraparound) {
                    //cout << "l3 right ww" << endl;
                    *ds += three16( dd(-stride, 0) + dx(adjacent, 0) );
                } else {
                    *ds += three8( dd(-stride, 0) );
                }
            }

            // First pixel
            *dx = firstPixel;
        }

        // Scaling
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            int i;
            for (i = 0; i < dCount; i++) {
                *ds <<= 1;
                *ds -= *dd & DestPixelType(1);
                *dd += *dd & DestPixelType(1);
                *dd >>= 1;
                ds.x += stride;
                dd.x += stride;
            }

            if (i < sCount) {
                //cout << "scale last s" << endl;
                *ds <<= 1;
            }
        }

    }
    } // end of all rows namespace

    sCount = src_h;
    for (unsigned int i = 0; i <= srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }

    // Now do the columns
    {
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dx.x != dend.x; ++dx.x) {

        DestImageIterator dy = dx;

        // Lifting 1:
        // ds = ds - 1/4 * (dd_prev + dd)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            // Top border - no wraparound
            *ds -= *dd >> 1;
            ds.y += stride;
            dd.y += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds -= (dd(0, -stride) + *dd) >> 2;
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < sCount) {
                *ds -= dd(0, -stride) >> 1;
            }
        }

        // Lifting 2:
        // dd = dd - (ds + ds_next)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            int i;
            for (i = 0; i < (sCount - 1); i++) {
                *dd -= *ds + ds(0, stride);
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < dCount) {
                *dd -= *ds << 1;
            }
        }

        // Lifting 3:
        // ds = ds + 3/16 * (dd_prev + dd)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            // Top border
            *ds += three8(*dd);
            ds.y += stride;
            dd.y += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds += three16( dd(0, -stride) + *dd );
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < sCount) {
                *ds += three8( dd(0, -stride) );
            }
        }

        // Scaling
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            int i;
            for (i = 0; i < dCount; i++) {
                *ds <<= 1;
                *ds -= *dd & DestPixelType(1);
                *dd += *dd & DestPixelType(1);
                *dd >>= 1;
                ds.y += stride;
                dd.y += stride;
            }

            if (i < sCount) {
                *ds <<= 1;
            }
        }

    }
    } // end of all columns namespace

}

// Inverse cubic B-spline wavelet transform with integer lifting
template <typename DestImageIterator, typename DestAccessor>
inline void _iwavelet(unsigned int srcLevel, bool wraparound,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;

    if (srcLevel == 0) {
        vigra_fail("Cannot run IDWT on level 0 image.");
    }

    // Size of input image
    int dest_w = dest_lowerright.x - dest_upperleft.x;
    int dest_h = dest_lowerright.y - dest_upperleft.y;

    // Distance between s pixels
    int stride = 1 << srcLevel;

    // Distance between s-d pair
    int adjacent = 1 << (srcLevel - 1);

    cout << "stride=" << stride << " adjacent=" << adjacent << endl;

    // Number of s and d pixels
    int sCount = dest_h;
    int dCount = 0;
    for (unsigned int i = 0; i < srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }
    cout << "dest_h=" << dest_h << " sCount=" << sCount << " dCount=" << dCount << endl;

    // First do the columns
    {
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dx.x != dend.x; ++dx.x) {

        DestImageIterator dy = dx;

        // Inverse scaling
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            int i;
            for (i = 0; i < dCount; i++) {
                *dd <<= 1;
                *dd -= *ds & DestPixelType(1);
                *ds += *ds & DestPixelType(1);
                *ds >>= 1;
                ds.y += stride;
                dd.y += stride;
            }

            if (i < sCount) {
                *ds = (*ds + DestPixelType(1)) >> 1;
            }
        }

        // Inverse lifting 3:
        // ds = ds - 3/16 * (dd_prev + dd)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            // Top border
            *ds -= three8(*dd);
            ds.y += stride;
            dd.y += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds -= three16( dd(0, -stride) + *dd );
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < sCount) {
                *ds -= three8( dd(0, -stride) );
            }
        }

        // Inverse lifting 2:
        // dd = dd + (ds + ds_next)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            int i;
            for (i = 0; i < (sCount - 1); i++) {
                *dd += *ds + ds(0, stride);
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < dCount) {
                *dd += *ds << 1;
            }
        }

        // Inverse lifting 1:
        // ds = ds + 1/4 * (dd_prev + dd)
        {
            DestImageIterator ds = dy;
            DestImageIterator dd = ds + Diff2D(0, adjacent);

            // Top border
            *ds += *dd >> 1;
            ds.y += stride;
            dd.y += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds += (dd(0, -stride) + *dd) >> 2;
                ds.y += stride;
                dd.y += stride;
            }

            // Bottom border
            if (i < sCount) {
                *ds += dd(0, -stride) >> 1;
            }
        }

    }
    } // end of all columns namespace

    sCount = dest_w;
    for (unsigned int i = 0; i < srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }
    cout << "dest_w=" << dest_w << " sCount=" << sCount << " dCount=" << dCount << endl;

    // Next do the rows
    {
    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dy.y != dend.y; ++dy.y) {

        DestImageIterator dx = dy;

        // Inverse scaling
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            int i;
            for (i = 0; i < dCount; i++) {
                *dd <<= 1;
                *dd -= *ds & DestPixelType(1);
                *ds += *ds & DestPixelType(1);
                *ds >>= 1;
                ds.x += stride;
                dd.x += stride;
            }

            if (i < sCount) {
                //cout << "iscale last s" << endl;
                *ds = (*ds + DestPixelType(1)) >> 1;
            }
        }

        // Inverse lifting 3:
        // ds = ds - 3/16 * (dd_prev + dd)
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            // In case sCount > dCount and wraparound:
            // last s will depend on first s
            DestPixelType firstPixel;

            // Left border
            if (wraparound) {
                //cout << "il3 left ww" << endl;
                firstPixel = *ds - three16( dx(adjacent + stride * (dCount - 1), 0) + *dd );
            } else {
                firstPixel = *ds - three8( *dd );
            }
            ds.x += stride;
            dd.x += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds -= three16( dd(-stride, 0) + *dd );
                ds.x += stride;
                dd.x += stride;
            }

            // Right border
            if (i < sCount) {
                if (wraparound) {
                    //cout << "il3 right ww" << endl;
                    *ds -= three16( dd(-stride, 0) + dx(adjacent, 0) );
                } else {
                    *ds -= three8( dd(-stride, 0) );
                }
            }

            // First pixel
            *dx = firstPixel;
        }

        // Inverse lifting 2:
        // dd = dd + (ds + ds_next)
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            int i;
            for (i = 0; i < (sCount - 1); i++) {
                *dd += *ds + ds(stride, 0);
                ds.x += stride;
                dd.x += stride;
            }

            // Right border
            if (i < dCount) {
                if (wraparound) {
                    //cout << "il2 right ww" << endl;
                    *dd += *ds + *dx;
                } else {
                    *dd += *ds << 1;
                }
            }
        }

        // Inverse lifting 1:
        // ds = ds + 1/4 * (dd_prev + dd)
        {
            DestImageIterator ds = dx;
            DestImageIterator dd = ds + Diff2D(adjacent, 0);

            // In case sCount > dCount and wraparound:
            // last s will depend on first s
            DestPixelType firstPixel;

            // Left border
            if (wraparound) {
                //cout << "il1 left ww" << endl;
                firstPixel = *ds + ( (dx(adjacent + stride * (dCount - 1), 0) + *dd) >> 2 );
            } else {
                firstPixel = *ds + ( *dd >> 1 );
            }
            ds.x += stride;
            dd.x += stride;

            int i;
            for (i = 1; i < dCount; i++) {
                *ds += (dd(-stride, 0) + *dd) >> 2;
                ds.x += stride;
                dd.x += stride;
            }

            // Right border
            if (i < sCount) {
                if (wraparound) {
                    //cout << "il1 right ww" << endl;
                    *ds += ( (dd(-stride, 0) + dx(adjacent, 0)) >> 2 );
                } else {
                    *ds += ( dd(-stride, 0) >> 1 );
                }
            }

            // First pixel
            *dx = firstPixel;
        }
    }
    } // end of all rows namespace

};

// Clear all detail coefficients at a certain level.
template <typename DestImageIterator, typename DestAccessor>
inline void _zeroDetailCoefficients(unsigned int srcLevel,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;

    if (srcLevel == 0) {
        vigra_fail("Cannot run zeroDetailCoefficients on level 0 image.");
    }

    // Size of input image
    int dest_w = dest_lowerright.x - dest_upperleft.x;
    int dest_h = dest_lowerright.y - dest_upperleft.y;

    // Distance between s pixels
    int stride = 1 << srcLevel;

    // Distance between s-d pair
    int adjacent = 1 << (srcLevel - 1);

    cout << "stride=" << stride << " adjacent=" << adjacent << endl;

    // Number of s and d pixels
    int sCount = dest_h;
    int dCount = 0;
    for (unsigned int i = 0; i < srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }
    cout << "dest_h=" << dest_h << " sCount=" << sCount << " dCount=" << dCount << endl;

    // First do the columns
    {
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dx.x != dend.x; ++dx.x) {

        DestImageIterator dy = dx;

        DestImageIterator ds = dy;
        DestImageIterator dd = ds + Diff2D(0, adjacent);

        for (int i = 0; i < dCount; i++) {
            *dd = DestPixelType(0);
            dd.y += stride;
        }
    }
    }

    sCount = dest_w;
    for (unsigned int i = 0; i < srcLevel; i++) {
        dCount = sCount >> 1;
        sCount = (sCount + 1) >> 1;
    }
    cout << "dest_w=" << dest_w << " sCount=" << sCount << " dCount=" << dCount << endl;

    // Now do the rows
    {
    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dy.y != dend.y; ++dy.y) {

        DestImageIterator dx = dy;

        DestImageIterator ds = dx;
        DestImageIterator dd = ds + Diff2D(adjacent, 0);

        for (int i = 0; i < dCount; i++) {
            *dd = DestPixelType(0);
            dd.x += stride;
        }
    }
    }

}

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

    _wavelet(0, wraparound,
            src_upperleft, src_lowerright, sa,
            dest_upperleft, dest_lowerright, da);

    for (unsigned int i = 1; i < levels; i++) {
        _wavelet(i, wraparound,
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

template <typename DestImageIterator, typename DestAccessor>
inline void iwavelet(unsigned int levels, bool wraparound,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    if (levels == 0) return;

    for (unsigned int i = levels; i > 0; i--) {
        _iwavelet(i, wraparound, dest_upperleft, dest_lowerright, da);
    }

}

// Version using argument object factories
template <typename DestImageType>
inline void iwavelet(unsigned int levels, bool wraparound,
        triple<typename DestImageType::traverser, typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    iwavelet(levels, wraparound, dest.first, dest.second, dest.third);
};

template <typename DestImageIterator, typename DestAccessor>
inline void zeroDetailCoefficients(unsigned int levels,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    if (levels == 0) return;

    for (unsigned int i = 1; i <= levels; i++) {
        _zeroDetailCoefficients(i, dest_upperleft, dest_lowerright, da);
    }

};

// Version using argument object factories
template<typename DestImageType>
inline void zeroDetailCoefficients(unsigned int levels,
        triple<typename DestImageType::traverser, typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    zeroDetailCoefficients(levels, dest.first, dest.second, dest.third);
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
