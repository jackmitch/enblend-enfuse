/*
 * Copyright (C) 2004-2005 Andrew Mihal
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
#ifndef __PYRAMID_H__
#define __PYRAMID_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>
#include <vector>

#include "vigra/convolution.hxx"
#include "vigra/error.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/transformimage.hxx"

#include "fixmath.h"

using std::cout;
using std::vector;

using vigra::linearRangeMapping;
using vigra::NumericTraits;
using vigra::transformImage;
using vigra::triple;
using vigra::Int8;
using vigra::Int16;
using vigra::Int32;
using vigra::Int64;
using vigra::RGBValue;
using vigra::UInt8;
using vigra::UInt16;
using vigra::UInt32;
using vigra::UInt64;
using vigra::UInt16Image;
using vigra::UInt16RGBImage;

namespace enblend {

struct Error_PyramidPromoteTraits_not_specialized_for_this_case { };

template<class A>
struct PyramidPromoteTraits {
    typedef Error_PyramidPromoteTraits_not_specialized_for_this_case Type;
    typedef Error_PyramidPromoteTraits_not_specialized_for_this_case Promote;
};

#define DEFINE_PYRAMIDPROMOTETRAITS(A, B) \
template<> \
struct PyramidPromoteTraits<A> { \
    typedef A Type; \
    typedef B Promote; \
};

// SKIPSM 5x5 math requires 6 more bits on top of base type
DEFINE_PYRAMIDPROMOTETRAITS(Int8, Int16);
DEFINE_PYRAMIDPROMOTETRAITS(Int16, Int32);
DEFINE_PYRAMIDPROMOTETRAITS(Int32, Int64);
DEFINE_PYRAMIDPROMOTETRAITS(Int64, double);
DEFINE_PYRAMIDPROMOTETRAITS(UInt8, UInt16);
DEFINE_PYRAMIDPROMOTETRAITS(UInt16, UInt32);
DEFINE_PYRAMIDPROMOTETRAITS(UInt32, UInt64);
DEFINE_PYRAMIDPROMOTETRAITS(UInt64, double);
DEFINE_PYRAMIDPROMOTETRAITS(float, float);
DEFINE_PYRAMIDPROMOTETRAITS(double, double);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<Int8>, RGBValue<Int16>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<Int16>, RGBValue<Int32>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<Int32>, RGBValue<Int64>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<Int64>, RGBValue<double>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<UInt8>, RGBValue<UInt16>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<UInt16>, RGBValue<UInt32>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<UInt32>, RGBValue<UInt64>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<UInt64>, RGBValue<double>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<float>, RGBValue<float>);
DEFINE_PYRAMIDPROMOTETRAITS(RGBValue<double>, RGBValue<double>);

// time with multiplies: 77.23
// time with shifts: 80.78
#define IMUL6(A) (A * SKIPSMImagePixelType(6))
//#define IMUL6(A) ((A + (A << 1)) << 1)
#define MMUL6(A) (A * SKIPSMMaskPixelType(6))
//#define MMUL6(A) ((A + (A << 1)) << 1)

// Pyramid filter coefficients.
static const double A = 0.4;
static const double W[] = {0.25 - A / 2.0, 0.25, A, 0.25, 0.25 - A / 2.0};
static const unsigned int A100 = 40;
static const unsigned int W100[] = {25 - A100 / 2, 25, A100, 25, 25 - A100 / 2};

/** Calculate the half-width of a n-level filter.
 *  Assumes that the input function is a left-handed function,
 *  and the last non-zero input is at location 0.
 *  Returns the location of the last non-zero output.
 */
template <typename PixelType>
unsigned int filterHalfWidth(const unsigned int levels) {
    // For levels >= 30, the full width will just barely fit in int32.
    // When this range is added to a bounding box it will certainly
    // overflow the Diff2D.
    vigra_precondition((levels >= 1 && levels <= 29),
            "filterHalfWidth: levels outside of range [1,29]");

    // This is the arithmetic half width.
    int halfWidth = (levels == 1) ? 0 : ((1 << (levels+1)) - 4);

    return halfWidth;
}

/** The Burt & Adelson Reduce operation.
 *  This version is for images with alpha channels.
 *  Gaussian blur, downsampling, and extrapolation in one pass over the input images using SKIPSM-based algorithm.
 *
 *  Frederick M. Waltz and John W.V. Miller. An efficient algorithm for Gaussian blur using finite-state machines.
 *  SPIE Conf. on Machine Vision Systems for Inspection and Metrology VII. November 1998.
 */
template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void reduce(bool wraparound,
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

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename DestAccessor::value_type DestPixelType;
    typedef typename NumericTraits<PixelType>::RealPromote RealPixelType;
    typedef typename MaskAccessor::value_type MaskPixelType;
    typedef typename DestMaskAccessor::value_type DestMaskPixelType;
    typedef typename PyramidPromoteTraits<PixelType>::Promote SKIPSMImagePixelType;
    typedef UInt16 SKIPSMMaskPixelType;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    int dst_w = dest_lowerright.x - dest_upperleft.x;
    //int dst_h = dest_lowerright.y - dest_upperleft.y;

    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in reduce");

#if 1
    SKIPSMImagePixelType isr0, isr1, isrp;
    SKIPSMImagePixelType *isc0 = new SKIPSMImagePixelType[dst_w + 1];
    SKIPSMImagePixelType *isc1 = new SKIPSMImagePixelType[dst_w + 1];
    SKIPSMImagePixelType *iscp = new SKIPSMImagePixelType[dst_w + 1];

    SKIPSMMaskPixelType msr0, msr1, msrp;
    SKIPSMMaskPixelType *msc0 = new SKIPSMMaskPixelType[dst_w + 1];
    SKIPSMMaskPixelType *msc1 = new SKIPSMMaskPixelType[dst_w + 1];
    SKIPSMMaskPixelType *mscp = new SKIPSMMaskPixelType[dst_w + 1];

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dx = dy;
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator sx = sy;
    MaskIterator my = mask_upperleft;
    MaskIterator mx = my;
    DestMaskIterator dmy = dest_mask_upperleft;
    DestMaskIterator dmx = dmy;

    bool evenY = true;
    bool evenX = true;
    int srcy = 0;
    int srcx = 0;
    int dsty = 0;
    int dstx = 0;

// without extrapolation = 2.63
// with extrapolation = 4.39
// with wraparound (no -w) = 4.11

    // First row
    {
        if (wraparound) {
            msr0 = SKIPSMMaskPixelType(ma(my, Diff2D(src_w-2, 0)) ?
                    NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
            msr1 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
            msrp = SKIPSMMaskPixelType(ma(my, Diff2D(src_w-1, 0)) ? 
                    NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero()) << 2;
            isr0 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-2, 0)));
            isr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            isrp = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1, 0))) << 2;
        } else {
            msr0 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
            msr1 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
            msrp = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
            isr0 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            isr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            isrp = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        }

        for (sx = sy, mx = my, evenX = true, srcx = 0, dstx = 0;  srcx < src_w; ++srcx, ++sx.x, ++mx.x) {
            SKIPSMImagePixelType icurrent(sa(sx));
            SKIPSMMaskPixelType mcurrent((ma(mx)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
            if (evenX) {
                msc1[dstx] = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msc0[dstx] = msr1 + MMUL6(msr0) + msrp + mcurrent;
                msr1 = msr0 + msrp;
                msr0 = mcurrent;
                isc1[dstx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isc0[dstx] = isr1 + IMUL6(isr0) + isrp + icurrent;
                isr1 = isr0 + isrp;
                isr0 = icurrent;
            }
            else {
                msrp = mcurrent << 2;
                isrp = icurrent << 2;
                ++dstx;
            }
            evenX = !evenX;
        }
        // Last entries in first row
        if (!evenX) {
            // previous srcx was even
            ++dstx;
            if (wraparound) {
                msc1[dstx] = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msc0[dstx] = msr1 + MMUL6(msr0)
                                  + (SKIPSMMaskPixelType(ma(my) ?
                                            NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero()) << 2)
                                  + SKIPSMMaskPixelType(ma(my, Diff2D(1,0)) ?
                                            NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                isc1[dstx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isc0[dstx] = isr1 + IMUL6(isr0) + (SKIPSMImagePixelType(sa(sy)) << 2) + SKIPSMImagePixelType(sa(sy, Diff2D(1,0)));
            } else {
                msc1[dstx] = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msc0[dstx] = msr1 + MMUL6(msr0);
                isc1[dstx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isc0[dstx] = isr1 + IMUL6(isr0);
            }
        }
        else {
            // previous srcx was odd
            if (wraparound) {
                msc1[dstx] = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msc0[dstx] = msr1 + MMUL6(msr0) + msrp
                                  + SKIPSMMaskPixelType(ma(my) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                isc1[dstx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isc0[dstx] = isr1 + IMUL6(isr0) + isrp + SKIPSMImagePixelType(sa(sy));
            } else {
                msc1[dstx] = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msc0[dstx] = msr1 + MMUL6(msr0) + msrp;
                isc1[dstx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isc0[dstx] = isr1 + IMUL6(isr0) + isrp;
            }
        }
    }
    ++sy.y;
    ++my.y;

    // Main Rows
    {
        for (evenY = false, srcy = 1, dsty = 0; srcy < src_h; ++srcy, ++sy.y, ++my.y) {

            if (wraparound) {
                msr0 = SKIPSMMaskPixelType(ma(my, Diff2D(src_w-2, 0)) ?
                        NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                msr1 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msrp = SKIPSMMaskPixelType(ma(my, Diff2D(src_w-1, 0)) ? 
                        NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero()) << 2;
                isr0 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-2, 0)));
                isr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isrp = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1, 0))) << 2;
            } else {
                msr0 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msr1 = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                msrp = SKIPSMMaskPixelType(NumericTraits<SKIPSMMaskPixelType>::zero());
                isr0 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
                isrp = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            }

            if (evenY) {
                // Even-numbered row

                // First entry in row
                sx = sy;
                mx = my;
                if (wraparound) {
                    msr1 = msr0 + msrp;
                    isr1 = isr0 + isrp;
                }
                msr0 = SKIPSMMaskPixelType((ma(mx)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                isr0 = SKIPSMImagePixelType(sa(sx));
                //isc1[0] = isc0[0] + (iscp[0] << 2);
                //isc0[0] = sa(sx);
                ++sx.x;
                ++mx.x;
                dx = dy;
                dmx = dmy;

                // Main entries in row
                for (evenX = false, srcx = 1, dstx = 0; srcx < src_w; ++srcx, ++sx.x, ++mx.x) {
                    SKIPSMMaskPixelType mcurrent((ma(mx)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                    SKIPSMImagePixelType icurrent(sa(sx));
                    if (evenX) {
                        SKIPSMMaskPixelType mp = msc1[dstx] + MMUL6(msc0[dstx]) + mscp[dstx];
                        msc1[dstx] = msc0[dstx] + mscp[dstx];
                        msc0[dstx] = msr1 + MMUL6(msr0) + msrp + mcurrent;
                        msr1 = msr0 + msrp;
                        msr0 = mcurrent;
                        mp += msc0[dstx];

                        SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
                        isc1[dstx] = isc0[dstx] + iscp[dstx];
                        isc0[dstx] = isr1 + IMUL6(isr0) + isrp + icurrent;
                        isr1 = isr0 + isrp;
                        isr0 = icurrent;
                        if (mp) {
                            ip += isc0[dstx];
                            ip /= SKIPSMImagePixelType(mp);
                            da.set(ip, dx);
                        } else {
                            da.set(DestPixelType(NumericTraits<DestPixelType>::zero()), dx);
                        }
                        dma.set(mp ? NumericTraits<DestMaskPixelType>::max() : NumericTraits<DestMaskPixelType>::zero(), dmx);

                        ++dx.x;
                        ++dmx.x;
                    }
                    else {
                        msrp = mcurrent << 2;
                        isrp = icurrent << 2;
                        ++dstx;
                    }
                    evenX = !evenX;
                }

                // Last entries in row
                if (!evenX) {
                    // previous srcx was even
                    ++dstx;

                    SKIPSMMaskPixelType mp = msc1[dstx] + MMUL6(msc0[dstx]) + mscp[dstx];
                    msc1[dstx] = msc0[dstx] + mscp[dstx];
                    if (wraparound) {
                        msc0[dstx] = msr1 + MMUL6(msr0)
                                          + (SKIPSMMaskPixelType(ma(my) ?
                                                    NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero()) << 2)
                                          + SKIPSMMaskPixelType(ma(my, Diff2D(1,0)) ?
                                                    NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                    } else {
                        msc0[dstx] = msr1 + MMUL6(msr0);
                    }
                    mp += msc0[dstx];

                    SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
                    isc1[dstx] = isc0[dstx] + iscp[dstx];
                    if (wraparound) {
                        isc0[dstx] = isr1 + IMUL6(isr0) + (SKIPSMImagePixelType(sa(sy)) << 2) + SKIPSMImagePixelType(sa(sy, Diff2D(1,0)));
                    } else {
                        isc0[dstx] = isr1 + IMUL6(isr0);
                    }
                    if (mp) {
                        ip += isc0[dstx];
                        ip /= SKIPSMImagePixelType(mp);
                        da.set(ip, dx);
                    } else {
                        da.set(DestPixelType(NumericTraits<DestPixelType>::zero()), dx);
                    }
                    dma.set(mp ? NumericTraits<DestMaskPixelType>::max() : NumericTraits<DestMaskPixelType>::zero(), dmx);
                }
                else {
                    // Previous srcx was odd
                    SKIPSMMaskPixelType mp = msc1[dstx] + MMUL6(msc0[dstx]) + mscp[dstx];
                    msc1[dstx] = msc0[dstx] + mscp[dstx];
                    if (wraparound) {
                        msc0[dstx] = msr1 + MMUL6(msr0) + msrp
                                          + SKIPSMMaskPixelType(ma(my) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                    } else {
                        msc0[dstx] = msr1 + MMUL6(msr0) + msrp;
                    }
                    mp += msc0[dstx];

                    SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
                    isc1[dstx] = isc0[dstx] + iscp[dstx];
                    if (wraparound) {
                        isc0[dstx] = isr1 + IMUL6(isr0) + isrp + SKIPSMImagePixelType(sa(sy));
                    } else {
                        isc0[dstx] = isr1 + IMUL6(isr0) + isrp;
                    }
                    if (mp) {
                        ip += isc0[dstx];
                        ip /= SKIPSMImagePixelType(mp);
                        da.set(ip, dx);
                    } else {
                        da.set(DestPixelType(NumericTraits<DestPixelType>::zero()), dx);
                    }
                    dma.set(mp ? NumericTraits<DestMaskPixelType>::max() : NumericTraits<DestMaskPixelType>::zero(), dmx);
                }

                ++dsty;
                ++dy.y;
                ++dmy.y;

            }
            else {
                // First entry in odd-numbered row
                sx = sy;
                mx = my;
                if (wraparound) {
                    msr1 = msr0 + msrp;
                    isr1 = isr0 + isrp;
                }
                msr0 = SKIPSMMaskPixelType((ma(mx)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                isr0 = SKIPSMImagePixelType(sa(sx));
                //isc1[0] = isc0[0] + (iscp[0] << 2);
                //isc0[0] = sa(sx);
                ++sx.x;
                ++mx.x;

                // Main entries in odd-numbered row
                for (evenX = false, srcx = 1, dstx = 0; srcx < src_w; ++srcx, ++sx.x, ++mx.x) {
                    SKIPSMMaskPixelType mcurrent((ma(mx)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero());
                    SKIPSMImagePixelType icurrent(sa(sx));
                    if (evenX) {
                        mscp[dstx] = (msr1 + MMUL6(msr0) + msrp + mcurrent) << 2;
                        msr1 = msr0 + msrp;
                        msr0 = mcurrent;
                        iscp[dstx] = (isr1 + IMUL6(isr0) + isrp + icurrent) << 2;
                        isr1 = isr0 + isrp;
                        isr0 = icurrent;
                    }
                    else {
                        msrp = mcurrent << 2;
                        isrp = icurrent << 2;
                        ++dstx;
                    }
                    evenX = !evenX;
                }
                // Last entries in row
                if (!evenX) {
                    // previous srcx was even
                    ++dstx;
                    if (wraparound) {
                        mscp[dstx] = (msr1 + MMUL6(msr0)
                                           + (SKIPSMMaskPixelType(ma(my) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero()) << 2)
                                           + SKIPSMMaskPixelType(ma(my, Diff2D(1,0)) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero())
                                     ) << 2;
                        iscp[dstx] = (isr1 + IMUL6(isr0) + (SKIPSMImagePixelType(sa(sy)) << 2) + SKIPSMImagePixelType(sa(sy, Diff2D(1,0)))) << 2;
                    } else {
                        mscp[dstx] = (msr1 + MMUL6(msr0)) << 2;
                        iscp[dstx] = (isr1 + IMUL6(isr0)) << 2;
                    }
                }
                else {
                    // previous srcx was odd
                    if (wraparound) {
                        mscp[dstx] = (msr1 + MMUL6(msr0) + msrp
                                           + SKIPSMMaskPixelType(ma(my) ? NumericTraits<SKIPSMMaskPixelType>::one() : NumericTraits<SKIPSMMaskPixelType>::zero())
                                     ) << 2;
                        iscp[dstx] = (isr1 + IMUL6(isr0) + isrp + SKIPSMImagePixelType(sa(sy))) << 2;
                    } else {
                        mscp[dstx] = (msr1 + MMUL6(msr0) + msrp) << 2;
                        iscp[dstx] = (isr1 + IMUL6(isr0) + isrp) << 2;
                    }
                }
            }
            evenY = !evenY;
        }
    }

    // Last Rows
    {
        if (!evenY) {
            // Last srcy was even
            // odd row will set all iscp[] to zero
            // even row will do:
            //isc0[dstx] = 0;
            //isc1[dstx] = isc0[dstx] + 4*iscp[dstx]
            //out = isc1[dstx] + 6*isc0[dstx] + 4*iscp[dstx] + newisc0[dstx]
            for (dstx = 1, dx = dy, dmx = dmy; dstx < (dst_w + 1); ++dstx, ++dx.x, ++dmx.x) {
                SKIPSMMaskPixelType mp = msc1[dstx] + MMUL6(msc0[dstx]);
                if (mp) {
                    SKIPSMImagePixelType ip = (isc1[dstx] + IMUL6(isc0[dstx])) / SKIPSMImagePixelType(mp);
                    da.set(ip, dx);
                } else {
                    da.set(DestPixelType(NumericTraits<DestPixelType>::zero()), dx);
                }
                dma.set(mp ? NumericTraits<DestMaskPixelType>::max() : NumericTraits<DestMaskPixelType>::zero(), dmx);
            }
        }
        else {
            // Last srcy was odd
            // even row will do:
            // isc0[dstx] = 0;
            // isc1[dstx] = isc0[dstx] + 4*iscp[dstx]
            // out = isc1[dstx] + 6*isc0[dstx] + 4*iscp[dstx] + newisc0[dstx]
            for (dstx = 1, dx = dy, dmx = dmy; dstx < (dst_w + 1); ++dstx, ++dx.x, ++dmx.x) {
                SKIPSMMaskPixelType mp = msc1[dstx] + MMUL6(msc0[dstx]) + mscp[dstx];
                if (mp) {
                    SKIPSMImagePixelType ip = (isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx]) / SKIPSMImagePixelType(mp);
                    da.set(ip, dx);
                } else {
                    da.set(DestPixelType(NumericTraits<DestPixelType>::zero()), dx);
                }
                dma.set(mp ? NumericTraits<DestMaskPixelType>::max() : NumericTraits<DestMaskPixelType>::zero(), dmx);
            }
        }
    }

    delete[] isc0;
    delete[] isc1;
    delete[] iscp;

    delete[] msc0;
    delete[] msc1;
    delete[] mscp;

#else

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    MaskIterator my = mask_upperleft;
    DestMaskIterator dmy = dest_mask_upperleft;
    for (int srcy = 0; dy.y != dend.y; ++dy.y, ++dmy.y, sy.y+=2, my.y+=2, srcy+=2) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        MaskIterator mx = my;
        DestMaskIterator dmx = dmy;
        for (int srcx = 0; dx.x != dend.x; ++dx.x, ++dmx.x, sx.x+=2, mx.x+=2, srcx+=2) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());
            unsigned int noContrib = 10000;

            for (int kx = -2; kx <= 2; kx++) {
                int bounded_kx = kx;

                if (wraparound) {
                    // Boundary condition: wrap around the image.
                    if (srcx + kx < 0) bounded_kx += src_w;
                    if (srcx + kx >= src_w) bounded_kx -= src_w;
                } else {
                    // Boundary condition: replicate first and last column.
                    if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                    if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                }

                for (int ky = -2; ky <= 2; ky++) {
                    int bounded_ky = ky;

                    // Boundary condition: replicate top and bottom rows.
                    if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                    if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                    if (mx(bounded_kx, bounded_ky)) {
                        p += (W[kx+2] * W[ky+2]) * sx(bounded_kx, bounded_ky);
                    } else {
                        // Transparent pixels don't count.
                        noContrib -= W100[kx+2] * W100[ky+2];
                    }

                }
            }

            // Adjust filter for any ignored transparent pixels.
            if (noContrib != 0) p /= ((double)noContrib / 10000.0);

            da.set(NumericTraits<PixelType>::fromRealPromote(p), dx);
            dma.set((noContrib == 0) ? NumericTraits<DestMaskPixelType>::zero()
                                     : NumericTraits<DestMaskPixelType>::nonZero(),
                    dmx);

        }
    }
#endif

};

// Version using argument object factories.
template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void reduce(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<MaskIterator, MaskAccessor> mask,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest,
        triple<DestMaskIterator, DestMaskIterator, DestMaskAccessor> destMask) {
    reduce(wraparound,
            src.first, src.second, src.third,
            mask.first, mask.second,
            dest.first, dest.second, dest.third,
            destMask.first, destMask.second, destMask.third);
};

/** The Burt & Adelson Reduce operation.
 *  This version is for images that do not have alpha channels.
 */
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void reduce(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename NumericTraits<PixelType>::RealPromote RealPixelType;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in reduce");

/*    RealPixelType isr0, isr1, isr2, isr3;
    RealPixelType *isc0 = new RealPixelType[src_w + 2];
    RealPixelType *isc1 = new RealPixelType[src_w + 2];
    RealPixelType *isc2 = new RealPixelType[src_w + 2];
    RealPixelType *isc3 = new RealPixelType[src_w + 2];

    UInt16 msr0, msr1, msr2, msr3;
    UInt16 *msc0 = new UInt16[src_w + 2];
    UInt16 *msc1 = new UInt16[src_w + 2];
    UInt16 *msc2 = new UInt16[src_w + 2];
    UInt16 *msc3 = new UInt16[src_w + 2];

    for (int i = 0; i < src_w + 2; i++) {
        isc0[i] = RealPixelType(NumericTraits<RealPixelType>::zero());
        isc1[i] = RealPixelType(NumericTraits<RealPixelType>::zero());
        isc2[i] = RealPixelType(NumericTraits<RealPixelType>::zero());
        isc3[i] = RealPixelType(NumericTraits<RealPixelType>::zero());
        msc0[i] = NumericTraits<UInt16>::zero();
        msc1[i] = NumericTraits<UInt16>::zero();
        msc2[i] = NumericTraits<UInt16>::zero();
        msc3[i] = NumericTraits<UInt16>::zero();
    }

    DestImageIterator dy = dest_upperleft;
    SrcImageIterator sy = src_upperleft;

    bool eveny = true;
    int srcy = 0;
    for (srcy = 0; srcy < src_h + 2; ++sy.y, ++srcy) {
        SrcImageIterator sx = sy;
        DestImageIterator dx = dy;

        isr0 = RealPixelType(NumericTraits<RealPixelType>::zero());
        isr1 = RealPixelType(NumericTraits<RealPixelType>::zero());
        isr2 = RealPixelType(NumericTraits<RealPixelType>::zero());
        isr3 = RealPixelType(NumericTraits<RealPixelType>::zero());
        msr0 = NumericTraits<UInt16>::zero();
        msr1 = NumericTraits<UInt16>::zero();
        msr2 = NumericTraits<UInt16>::zero();
        msr3 = NumericTraits<UInt16>::zero();

        bool evenx = true;
        int srcx = 0;
        for (srcx = 0; srcx < src_w + 2; ++sx.x, ++srcx) {

            SKIPSM5X5((((srcy < src_h) && (srcx < src_w)) ? 1 : 0),
                      (((srcy < src_h) && (srcx < src_w)) ? NumericTraits<PixelType>::toRealPromote(sa(sx)) : RealPixelType(NumericTraits<RealPixelType>::zero())))

            if (eveny && evenx) {
                if (srcy && srcx && mtmp1) da.set(NumericTraits<PixelType>::fromRealPromote(itmp1 / mtmp1), dx, Diff2D(-1, -1));
                ++dx.x;
            }

            evenx = !evenx;
        }

        if (eveny) ++dy.y;
        eveny = !eveny;
    }

    delete[] isc0;
    delete[] isc1;
    delete[] isc2;
    delete[] isc3;
    delete[] msc0;
    delete[] msc1;
    delete[] msc2;
    delete[] msc3;
*/

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    for (int srcy = 0; dy.y != dend.y; ++dy.y, sy.y+=2, srcy+=2) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        for (int srcx = 0; dx.x != dend.x; ++dx.x, sx.x+=2, srcx+=2) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());

            for (int kx = -2; kx <= 2; kx++) {
                int bounded_kx = kx;

                if (wraparound) {
                    // Boundary condition: wrap around the image.
                    if (srcx + kx < 0) bounded_kx += src_w;
                    if (srcx + kx >= src_w) bounded_kx -= src_w;
                } else {
                    // Boundary condition: replicate first and last column.
                    if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                    if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                }

                for (int ky = -2; ky <= 2; ky++) {
                    int bounded_ky = ky;

                    // Boundary condition: replicate top and bottom rows.
                    if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                    if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                    p += (W[kx+2] * W[ky+2]) * sx(bounded_kx, bounded_ky);

                }
            }

            da.set(NumericTraits<PixelType>::fromRealPromote(p), dx);
        }
    }

};

// Version using argument object factories.
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void reduce(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest) {
    reduce(wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

#define SKIPSM_EXPAND(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
    current = SKIPSMImagePixelType(sa(sx));                         \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                         \
    out01 = sc0a[srcx];                                             \
    out11 = sc0b[srcx];                                             \
    sc1a[srcx] = sc0a[srcx];                                        \
    sc1b[srcx] = sc0b[srcx];                                        \
    sc0a[srcx] = sr1 + IMUL6(sr0) + current;                        \
    sc0b[srcx] = (sr0 + current) << 2;                              \
    sr1 = sr0;                                                      \
    sr0 = current;                                                  \
    out00 += sc0a[srcx];                                            \
    out10 += sc0b[srcx];                                            \
    out01 += sc0a[srcx];                                            \
    out11 += sc0b[srcx];                                            \
    out00 /= SKIPSMImagePixelType(SCALE_OUT00);                     \
    out10 /= SKIPSMImagePixelType(SCALE_OUT10);                     \
    out01 /= SKIPSMImagePixelType(SCALE_OUT01);                     \
    out11 /= SKIPSMImagePixelType(SCALE_OUT11);                     \
    da.set(cf(SKIPSMImagePixelType(da(dx)), out00), dx);            \
    ++dx.x;                                                         \
    da.set(cf(SKIPSMImagePixelType(da(dx)), out10), dx);            \
    ++dx.x;                                                         \
    da.set(cf(SKIPSMImagePixelType(da(dxx)), out01), dxx);          \
    ++dxx.x;                                                        \
    da.set(cf(SKIPSMImagePixelType(da(dxx)), out11), dxx);          \
    ++dxx.x;

#define SKIPSM_EXPAND_SHIFT                                         \
    current = SKIPSMImagePixelType(sa(sx));                         \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                         \
    out01 = sc0a[srcx];                                             \
    out11 = sc0b[srcx];                                             \
    sc1a[srcx] = sc0a[srcx];                                        \
    sc1b[srcx] = sc0b[srcx];                                        \
    sc0a[srcx] = sr1 + IMUL6(sr0) + current;                        \
    sc0b[srcx] = (sr0 + current) << 2;                              \
    sr1 = sr0;                                                      \
    sr0 = current;                                                  \
    out00 += sc0a[srcx];                                            \
    out10 += sc0b[srcx];                                            \
    out01 += sc0a[srcx];                                            \
    out11 += sc0b[srcx];                                            \
    out00 >>= 6;                                                    \
    out10 >>= 6;                                                    \
    out01 >>= 4;                                                    \
    out11 >>= 4;                                                    \
    da.set(cf(SKIPSMImagePixelType(da(dx)), out00), dx);            \
    ++dx.x;                                                         \
    da.set(cf(SKIPSMImagePixelType(da(dx)), out10), dx);            \
    ++dx.x;                                                         \
    da.set(cf(SKIPSMImagePixelType(da(dxx)), out01), dxx);          \
    ++dxx.x;                                                        \
    da.set(cf(SKIPSMImagePixelType(da(dxx)), out11), dxx);          \
    ++dxx.x;

#define SKIPSM_EXPAND_ROW_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                         \
    out00 /= SKIPSMImagePixelType(SCALE_OUT00);                     \
    out10 /= SKIPSMImagePixelType(SCALE_OUT10);                     \
    da.set(cf(da(dx), out00), dx);                                  \
    ++dx.x;                                                         \
    da.set(cf(da(dx), out10), dx);                                  \
    ++dx.x;                                                         \
    if ((dst_h & 1) == 0) {                                         \
        out01 = sc0a[srcx];                                         \
        out11 = sc0b[srcx];                                         \
        out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
        out11 /= SKIPSMImagePixelType(SCALE_OUT11);                 \
        da.set(cf(da(dxx), out01), dxx);                            \
        ++dxx.x;                                                    \
        da.set(cf(da(dxx), out11), dxx);                            \
        ++dxx.x;                                                    \
    }

#define SKIPSM_EXPAND_COLUMN_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out01 = sc0a[srcx];                                             \
    out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                         \
    out11 = sc0b[srcx];                                             \
    sc1a[srcx] = sc0a[srcx];                                        \
    sc1b[srcx] = sc0b[srcx];                                        \
    sc0a[srcx] = sr1 + IMUL6(sr0);                                  \
    sc0b[srcx] = sr0 << 2;                                          \
    out00 += sc0a[srcx];                                            \
    out01 += sc0a[srcx];                                            \
    out00 /= SKIPSMImagePixelType(SCALE_OUT00);                     \
    out01 /= SKIPSMImagePixelType(SCALE_OUT01);                     \
    da.set(cf(da(dx), out00), dx);                                  \
    da.set(cf(da(dxx), out01), dxx);                                \
    if ((dst_w & 1) == 0) {                                         \
        ++dx.x;                                                     \
        ++dxx.x;                                                    \
        out10 += sc0b[srcx];                                        \
        out11 += sc0b[srcx];                                        \
        out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
        out11 /= SKIPSMImagePixelType(SCALE_OUT11);                 \
        da.set(cf(da(dx), out10), dx);                              \
        da.set(cf(da(dxx), out11), dxx);                            \
    }

#define SKIPSM_EXPAND_COLUMN_END_WRAPAROUND(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out01 = sc0a[srcx];                                             \
    out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                         \
    out11 = sc0b[srcx];                                             \
    sc1a[srcx] = sc0a[srcx];                                        \
    sc1b[srcx] = sc0b[srcx];                                        \
    sc0a[srcx] = sr1 + IMUL6(sr0) + SKIPSMImagePixelType(sa(sy));   \
    sc0b[srcx] = (sr0 + SKIPSMImagePixelType(sa(sy))) << 2;         \
    out00 += sc0a[srcx];                                            \
    out01 += sc0a[srcx];                                            \
    out00 /= SKIPSMImagePixelType(SCALE_OUT00);                     \
    out01 /= SKIPSMImagePixelType(SCALE_OUT01);                     \
    da.set(cf(da(dx), out00), dx);                                  \
    da.set(cf(da(dxx), out01), dxx);                                \
    if ((dst_w & 1) == 0) {                                         \
        ++dx.x;                                                     \
        ++dxx.x;                                                    \
        out10 += sc0b[srcx];                                        \
        out11 += sc0b[srcx];                                        \
        out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
        out11 /= SKIPSMImagePixelType(SCALE_OUT11);                 \
        da.set(cf(da(dx), out10), dx);                              \
        da.set(cf(da(dxx), out11), dxx);                            \
    }

#define SKIPSM_EXPAND_ROW_COLUMN_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
    out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                         \
    out00 /= SKIPSMImagePixelType(SCALE_OUT00);                     \
    da.set(cf(da(dx), out00), dx);                                  \
    if ((dst_w & 1) == 0) {                                         \
        out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
        out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
        ++dx.x;                                                     \
        da.set(cf(da(dx), out10), dx);                              \
    }                                                               \
    if ((dst_h & 1) == 0) {                                         \
        out01 = sc0a[srcx];                                         \
        out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
        da.set(cf(da(dxx), out01), dxx);                            \
        if ((dst_w & 1) == 0) {                                     \
            out11 = sc0b[srcx];                                     \
            out11 /= SKIPSMImagePixelType(SCALE_OUT11);             \
            ++dxx.x;                                                \
            da.set(cf(da(dxx), out11), dxx);                        \
        }                                                           \
    }

// with old expand: 27.42
// with new expand: 3.54
/** The Burt & Adelson Expand operation. */
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename CombineFunctor>
void expand(bool add, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da,
        CombineFunctor cf) {

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename NumericTraits<PixelType>::RealPromote RealPixelType;
    typedef typename DestAccessor::value_type DestPixelType;
    typedef typename PyramidPromoteTraits<PixelType>::Promote SKIPSMImagePixelType;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    int dst_w = dest_lowerright.x - dest_upperleft.x;
    int dst_h = dest_lowerright.y - dest_upperleft.y;

#if 1
    // Without templatized plus/minus functor: 2.15minus 1.15plus
    // with functor: 1.94minus 1.14plus
    // with functor and fromPromote 2.08minus 1.46plus
    // with functor and fromPromote only for add: 2.16 minus 1.39 plus
    // with wraparound, not used: 1.86minus 1.36plus

    SKIPSMImagePixelType current;
    SKIPSMImagePixelType out00, out10, out01, out11;
    SKIPSMImagePixelType sr0, sr1;
    SKIPSMImagePixelType *sc0a = new SKIPSMImagePixelType[src_w + 1];
    SKIPSMImagePixelType *sc0b = new SKIPSMImagePixelType[src_w + 1];
    SKIPSMImagePixelType *sc1a = new SKIPSMImagePixelType[src_w + 1];
    SKIPSMImagePixelType *sc1b = new SKIPSMImagePixelType[src_w + 1];

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dyy = dest_upperleft;
    DestImageIterator dx = dy;
    DestImageIterator dxx = dyy;
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator sx = sy;

    int srcy = 0;
    int srcx = 0;
    //int dsty = 0;
    //int dstx = 0;

    // First row
    {
        if (wraparound) {
            sr0 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1,0)));
            sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-2,0)));
        } else {
            sr0 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            sr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        }
        for (sx = sy, srcx = 0; srcx < src_w; ++srcx, ++sx.x) {
            current = SKIPSMImagePixelType(sa(sx));
            sc0a[srcx] = sr1 + IMUL6(sr0) + current;
            sc0b[srcx] = (sr0 + current) << 2;
            sc1a[srcx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            sc1b[srcx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
            sr1 = sr0;
            sr0 = current;
        }
        // extra column at end of first row
        if (wraparound) {
            sc0a[srcx] = sr1 + IMUL6(sr0) + SKIPSMImagePixelType(sa(sy));
            sc0b[srcx] = (sr0 + SKIPSMImagePixelType(sa(sy))) << 2;
        } else {
            sc0a[srcx] = sr1 + IMUL6(sr0);
            sc0b[srcx] = sr0 << 2;
        }
        sc1a[srcx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        sc1b[srcx] = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
    }

    // dy  = row 0
    // dyy = row 1
    ++dyy.y;
    // sy = row 1
    srcy = 1;
    ++sy.y;

    // Second row
    if (src_h > 1) {
        // First column
        srcx = 0;
        sx = sy;
        sr0 = SKIPSMImagePixelType(sa(sx));
        if (wraparound) {
            sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1,0)));
        } else {
            sr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        }
        // sc*[0] are irrelvant

        srcx = 1;
        ++sx.x;
        dx = dy;
        dxx = dyy;

        // Second column
        if (src_w > 1) {
            if (wraparound) {
                SKIPSM_EXPAND(56, 56, 16, 16)
            } else {
                SKIPSM_EXPAND(49, 56, 14, 16)
            }

            // Main columns
            for (srcx = 2, ++sx.x; srcx < src_w; ++srcx, ++sx.x) {
                SKIPSM_EXPAND(56, 56, 16, 16)
            }

            // extra column at end of second row
            if (wraparound) {
                SKIPSM_EXPAND_COLUMN_END_WRAPAROUND(56, 56, 16, 16)
            } else {
                SKIPSM_EXPAND_COLUMN_END(49, 28, 14, 8)
            }
        }
        else {
            // Math works out exactly the same for wraparound and no wraparound when src_w ==1
            SKIPSM_EXPAND_COLUMN_END(42, 28, 12, 8)
        }

    }
    else {
        // No Second Row
        // First Column
        srcx = 0;
        sr0 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        sr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());

        dx = dy;
        dxx = dyy;

        if (src_w > 1) {
            // Second Column
            srcx = 1;
            if (wraparound) {
                SKIPSM_EXPAND_ROW_END(48, 48, 8, 8)
            } else {
                SKIPSM_EXPAND_ROW_END(42, 48, 7, 8)
            }

            // Main columns
            for (srcx = 2; srcx < src_w; ++srcx) {
                SKIPSM_EXPAND_ROW_END(48, 48, 8, 8)
            }

            // extra column at end of row
            if (wraparound) {
                SKIPSM_EXPAND_ROW_COLUMN_END(48, 48, 8, 8)
            } else {
                SKIPSM_EXPAND_ROW_COLUMN_END(42, 24, 7, 4)
            }
        }
        else {
            // No Second Column
            // dst_w, dst_h must be at least 2
            SKIPSM_EXPAND_ROW_COLUMN_END(36, 24, 6, 4)
        }

        delete[] sc0a;
        delete[] sc0b;
        delete[] sc1a;
        delete[] sc1b;

        return;
    }

    // dy = row 2
    // dyy = row 3
    dy.y += 2;
    dyy.y += 2;
    // sy = row 2
    srcy = 2;
    ++sy.y;

    // Main Rows
    for (srcy = 2, sx = sy; srcy < src_h; ++srcy, ++sy.y, dy.y += 2, dyy.y += 2) {
        // First column
        srcx = 0;
        sx = sy;
        sr0 = SKIPSMImagePixelType(sa(sx));
        if (wraparound) {
            sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1,0)));
        } else {
            sr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        }
        // sc*[0] are irrelvant

        srcx = 1;
        ++sx.x;
        dx = dy;
        dxx = dyy;

        // Second column
        if (src_w > 1) {
            if (wraparound) {
                SKIPSM_EXPAND_SHIFT
            } else {
                SKIPSM_EXPAND(56, 64, 14, 16)
            }

            // Main columns
            for (srcx = 2, ++sx.x; srcx < src_w; ++srcx, ++sx.x) {
                //SKIPSM_EXPAND(64, 64, 16, 16)
                SKIPSM_EXPAND_SHIFT
            }

            // extra column at end of row
            if (wraparound) {
                SKIPSM_EXPAND_COLUMN_END_WRAPAROUND(64, 64, 16, 16)
            } else {
                SKIPSM_EXPAND_COLUMN_END(56, 32, 14, 8)
            }
        }
        else {
            // No second column
            // dst_w must be at least 2
            // Math works out exactly the same for wraparound and no wraparound when src_w ==1
            SKIPSM_EXPAND_COLUMN_END(48, 32, 12, 8)
        }
    }

    // Extra row at end
    {
        srcx = 0;
        sr0 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());
        sr1 = SKIPSMImagePixelType(NumericTraits<SKIPSMImagePixelType>::zero());

        dx = dy;
        dxx = dyy;

        if (src_w > 1) {
            // Second Column
            srcx = 1;
            if (wraparound) {
                SKIPSM_EXPAND_ROW_END(56, 56, 8, 8)
            } else {
                SKIPSM_EXPAND_ROW_END(49, 56, 7, 8)
            }

            // Main columns
            for (srcx = 2; srcx < src_w; ++srcx) {
                SKIPSM_EXPAND_ROW_END(56, 56, 8, 8)
            }

            // extra column at end of row
            if (wraparound) {
                SKIPSM_EXPAND_ROW_COLUMN_END(56, 56, 8, 8)
            } else {
                SKIPSM_EXPAND_ROW_COLUMN_END(49, 28, 7, 4)
            }
        }
        else {
            // No Second Column
            // dst_w, dst_h must be at least 2
            SKIPSM_EXPAND_ROW_COLUMN_END(42, 28, 6, 4)
        }
    }

    delete[] sc0a;
    delete[] sc0b;
    delete[] sc1a;
    delete[] sc1b;

#else

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    for (int srcy = 0, desty = 0; dy.y != dend.y; ++dy.y, desty++) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        for (int srcx = 0, destx = 0; dx.x != dend.x; ++dx.x, destx++) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());
            unsigned int totalContrib = 0;

            if ((destx & 1) == 1) {
                // This is an odd-numbered column.
                for (int kx = 0; kx <= 1; kx++) {
                    // kx=0 -> W[-1]
                    // kx=1 -> W[ 1]
                    int wIndexA = (kx==0) ? 1 : 3;
                    int bounded_kx = kx;

                    if (wraparound) {
                        // Boundary condition - wrap around the image.
                        // First boundary case cannot happen when destx is odd.
                        //if (srcx + kx < 0) bounded_kx += src_w;
                        if (srcx + kx >= src_w) bounded_kx -= src_w;
                    } else {
                        // Boundary condition - replicate first and last column.
                        // First boundary case cannot happen when destx is odd.
                        //if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                        if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                    }

                    if ((desty & 1) == 1) {
                        // This is an odd-numbered row.
                        for (int ky = 0; ky <= 1; ky++) {
                            // ky=0 -> W[-1]
                            // ky=1 -> W[ 1]
                            int wIndexB = (ky==0) ? 1 : 3;
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            // First boundary condition cannot happen when desty is odd.
                            //if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                    else {
                        // This is an even-numbered row.
                        for (int ky = -1; ky <= 1; ky++) {
                            // ky=-1 -> W[-2]
                            // ky= 0 -> W[ 0]
                            // ky= 1 -> W[ 2]
                            int wIndexB = (ky==-1) ? 0 : ((ky==0) ? 2 : 4);
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                }
            }
            else {
                // This is an even-numbered column.
                for (int kx = -1; kx <= 1; kx++) {
                    // kx=-1 -> W[-2]
                    // kx= 0 -> W[ 0]
                    // kx= 1 -> W[ 2]
                    int wIndexA = (kx==-1) ? 0 : ((kx==0) ? 2 : 4); 
                    int bounded_kx = kx;

                    if (wraparound) {
                        // Boundary condition - wrap around the image.
                        if (srcx + kx < 0) bounded_kx += src_w;
                        if (srcx + kx >= src_w) bounded_kx -= src_w;
                    } else {
                        // Boundary condition - replicate first and last column.
                        if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                        if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                    }

                    if ((desty & 1) == 1) {
                        // This is an odd-numbered row.
                        for (int ky = 0; ky <= 1; ky++) {
                            // ky=0 -> W[-1]
                            // ky=1 -> W[ 1]
                            int wIndexB = (ky==0) ? 1 : 3;
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            // First boundary condition cannot happen when desty is odd.
                            //if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                    else {
                        // This is an even-numbered row.
                        for (int ky = -1; ky <= 1; ky++) {
                            // ky=-1 -> W[-2]
                            // ky= 0 -> W[ 0]
                            // ky= 1 -> W[ 2]
                            int wIndexB = (ky==-1) ? 0 : ((ky==0) ? 2 : 4);
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                }
            }

            p /= (totalContrib / 10000.0);

            // Get dest pixel at dx.
            RealPixelType pOrig = NumericTraits<PixelType>::toRealPromote(da(dx));
            RealPixelType pMod;

            // Add or subtract p.
            if (add) {
                pMod = pOrig + p;
            }
            else {
                pMod = pOrig - p;
            }

            // Write back the result.
            da.set(NumericTraits<PixelType>::fromRealPromote(pMod), dx);

            if ((destx & 1) == 1) {
                sx.x++;
                srcx++;
            }
        }

        if ((desty & 1) == 1) {
            sy.y++;
            srcy++;
        }
    }
#endif

};

template<typename T1, typename T2, typename T3>
struct FromPromotePlusFunctorWrapper : public std::binary_function<T1, T2, T3> {
    inline T3 operator()(const T1 &a, const T2 &b) const {
        return NumericTraits<T3>::fromPromote(a + b);
    }
};

// Version using argument object factories.
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void expand(bool add, bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest) {

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename PyramidPromoteTraits<PixelType>::Promote SKIPSMImagePixelType;
    typedef typename DestAccessor::value_type DestPixelType;

    if (add) {
        expand(add, wraparound, src.first, src.second, src.third, dest.first, dest.second, dest.third,
                FromPromotePlusFunctorWrapper<DestPixelType, SKIPSMImagePixelType, DestPixelType>());
    }
    else {
        expand(add, wraparound, src.first, src.second, src.third, dest.first, dest.second, dest.third,
                std::minus<SKIPSMImagePixelType>());
    }

};

/** Calculate the Gaussian pyramid for the given SrcImage/AlphaImage pair. */
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename AlphaImageType::const_traverser alpha_upperleft,
        typename AlphaImageType::ConstAccessor aa) {

    vector<PyramidImageType*> *gp = new vector<PyramidImageType*>();

    // Size of pyramid level 0
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    // Pyramid level 0
    PyramidImageType *gp0 = new PyramidImageType(w, h);

    // Copy src image into gp0, using fixed-point conversions.
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft, src_lowerright, sa,
            gp0->upperLeft(), gp0->accessor());

    gp->push_back(gp0);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Gaussian pyramid: g0";
    }

    // Make remaining levels.
    PyramidImageType *lastGP = gp0;
    AlphaImageType *lastA = NULL;
    for (unsigned int l = 1; l < numLevels; l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " g" << l;
            cout.flush();
        }

        // Size of next level
        w = (w + 1) >> 1;
        h = (h + 1) >> 1;

        // Next pyramid level
        PyramidImageType *gpn = new PyramidImageType(w, h);
        AlphaImageType *nextA = new AlphaImageType(w, h);

        if (lastA == NULL) {
            reduce(wraparound,
                    srcImageRange(*lastGP), maskIter(alpha_upperleft, aa),
                    destImageRange(*gpn), destImageRange(*nextA));
        } else {
            reduce(wraparound,
                    srcImageRange(*lastGP), maskImage(*lastA),
                    destImageRange(*gpn), destImageRange(*nextA));
        }

        gp->push_back(gpn);
        lastGP = gpn;
        delete lastA;
        lastA = nextA;
    }

    delete lastA;

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

    return gp;

};

// Version using argument object factories.
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha) {
    return gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
            numLevels, wraparound,
            src.first, src.second, src.third,
            alpha.first, alpha.second);
};

/** Calculate the Gaussian pyramid for the given image (without an alpha channel). */
template <typename SrcImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa) {

    vector<PyramidImageType*> *gp = new vector<PyramidImageType*>();

    // Size of pyramid level 0
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    // Pyramid level 0
    PyramidImageType *gp0 = new PyramidImageType(w, h);

    // Copy src image into gp0, using fixed-point conversions.
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft, src_lowerright, sa,
            gp0->upperLeft(), gp0->accessor());

    gp->push_back(gp0);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Gaussian pyramid: g0";
    }

    // Make remaining levels.
    PyramidImageType *lastGP = gp0;
    for (unsigned int l = 1; l < numLevels; l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " g" << l;
            cout.flush();
        }

        // Size of next level
        w = (w + 1) >> 1;
        h = (h + 1) >> 1;

        // Next pyramid level
        PyramidImageType *gpn = new PyramidImageType(w, h);

        reduce(wraparound, srcImageRange(*lastGP), destImageRange(*gpn));

        gp->push_back(gpn);
        lastGP = gpn;
    }

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

    return gp;
};

// Version using argument object factories.
template <typename SrcImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src) {
    return gaussianPyramid<SrcImageType, PyramidImageType>(numLevels,
            wraparound,
            src.first, src.second, src.third);
};

/** Calculate the Laplacian pyramid of the given SrcImage/AlphaImage pair. */
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
vector<PyramidImageType*> *laplacianPyramid(const char* exportName, unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename AlphaImageType::const_traverser alpha_upperleft,
        typename AlphaImageType::ConstAccessor aa) {

    // First create a Gaussian pyramid.
    vector <PyramidImageType*> *gp =
            gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
                    numLevels, wraparound,
                    src_upperleft, src_lowerright, sa,
                    alpha_upperleft, aa);

    //exportPyramid(gp, exportName);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Laplacian pyramid:";
        cout.flush();
    }

    // For each level, subtract the expansion of the next level.
    // Stop if there is no next level.
    for (unsigned int l = 0; l < (numLevels-1); l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " l" << l;
            cout.flush();
        }

        expand(false, wraparound,
                srcImageRange(*((*gp)[l+1])),
                destImageRange(*((*gp)[l])));
    }

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << " l" << (numLevels-1) << endl;
    }

    //exportPyramid(gp, exportName);

    return gp;
};

// Version using argument object factories.
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *laplacianPyramid(const char* exportName, unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha) {
    return laplacianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
            exportName,
            numLevels, wraparound,
            src.first, src.second, src.third,
            alpha.first, alpha.second);
};

/** Collapse the given Laplacian pyramid. */
template <typename PyramidImageType>
void collapsePyramid(bool wraparound, vector<PyramidImageType*> *p) {

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Collapsing Laplacian pyramid: "
             << "l" << p->size()-1;
        cout.flush();
    }

    // For each level, add the expansion of the next level.
    // Work backwards from the smallest level to the largest.
    for (int l = (p->size()-2); l >= 0; l--) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " l" << l;
            cout.flush();
        }

        expand(true, wraparound,
                srcImageRange(*((*p)[l+1])),
                destImageRange(*((*p)[l])));
    }

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

};

// Export a scalar pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, const char *prefix, VigraTrueType) {
    typedef typename PyramidImageType::value_type PyramidValueType;

    //for (unsigned int i = 0; i < (v->size() - 1); i++) {
    //    // Clear all levels except last.
    //    initImage(destImageRange(*((*v)[i])), NumericTraits<PyramidValueType>::zero());
    //}
    //collapsePyramid(false, v);

    for (unsigned int i = 0; i < v->size(); i++) {
        char filenameBuf[512];
        snprintf(filenameBuf, 512, "%s%04u.tif", prefix, i);

        // Rescale the pyramid values to fit in UINT16.
        UInt16Image usPyramid((*v)[i]->width(), (*v)[i]->height());
        transformImage(srcImageRange(*((*v)[i])), destImage(usPyramid),
                linearRangeMapping(NumericTraits<PyramidValueType>::min(),
                                     NumericTraits<PyramidValueType>::max(),
                                     NumericTraits<unsigned short>::min(),
                                     NumericTraits<unsigned short>::max()));

        ImageExportInfo info(filenameBuf);
        exportImage(srcImageRange(usPyramid), info);
    }
};

// Export a vector pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, const char *prefix, VigraFalseType) {
    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidValueType;

    //for (unsigned int i = 0; i < (v->size() - 1); i++) {
    //    // Clear all levels except last.
    //    initImage(destImageRange(*((*v)[i])), NumericTraits<PyramidValueType>::zero());
    //}
    //collapsePyramid(false, v);

    for (unsigned int i = 0; i < v->size(); i++) {
        char filenameBuf[512];
        snprintf(filenameBuf, 512, "%s%04u.tif", prefix, i);

        // Rescale the pyramid values to fit in UINT16.
        UInt16RGBImage usPyramid((*v)[i]->width(), (*v)[i]->height());
        transformImage(srcImageRange(*((*v)[i])), destImage(usPyramid),
                linearRangeMapping(PyramidVectorType(NumericTraits<PyramidValueType>::min()),
                                   PyramidVectorType(NumericTraits<PyramidValueType>::max()),
                                   typename UInt16RGBImage::value_type(NumericTraits<unsigned short>::min()),
                                   typename UInt16RGBImage::value_type(NumericTraits<unsigned short>::max())));

        ImageExportInfo info(filenameBuf);
        exportImage(srcImageRange(usPyramid), info);
    }
};

// Export a pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, const char *prefix) {
    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar pyramid_is_scalar;
    exportPyramid(v, prefix, pyramid_is_scalar());
};


} // namespace enblend

#endif /* __PYRAMID_H__ */
