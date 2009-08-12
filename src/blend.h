/*
 * Copyright (C) 2004-2007 Andrew Mihal
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
#ifndef __BLEND_H__
#define __BLEND_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fixmath.h"

#include "vigra/combineimages.hxx"
#include "vigra/numerictraits.hxx"

//#include <sh/sh.hpp>

using std::cout;
using std::vector;

using vigra::combineThreeImages;
using vigra::NumericTraits;

//using namespace SH;

/*
#include <brook/brook.hpp>

void gpuBlendKernel(::brook::stream black,
                    ::brook::stream white,
                    ::brook::stream mask,
                    const float scale,
                    const float clamp_min,
                    const float clamp_max,
                    ::brook::stream o);
*/

namespace enblend {

/** Functor for blending a black and white pyramid level using a mask
 *  pyramid level.
 */
template <typename MaskPixelType>
class CartesianBlendFunctor {
public:
    CartesianBlendFunctor(MaskPixelType w) : white(NumericTraits<MaskPixelType>::toRealPromote(w)) {}

    template <typename ImagePixelType>
    ImagePixelType operator()(const MaskPixelType &maskP, const ImagePixelType &wP, const ImagePixelType &bP) const {

        typedef typename NumericTraits<ImagePixelType>::RealPromote RealImagePixelType;

        // Convert mask pixel to blend coefficient in range [0.0, 1.0].
        double whiteCoeff =
                NumericTraits<MaskPixelType>::toRealPromote(maskP) / white;
        // Sometimes masked image data is invalid.  For floating point samples
        // this includes possible NaN's in the data.   In that case, computing
        // the output sample will result in a NaN output if the weight on that
        // pixel is 0 (since 0*NaN = NaN )
        // Handle this by explicitly ignoring fully masked pixels
        if (whiteCoeff>=1.0)
            return wP;
        if (whiteCoeff<=0.0)
            return bP;
        double blackCoeff = 1.0 - whiteCoeff;

        RealImagePixelType rwP = NumericTraits<ImagePixelType>::toRealPromote(wP);
        RealImagePixelType rbP = NumericTraits<ImagePixelType>::toRealPromote(bP);

        RealImagePixelType blendP = (whiteCoeff * rwP) + (blackCoeff * rbP);

        return NumericTraits<ImagePixelType>::fromRealPromote(blendP);
    }

protected:
    double white;
};

/*
template <typename PROGRAM, typename CHANNEL>
void profile_run_program(PROGRAM & blend_prg, CHANNEL & result) {
    result = blend_prg;
};

template <typename MX, typename WX, typename BX>
void profile_to_float(unsigned int width, MX & mx, WX & wx, BX & bx, float *mask_ptr, float *white_ptr, float *black_ptr) {
    for (unsigned int i = 0, j = 0; i < width; ++i, ++mx.x, ++wx.x, ++bx.x) {
        mask_ptr[i] = static_cast<float>(*mx);
        white_ptr[j] = static_cast<float>(wx->red());
        black_ptr[j++] = static_cast<float>(bx->red());
        white_ptr[j] = static_cast<float>(wx->green());
        black_ptr[j++] = static_cast<float>(bx->green());
        white_ptr[j] = static_cast<float>(wx->blue());
        black_ptr[j++] = static_cast<float>(bx->blue());
    }
};

template <typename ImagePixelComponentType, typename BX>
void profile_from_float(BX & bx, float *results, unsigned int width) {
    for (unsigned int i = 0, j = 0; i < width; ++i, ++bx.x) {
        *bx = RGBValue<ImagePixelComponentType>(
                static_cast<ImagePixelComponentType>(lrintf(results[j++])),
                static_cast<ImagePixelComponentType>(lrintf(results[j++])),
                static_cast<ImagePixelComponentType>(lrintf(results[j++])));
    }
};
*/

/** Blend black and white pyramids using mask pyramid.
 */
template <typename MaskPyramidType, typename ImagePyramidType>
void blend(vector<MaskPyramidType*> *maskGP,
        vector<ImagePyramidType*> *whiteLP,
        vector<ImagePyramidType*> *blackLP,
        typename MaskPyramidType::value_type maskPyramidWhiteValue) {

    //typedef typename ImagePyramidType::value_type::value_type ImagePixelComponentType;

    if (Verbose > VERBOSE_BLEND_MESSAGES) {
        cout << "Blending layers:             ";
        cout.flush();
    }

/*
    float scale = static_cast<float>(maskPyramidWhiteValue);
    float clamp_min = static_cast<float>(NumericTraits<ImagePixelComponentType>::min());
    float clamp_max = static_cast<float>(NumericTraits<ImagePixelComponentType>::max());
*/

/*
    shInit();

    ShAttrib1f scale(static_cast<float>(maskPyramidWhiteValue));
    ShAttrib3f pyramidMin(static_cast<float>(NumericTraits<ImagePixelComponentType>::min()),
                          static_cast<float>(NumericTraits<ImagePixelComponentType>::min()),
                          static_cast<float>(NumericTraits<ImagePixelComponentType>::min()));
    ShAttrib3f pyramidMax(static_cast<float>(NumericTraits<ImagePixelComponentType>::max()),
                          static_cast<float>(NumericTraits<ImagePixelComponentType>::max()),
                          static_cast<float>(NumericTraits<ImagePixelComponentType>::max()));

    ShProgram prg = SH_BEGIN_PROGRAM("stream") {
        ShInputAttrib1f mask;
        ShInputAttrib3f white;
        ShInOutAttrib3f black;
        ShAttrib1f whiteCoeff = mask / scale;
        ShAttrib1f blackCoeff = ShAttrib1f(1.0) - whiteCoeff;
        black = SH::min(SH::max((whiteCoeff * white) + (blackCoeff * black), pyramidMin), pyramidMax);
    } SH_END;
*/

    for (unsigned int layer = 0; layer < maskGP->size(); layer++) {

        if (Verbose > VERBOSE_BLEND_MESSAGES) {
            cout << " l" << layer;
            cout.flush();
        }

        //if (!UseGPU) {
        if (1) {
            combineThreeImages(srcImageRange(*((*maskGP)[layer])),
                    srcImage(*((*whiteLP)[layer])),
                    srcImage(*((*blackLP)[layer])),
                    destImage(*((*blackLP)[layer])),
                    CartesianBlendFunctor<typename MaskPyramidType::value_type>(maskPyramidWhiteValue));
            continue;
        }

/*
        MaskPyramidType &maskImage = *((*maskGP)[layer]);
        ImagePyramidType &whiteImage = *((*whiteLP)[layer]);
        ImagePyramidType &blackImage = *((*blackLP)[layer]);

        typename MaskPyramidType::traverser my = maskImage.upperLeft();
        typename MaskPyramidType::traverser mend = maskImage.lowerRight();
        typename ImagePyramidType::traverser wy = whiteImage.upperLeft();
        typename ImagePyramidType::traverser by = blackImage.upperLeft();

        unsigned int width = mend.x - my.x;

        float *mask_data = new float[width];
        float *white_data = new float[width * 3];
        float *black_data = new float[width * 3];

        ::brook::stream black_stream(::brook::getStreamType((float3*)0), width, -1);
        ::brook::stream white_stream(::brook::getStreamType((float3*)0), width, -1);
        ::brook::stream output_stream(::brook::getStreamType((float3*)0), width, -1);
        ::brook::stream mask_stream(::brook::getStreamType((float*)0), width, -1);

        for (; my.y != mend.y; ++my.y, ++wy.y, ++by.y) {
            typename MaskPyramidType::traverser mx = my;
            typename ImagePyramidType::traverser wx = wy;
            typename ImagePyramidType::traverser bx = by;

            profile_to_float(width, mx, wx, bx, mask_data, white_data, black_data);

            streamRead(black_stream, black_data);
            streamRead(white_stream, white_data);
            streamRead(mask_stream, mask_data);

            gpuBlendKernel(black_stream, white_stream, mask_stream, scale, clamp_min, clamp_max, output_stream);

            streamWrite(output_stream, black_data);

            bx = by;
            profile_from_float<ImagePixelComponentType>(bx, black_data, width);
        }

        delete[] mask_data;
        delete[] white_data;
        delete[] black_data;
*/
/*
        ShHostMemoryPtr mask_data = new ShHostMemory(sizeof(float) * width, SH_FLOAT);
        ShHostMemoryPtr white_data = new ShHostMemory(sizeof(float) * width * 3, SH_FLOAT);
        ShHostMemoryPtr black_data = new ShHostMemory(sizeof(float) * width * 3, SH_FLOAT);

        ShChannel<ShAttrib1f> mask_channel(mask_data, width);
        ShChannel<ShAttrib3f> white_channel(white_data, width);
        ShChannel<ShAttrib3f> black_channel(black_data, width);

        ShHostStoragePtr mask_data_storage = shref_dynamic_cast<ShHostStorage>(mask_channel.memory()->findStorage("host"));
        ShHostStoragePtr white_data_storage = shref_dynamic_cast<ShHostStorage>(white_channel.memory()->findStorage("host"));
        ShHostStoragePtr black_data_storage = shref_dynamic_cast<ShHostStorage>(black_channel.memory()->findStorage("host"));

        ShProgram blend_prg = prg << mask_channel << white_channel << black_channel;

        for (; my.y != mend.y; ++my.y, ++wy.y, ++by.y) {
            typename MaskPyramidType::traverser mx = my;
            typename ImagePyramidType::traverser wx = wy;
            typename ImagePyramidType::traverser bx = by;

            mask_data_storage->dirty();
            white_data_storage->dirty();
            black_data_storage->dirty();

            float* mask_ptr = (float*)mask_data_storage->data();
            float* white_ptr = (float*)white_data_storage->data();
            float* black_ptr = (float*)black_data_storage->data();

            profile_to_float(mx, width, wx, bx, mask_ptr, white_ptr, black_ptr);

            profile_run_program(blend_prg, black_channel);

            black_data_storage->sync();
            float *results = (float*)black_data_storage->data();

            bx = by;
            profile_from_float<ImagePixelComponentType>(bx, results, width);
        }
*/

    }

    if (Verbose > VERBOSE_BLEND_MESSAGES) {
        cout << endl;
    }

};

} // namespace enblend

#endif /* __BLEND_H__ */
