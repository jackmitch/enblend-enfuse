/*
 * Copyright (C) 2004 Andrew Mihal
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
 * along with Foobar; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <tiffio.h>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

// Region of interest for this operation.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

void blend(std::vector<uint32*> *mask,
        std::vector<uint32*> *inputLP,
        std::vector<uint32*> *outputLP) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    for (uint32 layer = 0; layer < inputLP->size(); layer++) {
        uint32 layerWidth = roiWidth >> layer;
        uint32 layerHeight = roiHeight >> layer;
        for (uint32 index = 0; index < (layerWidth * layerHeight); index++) {
            uint32 maskPixel = ((*mask)[layer])[index];
            uint32 inPixel = ((*inputLP)[layer])[index];
            uint32 *outPixel = &(((*outputLP)[layer])[index]);

            uint32 maskR = TIFFGetR(maskPixel);
            uint32 maskG = TIFFGetG(maskPixel);
            uint32 maskB = TIFFGetB(maskPixel);
            uint32 maskA = TIFFGetA(maskPixel);

            uint32 inR = TIFFGetR(inPixel);
            uint32 inG = TIFFGetG(inPixel);
            uint32 inB = TIFFGetB(inPixel);
            uint32 inA = TIFFGetA(inPixel);

            uint32 outR = TIFFGetR(*outPixel);
            uint32 outG = TIFFGetG(*outPixel);
            uint32 outB = TIFFGetB(*outPixel);
            uint32 outA = TIFFGetA(*outPixel);

            double outRCoeff = maskR / 255.0;
            double outGCoeff = maskG / 255.0;
            double outBCoeff = maskB / 255.0;
            double outACoeff = maskA / 255.0;

            double outRD = outR * outRCoeff + inR * (1.0 - outRCoeff);
            double outGD = outG * outGCoeff + inG * (1.0 - outGCoeff);
            double outBD = outB * outBCoeff + inB * (1.0 - outBCoeff);
            double outAD = outA * outACoeff + inA * (1.0 - outACoeff);

            long outRL = lrint(outRD);
            long outGL = lrint(outGD);
            long outBL = lrint(outBD);
            long outAL = lrint(outAD);

            *outPixel = (outRL & 0xFF)
                    | ((outGL & 0xFF) << 8)
                    | ((outBL & 0xFF) << 16)
                    | ((outAL & 0xFF) << 24);

        }
    }

    return;
}

