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

void blend(std::vector<LPPixel*> *maskGP,
        std::vector<LPPixel*> *inputLP,
        std::vector<LPPixel*> *outputLP) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    for (uint32 layer = 0; layer < inputLP->size(); layer++) {
        uint32 layerWidth = roiWidth >> layer;
        uint32 layerHeight = roiHeight >> layer;

        LPPixel *maskPixel = (*maskGP)[layer];
        LPPixel *inPixel = (*inputLP)[layer];
        LPPixel *outPixel = (*outputLP)[layer];

        for (uint32 index = 0; index < (layerWidth * layerHeight); index++) {
            double outRCoeff = maskPixel->r / 255.0;
            double outGCoeff = maskPixel->g / 255.0;
            double outBCoeff = maskPixel->b / 255.0;
            double outACoeff = maskPixel->a / 255.0;

            double outRD = outPixel->r * outRCoeff
                    + inPixel->r * (1.0 - outRCoeff);
            double outGD = outPixel->g * outGCoeff
                    + inPixel->g * (1.0 - outGCoeff);
            double outBD = outPixel->b * outBCoeff
                    + inPixel->b * (1.0 - outBCoeff);
            double outAD = outPixel->a * outACoeff
                    + inPixel->a * (1.0 - outACoeff);

            outPixel->r = (int16)lrint(outRD);
            outPixel->g = (int16)lrint(outGD);
            outPixel->b = (int16)lrint(outBD);
            outPixel->a = (int16)lrint(outAD);

            maskPixel++;
            inPixel++;
            outPixel++;

        }
    }

    return;
}

