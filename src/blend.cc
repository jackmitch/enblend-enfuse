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
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <math.h>

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

/** Blend Laplacian pyramids whiteLP and blackLP using
 *  Gaussian pyramid maskGP.
 *  The result is written back to whiteLP.
 *  This function only operates within the region-of-interest.
 */
void blend(std::vector<LPPixel*> &whiteLP,
        std::vector<LPPixel*> &blackLP,
        std::vector<LPPixel*> &maskGP) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Do each layer individually.
    for (uint32 layer = 0; layer < blackLP.size(); layer++) {

        // Calculate the size of the layer.
        uint32 layerWidth = roiWidth >> layer;
        uint32 layerHeight = roiHeight >> layer;

        // Iterate over each pixel in the layer.
        LPPixel *maskPixel = maskGP[layer];
        LPPixel *blackPixel = blackLP[layer];
        LPPixel *whitePixel = whiteLP[layer];
        for (uint32 i = 0; i < (layerWidth * layerHeight); i++) {

            // whiteCoeff is the weight of whitePixel.
            double whiteCoeff = maskPixel->r / 255.0;
            // (1.0 - whiteCoeff) is the weight of blackPixel;
            double blackCoeff = 1.0 - whiteCoeff;

            double r = whitePixel->r * whiteCoeff + blackPixel->r * blackCoeff;
            double g = whitePixel->g * whiteCoeff + blackPixel->g * blackCoeff;
            double b = whitePixel->b * whiteCoeff + blackPixel->b * blackCoeff;
            double a = whitePixel->a * whiteCoeff + blackPixel->a * blackCoeff;

            // Convert back to int16
            whitePixel->r = (int16)lrint(r);
            whitePixel->g = (int16)lrint(g);
            whitePixel->b = (int16)lrint(b);
            whitePixel->a = (int16)lrint(a);

            maskPixel++;
            blackPixel++;
            whitePixel++;
        }

    }

    return;
}

