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

// Region of interest for pyramids.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

/** Blend Laplacian pyramids whiteLP and blackLP using
 *  Gaussian pyramid maskGP.
 *  blackLPFile and maskGPFile are roiWidth * roiHeight * LPPixel.
 *  The result is written back to whiteLP.
 */
void blend(vector<LPPixel*> &whiteLP, FILE *blackLPFile, FILE *maskGPFile) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Make sure we're at the beginning of the files.
    rewind(blackLPFile);
    rewind(maskGPFile);

    if (Verbose > 0) {
        cout << "Blending layers:";
        cout.flush();
    }

    uint32 layerWidth = roiWidth;
    uint32 layerHeight = roiHeight;

    // Do each layer individually.
    for (uint32 layer = 0; layer < whiteLP.size(); layer++) {

        if (Verbose > 0) {
            cout << " l" << layer;
            cout.flush();
        }

        // layerWidth-sized line for blackLPFile.
        LPPixel *blackLine = (LPPixel*)malloc(layerWidth * sizeof(LPPixel));
        if (blackLine == NULL) {
            cerr << endl
                 << "enblend: out of memory (in blend for blackLine)" << endl;
            exit(1);
        }

        // layerWidth-sized line for maskGPFile.
        LPPixel *maskLine = (LPPixel*)malloc(layerWidth * sizeof(LPPixel));
        if (maskLine == NULL) {
            cerr << endl
                 << "enblend: out of memory (in blend for maskLine)" << endl;
            exit(1);
        }

        // Iterate over each pixel in the layer.
        LPPixel *whitePixel = whiteLP[layer];
        for (uint32 y = 0; y < layerHeight; y++) {
            readFromTmpfile(blackLine, sizeof(LPPixel), layerWidth, blackLPFile);
            readFromTmpfile(maskLine, sizeof(LPPixel), layerWidth, maskGPFile);
            LPPixel *blackPixel = blackLine;
            LPPixel *maskPixel = maskLine;

            for (uint32 x = 0; x < layerWidth; x++) {
                // whiteCoeff is the weight of whitePixel.
                double whiteCoeff = maskPixel->r / 255.0;
                // (1.0 - whiteCoeff) is the weight of blackPixel;
                double blackCoeff = 1.0 - whiteCoeff;

                double r = whitePixel->r * whiteCoeff + blackPixel->r * blackCoeff;
                double g = whitePixel->g * whiteCoeff + blackPixel->g * blackCoeff;
                double b = whitePixel->b * whiteCoeff + blackPixel->b * blackCoeff;

                // Convert back to int16
                whitePixel->r = (int16)rint(r);
                whitePixel->g = (int16)rint(g);
                whitePixel->b = (int16)rint(b);

                maskPixel++;
                blackPixel++;
                whitePixel++;
            }
        }

        free(blackLine);
        free(maskLine);

        // Calculate the size of the next layer.
        layerWidth = (layerWidth + 1) >> 1;
        layerHeight = (layerHeight + 1) >> 1;
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return;
}

