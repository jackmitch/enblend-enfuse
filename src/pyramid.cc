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

const double A = 0.4;
const double W[] = {0.25 - A / 2.0, 0.25, A, 0.25, 0.25 - A / 2.0};

/** Expand in into out, either adding or subtracting from what is there.
 */
void expand(uint32 *in, uint32 inW, uint32 inH,
        uint32 *out, uint32 outW, uint32 outH,
        bool add) {
    cerr << "in = " << in << " out = " << out << endl;
    uint32 outIndex = 0;
    for (uint32 outY = 0; outY < outH; outY++) {
        for (uint32 outX = 0; outX < outW; outX++) {
            double outTmpR = 0.0;
            double outTmpG = 0.0;
            double outTmpB = 0.0;
            double outTmpA = 0.0;
            for (int m = 0; m < 5; m++) {
                if ((outX - m) & 1 == 1) continue;
                uint32 inX = abs((int32)outX - m) >> 1;
                if (inX >= inW) inX -= (2 * (inX - inW) + 2);
                for (int n = 0; n < 5; n++) {
                    if ((outY - n) & 1 == 1) continue;
                    uint32 inY = abs((int32)outY - n) >> 1;
                    if (inY >= inH) inY -= (2 * (inY - inH) + 2);

                    uint32 inPixel = in[inY * inW + inX];
                    outTmpR += W[m] * W[n] * TIFFGetR(inPixel);
                    outTmpG += W[m] * W[n] * TIFFGetG(inPixel);
                    outTmpB += W[m] * W[n] * TIFFGetB(inPixel);
                    outTmpA += W[m] * W[n] * TIFFGetA(inPixel);
                }
            }
            long outTmpRLong = lrint(outTmpR * 4.0);
            long outTmpGLong = lrint(outTmpG * 4.0);
            long outTmpBLong = lrint(outTmpB * 4.0);
            long outTmpALong = lrint(outTmpA * 4.0);

            if (add) {
                outTmpRLong = TIFFGetR(out[outIndex]) + outTmpRLong;
                outTmpGLong = TIFFGetG(out[outIndex]) + outTmpGLong;
                outTmpBLong = TIFFGetB(out[outIndex]) + outTmpBLong;
                outTmpALong = TIFFGetA(out[outIndex]) + outTmpALong;
            } else {
                outTmpRLong = TIFFGetR(out[outIndex]) - outTmpRLong;
                outTmpGLong = TIFFGetG(out[outIndex]) - outTmpGLong;
                outTmpBLong = TIFFGetB(out[outIndex]) - outTmpBLong;
                outTmpALong = TIFFGetA(out[outIndex]) - outTmpALong;
            }

            out[outIndex] = (outTmpRLong & 0xFF)
                    | ((outTmpGLong & 0xFF) << 8)
                    | ((outTmpBLong & 0xFF) << 16)
                    | ((outTmpALong & 0xFF) << 24);
            outIndex++;
        }
    }

    return;
}


uint32 *reduce(uint32 *in, uint32 w, uint32 h) {
    uint32 outW = w >> 1;
    uint32 outH = h >> 1;

    uint32* out = (uint32*)calloc(outW * outH, sizeof(uint32));
    if (out == NULL) {
        cerr << "enblend: malloc failed in reduce." << endl;
        exit(1);
    }

    for (uint32 outY = 0; outY < outH; outY++) {
        for (uint32 outX = 0; outX < outW; outX++) {
            double outTmpR = 0.0;
            double outTmpG = 0.0;
            double outTmpB = 0.0;
            double outTmpA = 0.0;
            for (int m = 0; m < 5; m++) {
                uint32 inX = abs(2 * (int)outX + m - 2);
                if (inX >= w) inX -= (2 * (inX - w) + 2);

                for (int n = 0; n < 5; n++) {
                    uint32 inY = abs(2 * (int)outY + n - 2);
                    if (inY >= h) inY -= (2 * (inY - h) + 2);

                    uint32 inPixel = in[inY * w + inX];
                    outTmpR += W[m] * W[n] * TIFFGetR(inPixel);
                    outTmpG += W[m] * W[n] * TIFFGetG(inPixel);
                    outTmpB += W[m] * W[n] * TIFFGetB(inPixel);
                    outTmpA += W[m] * W[n] * TIFFGetA(inPixel);
                }
            }
            long outTmpRLong = lrint(outTmpR);
            long outTmpGLong = lrint(outTmpG);
            long outTmpBLong = lrint(outTmpB);
            long outTmpALong = lrint(outTmpA);
            out[outY * outW + outX] = (outTmpRLong & 0xFF)
                    | ((outTmpGLong & 0xFF) << 8)
                    | ((outTmpBLong & 0xFF) << 16)
                    | ((outTmpALong & 0xFF) << 24);
        }
    }

    return out;
}

vector<uint32*> *gaussianPyramid(uint32 *mask, int32 levels) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    vector<uint32*> *v = new vector<uint32*>();

    if (Verbose > 0) {
        cout << "Generating gaussian pyramid g0" << endl;
    }

    // Build level 0
    uint32* g = (uint32*)malloc(roiWidth * roiHeight * sizeof(uint32));
    if (g == NULL) {
        cerr << "enblend: malloc failed for gaussian pyramid g0" << endl;
        exit(1);
    }

    // Copy mask ROI verbatim into g0.
    uint32 index = 0;
    for (uint32 maskY = ROIFirstY; maskY <= ROILastY; maskY++) {
        for (uint32 maskX = ROIFirstX; maskX <= ROILastX; maskX++) {
            g[index++] = mask[maskY * OutputWidth + maskX];
        }
    }

    levels--;
    cerr << "pushing g = " << g << endl;
    v->push_back(g);

    while (levels > 0) {
        if (Verbose > 0) {
            cout << "Generating gaussian pyramid g" << v->size() << endl;
        }

        g = reduce(g, roiWidth, roiHeight);
        v->push_back(g);
        cerr << "pushing g = " << g << endl;

        roiWidth = roiWidth >> 1;
        roiHeight = roiHeight >> 1;
        levels--;
    }

    return v;
}

vector<uint32*> *laplacianPyramid(uint32 *image, int32 levels) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    vector<uint32*> *gp = gaussianPyramid(image, levels);

    for (int32 i = 0; i < levels - 1; i++) {
        expand((*gp)[i + 1], roiWidth >> (i+1), roiHeight >> (i+1),
            (*gp)[i], roiWidth >> i, roiHeight >> i,
            false);
    }

    return gp;
}

vector<uint32*> *laplacianPyramid(TIFF *image, int32 levels) {

    // Allocate memory for the TIFF.
    uint32 *imageBuf = (uint32*)_TIFFmalloc(
            OutputWidth * OutputHeight * sizeof(uint32));
    if (imageBuf == NULL) {
        cerr << "enblend: malloc failed for imageBuf" << endl;
        exit(1);
    }

    for (uint32 i = 0; i < OutputHeight; i++) {
        TIFFReadScanline(image,
                &(imageBuf[i * OutputWidth]),
                i,
                8);
    }

    vector<uint32*> *lp = laplacianPyramid(imageBuf, levels);

    // lp contains the pyramid levels we need, the full-size input image
    // can be discarded now.
    _TIFFfree(imageBuf);

    return lp;

}

void collapsePyramid(vector<uint32*> *p, uint32 *dest) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    for (int i = p->size() - 2; i >= 0; i--) {
        expand((*p)[i + 1], roiWidth >> (i+1), roiHeight >> (i+1),
            (*p)[i], roiWidth >> i, roiHeight >> i,
            true);
    }

    cout << "p0 = " << (*p)[0] << endl;
    // Copy p0 into dest ROI, omitting partially transparent pixels.
    for (uint32 j = 0; j < roiHeight; j++) {
        for (uint32 i = 0; i < roiWidth; i++) {
            uint32 pPixel = ((*p)[0])[j * roiWidth + i];
            if (TIFFGetA(pPixel) != 0xFF) {
                dest[(j+ROIFirstY) * OutputWidth + (i+ROIFirstX)] = TRANS;
            } else {
                dest[(j+ROIFirstY) * OutputWidth + (i+ROIFirstX)] = pPixel;
            }
        }
    }

    return;
}
