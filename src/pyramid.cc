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
void expand(LPPixel *in, uint32 inW, uint32 inH,
        LPPixel *out, uint32 outW, uint32 outH,
        bool add) {

    cerr << "in = " << in << " out = " << out << endl;

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

                    LPPixel *inPixel = &(in[inY * inW + inX]);
                    outTmpR += W[m] * W[n] * inPixel->r;
                    outTmpG += W[m] * W[n] * inPixel->g;
                    outTmpB += W[m] * W[n] * inPixel->b;
                    outTmpA += W[m] * W[n] * inPixel->a;
                }
            }

            if (add) {
                out->r += (int16)lrint(outTmpR * 4.0);
                out->g += (int16)lrint(outTmpG * 4.0);
                out->b += (int16)lrint(outTmpB * 4.0);
                out->a += (int16)lrint(outTmpA * 4.0);
            } else {
                out->r -= (int16)lrint(outTmpR * 4.0);
                out->g -= (int16)lrint(outTmpG * 4.0);
                out->b -= (int16)lrint(outTmpB * 4.0);
                out->a -= (int16)lrint(outTmpA * 4.0);
            }

            out++;

            //long outTmpRLong = lrint(outTmpR * 4.0);
            //long outTmpGLong = lrint(outTmpG * 4.0);
            //long outTmpBLong = lrint(outTmpB * 4.0);
            //long outTmpALong = lrint(outTmpA * 4.0);
            ////if (outTmpRLong > 255 || outTmpGLong > 255 || outTmpBLong > 255 || outTmpALong > 255 || outTmpRLong < 0 || outTmpGLong < 0 || outTmpBLong < 0 || outTmpALong < 0) {
            ////    cerr << "expand overflow1 r=" << outTmpRLong
            ////         << " b=" << outTmpBLong
            ////         << " g=" << outTmpGLong
            ////         << " a=" << outTmpALong << endl;
            ////}

            //if (add) {
            //    outTmpRLong = TIFFGetR(out[outIndex]) + outTmpRLong;
            //    outTmpGLong = TIFFGetG(out[outIndex]) + outTmpGLong;
            //    outTmpBLong = TIFFGetB(out[outIndex]) + outTmpBLong;
            //    outTmpALong = TIFFGetA(out[outIndex]) + outTmpALong;
            //} else {
            //    outTmpRLong = TIFFGetR(out[outIndex]) - outTmpRLong;
            //    outTmpGLong = TIFFGetG(out[outIndex]) - outTmpGLong;
            //    outTmpBLong = TIFFGetB(out[outIndex]) - outTmpBLong;
            //    outTmpALong = TIFFGetA(out[outIndex]) - outTmpALong;
            //}
            ////if (outTmpRLong > 255 || outTmpGLong > 255 || outTmpBLong > 255 || outTmpALong > 255 || outTmpRLong < 0 || outTmpGLong < 0 || outTmpBLong < 0 || outTmpALong < 0) {
            ////    cerr << "expand overflow2 r=" << outTmpRLong
            ////         << " b=" << outTmpBLong
            ////         << " g=" << outTmpGLong
            ////         << " a=" << outTmpALong << endl;
            ////}

            //out[outIndex] = (outTmpRLong & 0xFF)
            //        | ((outTmpGLong & 0xFF) << 8)
            //        | ((outTmpBLong & 0xFF) << 16)
            //        | ((outTmpALong & 0xFF) << 24);
            //outIndex++;
        }
    }

    return;
}


LPPixel *reduce(LPPixel *in, uint32 w, uint32 h) {
    uint32 outW = w >> 1;
    uint32 outH = h >> 1;

    LPPixel *out = (LPPixel*)calloc(outW * outH, sizeof(LPPixel));
    if (out == NULL) {
        cerr << "enblend: malloc failed in reduce." << endl;
        exit(1);
    }

    LPPixel *outIndex = out;
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

                    LPPixel *inPixel = &(in[inY * w + inX]);
                    outTmpR += W[m] * W[n] * inPixel->r;
                    outTmpG += W[m] * W[n] * inPixel->g;
                    outTmpB += W[m] * W[n] * inPixel->b;
                    outTmpA += W[m] * W[n] * inPixel->a;
                }
            }
            outIndex->r = (int16)lrint(outTmpR);
            outIndex->g = (int16)lrint(outTmpG);
            outIndex->b = (int16)lrint(outTmpB);
            outIndex->a = (int16)lrint(outTmpA);

            outIndex++;

            //if (outTmpRLong > 255 || outTmpGLong > 255 || outTmpBLong > 255 || outTmpALong > 255 || outTmpRLong < 0 || outTmpGLong < 0 || outTmpBLong < 0 || outTmpALong < 0) {
            //    cerr << "reduce overflow r=" << outTmpRLong
            //         << " b=" << outTmpBLong
            //         << " g=" << outTmpGLong
            //         << " a=" << outTmpALong << endl;
            //}

            //out[outY * outW + outX] = (outTmpRLong & 0xFF)
            //        | ((outTmpGLong & 0xFF) << 8)
            //        | ((outTmpBLong & 0xFF) << 16)
            //        | ((outTmpALong & 0xFF) << 24);
        }
    }

    return out;
}

vector<LPPixel*> *gaussianPyramid(uint32 *mask, int32 levels) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    vector<LPPixel*> *v = new vector<LPPixel*>();

    if (Verbose > 0) {
        cout << "Generating gaussian pyramid g0" << endl;
    }

    // Build level 0
    LPPixel *g = (LPPixel*)malloc(roiWidth * roiHeight * sizeof(LPPixel));
    if (g == NULL) {
        cerr << "enblend: malloc failed for gaussian pyramid g0" << endl;
        exit(1);
    }

    // Copy mask ROI verbatim into g0.
    LPPixel *gIndex = g;
    for (uint32 maskY = ROIFirstY; maskY <= ROILastY; maskY++) {
        for (uint32 maskX = ROIFirstX; maskX <= ROILastX; maskX++) {
            uint32 maskPixel = mask[maskY * OutputWidth + maskX];
            gIndex->r = TIFFGetR(maskPixel);
            gIndex->g = TIFFGetG(maskPixel);
            gIndex->b = TIFFGetB(maskPixel);
            gIndex->a = TIFFGetA(maskPixel);
            gIndex++;
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

vector<LPPixel*> *laplacianPyramid(uint32 *image, int32 levels) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    vector<LPPixel*> *gp = gaussianPyramid(image, levels);

    for (int32 i = 0; i < levels - 1; i++) {
        expand((*gp)[i + 1], roiWidth >> (i+1), roiHeight >> (i+1),
            (*gp)[i], roiWidth >> i, roiHeight >> i,
            false);
    }

    return gp;
}

vector<LPPixel*> *laplacianPyramid(TIFF *image, int32 levels) {

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

    vector<LPPixel*> *lp = laplacianPyramid(imageBuf, levels);

    // lp contains the pyramid levels we need, the full-size input image
    // can be discarded now.
    _TIFFfree(imageBuf);

    return lp;

}

void collapsePyramid(vector<LPPixel*> *p, uint32 *dest) {
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    for (int i = p->size() - 2; i >= 0; i--) {
        expand((*p)[i + 1], roiWidth >> (i+1), roiHeight >> (i+1),
            (*p)[i], roiWidth >> i, roiHeight >> i,
            true);
    }

    cout << "p0 = " << (*p)[0] << endl;
    // Copy p0 into dest ROI, omitting partially transparent pixels.
    LPPixel *pixel = (*p)[0];
    for (uint32 j = 0; j < roiHeight; j++) {
        for (uint32 i = 0; i < roiWidth; i++) {
            pixel->r = min(255, max(0, (int)pixel->r));
            pixel->g = min(255, max(0, (int)pixel->g));
            pixel->b = min(255, max(0, (int)pixel->b));
            pixel->a = min(255, max(0, (int)pixel->a));
            uint32 p = (pixel->r & 0xFF)
                    | ((pixel->g & 0xFF) << 8)
                    | ((pixel->b & 0xFF) << 16)
                    | (0xFF << 24);
            dest[(j+ROIFirstY) * OutputWidth + (i+ROIFirstX)] = p;
            pixel++;
        }
    }

    return;
}
