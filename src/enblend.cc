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
#include <tiffio.h>

using namespace std;

int main(int argc, char** argv) {
    cout << "Left image = " << argv[1] << endl;
    cout << "Right image = " << argv[2] << endl;

    TIFF* leftTif = TIFFOpen(argv[1], "r");
    TIFF* rightTif = TIFFOpen(argv[2], "r");
    TIFF* outTif = TIFFOpen("out.tif", "w");

    uint32 w, h;
    uint16 planarConfig;
    uint16 photometric;

    TIFFGetField(leftTif, TIFFTAG_PLANARCONFIG, &planarConfig);
    TIFFGetField(leftTif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(leftTif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(leftTif, TIFFTAG_PHOTOMETRIC, &photometric);

    TIFFSetField(outTif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(outTif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(outTif, TIFFTAG_PHOTOMETRIC, photometric);
    TIFFSetField(outTif, TIFFTAG_PLANARCONFIG, planarConfig);
    TIFFSetField(outTif, TIFFTAG_ORIENTATION, 1);
    TIFFSetField(outTif, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(outTif, TIFFTAG_BITSPERSAMPLE, 8);

    uint32* leftScanline;
    uint32* rightScanline;
    uint32* outScanline;

    leftScanline = (uint32*)_TIFFmalloc(w * sizeof(uint32));
    rightScanline = (uint32*)_TIFFmalloc(w * sizeof(uint32));
    outScanline = (uint32*)_TIFFmalloc(w * sizeof(uint32));

    cout << hex;

    for (uint32 i = 0; i < h; i++) {
        TIFFReadScanline(leftTif, leftScanline, i, 8);
        TIFFReadScanline(rightTif, rightScanline, i, 8);
        for (uint32 j = 0; j < w; j++) {
            if ((TIFFGetA(leftScanline[j]) == 0xFF)
                    ^ (TIFFGetA(rightScanline[j]) == 0xFF)) {
                outScanline[j] = 0xFFFFFFFF; // white
            }
            else {
                outScanline[j] = 0xFF000000; // black
            }
        }
        TIFFWriteScanline(outTif, outScanline, i, 8);
    }

    TIFFClose(leftTif);
    TIFFClose(rightTif);
    TIFFClose(outTif);

    _TIFFfree(leftScanline);
    _TIFFfree(rightScanline);
    _TIFFfree(outScanline);

    return 0;
}
