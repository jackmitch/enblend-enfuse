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
#include <list>

#include <string.h>
#include <tiffio.h>

#include "enblend.h"

using namespace std;

// Global values from command line parameters.
extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of filenames.
 *  Returns an output-sized buffer created with _TIFFmalloc.
 */
FILE *assemble(std::list<char*> &filenames, bool pickOne) {

    if (filenames.empty()) return NULL;

    uint32 *image = (uint32*)calloc(OutputWidth * OutputHeight, sizeof(uint32));
    if (image == NULL) {
        cerr << "enblend: out of memory in assemble for image." << endl;
        exit(1);
    }

    uint32 *scanline = (uint32*)calloc(OutputWidth, sizeof(uint32));
    if (scanline == NULL) {
        cerr << "enblend: out of memory in assemble for scanline." << endl;
        exit(1);
    }

    // List of images we decide to assemble.
    list<list<char*>::iterator> toBeRemoved;

    list<char*>::iterator i;
    for (i = filenames.begin(); i != filenames.end(); i++) {
        char *filename = *i;

        if (Verbose > 0) {
            cout << "assemble: checking \"" << filename << "\"" << endl;
        }

        TIFF *tiff = TIFFOpen(filename, "r");
        if (tiff == NULL) {
            cerr << "enblend: error opening TIFF file \""
                 << filename
                 << "\""
                 << endl;
            exit(1);
        }

        // Check for overlap.
        bool overlapFound = false;
        for (uint32 y = 0; y < OutputHeight; y++) {
            TIFFReadScanline(tiff, scanline, y, 8);

            for (uint32 x = 0; x < OutputWidth; x++) {
                uint32 imagePixel = image[y * OutputWidth + x];
                uint32 tiffPixel = scanline[x];
                if (TIFFGetA(imagePixel) != 0
                        && TIFFGetA(tiffPixel) != 0) {
                    overlapFound = true;
                    break;
                }
            }

            if (overlapFound) break;
        }

        if (!overlapFound) {
            // No overlap. Copy tiff into image.
            for (uint32 y = 0; y < OutputHeight; y++) {
                TIFFReadScanline(tiff, scanline, y, 8);

                for (uint32 x = 0; x < OutputWidth; x++) {
                    if (TIFFGetA(scanline[x]) != 0) {
                        image[y * OutputWidth + x] = scanline[x];
                    }
                }
            }

            // Remove filename from list later.
            toBeRemoved.push_back(i);
        }

        TIFFClose(tiff);

        // If we only want one image at a time, break.
        if (!overlapFound && pickOne) break;
    }

    free(scanline);

    if (Verbose > 0) {
        cout << "assemble: combined";
    }

    // Erase the filenames we used.
    list<list<char*>::iterator>::iterator r;
    for (r = toBeRemoved.begin(); r != toBeRemoved.end(); r++) {
        if (Verbose > 0) {
            cout << " " << **r;
        }
        filenames.erase(*r);
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return dumpToTmpfile(image, sizeof(uint32), OutputWidth * OutputHeight);

}
