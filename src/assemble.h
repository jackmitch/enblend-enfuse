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
#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common.h"
#include "vigra/copyimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/impexalpha.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"

using namespace std;

/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of ImageImportInfos.
 *  Returns an ImageImportInfo for the temporary file.
 *  memory xsection = 2 * (ImageType*os + AlphaType*os)
 */
template <typename ImageType, typename AlphaType>
vigra::ImageImportInfo *assemble(list<vigra::ImageImportInfo*> &imageInfoList,
        EnblendROI &inputUnion,
        EnblendROI &bb) {

    typedef typename AlphaType::PixelType AlphaPixelType;
    typedef typename AlphaType::Iterator AlphaIteratorType;

    // No more images to assemble?
    if (imageInfoList.empty()) return NULL;

    // Create an image to assemble input images into.
    ImageType image(inputUnion.size());
    AlphaType imageA(inputUnion.size());

    if (Verbose > 0) {
        if (OneAtATime) {
            cout << "Loading next image: "
                 << imageInfoList.front()->getFileName()
                 << endl;
        } else {
            cout << "Combining non-overlapping images: "
                 << imageInfoList.front()->getFileName();
            cout.flush();
        }
    }

    // Load the first image into the destination.
    vigra::Diff2D imagePos = imageInfoList.front()->getPosition();
    vigra::importImageAlpha(*imageInfoList.front(),
            destIter(image.upperLeft() + imagePos - inputUnion.getUL()),
            destIter(imageA.upperLeft() + imagePos - inputUnion.getUL()));
    imageInfoList.erase(imageInfoList.begin());

    // Mask off pixels that are not totally opaque.
    vigra::transformImage(srcImageRange(imageA), destImage(imageA),
            vigra::Threshold<AlphaPixelType, AlphaPixelType>(
                    GetMaxAlpha<AlphaPixelType>(),
                    GetMaxAlpha<AlphaPixelType>(),
                    0,
                    GetMaxAlpha<AlphaPixelType>()
            )
    );

    if (!OneAtATime) {
        // Attempt to assemble additional non-overlapping images.

        // List of ImageImportInfos we decide to assemble.
        list<list<vigra::ImageImportInfo*>::iterator> toBeRemoved;

        list<vigra::ImageImportInfo*>::iterator i;
        for (i = imageInfoList.begin(); i != imageInfoList.end(); i++) {
            vigra::ImageImportInfo *info = *i;

            // Load the next image.
            ImageType src(info->size());
            AlphaType srcA(info->size());
            vigra::importImageAlpha(*info, destImage(src), destImage(srcA));

            // Mask off pixels that are not totally opaque.
            vigra::transformImage(srcImageRange(srcA), destImage(srcA),
                    vigra::Threshold<AlphaPixelType, AlphaPixelType>(
                            GetMaxAlpha<AlphaPixelType>(),
                            GetMaxAlpha<AlphaPixelType>(),
                            0,
                            GetMaxAlpha<AlphaPixelType>()
                    )
            );

            // Check for overlap.
            bool overlapFound = false;
            AlphaIteratorType dy =
                    imageA.upperLeft() - inputUnion.getUL() + info->getPosition();
            AlphaIteratorType sy = srcA.upperLeft();
            AlphaIteratorType send = srcA.lowerRight();
            for(; sy.y != send.y; ++sy.y, ++dy.y) {
                AlphaIteratorType sx = sy;
                AlphaIteratorType dx = dy;
                for(; sx.x != send.x; ++sx.x, ++dx.x) {
                    if (*sx == GetMaxAlpha<AlphaPixelType>()
                            && *dx == GetMaxAlpha<AlphaPixelType>()) {
                        overlapFound = true;
                        break;
                    }
                }
                if (overlapFound) break;
            }

            if (!overlapFound) {
                // Copy src and srcA into image and imageA.

                if (Verbose > 0) {
                    cout << info->getFileName() << " ";
                    cout.flush();
                }

                vigra::Diff2D srcPos = info->getPosition();
                vigra::copyImageIf(srcImageRange(src),
                        maskImage(srcA),
                        destIter(image.upperLeft() - inputUnion.getUL() + srcPos));
                vigra::copyImageIf(srcImageRange(srcA),
                        maskImage(srcA),
                        destIter(imageA.upperLeft() - inputUnion.getUL() + srcPos));

                // Remove info from list later.
                toBeRemoved.push_back(i);
            }

        }
        
        // Erase the ImageImportInfos we used.
        list<list<vigra::ImageImportInfo*>::iterator>::iterator r;
        for (r = toBeRemoved.begin(); r != toBeRemoved.end(); r++) {
            imageInfoList.erase(*r);
        }
    }

    if (Verbose > 0) cout << endl;

    // Calculate bounding box of image.
    vigra::FindBoundingRectangle unionRect;
    vigra::inspectImageIf(srcIterRange(vigra::Diff2D(), vigra::Diff2D() + image.size()),
            srcImage(imageA), unionRect);
    bb.setCorners(unionRect.upperLeft, unionRect.lowerRight);
    if (Verbose > 0) {
        cout << "Combined union bounding box: ("
             << unionRect.upperLeft.x
             << ", "
             << unionRect.upperLeft.y
             << ") -> ("
             << unionRect.lowerRight.x
             << ", "
             << unionRect.lowerRight.y
             << ")" << endl;
    }

    // Dump image+imageA to temp file.
    // Return ImageImportInfo for that file.
    char tmpFilename[] = ".enblend_assemble_XXXXXX";
    int tmpFD = mkstemp(tmpFilename);
    vigra::ImageExportInfo outputImageInfo(tmpFilename);
    outputImageInfo.setFileType("TIFF");
    vigra::exportImageAlpha(srcImageRange(image),
            srcImage(imageA),
            outputImageInfo);
    close(tmpFD);
    
    return new vigra::ImageImportInfo(tmpFilename);
};

#endif /* __ASSEMBLE_H__ */
