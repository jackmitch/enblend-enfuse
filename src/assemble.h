/*
 * Copyright (C) 2004-2005 Andrew Mihal
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
#include <lcms.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#include "common.h"
#include "vigra/copyimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"
#include "vigra_ext/FunctorAccessor.h"
#include "vigra_ext/impexalpha.hxx"

using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::copyImageIf;
using vigra::Diff2D;
using vigra::exportImageAlpha;
using vigra::FindBoundingRectangle;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::inspectImageIf;
using vigra::NumericTraits;
using vigra::Threshold;
using vigra::transformImage;

using vigra_ext::ReadFunctorAccessor;
using vigra_ext::WriteFunctorAccessor;

namespace enblend {

/** Write output images.
 */
template <typename ImageType, typename ImageComponentType, typename AlphaType>
void checkpoint(pair<ImageType*, AlphaType*> &p, ImageExportInfo &outputImageInfo) {

    typedef typename ImageType::Accessor ImageAccessor;
    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename AlphaType::Accessor AlphaAccessor;
    typedef typename AlphaType::PixelType AlphaPixelType;

    typedef ReadFunctorAccessor<
            Threshold<AlphaPixelType, ImageComponentType>, AlphaAccessor>
            ThresholdingAccessor;
 
    ThresholdingAccessor ata(
            Threshold<AlphaPixelType, ImageComponentType>(
                    NumericTraits<AlphaPixelType>::zero(),
                    NumericTraits<AlphaPixelType>::zero(),
                    NumericTraits<ImageComponentType>::max(),
                    NumericTraits<ImageComponentType>::zero()
            ),
            (p.second)->accessor());

    exportImageAlpha(srcImageRange(*(p.first)),
                     srcIter((p.second)->upperLeft(), ata),
                     outputImageInfo);

};

template <typename ImageComponentType,
          typename DestIterator, typename DestAccessor,
          typename AlphaIterator, typename AlphaAccessor>
void import(const ImageImportInfo &info, const pair<DestIterator, DestAccessor> &image, const pair<AlphaIterator, AlphaAccessor> &alpha) {

    typedef vigra::RGBValue<ImageComponentType> ImagePixelType;
    typedef typename AlphaIterator::PixelType AlphaPixelType;

    // Use a thresholding accessor to write to the alpha image.
    typedef WriteFunctorAccessor<
            Threshold<ImageComponentType, AlphaPixelType>, AlphaAccessor>
            ThresholdingAccessor;

    // Threshold the alpha mask so that all pixels are either contributing
    // or not contributing.
    ThresholdingAccessor ata(
            Threshold<ImageComponentType, AlphaPixelType>(
                    //NumericTraits<AlphaPixelType>::max() / 2,
                    NumericTraits<ImageComponentType>::max(),
                    NumericTraits<ImageComponentType>::max(),
                    NumericTraits<AlphaPixelType>::zero(),
                    NumericTraits<AlphaPixelType>::max()
            ),
            alpha.second);

    importImageAlpha(info, image, destIter(alpha.first, ata));


    if (UseCIECAM) {
        cmsHPROFILE sourceProfile;
        cmsHPROFILE destProfile = cmsCreateXYZProfile();

        ImageImportInfo::ICCProfile iccdata = info.getICCProfile();

        if (iccdata.empty()) {
            cout << "enblend: input image \"" << info.getFileName() << "\" does not have an ICC profile. Assuming sRGB." << endl;
            sourceProfile = cmsCreate_sRGBProfile();
        }
        else {
            sourceProfile = cmsOpenProfileFromMem(iccdata.data(), iccdata.size());
        }

        if (sourceProfile == NULL) {
            cerr << endl << "enblend: could not read ICC profile data from file \"" << info.getFileName() << "\"."
                 << endl << endl;
            exit(1);
        }

        cmsHTRANSFORM transform = cmsCreateTransform(sourceProfile, TYPE_RGB_DBL,
                                                     destProfile, TYPE_XYZ_DBL,
                                                     INTENT_PERCEPTUAL, 0);

        if (transform == NULL) {
            cerr << endl << "enblend: error building color transform from \"" << cmsTakeProductName(sourceProfile)
                 << "\" to XYZ." << endl << endl;
            exit(1);
        }

        cmsViewingConditions conditions;
        if (!cmsTakeMediaWhitePoint(&(conditions.whitePoint), sourceProfile)) {
            cerr << endl << "enblend: could not get media white point from \"" << cmsTakeProductName(sourceProfile)
                 << "\"." << endl << endl;
            exit(1);
        }
        conditions.whitePoint.X *= 100.0;
        conditions.whitePoint.Y *= 100.0;
        conditions.whitePoint.Z *= 100.0;
        conditions.Yb = 20.0;
        conditions.La = 20.0;
        conditions.surround = AVG_SURROUND;
        conditions.D_value = 1.0;
        cout << "source profile media white X=" << conditions.whitePoint.X << " Y=" << conditions.whitePoint.Y << " Z=" << conditions.whitePoint.Z << endl;

        LCMSHANDLE ciecamTransform = cmsCIECAM02Init(&conditions);
        if (!ciecamTransform) {
            cerr << endl << "enblend: error initializing CIECAM02 transform." << endl << endl;
            exit(1);
        }

        DestIterator sy = image.first;
        AlphaIterator my = alpha.first;
        for (int y = 0; y < info.height(); y++, ++sy.y, ++my.y) {
            DestIterator sx = sy;
            AlphaIterator mx = my;
            for (int x = 0; x < info.width(); x++, ++sx.x, ++mx.x) {
                if (alpha.second(mx)) {
                    double inputBuffer[3];
                    double outputBuffer[3];
                    double scaleFactor = 1.0 / NumericTraits<ImageComponentType>::max();
                    ImagePixelType pixel = image.second(sx);
                    inputBuffer[0] = scaleFactor * NumericTraits<ImageComponentType>::toRealPromote(pixel.red());
                    inputBuffer[1] = scaleFactor * NumericTraits<ImageComponentType>::toRealPromote(pixel.green());
                    inputBuffer[2] = scaleFactor * NumericTraits<ImageComponentType>::toRealPromote(pixel.blue());

                    cmsDoTransform(transform, inputBuffer, outputBuffer, 1);

                    cmsCIEXYZ xyz;
                    xyz.X = outputBuffer[0] * 100.0;
                    xyz.Y = outputBuffer[1] * 100.0;
                    xyz.Z = outputBuffer[2] * 100.0;

                    cmsJCh jch;
                    cmsCIECAM02Forward(ciecamTransform, &xyz, &jch);

                    cout << "pixel rgb value=(" << (int)pixel.red() << ", " << (int)pixel.green() << ", " << (int)pixel.blue() << ")" << endl;
                    cout << "pixel XYZ value=(" << xyz.X << ", " << xyz.Y << ", " << xyz.Z << ")" << endl;
                    cout << "pixel JCh value=(" << jch.J << ", " << jch.C << ", " << jch.h << ")" << endl;
                    //RGBValue<ImageComponentType> result;
                    //result.setRed(
                    //image.second.set(

                    goto DONE;
                }
            }
        }

DONE:
        cmsCIECAM02Done(ciecamTransform);
        cmsDeleteTransform(transform);
        cmsCloseProfile(sourceProfile);
        cmsCloseProfile(destProfile);
    }

};

/** Find images that don't overlap and assemble them into one image.
 *  Uses a greedy heuristic.
 *  Removes used images from given list of ImageImportInfos.
 *  Returns an ImageImportInfo for the temporary file.
 *  memory xsection = 2 * (ImageType*inputUnion + AlphaType*inputUnion)
 */
template <typename ImageType, typename ImageComponentType, typename AlphaType>
pair<ImageType*, AlphaType*> assemble(list<ImageImportInfo*> &imageInfoList,
        EnblendROI &inputUnion,
        EnblendROI &bb) {

    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename ImageType::Accessor ImageAccessor;
    typedef typename AlphaType::PixelType AlphaPixelType;
    typedef typename AlphaType::traverser AlphaIteratorType;
    typedef typename AlphaType::Accessor AlphaAccessor;

    // No more images to assemble?
    if (imageInfoList.empty()) return pair<ImageType*, AlphaType*>(NULL, NULL);

    // Create an image to assemble input images into.
    ImageType *image = new ImageType(inputUnion.size());
    AlphaType *imageA = new AlphaType(inputUnion.size());

    if (Verbose > VERBOSE_ASSEMBLE_MESSAGES) {
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

    //// Load the first image into the destination.
    //// Use a thresholding accessor to write to the alpha image.
    //typedef WriteFunctorAccessor<
    //        Threshold<ImageComponentType, AlphaPixelType>, AlphaAccessor>
    //        ThresholdingAccessor;
    //// Threshold the alpha mask so that all pixels are either contributing
    //// or not contributing.
    //ThresholdingAccessor imageATA(
    //        Threshold<ImageComponentType, AlphaPixelType>(
    //                //NumericTraits<AlphaPixelType>::max() / 2,
    //                NumericTraits<ImageComponentType>::max(),
    //                NumericTraits<ImageComponentType>::max(),
    //                NumericTraits<AlphaPixelType>::zero(),
    //                NumericTraits<AlphaPixelType>::max()
    //        ),
    //        imageA->accessor());

    Diff2D imagePos = imageInfoList.front()->getPosition();
    //importImageAlpha(*imageInfoList.front(),
    //        destIter(image->upperLeft() + imagePos - inputUnion.getUL()),
    //        destIter(imageA->upperLeft() + imagePos - inputUnion.getUL(), imageATA));
    import<ImageComponentType>(*imageInfoList.front(),
            destIter(image->upperLeft() + imagePos - inputUnion.getUL()),
            destIter(imageA->upperLeft() + imagePos - inputUnion.getUL()));
    imageInfoList.erase(imageInfoList.begin());

    if (!OneAtATime) {
        // Attempt to assemble additional non-overlapping images.

        // List of ImageImportInfos we decide to assemble.
        list<list<ImageImportInfo*>::iterator> toBeRemoved;

        list<ImageImportInfo*>::iterator i;
        for (i = imageInfoList.begin(); i != imageInfoList.end(); i++) {
            ImageImportInfo *info = *i;

            // Load the next image.
            ImageType *src = new ImageType(info->size());
            AlphaType *srcA = new AlphaType(info->size());

            // Use a thresholding accessor to write to the alpha image.
            //ThresholdingAccessor srcATA(
            //        Threshold<ImageComponentType, AlphaPixelType>(
            //                //NumericTraits<AlphaPixelType>::max() / 2,
            //                NumericTraits<ImageComponentType>::max(),
            //                NumericTraits<ImageComponentType>::max(),
            //                NumericTraits<AlphaPixelType>::zero(),
            //                NumericTraits<AlphaPixelType>::max()
            //        ),
            //        srcA->accessor());
            //importImageAlpha(*info, destImage(*src), destImage(*srcA, srcATA));
            import<ImageComponentType>(*info, destImage(*src), destImage(*srcA));

            // Check for overlap.
            bool overlapFound = false;
            AlphaIteratorType dy =
                    imageA->upperLeft() - inputUnion.getUL() + info->getPosition();
            AlphaAccessor da = imageA->accessor();
            AlphaIteratorType sy = srcA->upperLeft();
            AlphaIteratorType send = srcA->lowerRight();
            AlphaAccessor sa = srcA->accessor();
            for(; sy.y != send.y; ++sy.y, ++dy.y) {
                AlphaIteratorType sx = sy;
                AlphaIteratorType dx = dy;
                for(; sx.x != send.x; ++sx.x, ++dx.x) {
                    if (sa(sx) && da(dx)) {
                        overlapFound = true;
                        break;
                    }
                }
                if (overlapFound) break;
            }

            if (!overlapFound) {
                // Copy src and srcA into image and imageA.

                if (Verbose > VERBOSE_ASSEMBLE_MESSAGES) {
                    cout << " " << info->getFileName();
                    cout.flush();
                }

                Diff2D srcPos = info->getPosition();
                copyImageIf(srcImageRange(*src),
                        maskImage(*srcA),
                        destIter(image->upperLeft() - inputUnion.getUL() + srcPos));
                copyImageIf(srcImageRange(*srcA),
                        maskImage(*srcA),
                        destIter(imageA->upperLeft() - inputUnion.getUL() + srcPos));

                // Remove info from list later.
                toBeRemoved.push_back(i);
            }

            delete src;
            delete srcA;
        }
        
        // Erase the ImageImportInfos we used.
        list<list<ImageImportInfo*>::iterator>::iterator r;
        for (r = toBeRemoved.begin(); r != toBeRemoved.end(); r++) {
            imageInfoList.erase(*r);
        }
    }

    if (Verbose > VERBOSE_ASSEMBLE_MESSAGES && !OneAtATime) cout << endl;

    // Calculate bounding box of image.
    FindBoundingRectangle unionRect;
    inspectImageIf(srcIterRange(Diff2D(), Diff2D() + image->size()),
            srcImage(*imageA), unionRect);
    bb.setCorners(unionRect.upperLeft, unionRect.lowerRight);

    if (Verbose > VERBOSE_ABB_MESSAGES) {
        cout << "assembled images bounding box: ("
             << unionRect.upperLeft.x
             << ", "
             << unionRect.upperLeft.y
             << ") -> ("
             << unionRect.lowerRight.x
             << ", "
             << unionRect.lowerRight.y
             << ")" << endl;
    }

    return pair<ImageType*, AlphaType*>(image, imageA);

};

} // namespace enblend

#endif /* __ASSEMBLE_H__ */
