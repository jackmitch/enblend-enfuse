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
#ifndef __ENBLEND_H__
#define __ENBLEND_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>

#include "common.h"
#include "vigra/impex.hxx"

using namespace std;

template <typename ImageType, typename AlphaType, typename MaskType, typename PyramidType>
void enblend(list<vigra::ImageImportInfo*> &imageInfoList,
        vigra::ImageExportInfo &outputImageInfo,
        EnblendROI &inputUnion) {

    typedef typename ImageType::value_type image_value_type;
    typedef typename AlphaType::value_type alpha_value_type;
    typedef typename MaskType::value_type mask_value_type;
    typedef typename PyramidType::value_type pyramid_value_type;

    cout << "sizeof(image_value_type) = " << sizeof(image_value_type) << endl;
    cout << "sizeof(alpha_value_type) = " << sizeof(alpha_value_type) << endl;
    cout << "sizeof(mask_value_type) = " << sizeof(mask_value_type) << endl;
    cout << "sizeof(pyramid_value_type) = " << sizeof(pyramid_value_type) << endl;

    cout << "max_alpha=" << (int)GetMaxAlpha<alpha_value_type>() << endl;

    vigra::ImageImportInfo *whiteImageInfo =
            assemble<ImageType, AlphaType>(imageInfoList, inputUnion);

    delete whiteImageInfo;
}

#endif /* __ENBLEND_H__ */
