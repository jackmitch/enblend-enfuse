/*
 * Copyright (C) 2015-2017 Christoph L. Spiel
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
#ifndef OPENCL_EXPOSURE_WEIGHT_INCLUDED
#define OPENCL_EXPOSURE_WEIGHT_INCLUDED


#include <string>

#include "exposure_weight_base.h"


namespace opencl_exposure_weight
{
    bool is_opencl_file(const std::string& a_source_file_name);

    ExposureWeight* make_weight_function(const std::string& a_source_file_name,
                                         ExposureWeight::argument_const_iterator some_arguments_begin,
                                         ExposureWeight::argument_const_iterator some_arguments_end,
                                         double a_y_optimum, double a_width);
} // namespace opencl_exposure_weight


#endif // OPENCL_EXPOSURE_WEIGHT_INCLUDED


// Local Variables:
// mode: c++
// End:
