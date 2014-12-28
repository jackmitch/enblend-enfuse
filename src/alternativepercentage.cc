/*
 * Copyright (C) 2014 Christoph L. Spiel
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


#include <sstream>

#include "alternativepercentage.h"


AlternativePercentage::AlternativePercentage(double a_value, bool is_percentage) :
    value_(a_value), is_percentage_(is_percentage)
{}


void
AlternativePercentage::set_value(double a_value)
{
    value_ = a_value;
}


void
AlternativePercentage::set_percentage(bool is_percentage)
{
    is_percentage_ = is_percentage;
}


std::string
AlternativePercentage::str() const
{
    std::ostringstream oss;

    oss << value_;
    if (is_percentage_)
    {
        oss << "%";
    }

    return oss.str();
}


CompactifiedAlternativePercentage::CompactifiedAlternativePercentage(double a_value, bool is_percentage) :
    AlternativePercentage(a_value, is_percentage)
{}
