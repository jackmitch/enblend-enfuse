// Copyright (C) 2016 Christoph L. Spiel
//
// This file is part of Enblend.
//
// Enblend is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Enblend is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Enblend; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef OPTIONAL_TRANSITIONAL_HEADER_INCLUDED
#define OPTIONAL_TRANSITIONAL_HEADER_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#if defined(HAVE_OPTIONAL) || defined(HAVE_OPTIONAL_HPP)

#if defined(HAVE_OPTIONAL)
#include <optional>
#elif defined(HAVE_OPTIONAL_HPP)
#include <optional.hpp>
#endif

#if __cplusplus < 201701L
namespace std
{
    template <typename t> using optional = ::std::experimental::optional<t>;
    using experimental::nullopt;
}
#endif

#elif defined(HAVE_BOOST_OPTIONAL_HPP)

#include <boost/optional.hpp>
namespace std
{
    template <typename t> using optional = ::boost::optional<t>;
    constexpr ::boost::none_t nullopt;
}

#else

#error "no viable, C++-17-compatible <optional> header file defined"

#endif


#endif // OPTIONAL_TRANSITIONAL_HEADER_INCLUDED


// Local Variables:
// mode: c++
// End:
