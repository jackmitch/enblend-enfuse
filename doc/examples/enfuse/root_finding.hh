// Copyright (C) 2014-2017 Christoph L. Spiel
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


#ifndef ROOT_FINDING_INCLUDED
#define ROOT_FINDING_INCLUDED


#include <cassert>
#include <cmath>                  // std::abs(), std::sqrt()
#include <functional>             // std::unary_function<>
#include <limits>                 // std::numeric_limits<>


namespace root_finding
{
    // Rudimentary implementation of Dekker's algorithm

    template <class function>
    static double
    dekker(double x0, double x1, function f,
           double delta = std::sqrt(std::numeric_limits<double>::epsilon()),
           double epsilon = std::numeric_limits<double>::epsilon())
    {
        assert(x0 <= x1 && "degenerate interval");
        assert(f(x0) * f(x1) < 0.0 && "root not bracketed");

        while (x1 - x0 >= delta)
        {
            const double f0 = f(x0);
            const double f1 = f(x1);

            const double root = x1 - f1 * (x1 - x0) / (f1 - f0);
            const double f_root = f(root);

            if (std::fabs(f_root) <= epsilon)
            {
                return root;
            }

            const double middle = (x0 + x1) / 2.0;
            const double f_middle = f(middle);

            if (std::fabs(f_root) <= std::abs(f_middle))
            {
                if (f_root * f0 < 0.0)
                {
                    x1 = root;
                }
                else
                {
                    x0 = root;
                }
            }
            else
            {
                if (f_middle * f0 < 0.0)
                {
                    x1 = middle;
                }
                else
                {
                    x0 = middle;
                }
            }
        }

        return (x0 + x1) / 2.0;
    }
} // namespace root_finding


#endif // ROOT_FINDING_INCLUDED
