/*
 * Copyright (C) 2015-2017 Christoph Spiel
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
#ifndef PARAMETER_H_INCLUDED
#define PARAMETER_H_INCLUDED


#include <stdexcept>
#include <string>


namespace parameter
{
    // Identifier: [A-Za-z][A-Za-z0-9_-]*
    bool is_valid_identifier(const std::string& an_identifier);

    void insert(const std::string& a_key, const std::string& a_value);
    void erase(const std::string& a_key);
    void erase_all();


    struct not_found : public std::runtime_error
    {
        explicit not_found(const std::string& a_message) : std::runtime_error(a_message) {}
    };


    struct conversion_error : public std::runtime_error
    {
        explicit conversion_error(const std::string& a_message) : std::runtime_error(a_message) {}
    };


    // NOTES
    //
    // * The access of parameters through parameter::as_* is
    //   reasonably fast.  For time-critical parts of the code, the
    //   parameter's value can always be copied into a local variable.
    //
    // * The map from parameter keys to values is meant to be constant
    //   after the command line was parsed, i.e. neither the map
    //   itself, nor one of its values should be modified.  Following
    //   this convention makes all parameter::as_* functions thread
    //   safe.
    //
    // Some examples how to use parameters:
    //
    // (1) Check whether a parameter has been set.
    //         if (parameter::exists("foobar")) {...}
    //         else {...}
    //
    // (2) Use a parameter that is known to exist.
    //         std::string s = parameter::as_string("foobar");
    //         int i = parameter::as_integer("foobar");
    //         unsigned u = parameter::as_unsigned("foobar");
    //         double x = parameter::as_floating_point("foobar");
    //         bool b = parameter::as_boolean("foobar");
    //
    // (3) Substitute parameter value if it exists; otherwise go with
    //     the default.
    //         std::string s = parameter::as_string("foobar", "baz");
    //         int i = parameter::as_integer("foobar", 123);
    //         unsigned u = parameter::as_unsigned("foobar", 42U);
    //         double x = parameter::as_floating_point("foobar", 0.577215665);
    //         bool b = parameter::as_boolean("foobar", true);
    //
    // (4) React on parameter with a non-local change of control flow
    //         int i;
    //         try {i = parameter::as_integer("foobar");}
    //         catch (parameter::not_found&) {...}
    //
    // A parameter always can be retrieved as string with function
    // as_string().  All other as_* functions throw the exception
    // conversion_error, if the parameter's value cannot be
    // represented.


    bool exists(const std::string& a_key);

    std::string as_string(const std::string& a_key);
    std::string as_string(const std::string& a_key, const std::string& a_default_value);

    int as_integer(const std::string& a_key);
    int as_integer(const std::string& a_key, int a_default_value);

    unsigned as_unsigned(const std::string& a_key);
    unsigned as_unsigned(const std::string& a_key, unsigned a_default_value);

    double as_double(const std::string& a_key);
    double as_double(const std::string& a_key, double a_default_value);

    bool as_boolean(const std::string& a_key);
    bool as_boolean(const std::string& a_key, bool a_default_value);
} // namespace parameter


#endif // PARAMETER_H_INCLUDED

// Local Variables:
// mode: c++
// End:
