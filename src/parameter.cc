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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cctype>       // isalnum(), isalpha()
#include <cerrno>       // errno
#include <cstdlib>      // strtod(), strtol(), strtoul()

#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#else
#include <map>
#endif

#include "parameter.h"

namespace parameter
{
    bool
    is_valid_identifier(const std::string& an_identifier)
    {
        if (an_identifier.size() == 0)
        {
            return false;
        }
        else if (!isalpha(an_identifier[0]))
        {
            return false;
        }
        else
        {
            for (std::string::const_iterator x = an_identifier.begin(); x != an_identifier.end(); ++x)
            {
                if (!(isalnum(*x) || *x == '_' || *x == '-'))
                {
                    return false;
                }
            }
        }

        return true;
    }


    class ParameterValue
    {
    public:
        ParameterValue();
        explicit ParameterValue(const std::string& a_string);
        ParameterValue(const ParameterValue& another_parameter_value);
        ParameterValue& operator=(const ParameterValue& another_parameter_value);
        virtual ~ParameterValue() {release_memory();}

        std::string as_string() const {return value_as_string_;}
        const char* as_c_string() const {return value_as_string_.c_str();}

        int as_integer() const;
        unsigned as_unsigned() const;
        double as_double() const;
        bool as_boolean() const;

    private:
        void initialize();
        void initialize_integer();
        void initialize_unsigned_integer();
        void initialize_floating_point();
        void initialize_boolean();

        void copy_cached_values(const ParameterValue& a_parameter_value);

        void release_memory();

        std::string value_as_string_;
        int* integer_;
        unsigned* unsigned_integer_;
        double* floating_point_;
        bool* boolean_;
    }; // end class ParameterValue


    ParameterValue::ParameterValue() :
        value_as_string_(std::string()),
        integer_(nullptr), unsigned_integer_(nullptr), floating_point_(nullptr), boolean_(nullptr)
    {
        initialize();
    }


    ParameterValue::ParameterValue(const std::string& a_string) :
        value_as_string_(a_string),
        integer_(nullptr), unsigned_integer_(nullptr), floating_point_(nullptr), boolean_(nullptr)
    {
        initialize();
    }


    ParameterValue::ParameterValue(const ParameterValue& another_parameter_value) :
        value_as_string_(another_parameter_value.value_as_string_),
        integer_(nullptr), unsigned_integer_(nullptr), floating_point_(nullptr), boolean_(nullptr)
    {
        copy_cached_values(another_parameter_value);
    }


    ParameterValue&
    ParameterValue::operator=(const ParameterValue& another_parameter_value)
    {
        if (this != &another_parameter_value)
        {
            value_as_string_ = another_parameter_value.value_as_string_;
            release_memory();
            copy_cached_values(another_parameter_value);
        }

        return *this;
    }


    int
    ParameterValue::as_integer() const
    {
        if (integer_)
        {
            return *integer_;
        }
        else
        {
            throw conversion_error("cannot convert \"" + value_as_string_ + "\" to an integer");
        }
    }


    unsigned
    ParameterValue::as_unsigned() const
    {
        if (unsigned_integer_)
        {
            return *unsigned_integer_;
        }
        else
        {
            throw conversion_error("cannot convert \"" + value_as_string_ + "\" to an unsigned integer");
        }
    }


    double
    ParameterValue::as_double() const
    {
        if (floating_point_)
        {
            return *floating_point_;
        }
        else
        {
            throw conversion_error("cannot convert \"" + value_as_string_ + "\" to a floating-point number");
        }
    }


    bool
    ParameterValue::as_boolean() const
    {
        if (boolean_)
        {
            return *boolean_;
        }
        else
        {
            throw conversion_error("cannot convert \"" + value_as_string_ + "\" to a boolean");
        }
    }


    void
    ParameterValue::initialize()
    {
        initialize_integer();
        initialize_unsigned_integer();
        initialize_floating_point();
        initialize_boolean();
    }


    void
    ParameterValue::initialize_integer()
    {
        char* end;
        errno = 0;
        const int i = static_cast<int>(strtol(value_as_string_.c_str(), &end, 10));
        if (errno == 0 && *end == 0)
        {
            integer_ = new int;
            *integer_ = i;
        }
    }


    void
    ParameterValue::initialize_unsigned_integer()
    {
        char* end;
        errno = 0;
        const unsigned u = static_cast<unsigned>(strtoul(value_as_string_.c_str(), &end, 10));
        if (errno == 0 && *end == 0)
        {
            unsigned_integer_ = new unsigned;
            *unsigned_integer_ = u;
        }
    }


    void
    ParameterValue::initialize_floating_point()
    {
        char* end;
        errno = 0;
        const double x = strtod(value_as_string_.c_str(), &end);
        if (errno == 0 && *end == 0)
        {
            floating_point_ = new double;
            *floating_point_ = x;
        }
    }


    void
    ParameterValue::initialize_boolean()
    {
        std::string s(value_as_string_);

        bool b;
        if (s.empty() || s == "0" || s == "f" || s == "false")
        {
            b = false;
        }
        else if (s == "1" || s == "t" || s == "true")
        {
            b = true;
        }
        else
        {
            char* end;
            errno = 0;
            b = strtol(value_as_string_.c_str(), &end, 10) != 0;
            if (errno != 0 || *end != 0)
            {
                return;
            }
        }

        boolean_ = new bool;
        *boolean_ = b;
    }


    void
    ParameterValue::copy_cached_values(const ParameterValue& a_parameter_value)
    {
        if (a_parameter_value.integer_)
        {
            integer_ = new int;
            *integer_ = *a_parameter_value.integer_;
        }

        if (a_parameter_value.unsigned_integer_)
        {
            unsigned_integer_ = new unsigned;
            *unsigned_integer_ = *a_parameter_value.unsigned_integer_;
        }

        if (a_parameter_value.floating_point_)
        {
            floating_point_ = new double;
            *floating_point_ = *a_parameter_value.floating_point_;
        }

        if (a_parameter_value.boolean_)
        {
            boolean_ = new bool;
            *boolean_ = *a_parameter_value.boolean_;
        }
    }


    void
    ParameterValue::release_memory()
    {
        delete integer_;
        delete unsigned_integer_;
        delete floating_point_;
        delete boolean_;
    }


    //
    // Parameter Map
    //

#ifdef HAVE_UNORDERED_MAP
    typedef std::unordered_map<std::string, ParameterValue> parameter_map_t;
#else
    typedef std::map<std::string, ParameterValue> parameter_map_t;
#endif


    static parameter_map_t map;


    // Parameter Map - Modification

    void
    insert(const std::string& a_key, const std::string& a_value)
    {
        map.insert(parameter_map_t::value_type(a_key, ParameterValue(a_value)));
    }


    void
    erase(const std::string& a_key)
    {
        map.erase(a_key);
    }


    void
    erase_all()
    {
        map.clear();
    }


    // Parameter Map - Query

    bool
    exists(const std::string& a_key)
    {
        return map.find(a_key) != map.end();
    }


    std::string
    as_string(const std::string& a_key)
    {
        parameter_map_t::const_iterator x = map.find(a_key);
        if (x == map.end())
        {
            throw not_found(a_key);
        }
        else
        {
            return x->second.as_string();
        }
    }


    std::string
    as_string(const std::string& a_key, const std::string& a_default_value)
    {
        parameter_map_t::const_iterator x = map.find(a_key);
        if (x == map.end())
        {
            return a_default_value;
        }
        else
        {
            return x->second.as_string();
        }
    }


    int
    as_integer(const std::string& a_key)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            throw not_found(a_key);
        }
        else
        {
            return x->second.as_integer();
        }
    }


    int
    as_integer(const std::string& a_key, int a_default_value)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            return a_default_value;
        }
        else
        {
            return x->second.as_integer();
        }
    }


    unsigned
    as_unsigned(const std::string& a_key)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            throw not_found(a_key);
        }
        else
        {
            return x->second.as_unsigned();
        }
    }


    unsigned
    as_unsigned(const std::string& a_key, unsigned a_default_value)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            return a_default_value;
        }
        else
        {
            return x->second.as_unsigned();
        }
    }


    double
    as_double(const std::string& a_key)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            throw not_found(a_key);
        }
        else
        {
            return x->second.as_double();
        }
    }


    double
    as_double(const std::string& a_key, double a_default_value)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            return a_default_value;
        }
        else
        {
            return x->second.as_double();
        }
    }


    bool
    as_boolean(const std::string& a_key)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            throw not_found(a_key);
        }
        else
        {
            return x->second.as_boolean();
        }
    }


    bool
    as_boolean(const std::string& a_key, bool a_default_value)
    {
        parameter_map_t::iterator x = map.find(a_key);
        if (x == map.end())
        {
            return a_default_value;
        }
        else
        {
            return x->second.as_boolean();
        }
    }
} // end namespace parameter
