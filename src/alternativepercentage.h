/*
 * Copyright (C) 2015-2017  Christoph L. Spiel
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

#ifndef ALTERNATIVEPERCENTAGE_H_INCLUDED_
#define ALTERNATIVEPERCENTAGE_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vigra/numerictraits.hxx>


class AlternativePercentage
{
public:
    AlternativePercentage() = delete;
    AlternativePercentage(double a_value, bool is_percentage);
    virtual ~AlternativePercentage() {}

    virtual double value() const {return value_;}
    bool is_percentage() const {return is_percentage_;}

    virtual void set_value(double a_value);
    void set_percentage(bool is_percentage);

    virtual std::string str() const;

    template <class t>
    bool is_effective() const
    {
        return
            value_ > 0.0 &&
            ((is_percentage() && value() < 100.0) ||
             (!is_percentage() && value() < vigra::NumericTraits<t>::max()));
    }

    template <class t>
    t instantiate() const
    {
        const double max {static_cast<double>(vigra::NumericTraits<t>::max())};

        return is_percentage() ? value() * max / 100.0 : value();
    }

private:
    double value_;
    bool is_percentage_;
};


class CompactifiedAlternativePercentage : public AlternativePercentage
{
public:
    CompactifiedAlternativePercentage() = delete;
    CompactifiedAlternativePercentage(double a_value, bool is_percentage);

    template <class t>
    bool is_effective() const
    {
        return
            value() != 0.0 &&
            ((is_percentage() && std::abs(value()) != 100.0) ||
             (!is_percentage() && value() != vigra::NumericTraits<t>::max()));
    }

    template <class t>
    t instantiate() const
    {
        const double max {static_cast<double>(vigra::NumericTraits<t>::max())};

        if (is_percentage())
        {
            return (value() >= 0.0 ? value() : 100.0 + value()) * max / 100.0;
        }
        else
        {
            return value() >= 0.0 ? value() : max + value();
        }
    }
};


#endif // ALTERNATIVEPERCENTAGE_H_INCLUDED_


// Local Variables:
// mode: c++
// End:
