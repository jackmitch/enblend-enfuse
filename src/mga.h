/*
 * Copyright (C) 2004-2007 Andrew Mihal
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
#ifndef __MGA_H__
#define __MGA_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "common.h"

#include "vigra/colorconversions.hxx"

using std::cerr;
using std::endl;

using vigra::NumericTraits;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

namespace enblend {

template <typename InputType, typename ResultType>
class MultiGrayscaleAccessor
{
public:
    typedef ResultType value_type;

    MultiGrayscaleAccessor(const char* accessorName) {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        initializeTypeSpecific(srcIsScalar());
        initialize(accessorName);
    }

    template <class Iterator>
    ResultType operator()(const Iterator& i) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(i, srcIsScalar());
    }

    template <class Iterator, class Difference>
    ResultType operator()(const Iterator& i, Difference d) const {
        typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
        return f(i, d, srcIsScalar());
    }

private:
#define CHANNEL_MIXER "channel-mixer"
    void initialize(const char* accessorName) {
        if (accessorName == NULL) kind = AVERAGE;
        else
        {
            char dummy;
            double red, green, blue;
            std::string name(accessorName);

            std::transform(name.begin(), name.end(), name.begin(), tolower);
            if (name == "average") kind = AVERAGE;
            else if (name == "l-star") kind = LSTAR;
            else if (name == "lightness") kind = LIGHTNESS;
            else if (name == "value") kind = VALUE;
            else if (name == "luminance") kind = LUMINANCE;
            else if (name == CHANNEL_MIXER)
            {
                cerr << "enfuse: \"" CHANNEL_MIXER "\" is a grayscale projector requiring" << endl
                     << "enfuse:      arguments like e.g. \"channel-mixer:0.30:0.59:0.11\"." << endl;
                exit(1);
            }
            else if (sscanf(name.c_str(),
                            CHANNEL_MIXER "%["
                            OPTION_DELIMITERS "]%lf%["
                            OPTION_DELIMITERS "]%lf%["
                            OPTION_DELIMITERS "]%lf",
                            &dummy, &red, &dummy, &green, &dummy, &blue) == 6)
            {
                const double sum = red + green + blue;

                redWeight = red / sum;
                greenWeight = green / sum;
                blueWeight = blue / sum;
                kind = MIXER;
            }
            else
            {
                cerr << "enfuse: unknown grayscale projector \"" << name << "\"." << endl;
                exit(1);
            }
        }
    }

    void initializeTypeSpecific(VigraTrueType) {}

    void initializeTypeSpecific(VigraFalseType) {
        typedef typename InputType::value_type ValueType;
        labfun = vigra::RGB2LabFunctor<double>(NumericTraits<ValueType>::max());
    }

    template <class Iterator>
    ResultType f(const Iterator& i, VigraFalseType) const {
        typedef typename InputType::value_type ValueType;
        switch (kind)
        {
        case AVERAGE:
            return NumericTraits<ResultType>::fromRealPromote
                ((NumericTraits<ValueType>::toRealPromote((*i).red()) +
                  NumericTraits<ValueType>::toRealPromote((*i).green()) +
                  NumericTraits<ValueType>::toRealPromote((*i).blue())) /
                 3.0);
        case LSTAR:
        {
            typedef typename vigra::RGB2LabFunctor<double>::result_type LABResultType;
            const LABResultType y = labfun.operator()(*i) / 100.0;
            return NumericTraits<ResultType>::fromRealPromote(NumericTraits<ValueType>::max() * y[0]);
        }
        case LIGHTNESS:
            return NumericTraits<ResultType>::fromRealPromote
                ((std::min((*i).red(), std::min((*i).green(), (*i).blue())) +
                  std::max((*i).red(), std::max((*i).green(), (*i).blue()))) /
                 2.0);
        case VALUE:
            return std::max((*i).red(), std::max((*i).green(), (*i).blue()));
        case LUMINANCE:
            return NumericTraits<ResultType>::fromRealPromote((*i).luminance());
        case MIXER:
            return NumericTraits<ResultType>::fromRealPromote
                (redWeight * NumericTraits<ValueType>::toRealPromote((*i).red()) +
                 greenWeight * NumericTraits<ValueType>::toRealPromote((*i).green()) +
                 blueWeight * NumericTraits<ValueType>::toRealPromote((*i).blue()));
        }

        // never reached
        return ResultType();
    }

    template <class Iterator, class Difference>
    ResultType f(const Iterator& i, Difference d, VigraFalseType) const {
        typedef typename InputType::value_type ValueType;
        switch (kind)
        {
        case AVERAGE:
            return NumericTraits<ResultType>::fromRealPromote
                ((NumericTraits<ValueType>::toRealPromote(i[d].red()) +
                  NumericTraits<ValueType>::toRealPromote(i[d].green()) +
                  NumericTraits<ValueType>::toRealPromote(i[d].blue())) /
                 3.0);
        case LSTAR:
        {
            typedef typename vigra::RGB2LabFunctor<double>::result_type LABResultType;
            const LABResultType y = labfun.operator()(i[d]) / 100.0;
            return NumericTraits<ResultType>::fromRealPromote(NumericTraits<ValueType>::max() * y[0]);
        }
        case LIGHTNESS:
            return NumericTraits<ResultType>::fromRealPromote
                ((std::min(i[d].red(), std::min(i[d].green(), i[d].blue())) +
                  std::max(i[d].red(), std::max(i[d].green(), i[d].blue()))) /
                 2.0);
        case VALUE:
            return std::max(i[d].red(), std::max(i[d].green(), i[d].blue()));
        case LUMINANCE:
            return NumericTraits<ResultType>::fromRealPromote(i[d].luminance());
        case MIXER:
            return NumericTraits<ResultType>::fromRealPromote
                (redWeight * NumericTraits<ValueType>::toRealPromote(i[d].red()) +
                 greenWeight * NumericTraits<ValueType>::toRealPromote(i[d].green()) +
                 blueWeight * NumericTraits<ValueType>::toRealPromote(i[d].blue()));
        }

        // never reached
        return ResultType();
    }

    template <class Iterator>
    ResultType f(const Iterator& i, VigraTrueType) const {return *i;}

    template <class Iterator, class Difference>
    ResultType f(const Iterator& i, Difference d, VigraTrueType) const {return i[d];}

    enum AccKind {AVERAGE, LSTAR, LIGHTNESS, VALUE, LUMINANCE, MIXER} kind;
    double redWeight, greenWeight, blueWeight;
    vigra::RGB2LabFunctor<double> labfun;
};

} // namespace enblend

#endif /* __ENMIX_H__ */
