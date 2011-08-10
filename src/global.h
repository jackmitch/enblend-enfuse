/*
 * Copyright (C) 2009-2011 Dr. Christoph L. Spiel
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
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

// Here we define macros and types that we already need in the
// definitions of global variables.


#include <sstream>

#include <vigra/numerictraits.hxx>


// Defines to control how many -v flags are required for each type
// of message to be produced on stdout.
#define VERBOSE_ASSEMBLE_MESSAGES           1
#define VERBOSE_CHECKPOINTING_MESSAGES      1

#define VERBOSE_BLEND_MESSAGES              2
#define VERBOSE_GPU_MESSAGES                2
#define VERBOSE_MASK_MESSAGES               2
#define VERBOSE_NFT_MESSAGES                2
#define VERBOSE_PYRAMID_MESSAGES            2
#define VERBOSE_SIGNATURE_REPORTING         2
#define VERBOSE_VERSION_REPORTING           2


#define VERBOSE_COLOR_CONVERSION_MESSAGES   3
#define VERBOSE_LAYER_SELECTION             3
#define VERBOSE_RESPONSE_FILES              3

#define VERBOSE_ABB_MESSAGES                4
#define VERBOSE_IBB_MESSAGES                4
#define VERBOSE_INPUT_IMAGE_INFO_MESSAGES   4
#define VERBOSE_INPUT_UNION_SIZE_MESSAGES   4
#define VERBOSE_ROIBB_SIZE_MESSAGES         4
#define VERBOSE_UBB_MESSAGES                4

#define VERBOSE_CFI_MESSAGES                5
#define VERBOSE_GDA_MESSAGES                5

#define VERBOSE_MEMORY_ESTIMATION_MESSAGES  6


//< src::default-output-filename a.tif
#define DEFAULT_OUTPUT_FILENAME "a.tif"


// Safely retrieve the string associated with m_name from OpenGL.
#define GLGETSTRING(m_name) \
    (glGetString(m_name) == NULL ? \
     "<cannot retrieve " #m_name ">" : \
     (const char*) (glGetString(m_name)))


class AlternativePercentage
{
public:
    AlternativePercentage(double value, bool is_percentage) :
        value_(value), is_percentage_(is_percentage) {}

    double value() const {return value_;}
    double is_percentage() const {return is_percentage_;}

    void set_value(double value) {value_ = value;}
    void set_percentage(bool is_percentage) {is_percentage_ = is_percentage;}

    std::string str() const
    {
        std::ostringstream oss;
        oss << value_;
        if (is_percentage_)
        {
            oss << "%";
        }
        return oss.str();
    }

    template <class T>
    bool is_effective() const
    {
        return
            value_ > 0.0 &&
            ((is_percentage_ && value_ < 100.0) ||
             (!is_percentage_ && value_ < vigra::NumericTraits<T>::max()));
    }

    template <class T>
    T instantiate() const
    {
        return
            is_percentage_ ?
            value_ * static_cast<double>(vigra::NumericTraits<T>::max()) / 100.0 :
            value_;
    }

private:
    double value_;
    bool is_percentage_;
};


/** The different kinds of boundary conditions we can impose upon an
 *  image. */
typedef enum BoundaryKind
{
    UnknownWrapAround,          // unknown kind
    OpenBoundaries,             // contractible
    HorizontalStrip,            // contractible along 2nd axis
    VerticalStrip,              // contractible along 1st axis
    DoubleStrip                 // non-contractible
} boundary_t;

enum MainAlgo{
    NFT, GraphCut
};



//< src::default-tiff-resolution 300@dmn{dpi}
#define DEFAULT_TIFF_RESOLUTION 300.0f


struct TiffResolution {
    TiffResolution() : x(0.0f), y(0.0f) {}

    TiffResolution(float anXresolution, float aYresolution) :
        x(anXresolution), y(aYresolution) {}

    bool operator==(const TiffResolution& anOther) const {
        return this->x == anOther.x && this->y == anOther.y;
    }

    bool operator!=(const TiffResolution& anOther) const {
        return !operator==(anOther);
    }

    float x;
    float y;
};

#endif /* __GLOBAL_H__ */

// Local Variables:
// mode: c++
// End:
