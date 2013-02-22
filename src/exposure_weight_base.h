/*
 * Copyright (C) 2013 Dr. Christoph L. Spiel
 *
 * This file is part of Enblend.
 */
#ifndef EXPOSURE_WEIGHT_BASE_INCLUDED
#define EXPOSURE_WEIGHT_BASE_INCLUDED


#include <stdexcept>


class ExposureWeight
{
public:
    ExposureWeight() : y_optimum_(0.5), width_(0.25) {}
    ExposureWeight(double y_optimum, double width) : y_optimum_(y_optimum), width_(width)
    {
        check_invariant();
    }

    void initialize(double y_optimum, double width)
    {
        y_optimum_ = y_optimum;
        width_ = width;

        check_invariant();
    }

    double optimum() const {return y_optimum_;}
    double width() const {return width_;}

    virtual double normalize(double y) const {return (y - optimum()) / width();}
    virtual double weight(double) const = 0;

    virtual ~ExposureWeight() {}

private:
    void check_invariant() const
    {
        if (y_optimum_ < 0.0 || y_optimum_ > 1.0)
        {
            throw std::invalid_argument("y_optimum");
        }
        if (width_ <= 0.0)
        {
            throw std::invalid_argument("width");
        }
    }

    double y_optimum_;
    double width_;
};


#endif // EXPOSURE_WEIGHT_BASE_INCLUDED


// Local Variables:
// mode: c++
// End:
