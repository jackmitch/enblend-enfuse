/*
 * Copyright (C) 2013-2017 Dr. Christoph L. Spiel
 *
 * This file is part of Enblend.
 */
#ifndef EXPOSURE_WEIGHT_BASE_INCLUDED
#define EXPOSURE_WEIGHT_BASE_INCLUDED


#include <stdexcept>
#include <string>
#include <vector>
#include <iterator>


// Version number of the interface-to-user-weight functions.
#define EXPOSURE_WEIGHT_INTERFACE_VERSION 2


// The full width at half of the maximum of the Gauss-curve we use for
// exposure weighting is
//         FWHM = 2 * sqrt(2 * log(2))
// For compatibility in the sake of least surprise of the user, all
// other exposure weight functions rescale their native FWHM to match
// the Gauss-curve.
#define FWHM_GAUSSIAN 2.3548200450309493820231386529193992755


class ExposureWeight
{
public:
    typedef std::vector<std::string> argument_list_t;
    typedef argument_list_t::const_iterator argument_const_iterator;

    ExposureWeight() : y_optimum_(0.5), width_(0.2) {}

    ExposureWeight(double y_optimum, double width_parameter) :
        y_optimum_(y_optimum), width_(width_parameter)
    {
        check_invariant();
    }

    virtual void initialize(double y_optimum, double width_parameter,
                            argument_const_iterator arguments_begin,
                            argument_const_iterator arguments_end)
    {
        y_optimum_ = y_optimum;
        width_ = width_parameter;

        arguments_.clear();
        std::copy(arguments_begin, arguments_end, std::back_inserter(arguments_));

        check_invariant();
    }

    double optimum() const {return y_optimum_;}
    double width() const {return width_;}
    const argument_list_t& arguments() const {return arguments_;}

    virtual double normalize(double y) const {return (y - optimum()) / width();}

    virtual double weight(double) = 0;

    virtual ~ExposureWeight() {}

    struct error : public std::runtime_error
    {
        error(const std::string& message) : std::runtime_error(message) {}
    };

    virtual int interface_version() const {return EXPOSURE_WEIGHT_INTERFACE_VERSION;}

private:
    void check_invariant() const
    {
        if (y_optimum_ < 0.0 || y_optimum_ > 1.0)
        {
            throw std::invalid_argument("y_optimum");
        }
        if (width_ <= 0.0)
        {
            throw std::invalid_argument("width_parameter");
        }
    }

    double y_optimum_;
    double width_;
    argument_list_t arguments_;
};


#endif // EXPOSURE_WEIGHT_BASE_INCLUDED


// Local Variables:
// mode: c++
// End:
