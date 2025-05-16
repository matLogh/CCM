#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include <TH2.h>
#pragma GCC diagnostic pop

#include <iostream>
#include <stdexcept>

#include "CCMInterpolator.h"

namespace TEC
{

class RegionOfInterest
{
  private:
    friend class CCM;
    friend class CrossCorrel;

    double energy_window_low;
    double energy_window_high;
    double energy_displacement_low;
    double energy_displacement_high;
    double desired_energy;

    /// @brief  bin_window_low: bin number of the energy_window_low
    int bin_window_low;
    /// @brief  bin_window_high: bin number of the energy_window_high
    int bin_window_high;
    /// @brief how many bins below window_low are we shifting
    int bin_displacement_low;
    /// @brief how many bins above window_high are we shifting
    int bin_displacement_high;

    /// @brief dimension of vectors that is used for cross-correlation
    int vector_dimension;
    /// @brief how many bins are in the ROI it total - number of bins between
    /// bin_window_low minus displacement_low and bin_window_high plus displacement_high
    int displacement_range;
    /// @brief how many shifts of the vector are to be performed
    int displacement_steps;

    /// @brief just a helper variable, used to calculate initial shift of the vector
    int base_shift_value;

    TEC::CCMInterpolator interpolator;
    bool                 is_interpolator_initialized{false};

  public:
    RegionOfInterest() = delete;
    RegionOfInterest(const std::shared_ptr<TH2> matrix,
                     const double               _energy_window_low,
                     const double               _energy_window_high,
                     const double               _energy_displacement_low,
                     const double               _energy_displacement_high,
                     const double               _desired_energy);
    // Parameterized constructor
    RegionOfInterest(double _energy_window_low,
                     double _energy_window_high,
                     double _energy_displacement_low,
                     double _energy_displacement_high,
                     double _desired_energy,
                     int    _bin_window_low,
                     int    _bin_window_high,
                     int    _bin_displacement_low,
                     int    _bin_displacement_high,
                     int    _vector_dimension,
                     int    _displacement_range,
                     int    _displacement_steps,
                     int    _base_shift_value)
        : energy_window_low(_energy_window_low), energy_window_high(_energy_window_high),
          energy_displacement_low(_energy_displacement_low),
          energy_displacement_high(_energy_displacement_high),
          desired_energy(_desired_energy), bin_window_low(_bin_window_low),
          bin_window_high(_bin_window_high), bin_displacement_low(_bin_displacement_low),
          bin_displacement_high(_bin_displacement_high),
          vector_dimension(_vector_dimension), displacement_range(_displacement_range),
          displacement_steps(_displacement_steps), base_shift_value(_base_shift_value)
    {
        if (_energy_window_low > _energy_window_high)
            throw std::runtime_error("ERROR in constructor of RegionOfInterest: "
                                     "energy_window_low > energy_window_high");
        if (_energy_displacement_low > _energy_displacement_high)
            throw std::runtime_error(
                "ERROR in constructor of RegionOfInterest: "
                "energy_displacement_low > energy_displacement_high");

        if (bin_window_low - bin_displacement_low < 1)
            throw std::runtime_error("ERROR in constructor of RegionOfInterest: Energy "
                                     "window - displacement is "
                                     "smaller then range of matrix!");
    }

    // Copy constructor
    RegionOfInterest(const RegionOfInterest &other) = default;

    // Copy assignment operator
    RegionOfInterest &operator=(const RegionOfInterest &other) = default;

    // Destructor
    ~RegionOfInterest() = default;

    void Print() const noexcept;

    RegionOfInterest Clone(const std::shared_ptr<TH2> matrix) const noexcept;
};

} // End of namespace TEC