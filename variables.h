#include "TH2D.h"
#include <atomic>
#include <iostream>
#include <vector>

#ifndef HEADER_H
#define HEADER_H

namespace TEC
{

extern std::atomic<int> thread_task;
// extern std::atomic<bool> wait;

struct Region_of_interest
{
    const double energy_window_low;
    const double energy_window_high;
    const double energy_displacement_low;
    const double energy_displacement_high;
    const double desired_energy;

    /// @brief  bin_window_low: bin number of the energy_window_low
    const int bin_window_low;
    /// @brief  bin_window_high: bin number of the energy_window_high
    const int bin_window_high;
    /// @brief how many bins below window_low are we moving
    const int bin_displacement_low;
    /// @brief how many bins above window_high are we moving
    const int bin_displacement_high;

    /// @brief dimension of vectors that is used for cross-correlation
    const int vector_dimension;
    /// @brief how many bins are in the ROI it total - number of bins between bin_window_low minus displacement_low and
    /// bin_window_high plus displacement_high
    const int displacement_range;
    /// @brief how many shifts of the vector are to be performed
    const int displacement_steps;

    const int base_shift_value;

    Region_of_interest() = delete;
    Region_of_interest(const TH2D &matrix, double _energy_window_low, double _energy_window_high,
                       double _energy_displacement_low, double _energy_displacement_high, double _desired_energy)
        : energy_window_low(_energy_window_low), energy_window_high(_energy_window_high),
          energy_displacement_low(_energy_displacement_low), energy_displacement_high(_energy_displacement_high),
          desired_energy(_desired_energy), bin_window_low(matrix.GetYaxis()->FindBin(_energy_window_low)),
          bin_window_high(matrix.GetYaxis()->FindBin(_energy_window_high)),
          bin_displacement_low(static_cast<int>(energy_displacement_low / matrix.GetYaxis()->GetBinWidth(1))),
          bin_displacement_high(static_cast<int>(energy_displacement_high / matrix.GetYaxis()->GetBinWidth(1))),
          vector_dimension(bin_window_high - bin_window_low),
          displacement_range(bin_displacement_high - bin_displacement_low + vector_dimension),
          displacement_steps(displacement_range - vector_dimension),
          //   base_shift_value(bin_window_low - bin_displacement_low)
          base_shift_value(bin_displacement_low)
    {
        if (_energy_window_low > _energy_window_high)
            throw std::runtime_error("energy_window_low > energy_window_high");
        if (_energy_displacement_low > _energy_displacement_high)
            throw std::runtime_error("energy_displacement_low > energy_displacement_high");
        if (bin_window_high + bin_displacement_high > matrix.GetYaxis()->GetNbins())
            throw std::runtime_error("Energy window + displacement is larger then range of matrix!");
        if (bin_window_low - bin_displacement_low < 1)
            throw std::runtime_error("Energy window - displacement is smaller then range of matrix!");
    }
    void print() const
    {
        std::cout << "energy_window_low: " << energy_window_low << '\n';
        std::cout << "energy_window_high: " << energy_window_high << '\n';
        std::cout << "energy_displacement_low: " << energy_displacement_low << '\n';
        std::cout << "energy_displacement_high: " << energy_displacement_high << '\n';
        std::cout << "desired_energy: " << desired_energy << '\n';
        std::cout << "bin_window_low: " << bin_window_low << '\n';
        std::cout << "bin_window_high: " << bin_window_high << '\n';
        std::cout << "bin_displacement_low: " << bin_displacement_low << '\n';
        std::cout << "bin_displacement_high: " << bin_displacement_high << '\n';
        std::cout << "vector_dimension: " << vector_dimension << '\n';
        std::cout << "displacement_range: " << displacement_range << '\n';
        std::cout << "displacement_steps: " << displacement_steps << '\n';
        std::cout << "base_shift_value: " << base_shift_value << '\n';
    }
};

struct VarManager
{
    // INPUT VARIABLES
    TH2D *TEMAT;
    double ***TEMATarr; // TH2D loaded as a data cube [ROI][x][y]
                        // keep in mind that coord x and y starts with 0 for every ROI

    std::vector<Region_of_interest> ROIs;
    std::vector<std::vector<double>> sample_vector;

    int number_of_ROIs;
    int total_tasks;
    int time_bins;
};
// definition of result container
struct ResCont
{
    bool isValid{true};
    double bin_shift{0.};
    double energy_shift{0.};
    double dp{0.};
    double poly_shift{0.};
    double gfit_sigma{0.};
    double gfit_mu{0.};
    double gfit_chi2{0.};
    std::vector<double> dp_vec;

    void print()
    {
        std::cout << "isValid: " << (isValid ? "true" : "false") << "\n";
        std::cout << "bin_shift: " << bin_shift << "\n";
        std::cout << "energy_shift: " << energy_shift << "\n";
        std::cout << "dp: " << dp << "\n";
        std::cout << "gfit_sigma: " << gfit_sigma << "\n";
        std::cout << "gfit_mu: " << gfit_mu << "\n";
        std::cout << "gfit_chi2: " << gfit_chi2 << "\n";
        std::cout << "dp_vec: ";
        for (const auto &val : dp_vec)
        {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
};

struct FitCont
{
    std::string functionUsed;
    std::vector<double> coef;
};

struct ShiftInterpolator
{
    TGraph valid_points;
    // TGraph invalid
};

} // namespace TEC

#endif