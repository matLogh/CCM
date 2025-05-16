
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include <TH2.h>
#pragma GCC diagnostic pop

#include <atomic>
#include <iostream>
#include <vector>

#ifndef HEADER_H
#define HEADER_H

namespace TEC
{

extern std::atomic<int> thread_task;

/// @brief container to store/reference variables used across the classes. This container
/// is managed by CCM object
struct VarManager
{
    // INPUT VARIABLES
    std::shared_ptr<TH2> TEMAT;
    float             ***TEMATarr; // TH2 loaded as a data cube [ROI][x][y]
                       // keep in mind that coord x and y starts with 0 for every ROI
    std::vector<RegionOfInterest> ROIs;

    std::vector<std::vector<float>> sample_vector;

    size_t number_of_ROIs;
    size_t total_tasks;
    size_t time_bins;
};
// definition of result container
struct ResCont
{
    bool               isValid{true};
    double             bin_shift{0.};
    double             energy_shift{0.};
    double             dp{0.};
    double             poly_shift{0.};
    double             gfit_sigma{0.};
    double             gfit_mu{0.};
    double             gfit_chi2{0.};
    std::vector<float> dp_vec;

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
        for (const auto &val : dp_vec) { std::cout << val << " "; }
        std::cout << "\n";
    }
};

struct FitCont
{
    std::string         functionUsed;
    std::vector<double> coef;
};

} // namespace TEC

#endif