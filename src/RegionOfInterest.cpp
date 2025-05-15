#include "RegionOfInterest.h"

TEC::RegionOfInterest::RegionOfInterest(const std::shared_ptr<TH2> matrix,
                                        const double               _energy_window_low,
                                        const double               _energy_window_high,
                                        const double _energy_displacement_low,
                                        const double _energy_displacement_high,
                                        const double _desired_energy)
    : energy_window_low(_energy_window_low), energy_window_high(_energy_window_high),
      energy_displacement_low(_energy_displacement_low),
      energy_displacement_high(_energy_displacement_high),
      desired_energy(_desired_energy),
      bin_window_low(matrix->GetYaxis()->FindBin(_energy_window_low)),
      bin_window_high(matrix->GetYaxis()->FindBin(_energy_window_high)),
      bin_displacement_low(
          static_cast<int>(energy_displacement_low / matrix->GetYaxis()->GetBinWidth(1))),
      bin_displacement_high(static_cast<int>(energy_displacement_high /
                                             matrix->GetYaxis()->GetBinWidth(1))),
      vector_dimension(bin_window_high - bin_window_low),
      displacement_range(bin_displacement_high - bin_displacement_low + vector_dimension),
      displacement_steps(displacement_range - vector_dimension),
      //   base_shift_value(bin_window_low - bin_displacement_low)
      base_shift_value(bin_displacement_low)
{
    if (_energy_window_low > _energy_window_high)
        throw std::runtime_error("ERROR in constructor of RegionOfInterest: "
                                 "energy_window_low > energy_window_high");
    if (_energy_displacement_low > _energy_displacement_high)
        throw std::runtime_error("ERROR in constructor of RegionOfInterest: "
                                 "energy_displacement_low > energy_displacement_high");
    if (bin_window_high + bin_displacement_high > matrix->GetYaxis()->GetNbins())
        throw std::runtime_error("ERROR in constructor of RegionOfInterest: Energy "
                                 "window + displacement is larger then range of matrix!");
    if (bin_window_low - bin_displacement_low < 1)
        throw std::runtime_error(
            "ERROR in constructor of RegionOfInterest: Energy window - displacement is "
            "smaller then range of matrix!");
}

void TEC::RegionOfInterest::Print() const noexcept
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

TEC::RegionOfInterest TEC::RegionOfInterest::Clone(
    const std::shared_ptr<TH2> matrix) const noexcept
{
    return RegionOfInterest(matrix, this->energy_window_low, this->energy_window_high,
                            this->energy_displacement_low, this->energy_displacement_high,
                            this->desired_energy);
}