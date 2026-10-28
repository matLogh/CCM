#pragma once

#include <iomanip>
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include "Math/Interpolator.h"
#pragma GCC diagnostic pop

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace TEC
{

class CCMInterpolator
{
  public:
    struct Diagnostics
    {
        std::uint64_t initialization_failures{0};
        std::uint64_t non_finite_results{0};
        std::uint64_t linear_fallbacks{0};
        std::uint64_t nearest_fallbacks{0};
        std::uint64_t out_of_range_evaluations{0};

        bool HasProblems() const noexcept
        {
            return initialization_failures != 0 || non_finite_results != 0 || linear_fallbacks != 0 || nearest_fallbacks != 0 ||
                   out_of_range_evaluations != 0;
        }
    };

    CCMInterpolator();
    /// @brief Interpolator constructor, default type is kLINEAR. List of available types:
    // LINEAR interpolation between 2 neighbors;
    // CSPLINE cubic spline with natural boundary conditions;
    // CSPLINE_PERIODIC cubic spline with periodic boundary conditions;
    // AKIMA, Akima spline with natural boundary conditions ( requires a minimum of 5
    // points); AKIMA_PERIODIC, Akima spline with periodic boundaries ( requires a minimum
    // of 5 points);
    // Note: POLYNOMIAL interpolation is disabled because it introduces HUGE oscillations
    // for more than just a few points
    /// @param type
    CCMInterpolator(const std::string &type, bool valid_only = true);
    CCMInterpolator(ROOT::Math::Interpolation::Type type, bool valid_only = true);

    CCMInterpolator(const CCMInterpolator &other);
    CCMInterpolator(CCMInterpolator &&other) = default;
    CCMInterpolator &operator=(const CCMInterpolator &other);

    void AddPoint(const double x, const double y, const bool valid);

    void DisableInterpolation();
    void EnableInterpolation();

    ///@brief check if the point x is valid
    bool IsValueValid(const double x) const;

    /// @brief Evaluates the interpolation at point x. Function ALWAYS returns a value. If
    /// interpolation is marked as invalid, function returned the shift of the closest
    /// point. Interpolation is marked invalid if the closest point is invalid or when
    /// enough neighboring points are marked as invalid.
    /// @param x
    /// @return
    double Eval(const double x);

    double GetFirstX() const { return fX.front(); }

    double GetLastX() const { return fX.back(); }
    int    GetNPoints() const { return static_cast<int>(fX.size()); }

    ROOT::Math::Interpolation::Type GetType() const { return fType; }

    static Diagnostics GetDiagnostics() noexcept;
    static void        ResetDiagnostics() noexcept;

    void Reset();
    void ClearPoints();
    void SetType(std::string type, const bool valid_only = true);
    void SetType(const ROOT::Math::Interpolation::Type type, const bool valid_only = true);

  private:
    std::unique_ptr<ROOT::Math::Interpolator> fInterpolator{nullptr};

    ROOT::Math::Interpolation::Type fType;
    bool                            fValidOnly;
    bool                            fEnableInterpolation;
    bool                            fInterpolatorInitialized;

    std::vector<double> fX{};
    std::vector<double> fY{};
    std::vector<bool>   fValid{};

    int GetIndex(const double x) const;

    /// @brief Answers question: do we have enough valid points around x to interpolate?
    /// @param x
    /// @param index index representing x in the array, obtained by GetIndex(x). Used as
    /// an input to avoid duplicit calls of GetIndex
    /// @return 0 if not valid; 1 if valid; 2 if not valid but linear interpolation is
    /// possible
    int InterpolationValid(const double x, const int index) const;

    bool   InitializeInterpolator();
    double EvalLinearOrNearest(const double x, const uint index, bool &used_linear) const;
    size_t MinimumRequiredPoints() const noexcept;
    bool   HasValidInterpolationData() const noexcept;

  private:
    /// @brief Get the value of the closest point to x. If x is outside the range of the
    /// array, function returns 0.
    /// @param x
    /// @return
    double Eval_noInterpolation(const double x);
};

} // namespace TEC
