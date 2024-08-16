#pragma once

#include <iomanip>
#include <iostream>

#include "Math/Interpolator.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

namespace TEC
{

class CCMInterpolator
{
  public:
    /// @brief Interpolator constructor, default type is kLINEAR. List of available types:
    // LINEAR interpolation;
    // POLYNOMIAL interpolation, to be used for small number of points since introduces large oscillations;
    // CSPLINE cubic spline with natural boundary conditions;
    // CSPLINE_PERIODIC cubic spline with periodic boundary conditions;
    // AKIMA, Akima spline with natural boundary conditions ( requires a minimum of 5 points);
    // AKIMA_PERIODIC, Akima spline with periodic boundaries ( requires a minimum of 5 points);
    /// @param type
    CCMInterpolator(std::string type, bool valid_only = true);
    CCMInterpolator(ROOT::Math::Interpolation::Type type, bool valid_only = true)
        : fType(type), fValidOnly(valid_only){};

    // const ROOT::Math::Interpolation::Type type = ROOT::Math::Interpolation::kLINEAR);
    ~CCMInterpolator();

    void AddPoint(const double x, const double y, const bool valid);

    void DisableInterpolation()
    {
        fEnableInterpolation = false;
    };
    void EnableInterpolation()
    {
        fEnableInterpolation = true;
    };

    ///@brief check if the point x is valid
    const bool IsValueValid(const double x) const;

    /// @brief Evaluates the interpolation at point x. Function ALWAYS retunrs a value. If interpolation is marked as
    /// invalid, function returned the shift of the closest point. Interpolation is marked invalid if the closest point
    /// is invalid or when enough neighbouring points are marked as invalid.
    /// @param x
    /// @return
    const double Eval(const double x);

    /// @brief Get the value of the closest point to x. If x is outside the range of the array, function returns 0.
    /// @param x
    /// @return
    const double Eval_noInterpolation(const double x);

    const double GetFirstX() const
    {
        return fX.front();
    }
    const double GetLastX() const
    {
        return fX.back();
    }
    const int GetNPoints() const
    {
        return fX.size();
    }

    ROOT::Math::Interpolation::Type GetType() const
    {
        return fType;
    }

    void Reset();

  private:
    ROOT::Math::Interpolator *fInterpolator{nullptr};

    bool fValidOnly;
    bool fEnableInterpolation{true};
    ROOT::Math::Interpolation::Type fType;

    std::vector<double> fX{};
    std::vector<double> fY{};
    std::vector<bool> fValid{};

    const int GetIndex(const double x) const;
    const bool InterpolationValid(const double x) const;
};

} // namespace TEC