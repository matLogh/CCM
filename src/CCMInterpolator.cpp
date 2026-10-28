#include "CCMInterpolator.h"

#include <atomic>
#include <cmath>

namespace
{
std::atomic<std::uint64_t> gInitializationFailures{0};
std::atomic<std::uint64_t> gNonFiniteResults{0};
std::atomic<std::uint64_t> gLinearFallbacks{0};
std::atomic<std::uint64_t> gNearestFallbacks{0};
std::atomic<std::uint64_t> gOutOfRangeEvaluations{0};
} // namespace

TEC::CCMInterpolator::CCMInterpolator() { this->SetType("akima", true); }

TEC::CCMInterpolator::CCMInterpolator(const std::string &type, bool valid_only) { this->SetType(type, valid_only); }

TEC::CCMInterpolator::CCMInterpolator(ROOT::Math::Interpolation::Type type, bool valid_only) { this->SetType(type, valid_only); }

void TEC::CCMInterpolator::SetType(const ROOT::Math::Interpolation::Type type, const bool valid_only)
{
    if (type == ROOT::Math::Interpolation::kPOLYNOMIAL)
    {
        throw std::runtime_error("Polynomial interpolation is disabled on purpose, it produces unphysical "
                                 "results! Use AKIMA or CSPLINE instead.");
    }
    fType                    = type;
    fInterpolator            = nullptr;
    fInterpolatorInitialized = false;
    fValidOnly               = valid_only;
    fEnableInterpolation     = true;
}

TEC::CCMInterpolator::CCMInterpolator(const TEC::CCMInterpolator &other)
{
    this->fValidOnly           = other.fValidOnly;
    this->fEnableInterpolation = other.fEnableInterpolation;
    this->fType                = other.fType;

    this->fX     = other.fX;
    this->fY     = other.fY;
    this->fValid = other.fValid;

    // ROOT::Math::Interpolator cannot be copied. Rebuild it lazily if it is needed.
    this->fInterpolator            = nullptr;
    this->fInterpolatorInitialized = false;
}

TEC::CCMInterpolator &TEC::CCMInterpolator::operator=(const TEC::CCMInterpolator &other)
{
    if (this != &other)
    {

        this->fValidOnly           = other.fValidOnly;
        this->fEnableInterpolation = other.fEnableInterpolation;
        this->fType                = other.fType;

        this->fX     = other.fX;
        this->fY     = other.fY;
        this->fValid = other.fValid;

        // ROOT::Math::Interpolator cannot be copied. Rebuild it lazily if it is needed.
        this->fInterpolator            = nullptr;
        this->fInterpolatorInitialized = false;
    }
    return *this;
}

void TEC::CCMInterpolator::SetType(std::string type, const bool valid_only)
{
    // parse string
    ROOT::Math::Interpolation::Type _type;

    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    type.erase(std::remove(type.begin(), type.end(), ' '), type.end());
    if (type == "linear") { _type = ROOT::Math::Interpolation::kLINEAR; }
    // polynomial interpolation is not useful due to large oscillations
    // else if (type == "polynomial") { _type = ROOT::Math::Interpolation::kPOLYNOMIAL; }
    else if (type == "cspline") { _type = ROOT::Math::Interpolation::kCSPLINE; }
    else if (type == "cspline_periodic") { _type = ROOT::Math::Interpolation::kCSPLINE_PERIODIC; }
    else if (type == "akima") { _type = ROOT::Math::Interpolation::kAKIMA; }
    else if (type == "akima_periodic") { _type = ROOT::Math::Interpolation::kAKIMA_PERIODIC; }
    else
    {
        std::string err_msg = "Unknown interpolation type: " + type;
        err_msg += "\nPossible types are: linear, cspline, cspline_periodic, akima, "
                   "akima_periodic";
        throw std::runtime_error(err_msg);
    }
    this->SetType(_type, valid_only);
}

void TEC::CCMInterpolator::Reset()
{
    fInterpolator            = nullptr;
    fInterpolatorInitialized = false;
    this->ClearPoints();
}

void TEC::CCMInterpolator::ClearPoints()
{
    fX.clear();
    fY.clear();
    fValid.clear();
    fInterpolator            = nullptr;
    fInterpolatorInitialized = false;
}

void TEC::CCMInterpolator::DisableInterpolation()
{
    fEnableInterpolation = false;
    fInterpolator        = nullptr;
}

void TEC::CCMInterpolator::EnableInterpolation()
{
    fEnableInterpolation     = true;
    fInterpolator            = nullptr;
    fInterpolatorInitialized = false;
};

TEC::CCMInterpolator::Diagnostics TEC::CCMInterpolator::GetDiagnostics() noexcept
{
    return {gInitializationFailures.load(),
            gNonFiniteResults.load(),
            gLinearFallbacks.load(),
            gNearestFallbacks.load(),
            gOutOfRangeEvaluations.load()};
}

void TEC::CCMInterpolator::ResetDiagnostics() noexcept
{
    gInitializationFailures.store(0);
    gNonFiniteResults.store(0);
    gLinearFallbacks.store(0);
    gNearestFallbacks.store(0);
    gOutOfRangeEvaluations.store(0);
}

int TEC::CCMInterpolator::GetIndex(const double x) const
{
    auto lower = std::lower_bound(fX.begin(), fX.end(), x);
    if (lower == fX.end()) { return -1; }
    else if (lower == fX.begin()) { return 0; }

    auto index = static_cast<int>(std::distance(fX.begin(), lower));
    if (std::abs(fX[static_cast<uint>(index)] - x) < std::abs(fX[static_cast<uint>(index) - 1] - x)) { return index; }
    return index - 1;
}

int TEC::CCMInterpolator::InterpolationValid(const double x, const int index) const
{
    // int    index = this->GetIndex(x);
    if (index < 0 || index >= static_cast<int>(fX.size()))
    {
        throw std::runtime_error("TEC::CCMInterpolator::InterpolationValid: Error! Index out of bounds");
    }
    double diff = x - fX[static_cast<uint>(index)];

    // check if the closest point is valid
    if (!fValid[static_cast<uint>(index)]) return 0;

    if (diff < 0)
    {
        if (index - 1 < 0) return 2;
        if (!fValid[static_cast<uint>(index) - 1]) return 0;
    }
    if (diff > 0)
    {
        if (index + 1 >= static_cast<int>(fX.size())) return 2;
        if (!fValid[static_cast<uint>(index) + 1]) return 0;
    }

    // interpolating before first point in array
    // if (index == 0 && diff < 0) return false;
    // interpolating after last point in array
    // if (index == fX.size() - 1 && diff > 0) return false;

    // check validity of just closest neighbour
    switch (fType)
    {
    case ROOT::Math::Interpolation::kLINEAR: {
        // we need just our point and closest neighbour, check already performed above
        return 1;
    }
    case ROOT::Math::Interpolation::kCSPLINE:
    case ROOT::Math::Interpolation::kCSPLINE_PERIODIC: {
        // we need 3 points, check if we have them
        if (index - 1 < 0 || index + 1 >= static_cast<int>(fX.size())) return 2;
        return fValid[static_cast<uint>(index - 1)] && fValid[static_cast<uint>(index + 1)] ? 1 : 2;
    }
    default: {
        // we need 5 points, check if we have them
        if (index - 2 < 0 || index + 2 >= static_cast<int>(fX.size())) return 2;
        return fValid[static_cast<uint>(index - 2)] && fValid[static_cast<uint>(index - 1)] && fValid[static_cast<uint>(index + 1)] &&
                       fValid[static_cast<uint>(index + 2)]
                   ? 1
                   : 2;
    }
    }
}

bool TEC::CCMInterpolator::IsValueValid(const double x) const
{
    int index = this->GetIndex(x);
    if (index == -1) return false;
    return fValid[static_cast<uint>(index)];
}

double TEC::CCMInterpolator::Eval_noInterpolation(const double x)
{
    int index = this->GetIndex(x);
    if (index == -1) return 0;
    return fY[static_cast<uint>(index)];
}

double TEC::CCMInterpolator::Eval(const double x)
{
    if (!fEnableInterpolation) return Eval_noInterpolation(x);

    const int _index = this->GetIndex(x);
    // closest value is not valid, interpolation is invalid by default
    if (_index == -1)
    {
        ++gOutOfRangeEvaluations;
        return 0;
    }
    const uint index = static_cast<uint>(_index);
    // is interpolation possible? If no, get the closest value, if yes get interpolated
    // value
    int valid = this->InterpolationValid(x, _index);
    if (valid == 0)
    {
        ++gNearestFallbacks;
        return fY[index];
    }

    if (valid == 1 && this->InitializeInterpolator())
    {
        const double result = fInterpolator->Eval(x);
        if (std::isfinite(result)) { return result; }
        ++gNonFiniteResults;
    }

    bool         used_linear = false;
    const double fallback    = this->EvalLinearOrNearest(x, index, used_linear);
    if (used_linear) { ++gLinearFallbacks; }
    else
    {
        ++gNearestFallbacks;
    }
    return fallback;
}

bool TEC::CCMInterpolator::InitializeInterpolator()
{
    if (fInterpolator) { return true; }
    if (fInterpolatorInitialized) { return false; }

    // Mark the attempt even when it fails, otherwise every evaluation would emit the same
    // low-level GSL error again.
    fInterpolatorInitialized = true;
    if (!this->HasValidInterpolationData())
    {
        ++gInitializationFailures;
        return false;
    }

    auto interpolator = std::make_unique<ROOT::Math::Interpolator>(static_cast<unsigned int>(fX.size()), fType);
    if (!interpolator->SetData(fX, fY))
    {
        ++gInitializationFailures;
        return false;
    }

    fInterpolator = std::move(interpolator);
    return true;
}

double TEC::CCMInterpolator::EvalLinearOrNearest(const double x, const uint index, bool &used_linear) const
{
    used_linear       = false;
    const double diff = x - fX[index];
    if (diff == 0) { return fY[index]; }

    size_t left  = index;
    size_t right = index;
    if (diff < 0)
    {
        if (index == 0) { return fY[index]; }
        left  = static_cast<size_t>(index) - 1;
        right = index;
    }
    else
    {
        if (static_cast<size_t>(index) + 1 >= fY.size()) { return fY[index]; }
        left  = index;
        right = static_cast<size_t>(index) + 1;
    }

    const double dx = fX[right] - fX[left];
    if (!(dx > 0) || !std::isfinite(dx)) { return fY[index]; }

    used_linear = true;
    return fY[left] + (fY[right] - fY[left]) / dx * (x - fX[left]);
}

size_t TEC::CCMInterpolator::MinimumRequiredPoints() const noexcept
{
    switch (fType)
    {
    case ROOT::Math::Interpolation::kLINEAR: return 2;
    case ROOT::Math::Interpolation::kCSPLINE:
    case ROOT::Math::Interpolation::kCSPLINE_PERIODIC: return 3;
    case ROOT::Math::Interpolation::kAKIMA:
    case ROOT::Math::Interpolation::kAKIMA_PERIODIC: return 5;
    default: return 2;
    }
}

bool TEC::CCMInterpolator::HasValidInterpolationData() const noexcept
{
    if (fX.size() != fY.size() || fX.size() < this->MinimumRequiredPoints()) { return false; }

    for (size_t i = 0; i < fX.size(); ++i)
    {
        if (!std::isfinite(fX[i]) || !std::isfinite(fY[i])) { return false; }
        if (i > 0 && !(fX[i] > fX[i - 1])) { return false; }
    }
    return true;
}

void TEC::CCMInterpolator::AddPoint(const double x, const double y, const bool valid)
{
    if (fX.size() != 0 && x < fX.back())
    {
        int index = 0;
        for (uint i = 0; i < fX.size(); i++)
        {
            if (fX[i] > x)
            {
                index = static_cast<int>(i);
                break;
            }
        }
        fX.insert(fX.begin() + index, x);
        fY.insert(fY.begin() + index, y);
        fValid.insert(fValid.begin() + index, !fValidOnly || valid);
    }
    else
    {
        fX.push_back(x);
        fY.push_back(y);
        fValid.push_back(!fValidOnly || valid);
    }
    fInterpolator            = nullptr;
    fInterpolatorInitialized = false;
}
