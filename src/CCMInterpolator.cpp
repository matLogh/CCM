#include "CCMInterpolator.h"

TEC::CCMInterpolator::CCMInterpolator() { this->SetType("akima", true); }

TEC::CCMInterpolator::CCMInterpolator(const std::string &type, bool valid_only)
{
    this->SetType(type, valid_only);
}

TEC::CCMInterpolator::CCMInterpolator(ROOT::Math::Interpolation::Type type,
                                      bool                            valid_only)
{
    this->SetType(type, valid_only);
}

void TEC::CCMInterpolator::SetType(const ROOT::Math::Interpolation::Type type,
                                   const bool                            valid_only)
{
    if (type == ROOT::Math::Interpolation::kPOLYNOMIAL)
    {
        throw std::runtime_error(
            "Polynomial interpolation is disabled on purpose, it produces unphysical "
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
    this->fValidOnly               = other.fValidOnly;
    this->fEnableInterpolation     = other.fEnableInterpolation;
    this->fInterpolatorInitialized = other.fInterpolatorInitialized;
    this->fType                    = other.fType;

    this->fX     = other.fX;
    this->fY     = other.fY;
    this->fValid = other.fValid;

    // making copy of interpolator is not public method
    if (other.fInterpolator)
    {
        this->fInterpolator =
            std::make_unique<ROOT::Math::Interpolator>(this->fX, this->fY, this->fType);
    }
    else { this->fInterpolator = nullptr; }
}

TEC::CCMInterpolator &TEC::CCMInterpolator::operator=(const TEC::CCMInterpolator &other)
{
    if (this != &other)
    {

        this->fValidOnly               = other.fValidOnly;
        this->fEnableInterpolation     = other.fEnableInterpolation;
        this->fInterpolatorInitialized = other.fInterpolatorInitialized;
        this->fType                    = other.fType;

        this->fX     = other.fX;
        this->fY     = other.fY;
        this->fValid = other.fValid;

        // making copy of interpolator is not public method
        if (other.fInterpolator)
        {
            this->fInterpolator = std::make_unique<ROOT::Math::Interpolator>(
                this->fX, this->fY, this->fType);
        }
        else { this->fInterpolator = nullptr; }
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
    else if (type == "cspline_periodic")
    {
        _type = ROOT::Math::Interpolation::kCSPLINE_PERIODIC;
    }
    else if (type == "akima") { _type = ROOT::Math::Interpolation::kAKIMA; }
    else if (type == "akima_periodic")
    {
        _type = ROOT::Math::Interpolation::kAKIMA_PERIODIC;
    }
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
}

void TEC::CCMInterpolator::DisableInterpolation()
{
    fEnableInterpolation = false;
    fInterpolator        = nullptr;
}

void TEC::CCMInterpolator::EnableInterpolation()
{
    fEnableInterpolation = true;
    fInterpolator        = nullptr;
};

int TEC::CCMInterpolator::GetIndex(const double x) const
{
    double diff = std::numeric_limits<double>::max();

    auto lower = std::lower_bound(fX.begin(), fX.end(), x);
    if (lower == fX.end()) { return -1; }
    else if (lower == fX.begin()) { return 0; }

    auto index = std::distance(fX.begin(), lower);
    if (std::abs(fX[index] - x) < std::abs(fX[index - 1] - x)) { return index; }
    return index - 1;
}

int TEC::CCMInterpolator::InterpolationValid(const double x, const int index) const
{
    // int    index = this->GetIndex(x);
    if (index < 0 || index >= fX.size())
    {
        throw std::runtime_error(
            "TEC::CCMInterpolator::InterpolationValid: Error! Index out of bounds");
    }
    double diff = x - fX[index];

    // check if the closest point is valid
    if (!fValid[index]) return 0;

    if (diff < 0 && !fValid[index - 1]) return 0;
    if (diff > 0 && !fValid[index + 1]) return 0;

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
    case ROOT::Math::Interpolation::kCSPLINE ||
        ROOT::Math::Interpolation::kCSPLINE_PERIODIC: {
        // we need 3 points, check if we have them
        if (index - 1 < 0 || index + 1 >= fX.size()) return 2;
        return fValid[index - 1] && fValid[index + 1] ? 1 : 2;
    }
    default: {
        // we need 5 points, check if we have them
        if (index - 2 < 0 || index + 2 >= fX.size()) return 2;
        return fValid[index - 2] && fValid[index - 1] && fValid[index + 1] &&
                       fValid[index + 2]
                   ? 1
                   : 2;
    }
    }
}

bool TEC::CCMInterpolator::IsValueValid(const double x) const
{
    int index = this->GetIndex(x);
    if (index == -1) return false;
    return fValid[index];
}

double TEC::CCMInterpolator::Eval_noInterpolation(const double x)
{
    int index = this->GetIndex(x);
    if (index == -1) return 0;
    return fY[index];
}

double TEC::CCMInterpolator::Eval(const double x)
{
    if (!fEnableInterpolation) return Eval_noInterpolation(x);

    if (!fInterpolator)
    {
        fInterpolator = std::make_unique<ROOT::Math::Interpolator>(fX, fY, fType);
    }

    const int index = this->GetIndex(x);
    // closest value is not valid, interpolation is invalid by default
    if (index == -1) { return 0; }

    // is interpolation possible? If no, get the closest value, if yes get interpolated
    // value
    int valid = this->InterpolationValid(x, index);
    if (valid == 0) { return fY[index]; }
    if (valid == 1) { return fInterpolator->Eval(x); }

    // return val is 2, if we cannot do selected interpolation, but simple linear
    // interpolation is possible
    double diff = x - fX[index];
    if (diff < 0)
    {
        return fY[index] +
               (fY[index] - fY[index - 1]) / (fX[index] - fX[index - 1]) * diff;
    }
    return fY[index] + (fY[index + 1] - fY[index]) / (fX[index + 1] - fX[index]) * diff;
}

void TEC::CCMInterpolator::AddPoint(const double x, const double y, const bool valid)
{
    if (fX.size() != 0 && x < fX.back())
    {
        int index = 0;
        for (int i = 0; i < fX.size(); i++)
        {
            if (fX[i] > x)
            {
                index = i;
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
    fInterpolator = nullptr;
}