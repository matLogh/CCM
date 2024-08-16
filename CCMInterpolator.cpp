#include "CCMInterpolator.h"

TEC::CCMInterpolator::CCMInterpolator(std::string type, bool valid_only) : fValidOnly(valid_only)
{
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    type.erase(std::remove(type.begin(), type.end(), ' '), type.end());
    if (type == "linear")
    {
        fType = ROOT::Math::Interpolation::kLINEAR;
    }
    else if (type == "polynomial")
    {
        fType = ROOT::Math::Interpolation::kPOLYNOMIAL;
    }
    else if (type == "cspline")
    {
        fType = ROOT::Math::Interpolation::kCSPLINE;
    }
    else if (type == "cspline_periodic")
    {
        fType = ROOT::Math::Interpolation::kCSPLINE_PERIODIC;
    }
    else if (type == "akima")
    {
        fType = ROOT::Math::Interpolation::kAKIMA;
    }
    else if (type == "akima_periodic")
    {
        fType = ROOT::Math::Interpolation::kAKIMA_PERIODIC;
    }
    else
    {
        throw std::runtime_error("Unknown interpolation type");
    }
}

TEC::CCMInterpolator::~CCMInterpolator()
{
    if (fInterpolator != nullptr)
        delete fInterpolator;
}

void TEC::CCMInterpolator::Reset()
{
    if (fInterpolator != nullptr)
    {
        delete fInterpolator;
        fInterpolator = nullptr;
    }

    fX.clear();
    fY.clear();
    fValid.clear();
}

const int TEC::CCMInterpolator::GetIndex(const double x) const
{
    double diff = std::numeric_limits<double>::max();
    int index = -1;

    for (int i = 0; i < fX.size(); i++)
    {
        if (std::abs(fX[i] - x) < diff)
        {
            diff = std::abs(fX[i] - x);
            index = i;
        }
    }

    return index;
}

const bool TEC::CCMInterpolator::InterpolationValid(const double x) const
{
    int index = this->GetIndex(x);
    double diff = x - fX[index];

    // check if the closest point is valid
    if (!fValid[index])
        return false;

    // interpolating before first point in array
    if (index == 0 && diff < 0)
        return false;
    // interpolating after last point in array
    if (index == fX.size() - 1 && diff > 0)
        return false;

    // check validity of just closest neighbour
    if (fType == ROOT::Math::Interpolation::kLINEAR)
    {
        if (diff < 0)
            return fValid[index - 1];
        else
            return fValid[index + 1];
    }
    else
    {
        for (int i = index - 1; i <= index + 1; i++)
        {
            if (i < 0 || i >= fX.size())
                continue;
            if (!fValid[i])
                return false;
        }
    }

    return true;
}

const bool TEC::CCMInterpolator::IsValueValid(const double x) const
{
    int index = this->GetIndex(x);
    if (index == -1)
        return false;
    return fValid[index];
}

const double TEC::CCMInterpolator::Eval_noInterpolation(const double x)
{
    int index = this->GetIndex(x);
    if (index == -1)
        return 0;
    return fY[index];
}

const double TEC::CCMInterpolator::Eval(const double x)
{
    if (!fEnableInterpolation)
        return Eval_noInterpolation(x);

    if (!fInterpolator)
        fInterpolator = new ROOT::Math::Interpolator(fX, fY, fType);

    // if neighbours are both valid, return interpolated value

    int index = this->GetIndex(x);
    if (index == -1)
    {
        return 0;
    }

    if (!this->InterpolationValid(x))
    {
        return fY[index];
    }

    return fInterpolator->Eval(x);
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
}