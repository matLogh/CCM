#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#include "TCanvas.h"
#include "TCollection.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include "TVirtualFitter.h"

#include "Util.hpp"
#pragma once

class TheuerkaufPeak
{
    friend class TheuerkaufFitter;

  public:
    enum ParamState
    {
        FREE,
        FIXED,
        NONE,
        SAME
    };

  public:
    TheuerkaufPeak() : fFcn(nullptr){};
    TheuerkaufPeak(double min, double max, int id = 0, bool hasTL = false, bool hasTR = false, bool hasStep = false);
    TheuerkaufPeak(const TheuerkaufPeak &peak);
    ~TheuerkaufPeak();

    void Print() const;

    const int GetID() const noexcept
    {
        return fId;
    };

    const int GetIndex_Volume() const noexcept;
    const int GetIndex_Position() const noexcept;
    const int GetIndex_Sigma() const noexcept;
    const int GetIndex_TailLeft() const noexcept;
    const int GetIndex_TailRight() const noexcept;
    const int GetIndex_StepHeight() const noexcept;
    const int GetIndex_StepWidth() const noexcept;

    void Draw(std::string options = "")
    {
        fFcn->Draw(options.c_str());
    };

    void SetProperty_Volume(const ParamState prop) noexcept
    {
        fParamState[0] = prop;
    };
    void SetProperty_Position(const ParamState prop) noexcept
    {
        fParamState[1] = prop;
    };
    void SetProperty_Sigma(const ParamState prop) noexcept
    {
        fParamState[2] = prop;
    };
    void SetProperty_TailLeft(const ParamState prop) noexcept
    {
        fParamState[3] = prop;
    };
    void SetProperty_TailRight(const ParamState prop) noexcept
    {
        fParamState[4] = prop;
    };
    void SetProperty_StepHeight(const ParamState prop) noexcept
    {
        fParamState[5] = prop;
    };
    void SetProperty_StepWidth(const ParamState prop) noexcept
    {
        fParamState[6] = prop;
    };

    /// @brief It is expected to get full range of parameters in this order:
    ///     p[fParamIndex[0]] = volume
    ///     p[fParamIndex[1]] = position
    ///     p[fParamIndex[2]] = sigma
    ///     p[fParamIndex[3]] = left tail
    ///     p[fParamIndex[4]] = right tail
    ///     p[fParamIndex[5]] = step height
    ///     p[fParamIndex[6]] = step width
    /// @param x
    /// @param p
    /// @return
    double Eval(const double *x, const double *p) const;
    TF1 *GetFunction() const noexcept
    {
        return fFcn;
    };

    /// @brief These functions are relevant only with Theuerkauf fitter, or when peak function is coupled with another
    /// function and you want to share some parameters between them, such as sigma etc
    /// @param index
    /// @return
    TheuerkaufPeak *SetIndex_Volume(const int index);
    TheuerkaufPeak *SetIndex_Position(const int index);
    TheuerkaufPeak *SetIndex_Sigma(const int index);
    TheuerkaufPeak *SetIndex_TailLeft(const int index);
    TheuerkaufPeak *SetIndex_TailRight(const int index);
    TheuerkaufPeak *SetIndex_StepHeight(const int index);
    TheuerkaufPeak *SetIndex_StepWidth(const int index);

    TheuerkaufPeak *SetParameter_Volume(double val, ParamState prop = FREE, double min = -1e9, double max = 1e9);
    TheuerkaufPeak *SetParameter_Position(double val, ParamState prop = FREE, double min = -1e9, double max = 1e9);
    TheuerkaufPeak *SetParameter_Sigma(double val, ParamState prop = FREE,
                                       double min = std::numeric_limits<double>::epsilon(), double max = 1e9);
    TheuerkaufPeak *SetParameter_TailLeft(double val = 10., ParamState prop = NONE, double min = 0, double max = 1e12);
    TheuerkaufPeak *SetParameter_TailRight(double val = 10., ParamState prop = NONE, double min = 0, double max = 1e12);
    TheuerkaufPeak *SetParameter_StepHeight(double val, ParamState prop = NONE, double min = 0, double max = 1e9);
    TheuerkaufPeak *SetParameter_StepWidth(double val, ParamState prop = NONE, double min = 0, double max = 1e9);

    bool HasStep() const noexcept
    {
        return fHasStep;
    };
    bool HasTL() const noexcept
    {
        return fHasLeftTail;
    };
    bool HasTR() const noexcept
    {
        return fHasRightTail;
    };

    TheuerkaufPeak *SetRange(double min, double max);

    const ParamState GetState(int index) const;
    const ParamState GetState_Volume() const noexcept
    {
        return fParamState[0];
    };
    const ParamState GetState_Position() const noexcept
    {
        return fParamState[1];
    };
    const ParamState GetState_Sigma() const noexcept
    {
        return fParamState[2];
    };
    const ParamState GetState_TailLeft() const noexcept
    {
        return fParamState[3];
    };
    const ParamState GetState_TailRight() const noexcept
    {
        return fParamState[4];
    };
    const ParamState GetState_StepHeight() const noexcept
    {
        return fParamState[5];
    };
    const ParamState GetState_StepWidth() const noexcept
    {
        return fParamState[6];
    };

    const double GetVol() const noexcept
    {
        return fHistBinning_normalization <= 0 ? fFcn->GetParameter(0)
                                               : fFcn->GetParameter(0) * fHistBinning_normalization;
        // return fFcn->GetParameter(0);
    };
    const double GetPos() const noexcept
    {
        return fFcn->GetParameter(1);
    };
    const double GetSig() const noexcept
    {
        return fFcn->GetParameter(2);
    };
    const double GetTL() const noexcept
    {
        return fFcn->GetParameter(3);
    };
    const double GetTR() const noexcept
    {
        return fFcn->GetParameter(4);
    };
    const double GetSH() const noexcept
    {
        return fFcn->GetParameter(5);
    };
    const double GetSW() const noexcept
    {
        return fFcn->GetParameter(6);
    };

    const double GetFWHM() const noexcept
    {
        return sig_to_fwhm * this->GetSig();
    };

    // Filippo
    const double GetVolErr() const noexcept
    {
        return fHistBinning_normalization <= 0 ? fFcn->GetParError(0)
                                               : fFcn->GetParError(0) * fHistBinning_normalization;
        // return fFcn->GetParError(0);
    };
    const double GetPosErr() const noexcept
    {
        return fFcn->GetParError(1);
    };
    const double GetSigErr() const noexcept
    {
        return fFcn->GetParError(2);
    };
    const double GetTLErr() const noexcept
    {
        return fFcn->GetParError(3);
    };
    const double GetTRErr() const noexcept
    {
        return fFcn->GetParError(4);
    };
    const double GetSHErr() const noexcept
    {
        return fFcn->GetParError(5);
    };
    const double GetSWErr() const noexcept
    {
        return fFcn->GetParError(6);
    };

    const double GetFWHMErr() const noexcept
    {
        return sig_to_fwhm * this->GetSigErr();
    };

    void ResetIndexes(bool disregard_id = false) const;
    void SetBinning(const double &energy_per_bin)
    {
        fHistBinning_normalization = 1. / energy_per_bin;
    };
    const double GetBinning() const noexcept
    {
        return 1. / fHistBinning_normalization;
    };

  private:
    double EvalNoStep(const double *x, const double *p) const;
    double EvalStep(const double *x, const double *p) const;
    double GetNorm(const double sigma, const double tl, const double tr) const;
    // void SetMinMax(double min, double max);

    int fId;
    double fMin, fMax;
    bool fHasStep, fHasLeftTail, fHasRightTail;
    TF1 *fFcn;
    double fHistBinning_normalization{-1.}; // necessary to get the real volume, cleared of the

    static const double sig_to_fwhm;

    mutable double fCachedNorm{-1.}, fCachedSigma{-1.}, fCachedTL{-1.}, fCachedTR{-1.};
    mutable std::array<int, 7> fParamIndex;
    mutable std::array<ParamState, 7> fParamState;
};

const double TheuerkaufPeak::sig_to_fwhm = 2. * TMath::Sqrt(2. * TMath::Log(2.));

TheuerkaufPeak::TheuerkaufPeak(const TheuerkaufPeak &peak)
{
    fId = peak.fId;
    fMin = peak.fMin;
    fMax = peak.fMax;
    fHasStep = peak.fHasStep;
    fHasLeftTail = peak.fHasLeftTail;
    fHasRightTail = peak.fHasRightTail;
    fParamIndex = peak.fParamIndex;
    fParamState = peak.fParamState;

    // we need to make a new function to reference the ::Eval of THIS member
    std::string fcn_name = "theurekauf_" + std::to_string(fId);
    fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, fMin, fMax, 7,
                   "TheuerkaufFitter", "Eval");
    fFcn->SetNpx(10000);
    fFcn->SetLineColor(GetColor());
    for (int i = 0; i < 7; i++)
    {
        double min, max;
        peak.fFcn->GetParLimits(i, min, max);
        fFcn->SetParName(i, peak.fFcn->GetParName(i));
        fFcn->SetParLimits(i, min, max);
        fFcn->SetParameter(i, peak.fFcn->GetParameter(i));
    }
}

TheuerkaufPeak::TheuerkaufPeak(double min, double max, int id, bool hasTL, bool hasTR, bool hasStep)
    : fId(id), fHasLeftTail(hasTL), fHasRightTail(hasTR), fHasStep(hasStep), fMin(min), fMax(max)
{
    assert(id >= 0);
    std::string fcn_name = "theurekauf_" + std::to_string(id);
    fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, min, max, 7,
                   "TheuerkaufFitter", "Eval");
    fFcn->SetNpx(10000);
    fFcn->SetLineColor(GetColor());

    this->ResetIndexes();
    this->SetParameter_Volume(0., FREE)
        ->SetParameter_Position((min + max) / 2, FREE)
        ->SetParameter_Sigma(std::numeric_limits<double>::epsilon(), FREE)
        ->SetParameter_TailLeft(0., hasTL ? FREE : NONE)
        ->SetParameter_TailRight(0., hasTR ? FREE : NONE)
        ->SetParameter_StepHeight(0., hasStep ? FREE : NONE)
        ->SetParameter_StepWidth(0., hasStep ? FREE : NONE);
};

TheuerkaufPeak::~TheuerkaufPeak()
{
    if (fFcn)
    {
        delete fFcn;
    }
}

inline const int TheuerkaufPeak::GetIndex_Volume() const noexcept
{
    return fParamIndex[0];
};
inline const int TheuerkaufPeak::GetIndex_Position() const noexcept
{
    return fParamIndex[1];
};
inline const int TheuerkaufPeak::GetIndex_Sigma() const noexcept
{
    return fParamIndex[2];
};
inline const int TheuerkaufPeak::GetIndex_TailLeft() const noexcept
{
    return fParamIndex[3];
};
inline const int TheuerkaufPeak::GetIndex_TailRight() const noexcept
{
    return fParamIndex[4];
};
inline const int TheuerkaufPeak::GetIndex_StepHeight() const noexcept
{
    return fParamIndex[5];
};
inline const int TheuerkaufPeak::GetIndex_StepWidth() const noexcept
{
    return fParamIndex[6];
};

inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_Volume(const int index)
{
    fParamIndex[0] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_Position(const int index)
{
    fParamIndex[1] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_Sigma(const int index)
{
    fParamIndex[2] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_TailLeft(const int index)
{
    fParamIndex[3] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_TailRight(const int index)
{
    fParamIndex[4] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_StepHeight(const int index)
{
    fParamIndex[5] = index;
    return this;
};
inline TheuerkaufPeak *TheuerkaufPeak::SetIndex_StepWidth(const int index)
{
    fParamIndex[6] = index;
    return this;
};

inline TheuerkaufPeak *TheuerkaufPeak::SetRange(double min, double max)
{
    assert(min < max);

    fMin = min;
    fMax = max;
    fFcn->SetRange(fMin, fMax);

    double _min, _max;
    fFcn->GetParLimits(1, _min, _max);
    if (_min < min)
        _min = min;
    if (_max > max)
        _max = max;
    fFcn->SetParLimits(1, _min, _max);

    return this;
}

inline const TheuerkaufPeak::ParamState TheuerkaufPeak::GetState(int index) const
{
    assert(index >= 0 && index < 7);
    return fParamState[index];
};

// inline void TheuerkaufPeak::SetMinMax(double min, double max)
// {
//     assert(min < max);
//     fMin = min;
//     fMax = max;
//     fFcn->SetRange(min, max);
// };

inline void TheuerkaufPeak::Print() const
{
    std::cout << std::endl;
    std::cout << "TheuerkaufPeak model, peakID " << fId << std::endl;
    std::cout << "Function name " << fFcn->GetName() << std::endl << std::endl;
    std::cout << "PARAM   STATUS   VAL   MIN   MAX" << std::endl;

    for (int i = 0; i < 7; i++)
    {
        std::cout << fFcn->GetParName(i) << " ";
        if (GetState(i) == FREE)
            std::cout << "FREE  ";
        if (GetState(i) == FIXED)
            std::cout << "FIXED ";
        if (GetState(i) == NONE)
            std::cout << "NONE  ";
        if (GetState(i) == SAME)
            std::cout << "SAME  ";
        std::cout << fFcn->GetParameter(i) << " ";
        double min, max;
        fFcn->GetParLimits(i, min, max);
        std::cout << min << " " << max << std::endl;
    }
    std::cout << std::endl;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_Volume(double val, ParamState prop, double min, double max)
{
    assert(prop != NONE);
    assert(val >= min && val <= max);

    int index = 0;
    std::string name = std::string("vol_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());

    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    default:
        throw std::runtime_error("It is not allowed to have volume set as NONE or SAME!");
    }
    fFcn->SetParLimits(index, min, max);

    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_Position(double val, ParamState prop, double min, double max)
{
    assert(prop != NONE && prop != SAME);
    if (min < fMin)
        min = fMin;
    if (max > fMax)
        max = fMax;
    assert(val >= min && val <= max);

    int index = 1;
    std::string name = std::string("pos_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    default:
        throw std::runtime_error("It is not allowed to have volume set as NONE or SAME!");
    }
    fFcn->SetParLimits(index, min, max);

    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_Sigma(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 2;
    std::string name = std::string("sig_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    case NONE:
        std::cerr
            << "Setting sigma to NONE will turns peak to delta function, are you sure this is what you wanted to do?"
            << std::endl;
        fFcn->FixParameter(index, std::numeric_limits<double>::epsilon());
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_TailLeft(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 3;
    std::string name = std::string("TL__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasLeftTail = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasLeftTail = true;
        break;
    case NONE:
        fHasLeftTail = false;
        fFcn->FixParameter(index, std::numeric_limits<double>::max());
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasLeftTail = true;
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_TailRight(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 4;
    std::string name = std::string("TR__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasRightTail = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasRightTail = true;
        break;
    case NONE:
        fHasRightTail = false;
        fFcn->FixParameter(index, std::numeric_limits<double>::max());
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasRightTail = true;
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_StepHeight(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 5;
    std::string name = std::string("SH__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    bool hasStep = (prop != NONE);
    if (fHasStep != hasStep)
    {
        fHasStep = hasStep;
        this->SetParameter_StepWidth(1., prop);
    }

    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasStep = true;
        break;
    case NONE:
        fFcn->FixParameter(index, 0.);
        fHasStep = false;
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    }
    fFcn->SetParLimits(index, min, max);
    return this;
}

inline TheuerkaufPeak *TheuerkaufPeak::SetParameter_StepWidth(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 6;
    std::string name = std::string("SW__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    bool hasStep = (prop != NONE);
    if (fHasStep != hasStep)
    {
        fHasStep = hasStep;
        this->SetParameter_StepHeight(1., prop);
    }

    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasStep = true;
        break;
    case NONE:
        fHasStep = false;
        fFcn->FixParameter(index, 0.);
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    }
    fFcn->SetParLimits(index, min, max);
    return this;
}

double TheuerkaufPeak::Eval(const double *x, const double *p) const
{
    return EvalNoStep(x, p) + EvalStep(x, p);
}

double TheuerkaufPeak::EvalNoStep(const double *x, const double *p) const
{
    double vol = p[fParamIndex[0]];
    double dx = *x - p[fParamIndex[1]];
    double sigma = p[fParamIndex[2]];

    double tl = p[fParamIndex[3]];
    double tr = p[fParamIndex[4]];
    double norm = GetNorm(sigma, tl, tr);
    double _x;

    // Peak function
    if (dx < -tl && fHasLeftTail)
    {
        _x = tl / (sigma * sigma) * (dx + tl / 2.0);
    }
    else if (dx < tr || !fHasRightTail)
    {
        _x = -dx * dx / (2.0 * sigma * sigma);
    }
    else
    {
        _x = -tr / (sigma * sigma) * (dx - tr / 2.0);
    }

    return vol * norm * std::exp(_x);
}

double TheuerkaufPeak::EvalStep(const double *x, const double *p) const
{
    //! Step function

    if (fHasStep)
    {
        double vol = p[fParamIndex[0]];
        double dx = *x - p[fParamIndex[1]];
        double sigma = p[fParamIndex[2]];
        double tl = p[fParamIndex[3]];
        double tr = p[fParamIndex[4]];
        double sh = p[fParamIndex[5]];
        double sw = p[fParamIndex[6]];

        double norm = GetNorm(sigma, tl, tr);

        return vol * norm * sh * (M_PI / 2. + std::atan(sw * dx / (std::sqrt(2.) * sigma)));
    }
    else
    {
        return 0.0;
    }
}

double TheuerkaufPeak::GetNorm(const double sigma, const double tl, const double tr) const
{
    if (fCachedSigma == sigma && fCachedTL == tl && fCachedTR == tr)
    {
        return fCachedNorm;
    }

    double vol;

    // Contribution from left tail + left half of truncated gaussian
    if (fHasLeftTail)
    {
        vol = (sigma * sigma) / tl * std::exp(-(tl * tl) / (2.0 * sigma * sigma));
        vol += std::sqrt(M_PI / 2.0) * sigma * std::erf(tl / (std::sqrt(2.0) * sigma));
    }
    else
    {
        vol = std::sqrt(M_PI / 2.0) * sigma;
    }

    // Contribution from right tail + right half of truncated gaussian
    if (fHasRightTail)
    {
        vol += (sigma * sigma) / tr * std::exp(-(tr * tr) / (2.0 * sigma * sigma));
        vol += std::sqrt(M_PI / 2.0) * sigma * std::erf(tr / (std::sqrt(2.0) * sigma));
    }
    else
    {
        vol += std::sqrt(M_PI / 2.0) * sigma;
    }

    fCachedSigma = sigma;
    fCachedTL = tl;
    fCachedTR = tr;
    // fCachedNorm  = fHistBinning_normalization <= 0 ? 1. / vol : fHistBinning_normalization / vol;
    fCachedNorm = 1. / vol;

    return fCachedNorm;
}

void TheuerkaufPeak::ResetIndexes(bool disregard_id) const
{
    if (disregard_id)
        std::iota(fParamIndex.begin(), fParamIndex.end(), 0);
    else
        std::iota(fParamIndex.begin(), fParamIndex.end(), fId * 7);
};

//!  THEUERKAUF FITTER
//!  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TheuerkaufFitter
{
  public:
    TheuerkaufFitter(double min, double max);
    ~TheuerkaufFitter();

    /// @brief Add peak to fitter
    /// @param position rought position of the peak
    /// @param leftTail has left tail?
    /// @param rightTail has right tail?
    /// @param step has step?
    /// @return id of the peak
    const int AddPeak(const double position, const bool leftTail, const bool rightTail, const bool step);

    /// @brief This function will create a copy of the peak and set id to one of the internal structure
    /// @param peak
    /// @return id of added peak
    const int AddPeak(const TheuerkaufPeak &peak);
    TheuerkaufFitter *SetBacground(std::unique_ptr<TF1> &bcg_fcn) noexcept;

    /// @brief get peak based on its ID
    /// @param id
    /// @return
    TheuerkaufPeak *GetPeak(const int id);

    /// @brief Get number of peaks
    /// @return
    const int GetNPeaks();

    /// @brief Set internal polynomial background
    /// @param poly_order number of parameters (poly_order=0 no background, poly_order=3 quadratic etc)
    /// @return
    TheuerkaufFitter *SetBacgroundPoly(unsigned int poly_order) noexcept;

    /// @brief remove function set by the TheuerkaufFitter::SetBackground, internal polynomial background remains
    /// @return
    TheuerkaufFitter *RemoveBackground() noexcept;

    /// @brief If true, allows peaks to have negative volume
    /// @param allow_negative_peaks
    /// @return
    TheuerkaufFitter *AllowNegativePeaks(bool allow_negative_peaks) noexcept;

    /// @brief Change range of the fitter
    /// @param min X min
    /// @param max X max
    /// @return
    TheuerkaufFitter *SetRange(double min, double max);

    double Eval(const double *x, const double *p) const;

    /// @brief Fit histogram. Only here the total function is being computed.
    /// @param histFit histogram to be fitted
    /// @param options
    ///        "OUTPUT_NONE" - do not print fit results,
    ///        "OUTPUT_STANDARD" - print standard output,
    ///        "OUTPUT_ROOT" - print ROOT fit output,
    ///        "OUTPUT_PLAIN" - plain output without table (id, pos, pos_err, fwhm, fwhm_err, vol, vol_err)
    ///        "LIKELIHOOD" - use likelihood in fitting,
    ///        "POISSON" - use likelihood in fitting
    void Fit(TH1 *histFit, std::string options = "");
    void Analyze(TH1 *histAna, TPad *toPrint);
    TFitResultPtr GetFitResults()
    {
        return fFitResults;
    };

    std::shared_ptr<TF1> GetTotalFunction()
    {
        return fSumFunc;
    };

    // std::vector<std::pair<double,double>> GetPeaksIntegral(); //Filippo
    void PrintFitResults() const;
    void PrintFitResults_plain(std::ostream &os = std::cout) const;

    ///@brief Set output verbosity
    ///@param verbose 0 - do not print anything, 1 - standard print, 2 - print all
    void SetVerbosity(const int verbose)
    {
        fVerbose = verbose;
    };

  private:
    double EvalTotalBacground(const double *x, const double *p);
    const int GetNumParams() const noexcept
    {
        return static_cast<int>(fPeaks.size() * 7 + fPolyBcgDegree);
    };
    void HandleParameterStates();
    void DistributeParametersToPeaks();

  private:
    int fVerbose;
    double fMin, fMax;
    double fChiSquare;
    bool fOnlypositivepeaks;
    unsigned int fPolyBcgDegree;
    TFitResultPtr fFitResults;
    /// @brief points to the fitter following the main Fit
    TVirtualFitter *fFitter;

    std::shared_ptr<TF1> fSumFunc;
    std::unique_ptr<TF1> fBcgFunc;
    std::shared_ptr<TH1> fTempHist;

    std::vector<TheuerkaufPeak> fPeaks;
    std::vector<std::unique_ptr<TObject>> fTempObjects;
    std::vector<TheuerkaufPeak> fTempPeaks;
};

TheuerkaufFitter::~TheuerkaufFitter()
{
}

TheuerkaufFitter::TheuerkaufFitter(double min, double max) : fMin(min), fMax(max)
{
    fVerbose = 1;
    fSumFunc = nullptr;
    fBcgFunc = nullptr;
    fPolyBcgDegree = 2;
    fOnlypositivepeaks = true;
    fChiSquare = std::numeric_limits<double>::quiet_NaN();
}

TheuerkaufFitter *TheuerkaufFitter::SetBacground(std::unique_ptr<TF1> &bcg_fcn) noexcept
{
    fBcgFunc = std::move(bcg_fcn);
    return this;
};

void TheuerkaufFitter::PrintFitResults_plain(std::ostream &os) const
{
    for (auto &peak : fPeaks)
    {
        std::string _id = std::to_string(peak.GetID());
        std::string _pos = std::to_string(peak.GetPos());
        std::string _pos_e = std::to_string(peak.GetPosErr());
        std::string _fwhm = std::to_string(peak.GetFWHM());
        std::string _fwhm_err = std::to_string(peak.GetFWHMErr());
        std::string _vol = std::to_string(peak.GetVol());
        std::string _vol_err = std::to_string(peak.GetVolErr());
        std::string _tl = peak.HasTL() ? std::to_string(peak.GetTL()) : "--";
        std::string _tr = peak.HasTR() ? std::to_string(peak.GetTR()) : "--";
        std::string _sw = peak.HasStep() ? std::to_string(peak.GetSW()) : "--";
        std::string _sh = peak.HasStep() ? std::to_string(peak.GetSH()) : "--";

        os << _id << " " << _pos << " " << _pos_e << " " << _fwhm << " " << _fwhm_err << " " << _vol << " " << _vol_err
           << std::endl;
    }
}

void TheuerkaufFitter::PrintFitResults() const
{
    // using namespace tabulate;
    assert(fSumFunc && fFitResults && "Fit was likely not yet performed!");
    tabulate::Table fit_results;

    std::string _status;

    if (fFitResults->IsValid())
        _status = "Fit OK";
    else
        _status = "Fit not valid, minimizer status " + std::to_string(fFitResults->Status());
    std::string _chi2 = "chi2 " + std::to_string(fChiSquare);
    std::string _rchi2 = "rchi2 " + std::to_string(fChiSquare / fSumFunc.get()->GetNDF());
    std::string _errmat = "error matrix ";

    tabulate::Color _matrix_status;
    switch (fFitResults->CovMatrixStatus())
    {
    case 0:
        _errmat += "not calculated";
        _matrix_status = tabulate::Color::red;
        break;
    case 1:
        _errmat += "approximated";
        _matrix_status = tabulate::Color::yellow;
        break;
    case 2:
        _errmat += "pos def";
        _matrix_status = tabulate::Color::yellow;
        break;
    case 3:
        _errmat += "accurate";
        _matrix_status = tabulate::Color::green;
        break;
    default:
        _errmat += "unknown";
        _matrix_status = tabulate::Color::red;
    }

    std::string _hname = fTempHist != nullptr ? fTempHist->GetName() : "unknown hist";

    fit_results.add_row({_hname, _status, _chi2, _rchi2, _errmat});
    // fit_results.format()
    //            .border_top("*")
    //            .border_bottom("-")
    //            .border_left(" ")
    //            .border_right(" ");
    fit_results.row(0)
        .format()
        .background_color(tabulate::Color::white)
        .font_style({tabulate::FontStyle::bold})
        .font_align(tabulate::FontAlign::center)
        .background_color(tabulate::Color::grey)
        .font_color(tabulate::Color::white);
    if (fFitResults->IsValid())
        fit_results.row(0).format().background_color(tabulate::Color::green);
    else
        fit_results.row(0).format().background_color(tabulate::Color::red);
    fit_results[0][4].format().background_color(_matrix_status);
    fit_results[0][0].format().background_color(tabulate::Color::magenta);

    tabulate::Table peak_parameters;
    // peak_parameters.format()
    //            .border_top("*")
    //            .border_bottom("-")
    //            .border_left(" ")
    //            .border_right(" ");
    peak_parameters.add_row(
        {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Volume", "VolErr", "TL", "TR", "SW", "SH"});
    // check if we indeed have Volumes or just normalizations by checking the binwidth parameter of each peak
    bool binning_set = std::find_if(fPeaks.begin(), fPeaks.end(),
                                    [](const auto &peak) { return peak.GetBinning() <= 0; }) == fPeaks.end();
    if (binning_set)
    {
        peak_parameters.add_row(
            {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Volume", "VolErr", "TL", "TR", "SW", "SH"});
    }
    else
    {
        peak_parameters.add_row(
            {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Normalization", "NormErr", "TL", "TR", "SW", "SH"});
    }

    // bool has_tl = std::find(fPeaks.begin(), fPeaks.end(), [](auto& peak) {return peak.HasTL();}) == fPeaks.end() ?
    // false : true; bool has_tr = std::find(fPeaks.begin(), fPeaks.end(), [](auto& peak) {return peak.HasTR();}) ==
    // fPeaks.end() ? false : true; bool has_step = std::find(fPeaks.begin(), fPeaks.end(), [](auto& peak) {return
    // peak.HasStep();}) == fPeaks.end() ? false : true;

    peak_parameters.row(0)
        .format()
        .font_style({tabulate::FontStyle::bold})
        .font_align(tabulate::FontAlign::center)
        .font_color(tabulate::Color::yellow);

    for (auto &peak : fPeaks)
    {
        std::string _id = std::to_string(peak.GetID());
        std::string _pos = std::to_string(peak.GetPos());
        std::string _pos_e = std::to_string(peak.GetPosErr());
        std::string _fwhm = std::to_string(peak.GetFWHM());
        std::string _fwhm_err = std::to_string(peak.GetFWHMErr());
        std::string _vol = std::to_string(peak.GetVol());
        std::string _vol_err = std::to_string(peak.GetVolErr());
        std::string _tl = peak.HasTL() ? std::to_string(peak.GetTL()) : "--";
        std::string _tr = peak.HasTR() ? std::to_string(peak.GetTR()) : "--";
        std::string _sw = peak.HasStep() ? std::to_string(peak.GetSW()) : "--";
        std::string _sh = peak.HasStep() ? std::to_string(peak.GetSH()) : "--";

        peak_parameters.add_row({_id, _pos, _pos_e, _fwhm, _fwhm_err, _vol, _vol_err, _tl, _tr, _sw, _sh});
    }

    peak_parameters.column(0).format().font_style({tabulate::FontStyle::bold}).font_align(tabulate::FontAlign::center);

    std::cout << std::endl << fit_results << std::endl;
    std::cout << peak_parameters << std::endl << std::endl;
}

TheuerkaufPeak *TheuerkaufFitter::GetPeak(const int id)
{
    auto it = std::find_if(fPeaks.begin(), fPeaks.end(), [&](const TheuerkaufPeak &pk) { return pk.GetID() == id; });
    assert(it != fPeaks.end());

    return &(*it);
}

const int TheuerkaufFitter::GetNPeaks()
{
    return fPeaks.size();
}

inline TheuerkaufFitter *TheuerkaufFitter::SetRange(double min, double max)
{
    fMin = min;
    fMax = max;

    for (auto &peak : fPeaks)
    {
        peak.SetRange(min, max);
    }
    return this;
}

inline TheuerkaufFitter *TheuerkaufFitter::SetBacgroundPoly(unsigned int poly_order) noexcept
{
    fPolyBcgDegree = poly_order;
    return this;
};

inline TheuerkaufFitter *TheuerkaufFitter::RemoveBackground() noexcept
{
    fBcgFunc.reset(nullptr);
    return this;
}

inline TheuerkaufFitter *TheuerkaufFitter::AllowNegativePeaks(bool allow_negative_peaks) noexcept
{
    fOnlypositivepeaks = !allow_negative_peaks;
    return this;
};

const int TheuerkaufFitter::AddPeak(const double position, const bool leftTail, const bool rightTail, const bool step)
{
    auto &pk = fPeaks.emplace_back(TheuerkaufPeak(fMin, fMax, fPeaks.size(), leftTail, rightTail, step));
    pk.SetParameter_Sigma(1., TheuerkaufPeak::ParamState::SAME);
    pk.SetParameter_Position(position, TheuerkaufPeak::ParamState::FREE, fMin, fMax);
    if (leftTail)
        pk.SetParameter_TailLeft(10., TheuerkaufPeak::ParamState::FREE);
    if (rightTail)
        pk.SetParameter_TailRight(10., TheuerkaufPeak::ParamState::FREE);
    if (step)
    {
        pk.SetParameter_StepHeight(0., TheuerkaufPeak::ParamState::FREE);
        pk.SetParameter_StepWidth(0., TheuerkaufPeak::ParamState::FREE);
    }

    return pk.fId;
}

const int TheuerkaufFitter::AddPeak(const TheuerkaufPeak &peak)
{
    TheuerkaufPeak pk(peak);
    pk.fId = fPeaks.size();
    fPeaks.emplace_back(pk);
    return pk.fId;
}

double TheuerkaufFitter::Eval(const double *x, const double *p) const
{
    //! Private: evaluation function for fit
    auto num_params = GetNumParams();

    // Evaluate background function, if it has been given
    double sum = fBcgFunc ? fBcgFunc->Eval(*x) : 0.0;

    // Evaluate internal background
    sum += std::accumulate(std::reverse_iterator<const double *>(p + num_params),
                           std::reverse_iterator<const double *>(p + num_params - fPolyBcgDegree), 0.0,
                           [&x](double bg, double param) { return bg * *x + param; });

    // Evaluate peaks
    return std::accumulate(fPeaks.begin(), fPeaks.end(), sum,
                           [x, p](double _sum, const TheuerkaufPeak &peak) { return _sum + peak.Eval(x, p); });
}

double TheuerkaufFitter::EvalTotalBacground(const double *x, const double *p)
{
    int num_params = GetNumParams();

    // eval bcg fcn
    double sum = fBcgFunc ? fBcgFunc->Eval(*x) : 0.0;
    // eval polynomial bcg
    std::accumulate(std::reverse_iterator<const double *>(p + num_params),
                    std::reverse_iterator<const double *>(p + num_params - fPolyBcgDegree), sum,
                    [&x](double bg, double param) { return bg * *x + param; });

    // eval peak steps
    std::accumulate(fPeaks.begin(), fPeaks.end(), sum,
                    [x, p](double _sum, const TheuerkaufPeak &peak) { return _sum + peak.EvalStep(x, p); });

    return sum;
}

void TheuerkaufFitter::Analyze(TH1 *histAna, TPad *toPrint = nullptr)
{
    TH1 *histo = (TH1 *)histAna->Clone();
    if (!fSumFunc)
        return;
    fTempObjects.clear();
    fTempPeaks.clear();

    int num_params = GetNumParams();

    auto tot_bcg =
        std::make_unique<TF1>(GetFuncUniqueName("total_bcg", this).c_str(), this, &TheuerkaufFitter::EvalTotalBacground,
                              fMin, fMax, num_params, "TheuerkaufFitter", "EvalTotalBacground");
    fSumFunc->SetLineColor(kRed);
    fSumFunc->SetLineWidth(2);
    tot_bcg->SetLineColor(kGray);
    tot_bcg->SetNpx(10000);

    fTempPeaks = fPeaks;
    int j = 0;
    for (auto &temp_peak : fTempPeaks)
    {
        temp_peak.ResetIndexes(true);
        TF1 *fcn = temp_peak.GetFunction();

        fcn->SetParameter(5, 0);
        fcn->SetParameter(6, 0);

        // std::cout << "temp_peak " << fcn->GetParameter(0) << " " << fcn->GetParameter(1) << " "
        // <<fcn->GetParameter(2) << " " << std::endl; std::cout << temp_peak.GetIndex_Sigma() << " " <<
        // temp_peak.GetIndex_Position() << std::endl;

        // double position = temp_peak.GetPos();
        // double sigma = temp_peak.GetSig();
        // double range_min = position - 3 * sigma < fMin ? fMin : position - 3 * sigma;
        // double range_max = position + 3 * sigma > fMax ? fMax : position + 3 * sigma;

        // fcn->SetRange(range_min, range_max);
        fcn->SetLineColor(GetColor(j + 2));
        fcn->SetLineStyle(2);
        // fcn->SetNpx(10000. * (range_max - range_min) / (double)(fMax - fMin));
        fcn->SetNpx(10000.);
        j++;
    }

    //   std::unique_ptr<TCanvas> canvas(new TCanvas("canvas", "TheuerkaufFitter", 800, 800));

    TPad *can;

    if (toPrint == nullptr)
    {
        can = new TCanvas("canvas", "TheuerkaufFitter", 800, 800);
    }
    else
    {
        can = toPrint;
    }

    can->Divide(1, 2, 0, 0);
    can->cd(1);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.02);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.01);

    double yMin = std::numeric_limits<double>::max(), yMax = 0;
    double xMinRange, xMaxRange;

    for (auto p : fTempPeaks)
    {
        auto fP = p.GetFunction();
        if (yMin > fSumFunc->Eval(p.GetPos()))
            yMin = fSumFunc->Eval(p.GetPos()) * 1.1;
        if (yMax < fSumFunc->Eval(p.GetPos()))
            yMax = fSumFunc->Eval(p.GetPos()) * 1.1;
    }

    fSumFunc->GetRange(xMinRange, xMaxRange);

    if (yMin > tot_bcg->Eval(xMinRange))
        yMin = tot_bcg->Eval(xMinRange) * 1.1;
    if (yMin > tot_bcg->Eval(xMaxRange))
        yMin = tot_bcg->Eval(xMaxRange) * 1.1;
    if (yMax < tot_bcg->Eval(xMinRange))
        yMax = tot_bcg->Eval(xMinRange) * 1.1;
    if (yMax < tot_bcg->Eval(xMaxRange))
        yMax = tot_bcg->Eval(xMaxRange) * 1.1;

    histo->GetXaxis()->SetRangeUser(fTempPeaks.front().GetPos() - 6 * fTempPeaks.front().GetSig() - 10,
                                    fTempPeaks.back().GetPos() + 6 * fTempPeaks.back().GetSig() + 10);

    std::unique_ptr<TH1> fit_confidence_histo((TH1 *)histo->Clone());
    std::unique_ptr<TH1> substracted_histo((TH1 *)histo->Clone());
    std::unique_ptr<TH1> sigma_histo_plus((TH1 *)histo->Clone());
    std::unique_ptr<TH1> sigma_histo_minus((TH1 *)histo->Clone());

    histo->GetYaxis()->SetRangeUser(yMin, yMax);
    histo->SetTitle(";;Counts");
    histo->GetXaxis()->SetLabelOffset(1.5);
    histo->Draw("HIST");

    // histo->Draw("HIST");
    fSumFunc->Draw("SAME");
    tot_bcg->Draw("SAME");

    std::for_each(fTempPeaks.begin(), fTempPeaks.end(), [](auto &peak) { peak.GetFunction()->Draw("SAME"); });

    fit_confidence_histo->Reset();
    fFitter->GetConfidenceIntervals(fit_confidence_histo.get(), 0.95);
    histo->GetXaxis()->SetLabelOffset(1);
    substracted_histo->GetXaxis()->SetLabelSize(0.05);
    for (int i = 0; i < histo->GetNbinsX() + 2; i++)
    {
        double low_edge = substracted_histo->GetXaxis()->GetBinLowEdge(i);
        double high_edge = substracted_histo->GetXaxis()->GetBinUpEdge(i);
        double fcn_integral = fSumFunc->Integral(low_edge, high_edge) / (double)histo->GetXaxis()->GetBinWidth(i);
        substracted_histo->SetBinContent(i, histo->GetBinContent(i) - fcn_integral);
        fit_confidence_histo->SetBinContent(i, fit_confidence_histo->GetBinContent(i) - fcn_integral);
    }

    fit_confidence_histo->SetFillColor(kRed);
    fit_confidence_histo->SetFillStyle(3002);

    sigma_histo_plus->Reset();
    sigma_histo_minus->Reset();
    for (int i = 0; i < histo->GetNbinsX() + 2; i++)
    {
        sigma_histo_plus->SetBinContent(i, 2. * sqrt(histo->GetBinContent(i)));
        sigma_histo_minus->SetBinContent(i, -2. * sqrt(histo->GetBinContent(i)));
    }
    sigma_histo_plus->SetLineColor(kBlack);
    sigma_histo_minus->SetLineColor(kBlack);

    double maxY = substracted_histo->GetMaximum() > sigma_histo_plus->GetMaximum()
                      ? substracted_histo->GetMaximum() * 1.1
                      : sigma_histo_plus->GetMaximum() * 1.1;
    double minY = substracted_histo->GetMinimum() < sigma_histo_minus->GetMinimum()
                      ? substracted_histo->GetMinimum() * 1.1
                      : sigma_histo_minus->GetMinimum() * 1.1;
    double range = abs(maxY) > abs(minY) ? abs(maxY) : abs(minY);

    can->cd(2);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.01);
    substracted_histo->SetTitle(";Energy (keV); Counts - Fit");
    substracted_histo->Draw("hist");
    substracted_histo->GetYaxis()->SetRangeUser(minY, maxY);
    fit_confidence_histo->Draw("e3 same");
    sigma_histo_plus->Draw("SAME HIST");
    sigma_histo_minus->Draw("SAME HIST");

    /*
        can->cd(3);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.2);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.01);
        tot_bcg->GetYaxis()->SetRangeUser(tot_bcg->GetMinimum(xMinRange,xMaxRange),tot_bcg->GetMaximum(xMinRange,xMaxRange));
        tot_bcg->Draw("");
    */
    // fTempObjects.emplace_back(std::move(histo));
    // fTempObjects.emplace_back(std::move(fSumFunc));
    fTempObjects.emplace_back(std::move(tot_bcg));
    fTempObjects.emplace_back(std::move(can));
    fTempObjects.emplace_back(std::move(substracted_histo));
    fTempObjects.emplace_back(std::move(fit_confidence_histo));
    fTempObjects.emplace_back(std::move(sigma_histo_plus));
    fTempObjects.emplace_back(std::move(sigma_histo_minus));
}

void TheuerkaufFitter::HandleParameterStates()
{
    if (!fSumFunc)
        return;
    // The fSumFunct has complete set of parameters for every TheuerkaufPeak + polynomial background.
    // Each TheuerkaufPeak knows which parameters it should take (based on the peakID)
    // However, some of thems are not used by the TheuerkaufPeak function, e.i. when NONE is set,
    // but TheuerkaufFitter does not know about this. We need to fix them otherwise minimizer will vary parameters
    // that are not really used, which will slow down minimization process, produce warnings and who knows what else.
    // The parameters with NONE and FIXED are fixed to its values with fSumFunc->FixParameter() to solve this.
    // If SAME state is observed, all parameters of the same type (e.g. sigma) must point to the same index
    // of the TheuerkaufFitter parameter array, while first being set the FREE and rest fixed.
    // For example we have 3 peaks and linear background, the number of Fitter parameters is 21+2=23.
    // The TheuerkaufPeak internally knows which index to take: peak with ID=0 knows, that it in Eval() it
    // should take indexes no, 0-6 of the fitter, peak with ID=1 indexes 7-13 of Fitter etc.
    // Let's say that SAME is set to all sigma parameters. Sigma is third parameter of TheuerkaufPeak,
    // so in Fitter the sigma's will have indexes 2(peakID=0), 9(peakID=1) and 16(peakID=2).
    // So, to apply the "SAME" option, we tell TheuerkaufPeaks with ID 1 and 2, that they should use
    // index 2 for their sigma. The indexes 9 and 16 are set as Fix in the Fitter, thus only index 2 will remain free.
    // This index redirection is done using the SetIndex_Sigma() function of TheuerkaufPeak

    using PAR_STATE = TheuerkaufPeak::ParamState;
    int first_same_index;

    // set volume
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_Volume();
        auto param_index = peak.GetIndex_Volume();
        // NONE is invalid state for volume
        if (state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_Volume(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set POSITION
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_Position();
        auto param_index = peak.GetIndex_Position();
        // NONE and SAME are invalid states for position
        if (state == PAR_STATE::FIXED) // NONE is invalid state for position
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SIGMA
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_Sigma();
        auto param_index = peak.GetIndex_Sigma();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_Sigma(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set TL
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_TailLeft();
        auto param_index = peak.GetIndex_TailLeft();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_TailLeft(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set TR
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_TailRight();
        auto param_index = peak.GetIndex_TailRight();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_TailRight(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SH
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_StepHeight();
        auto param_index = peak.GetIndex_StepHeight();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_StepHeight(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SW
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak.GetState_StepWidth();
        auto param_index = peak.GetIndex_StepWidth();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to fixed
            peak.SetIndex_StepWidth(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
}

void TheuerkaufFitter::DistributeParametersToPeaks()
{
    if (!fSumFunc)
        return;
    for (auto &peak : fPeaks)
    {
        int i_vol = peak.GetIndex_Volume();
        int i_pos = peak.GetIndex_Position();
        int i_sig = peak.GetIndex_Sigma();
        int i_tl = peak.GetIndex_TailLeft();
        int i_tr = peak.GetIndex_TailRight();
        int i_sh = peak.GetIndex_StepWidth();
        int i_sw = peak.GetIndex_StepHeight();

        peak.GetFunction()->SetParameter(0, fSumFunc->GetParameter(i_vol));
        peak.GetFunction()->SetParameter(1, fSumFunc->GetParameter(i_pos));
        peak.GetFunction()->SetParameter(2, fSumFunc->GetParameter(i_sig));
        peak.GetFunction()->SetParameter(3, fSumFunc->GetParameter(i_tl));
        peak.GetFunction()->SetParameter(4, fSumFunc->GetParameter(i_tr));
        peak.GetFunction()->SetParameter(5, fSumFunc->GetParameter(i_sh));
        peak.GetFunction()->SetParameter(6, fSumFunc->GetParameter(i_sw));

        peak.GetFunction()->SetParError(0, fSumFunc->GetParError(i_vol));
        peak.GetFunction()->SetParError(1, fSumFunc->GetParError(i_pos));
        peak.GetFunction()->SetParError(2, fSumFunc->GetParError(i_sig));
        peak.GetFunction()->SetParError(3, fSumFunc->GetParError(i_tl));
        peak.GetFunction()->SetParError(4, fSumFunc->GetParError(i_tr));
        peak.GetFunction()->SetParError(5, fSumFunc->GetParError(i_sh));
        peak.GetFunction()->SetParError(6, fSumFunc->GetParError(i_sw));
    }

    // Filippo
    for (int i = 0; i < fPolyBcgDegree; ++i)
    {
        //        fBcgFunc->SetParameter(i,fSumFunc->GetParameter(fSumFunc->GetNpar()-fPolyBcgDegree-1+i));
        //        fBcgFunc->SetParError(i,fSumFunc->GetParError(fSumFunc->GetNpar()-fPolyBcgDegree-1+i));
    }
}

// std::vector<std::pair<double,double>> TheuerkaufFitter::GetPeaksIntegral(){
//     assert(fPeaks.size()>0);
//     std::vector<std::pair<double,double>> Integrals;
//     if(fPeaks.size()!=0){
//
//         *analyzedHisto
//     }
//
// }

void TheuerkaufFitter::Fit(TH1 *histFit, std::string options)
{

    std::shared_ptr<TH1> ptr((TH1 *)(histFit->Clone()));
    fTempHist = ptr;
    // TH1 *hist = (TH1 *)histFit->Clone();

    assert(fPeaks.size() > 0);

    int num_params = GetNumParams();

    // set binning, becasue peak volume depends on the histogram binning
    double hist_binning = fTempHist->GetXaxis()->GetBinWidth(1);
    if (fVerbose > 2)
        std::cout << "binning " << hist_binning << std::endl;
    std::for_each(fPeaks.begin(), fPeaks.end(), [&hist_binning](auto &peak) { peak.SetBinning(hist_binning); });

    // Create fit function
    fSumFunc = std::make_unique<TF1>(GetFuncUniqueName("f", this).c_str(), this, &TheuerkaufFitter::Eval, fMin, fMax,
                                     num_params, "TheuerkaufFitter", "Eval");
    int nbins = fTempHist->FindBin(fMax) - fTempHist->FindBin(fMin);
    fSumFunc->SetNpx(nbins * 10);
    // set bcg params names
    {
        int polynomial_order = 0;
        for (int i = fPeaks.size() * 7; i < num_params; i++)
        {
            std::string par_name = "bcg_p" + std::to_string(polynomial_order);
            fSumFunc->SetParName(i, par_name.c_str());
            polynomial_order++;
        }
    }

    // *** Initial parameter estimation ***
    int b1 = fTempHist->FindBin(fMin);
    int b2 = fTempHist->FindBin(fMax);
    // Check if any of the peaks contain steps
    bool steps = std::find_if(fPeaks.begin(), fPeaks.end(),
                              [](const TheuerkaufPeak &peak) { return peak.HasStep(); }) != fPeaks.end();
    // If there is internal background, we need to estimate it first. We will
    // estimate that the background is constant at the level of the bin with the
    // lowest content if there are no steps, or constant at the level of the
    // leftmost bin in the fit region otherwise (both after substraction of
    // possible external background).
    // NOTE: we generally assume that the step width is positive, so the step
    // function goes to zero at the far left side of the peak. This seems
    // reasonable, as the step width is usually fixed at 1.0.

    double intBg0 = 0.0;
    if (fPolyBcgDegree >= 1)
    {
        if (steps)
        {
            intBg0 = fTempHist->GetBinContent(b1);
            if (fBcgFunc != nullptr)
            {
                intBg0 -= fBcgFunc->Eval(fTempHist->GetBinCenter(b1));
            }
        }
        else
        {
            intBg0 = std::numeric_limits<double>::infinity();

            if (fBcgFunc != nullptr)
            {
                for (int b = b1; b <= b2; ++b)
                {
                    double bc = fTempHist->GetBinContent(b) - fBcgFunc->Eval(fTempHist->GetBinCenter(b));
                    if (bc < intBg0)
                    {
                        intBg0 = bc;
                    }
                }
            }
            else
            {
                for (int b = b1; b <= b2; ++b)
                {
                    double bc = fTempHist->GetBinContent(b);
                    if (bc < intBg0)
                    {
                        intBg0 = bc;
                    }
                }
            }
        }

        // Set background parameters of sum function
        fSumFunc->SetParameter(num_params - fPolyBcgDegree, intBg0);
        if (fPolyBcgDegree >= 2)
        {
            for (int i = num_params - fPolyBcgDegree + 1; i < num_params; ++i)
            {
                fSumFunc->SetParameter(i, 0.0);
            }
        }
    }
    // Next, we must estimate possible steps in the background. We estimate the
    // sum of all step heights as the difference between the first and the last
    // bin in the region (with background substracted). From this, we substract
    // all step heights that have been fixed, and evenly distribute the difference
    // among the others.
    // For the rest of the initial parameter estimation process, we will generally
    // assume the steps to be sharp (step width zero). This is because the step
    // width depends on the peak width, which will only be estimated in the end.
    double avgFreeStep = 0.0;
    if (steps)
    {
        struct Result
        {
            int nStepFree;
            double sumFixedStep;
        };

        auto result = std::accumulate(fPeaks.begin(), fPeaks.end(), Result{0, 0.0},
                                      [](Result _result, const TheuerkaufPeak &peak) {
                                          if (peak.HasStep())
                                          {
                                              if (peak.GetState_StepHeight() == TheuerkaufPeak::ParamState::FIXED)
                                              {
                                                  _result.sumFixedStep += peak.GetSH();
                                              }
                                              else
                                              {
                                                  ++_result.nStepFree;
                                              }
                                          }
                                          return _result;
                                      });

        double sumStep = fTempHist->GetBinContent(b2) - fTempHist->GetBinContent(b1);
        if (result.nStepFree != 0)
        {
            avgFreeStep = (sumStep - result.sumFixedStep) / result.nStepFree;
        }
    }

    // Estimate peak amplitudes:
    // We assume that the peak positions provided are already a good estimate of
    // the peak centers. The peak amplitude is then estimated as the bin content
    // at the center, with possible external and internal background, and a
    // possible step, substracted. Note that our estimate gets bad if peaks
    // overlap a lot, but it seems hard to do something about this, because we do
    // not know the peak width yet.
    std::vector<double> amps;
    amps.reserve(fPeaks.size());
    double sumAmp = 0.0;

    // First: no steps
    std::transform(fPeaks.begin(), fPeaks.end(), std::back_inserter(amps), [&](const TheuerkaufPeak &peak) {
        double pos = peak.GetPos();
        double amp = fTempHist->GetBinContent(fTempHist->FindBin(pos)) - intBg0;
        if (fBcgFunc)
        {
            amp -= fBcgFunc->Eval(pos);
        }
        sumAmp += amp;
        return amp;
    });

    // Second: include steps
    using PeakVector_t = std::vector<TheuerkaufPeak>;
    using PeakID_t = PeakVector_t::size_type;

    if (steps)
    {
        // Generate a list of peak IDs sorted by position
        std::vector<PeakID_t> sortedPeakIDs(fPeaks.size());
        std::iota(sortedPeakIDs.begin(), sortedPeakIDs.end(), 0);
        std::sort(sortedPeakIDs.begin(), sortedPeakIDs.end(), [&](const PeakID_t &lhs, const PeakID_t &rhs) {
            return fPeaks[lhs].GetPos() < fPeaks[rhs].GetPos();
        });

        struct Sums
        {
            double step, amp;
        };
        auto sums =
            std::accumulate(sortedPeakIDs.begin(), sortedPeakIDs.end(), Sums{0.0, 0.0}, [&](Sums _sums, PeakID_t id) {
                auto &peak = fPeaks[id];
                double curStep = 0.0;
                if (peak.HasStep())
                {
                    curStep = peak.HasStep() ? peak.GetSH() : avgFreeStep;
                }
                amps[id] -= _sums.step + curStep / 2.0;
                _sums.amp -= _sums.step + curStep / 2.0;
                _sums.step += curStep;
                return _sums;
            });
        sumAmp -= sums.amp;
    }

    // Estimate peak parameters
    //
    // Assuming that all peaks in the fit have the same width, their volume is
    // proportional to their amplitude. We thus calculate the total volume as a
    // sum over bin contents (with background substracted) and distribute it among
    // the peaks according to their amplitude. Last, we can estimate the common
    // width from the total volume and the sum of the amplitudes.
    //
    // NOTE: This assumes that the peaks are purely Gaussian, i.e. that there are
    // no tails.

    // First: calculate total volume
    double sumVol = 0.0;
    for (int b = b1; b <= b2; ++b)
    {
        sumVol += fTempHist->GetBinContent(b);
    }
    sumVol -= intBg0 * (b2 - b1 + 1.);
    if (fBcgFunc != nullptr)
    {
        for (int b = b1; b <= b2; ++b)
        {
            sumVol -= fBcgFunc->Eval(fTempHist->GetBinCenter(b));
        }
    }

    if (steps)
    {
        sumVol -= std::accumulate(fPeaks.begin(), fPeaks.end(), 0.0, [&](double sum, const TheuerkaufPeak &peak) {
            if (peak.HasStep())
            {
                double curStep = peak.HasStep() ? avgFreeStep : peak.GetSH();
                int b = fTempHist->FindBin(peak.GetPos());
                sum -= curStep * (b2 - std::min(b, b2) + 0.5);
            }
            return sum;
        });
    }

    // normalize on binning
    sumAmp = sumAmp / fTempHist->GetXaxis()->GetBinWidth(1);

    // Second: calculate average peak width (sigma)
    double avgSigma = std::abs(sumVol / (sumAmp * std::sqrt(2. * M_PI)));
    if (fVerbose > 2)
        std::cout << "sumvVol " << sumVol << " sumAmp " << sumAmp << " " << std::endl;
    if (fVerbose > 2)
        std::cout << "avgSigma " << avgSigma << std::endl;

    // Third: calculate sum of free volumes and amplitudes
    double sumFreeAmp = sumAmp;
    double sumFreeVol = sumVol;
    auto ampIter = amps.begin();
    for (auto &peak : fPeaks)
    {
        if (peak.GetState_Volume() != TheuerkaufPeak::ParamState::FREE)
        {
            sumFreeAmp -= *(ampIter++);
            sumFreeVol -= peak.GetVol();
        }
    }

    // Init fit parameters for peaks

    int par_n = 0;
    // first initialize all parameters to the their peak values
    for (auto &peak : fPeaks)
    {
        for (int i = 0; i < 7; i++)
        {
            double _min, _max;
            peak.GetFunction()->GetParLimits(i, _min, _max);
            fSumFunc->SetParLimits(par_n, _min, _max);

            std::string par_name = peak.GetFunction()->GetParName(i);
            if (peak.GetState(i) == TheuerkaufPeak::ParamState::SAME)
            {
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak.GetID() << " is set to \"SAME\""
                              << std::endl;
                par_name += "_SAME ";
                fSumFunc->FixParameter(par_n, peak.GetFunction()->GetParameter(i));
            }
            else if (peak.GetState(i) == TheuerkaufPeak::ParamState::FIXED)
            {
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak.GetID() << " is set to \"FIXED\""
                              << std::endl;
                par_name += "_FIXED";
                fSumFunc->FixParameter(par_n, peak.GetFunction()->GetParameter(i));
            }
            else if (peak.GetState(i) == TheuerkaufPeak::ParamState::NONE)
            {
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak.GetID() << " is set to \"NONE\""
                              << std::endl;
                par_name += "_NONE ";
                fSumFunc->FixParameter(par_n, peak.GetFunction()->GetParameter(i));
            }
            else
            {
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak.GetID() << " is set to \"FREE\""
                              << std::endl;
                fSumFunc->SetParameter(par_n, peak.GetFunction()->GetParameter(i));
            }
            fSumFunc->SetParName(par_n, par_name.c_str());
            par_n++;
        }
    }
    this->HandleParameterStates();

    ampIter = amps.begin();
    for (auto &peak : fPeaks)
    {
        int i_vol = peak.GetIndex_Volume();
        int i_pos = peak.GetIndex_Position();
        int i_sig = peak.GetIndex_Sigma();
        int i_tl = peak.GetIndex_TailLeft();
        int i_tr = peak.GetIndex_TailRight();
        int i_sh = peak.GetIndex_StepWidth();
        int i_sw = peak.GetIndex_StepHeight();

        if (fVerbose > 1)
            std::cout << "peak " << peak.GetID() << " i_vol " << i_vol << " i_pos " << i_pos << " i_sig " << i_sig
                      << " i_tl " << i_tl << " i_tr " << i_tr << " i_sw " << i_sw << " i_sh " << i_sh << std::endl;

        double amp = *(ampIter++);

        if (fOnlypositivepeaks)
        {
            // The root fit algorithm produces strange results when fitting with extreme limits/boundary conditions,
            // hence set some sane upper limits
            fSumFunc->SetParameter(i_vol, std::max(sumFreeVol * amp / sumFreeAmp, 1.));
            fSumFunc->SetParLimits(i_vol, 0., std::max(100. * sumVol, 0.) + 1e9);
        }
        else
        {
            fSumFunc->SetParameter(i_vol, std::max(sumFreeVol * amp / sumFreeAmp, 1.));
        }

        fSumFunc->SetParameter(i_sig, avgSigma);
        fSumFunc->SetParLimits(i_sig, 0., 10 * (fMax - fMin) + 1e3);

        if (peak.HasTL())
        {
            fSumFunc->ReleaseParameter(i_tl);
            fSumFunc->SetParameter(i_tl, avgSigma * 3);
            fSumFunc->SetParLimits(i_tl, 0, 1E9);
            // fSumFunc->SetParLimits(i_tl, 0, 1e9);
        }
        if (peak.HasTR())
        {
            fSumFunc->ReleaseParameter(i_tr);
            fSumFunc->SetParameter(i_tr, avgSigma * 3);
            fSumFunc->SetParLimits(i_tr, 0, 1E9);
            // fSumFunc->SetParameter(i_tr, peak.GetTR());
            // fSumFunc->SetParLimits(i_tr, 0, 100);
        }
        if (peak.HasStep())
        {
            fSumFunc->ReleaseParameter(i_sh);
            fSumFunc->ReleaseParameter(i_sw);
            fSumFunc->SetParameter(i_sh, avgFreeStep / (amp * M_PI));
            // fSumFunc->SetParLimits(i_sh, -10, 10);
            // fSumFunc->SetParLimits(i_sw, -10, 10);
        }
    }

    if (fVerbose > 1)
        std::cout << "Sum function has " << fSumFunc->GetNpar() << " from which " << fSumFunc->GetNumberFreeParameters()
                  << " are set free" << std::endl;

    // Now, do the fit
    std::string fit_options = "0NRMIS";
    std::transform(options.begin(), options.end(), options.begin(), ::toupper);
    if (options.find("OUTPUT_ROOT") == std::string::npos)
        fit_options += "Q";
    if ((options.find("LIKELIHOOD") != std::string::npos) || options.find("POISSON") != std::string::npos)
        fit_options += "L";

    fFitResults = fTempHist->Fit(fSumFunc.get(), fit_options.c_str());

    fFitter = TVirtualFitter::GetFitter();
    this->DistributeParametersToPeaks();
    // Store Chi^2
    fChiSquare = fSumFunc->GetChisquare();
    if (fVerbose != 0 && options.find("OUTPUT_NONE") == std::string::npos ||
        options.find("OUTPUT_STANDARD") != std::string::npos)
        this->PrintFitResults();
    if (fVerbose != 0 && options.find("OUTPUT_PLAIN") != std::string::npos)
        this->PrintFitResults_plain();
}
