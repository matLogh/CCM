/*
list of throw exceptions:

1 = vectors fed to the constructor are not the same size
2 = normalization of sample vector error = its field of zeroes
3 = vectors in dot product are not the same size
*/
#include <algorithm> // copy
#include <atomic>
#include <bits/stdc++.h>
#include <fstream> //open file
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept> // std::out_of_range
#include <stdio.h>   //strcat
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include <TFile.h>
#include <TGraphSmooth.h>
#include <TTree.h>
#pragma GCC diagnostic pop

#include "CCM.h"
#include "Cross_correlation.h"
#include "variables.h"

std::mutex        TEC::CCM::fMtx_ROOTfit;
const std::string TEC::CCM::EMPTY_FUNCTION_NAME{"EMPTY_FUNCTION"};

TEC::CCM::CCM(std::shared_ptr<TH2>                matrix,
              const std::vector<RegionOfInterest> _ROIs,
              const double                        reference_time_low,
              const double                        reference_time_high)
    : CCM(matrix, _ROIs)
{
    // create sample_vector
    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        this->CreateReferenceVector(ROI_index, reference_time_low, reference_time_high);
    }
}

TEC::CCM::CCM(std::shared_ptr<TH2> matrix, const std::vector<RegionOfInterest> _ROIs)
{

    // V.TEMAT = dynamic_cast<TH2D *>(matrix.Clone(Form("%s_clone", matrix.GetName())));
    V.TEMAT = matrix;
    if (!V.TEMAT) { throw std::runtime_error("Error! Input matrix is null!"); }
    fXbins = static_cast<size_t>(V.TEMAT->GetXaxis()->GetNbins());
    fYbins = static_cast<size_t>(V.TEMAT->GetYaxis()->GetNbins());
    V.ROIs = _ROIs;

    // V.ROIs = _ROIs;
    V.total_tasks = fXbins * (V.ROIs.size());
    V.time_bins   = fXbins;

    // reserve space for result container structure
    ResVec = new ResCont *[V.ROIs.size()];
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        ResVec[roi_index] = new ResCont[V.time_bins];
    }
    this->CopyMatrixContent();
}

TEC::CCM::CCM(const CCM &other)
    : fXbins(other.fXbins), fYbins(other.fYbins), fFixedTEMAT(other.fFixedTEMAT),
      fCorrectionFunctions(other.fCorrectionFunctions),
      fForceRebuildInterpolators(other.fForceRebuildInterpolators),
      fThreadTask(other.fThreadTask.load()), fNthreads(other.fNthreads), V(other.V),
      fFitDone(other.fFitDone), fCorrectionFits(other.fCorrectionFits)
{
    std::cout << "copy constructor called" << std::endl;
    if (other.ResVec)
    {
        ResVec = new ResCont *[V.ROIs.size()];
        for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
        {
            ResVec[roi_index] = new ResCont[V.time_bins];
        }
    }
    else { ResVec = nullptr; }
    this->CopyMatrixContent();
}

TEC::CCM::CCM(CCM &&other)
    : fXbins(other.fXbins), fYbins(other.fYbins),
      fFixedTEMAT(std::move(other.fFixedTEMAT)),
      fCorrectionFunctions(std::move(other.fCorrectionFunctions)),
      fForceRebuildInterpolators(other.fForceRebuildInterpolators),
      fThreadTask(other.fThreadTask.load()), fNthreads(other.fNthreads),
      V(std::move(other.V)), fFitDone(other.fFitDone),
      fCorrectionFits(std::move(other.fCorrectionFits)),
      ResVec(other.ResVec) // Transfer ownership of ResVec
{
    std::cout << "move constructor called" << std::endl;
    // Transfer ownership of ResVec
    other.ResVec     = nullptr;
    other.V.TEMATarr = nullptr;
}

TEC::CCM::~CCM()
{
    if (ResVec)
    {
        for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
        {
            delete[] ResVec[roi_index];
        }
        delete[] ResVec;
        ResVec == nullptr;
    }

    if (V.TEMATarr)
    {
        for (auto roi_index = 0; roi_index < static_cast<int>(V.ROIs.size()); roi_index++)
        {
            for (size_t j = 0; j < fXbins; j++) { delete[] V.TEMATarr[roi_index][j]; }
            delete[] V.TEMATarr[roi_index];
        }
        delete[] V.TEMATarr;
        V.TEMATarr = nullptr;
    }

    for (auto &f : fCorrectionFunctions)
    {
        delete f.first;
        f.first = nullptr;
    }
}

void TEC::CCM::SetCorrectionFunction(const TF1 &fcn, const std::string &fit_options)
{
    if (!fcn.GetFormula())
    {
        throw std::runtime_error("Error! Function must be created using "
                                 "TFormula, e.i. using the text string, "
                                 "otherwise it cannot be saved in a ROOT file "
                                 "and used in the CCM object!");
    }
    if (fCorrectionFunctions.size() > 0)
    {
        throw std::runtime_error("Error! Primary function needs to be set "
                                 "first, only one primary function can be set");
    }

    // add empty function for cases when no ROI is valid - this essentially does
    // nothing
    if (fCorrectionFunctions.empty())
    {
        TF1 *empty_function =
            new TF1("EMPTY_FUNCTION", "x", V.TEMAT->GetYaxis()->GetBinLowEdge(1),
                    V.TEMAT->GetYaxis()->GetBinUpEdge(V.TEMAT->GetYaxis()->GetNbins()));
        fCorrectionFunctions.emplace_back(std::make_pair(empty_function, "NQ"));
    }

    TF1 *fclone = (TF1 *)fcn.Clone();

    fclone->SetRange(V.TEMAT->GetYaxis()->GetBinLowEdge(1),
                     V.TEMAT->GetYaxis()->GetBinUpEdge(V.TEMAT->GetYaxis()->GetNbins()));
    fclone->SetNpx(V.TEMAT->GetYaxis()->GetNbins() * 100);
    fclone->Update();
    fCorrectionFunctions.emplace_back(std::make_pair(fclone, fit_options + "NQ"));

    std::sort(fCorrectionFunctions.begin(), fCorrectionFunctions.end(),
              [](const auto &a, const auto &b) {
                  return a.first->GetNpar() > b.first->GetNpar();
              });
}

std::unique_ptr<TGraph> TEC::CCM::GetROIShifts(const size_t roi_index,
                                               const bool   valid_only)
{
    std::unique_ptr<TGraph> gr = std::make_unique<TGraph>();
    gr->SetBit(TGraph::kIsSortedX);
    gr->SetName(Form("shift_ROI_%li", roi_index));
    gr->SetTitle(Form("shifts for ROI %li;Time;Energy shift", roi_index));

    for (size_t time = 0; time < V.time_bins; time++)
    {
        if (valid_only && !ResVec[roi_index][time].isValid) { continue; }
        gr->AddPoint(GetMatrixTime(time), ResVec[roi_index][time].energy_shift);
    }
    return gr;
}

std::unique_ptr<TGraph> TEC::CCM::GetShiftProfile(const int  time_bin,
                                                  const bool valid_only)
{
    std::unique_ptr<TGraph> gr = std::make_unique<TGraph>();
    gr->SetBit(TGraph::kIsSortedX);
    gr->SetTitle(Form("%i", time_bin));
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        if (ResVec[roi_index][time_bin].isValid)
        {
            if (valid_only && !ResVec[roi_index][time_bin].isValid) { continue; }
            gr->AddPoint(V.ROIs[roi_index].desired_energy,
                         ResVec[roi_index][time_bin].bin_shift);
        }
    }
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(20);
    return gr;
}

void TEC::CCM::CopyMatrixContent()
{
    // V.TEMAT = matrix;
    // fXbins = V.TEMAT->GetXaxis()->GetNbins();
    // V.time_bins = fXbins;

    V.TEMATarr = new float **[V.ROIs.size()];

    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        V.TEMATarr[ROI_index] = new float *[fXbins];
        for (size_t j = 0; j < fXbins; j++)
        {
            V.TEMATarr[ROI_index][j] = new float[V.ROIs[ROI_index].displacement_range];
        }
    }

    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        for (int x = 0; x < static_cast<int>(fXbins); x++)
        {
            for (int y = 0; y < V.ROIs[ROI_index].displacement_range; y++)
            {
                V.TEMATarr[ROI_index][x][y] = static_cast<float>(V.TEMAT->GetBinContent(
                    x + 1, y + V.ROIs[ROI_index].bin_window_low +
                               V.ROIs[ROI_index].bin_displacement_low)); // y+ 1?? - check
            }
        }
    }
}

void TEC::CCM::CreateReferenceVector(const uint   ROI_index,
                                     const double sample_time_low,
                                     const double sample_time_high)
{
    assert(static_cast<size_t>(ROI_index) < V.ROIs.size());
    const size_t       vec_size = static_cast<size_t>(V.ROIs[ROI_index].vector_dimension);
    std::vector<float> vec(vec_size, 0.0);

    size_t vector_iterator{0};
    // loop over energy

    int tbin_start = V.TEMAT->GetXaxis()->FindBin(sample_time_low);
    int tbin_end   = V.TEMAT->GetXaxis()->FindBin(sample_time_high);

    for (int e_bin = V.ROIs[ROI_index].bin_window_low;
         e_bin < V.ROIs[ROI_index].bin_window_high; e_bin++, vector_iterator++)
    {
        // loop over time, since sample vector should be formed of more than 1
        // time-bin width of TEMAT
        for (int t_bin = tbin_start; t_bin <= tbin_end; t_bin++)
        {
            vec[vector_iterator] +=
                static_cast<float>(V.TEMAT->GetBinContent(t_bin, e_bin));
        }
    }

    this->Normalize(vec);

    V.sample_vector.emplace_back(std::move(vec));
}

void TEC::CCM::Normalize(std::vector<float> &v)
{
    double norm = 0;
    for (uint i = 0; i < v.size(); i++) { norm += (v[i] * v[i]); }
    if (norm <= 0)
    {
        throw std::runtime_error(
            "Normalization ERROR, sample vector cannot consist of zeros!");
    }

    norm = 1. / (double)sqrt(norm);
    for (uint i = 0; i < v.size(); i++) { v[i] = v[i] * norm; }
    // double norm = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    // if (norm <= 0)
    // {
    //     throw std::runtime_error(
    //         "Normalization ERROR, sample vector cannot consist of zeros!");
    // }
    // norm = sqrt(norm);
    // std::transform(v.begin(), v.end(), v.begin(),
    //                [norm](const float &x) { return x / (float)norm; });
}

void TEC::CCM::CalculateEnergyShifts(const unsigned int threads)
{
    this->CheckReferenceVectors();
    fNthreads = threads;

    std::vector<std::thread> t;
    t.reserve(fNthreads);

    // std::cout << "***** Starting calculations of Cross-correlation corrections  *****"
    //           << std::endl;

    fThreadTask = 0; // so that the first job start at zero

    std::vector<std::unique_ptr<CrossCorrel>> correlation_object(fNthreads);
    for (auto &co : correlation_object)
    {
        // correlation_object[i] = new CrossCorrel(&V, ResVec);
        co = std::make_unique<CrossCorrel>(&V, ResVec);
    }

    std::mutex mtx_task;
    // std::mutex mtx_fit;
    // start threads
    // std::cout << "Total tasks: " << V.total_tasks << std::endl;
    if (fNthreads > static_cast<uint>(V.total_tasks))
    {
        for (size_t i = 0; i < V.total_tasks; i++)
        {
            // pass thread task, VarManger variable, mutex and ResVec as
            // reference
            t.push_back(std::thread(&CrossCorrel::Process, correlation_object[i].get(),
                                    i + 1, &fThreadTask, std::ref(mtx_task),
                                    std::ref(fMtx_ROOTfit)));
        }
    }
    else
    {
        for (size_t i = 0; i < fNthreads; i++)
        {
            // pass thread task, VarManger variable, mutex and ResVec as
            // reference
            t.push_back(std::thread(&CrossCorrel::Process, correlation_object[i].get(),
                                    i + 1, &fThreadTask, std::ref(mtx_task),
                                    std::ref(fMtx_ROOTfit)));
        }
    }
    for (auto &th : t) th.join();
    // std::cout << "Progress: 100%                " << std::endl;

    this->BuildInterpolators();
}

void TEC::CCM::SaveToRootFile(const std::string &outroot_file)
{
    TFile file(outroot_file.c_str(), "RECREATE");
    {
        TH2 *mat = dynamic_cast<TH2 *>(V.TEMAT->Clone("input_CCM_matrix"));
        mat->Write();
        if (fFixedTEMAT) fFixedTEMAT->Write();
    }
    for (uint roi = 0; roi < V.ROIs.size(); roi++)
    {
        ResCont results;
        double  time;

        // TTree *t = new TTree(Form("ROI_%i", roi), Form("CCM tree of ROI %i",
        // roi));
        TTree *t = new TTree(Form("ROI_%i", roi), Form("CCM tree of ROI %i", roi));
        t->Branch("fit_valid", &results.isValid);
        t->Branch("bin_shift", &results.bin_shift);
        t->Branch("energy_shift", &results.energy_shift);
        t->Branch("dot_product", &results.dp);
        t->Branch("dp_vector", &results.dp_vec);
        t->Branch("poly_shift", &results.poly_shift);
        t->Branch("gfit_chi2", &results.gfit_chi2);
        t->Branch("gfit_sigma", &results.gfit_sigma);
        t->Branch("gfit_mu", &results.gfit_sigma);
        t->Branch("time", &time);

        for (size_t time_bin = 0; time_bin < V.time_bins; time_bin++)
        {
            time    = this->GetMatrixTime(time_bin);
            results = ResVec[roi][time_bin];
            t->Fill();
        }
        file.cd();
        t->Write();
    }

    {
        std::string         fcn_name;
        std::vector<double> fcn_coefs;
        double              time;
        TTree              *t =
            new TTree("fit_coef", "Fit coefficients for given time, parameters = 0.");

        t->Branch("time", &time);
        t->Branch("fcn_name", &fcn_name);
        t->Branch("fcn_coef", &fcn_coefs);

        for (size_t time_bin = 0; time_bin < V.time_bins; time_bin++)
        {
            time     = this->GetMatrixTime(time_bin);
            auto fit = this->GetCorrectionFit(time);
            // first zero out params
            fcn_name  = fit.functionUsed;
            fcn_coefs = fit.coef;
            t->Fill();
        }
        file.cd();
        t->Write();
    }
    file.mkdir("fit_fcns");
    file.cd("fit_fcns");
    for (unsigned int i = 0; i < fCorrectionFunctions.size(); i++)
    {
        fCorrectionFunctions[i].first->Write();
    }
    // save coefficients to file
    file.Close();
}

void TEC::CCM::SaveShiftTable(const std::string &table_filename)
{

    std::ofstream output;
    output.open(table_filename);
    // table head
    output << "# number of ROIs (regions of interest) in the matrix \n";
    output << V.ROIs.size() << '\n';
    output << "# desired energies for each ROI \n";
    for (const auto &r : V.ROIs)
    {
        output << std::setprecision(15) << r.desired_energy << " ";
    }
    output << '\n';
    output << "################ \n";
    output << "# time_bin time_start time_end ";
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        output << "ROI" << roi_index << "_ok(bool) "
               << "ROI" << roi_index << "_shift(double) ";
    }
    output << "\n";

    for (int time_slice_index = 0; time_slice_index < (int)V.time_bins;
         time_slice_index++)
    {
        output << time_slice_index + 1 << std::setprecision(15) << '\t';
        output << std::setprecision(15)
               << V.TEMAT->GetXaxis()->GetBinLowEdge(time_slice_index + 1) << '\t';
        output << std::setprecision(15)
               << V.TEMAT->GetXaxis()->GetBinUpEdge(time_slice_index + 1) << '\t';
        for (size_t ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
        {
            // check if shift is converted to energy units
            output << '\t' << ResVec[ROI_index][time_slice_index].isValid << '\t'
                   << ResVec[ROI_index][time_slice_index].energy_shift;
        }
        output << '\n';
    }
}

void TEC::CCM::BuildInterpolator(const size_t ROI_index)
{
    fFitDone           = false;
    auto *interpolator = &V.ROIs.at(ROI_index).interpolator;
    interpolator->Reset();

    // fill with values
    for (size_t time_index = 0; time_index < V.time_bins; time_index++)
    {
        interpolator->AddPoint(this->GetMatrixTime(time_index),
                               ResVec[ROI_index][time_index].energy_shift,
                               ResVec[ROI_index][time_index].isValid);
    }
}

void TEC::CCM::BuildInterpolators()
{

    fFitDone = false;

    for (size_t ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        this->BuildInterpolator(ROI_index);
    }
    fForceRebuildInterpolators = false;
}

void TEC::CCM::ConfigureShiftInterpolator(const size_t      ROI_index,
                                          const std::string type,
                                          const bool        valid_only)
{
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("Error! ROI index out of range!");
    }
    V.ROIs.at(ROI_index).interpolator.SetType(type, valid_only);
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::ConfigureShiftInterpolator(const size_t                          ROI_index,
                                          const ROOT::Math::Interpolation::Type type,
                                          const bool valid_only)
{
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("Error! ROI index out of range!");
    }

    V.ROIs.at(ROI_index).interpolator.SetType(type, valid_only);
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::ConfigureShiftInterpolator(const std::string type, const bool valid_only)
{
    for (auto &roi : V.ROIs) { roi.interpolator.SetType(type, valid_only); }
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::ConfigureShiftInterpolator(const ROOT::Math::Interpolation::Type type,
                                          const bool valid_only)
{
    for (auto &roi : V.ROIs) { roi.interpolator.SetType(type, valid_only); }
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::DisableInterpolation()
{
    for (auto &roi : V.ROIs) { roi.interpolator.DisableInterpolation(); }
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::DisableInterpolation(const size_t ROI_index)
{
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("Error! ROI index out of range!");
    }
    V.ROIs.at(ROI_index).interpolator.DisableInterpolation();
    fFitDone                   = false;
    fForceRebuildInterpolators = true;
}

void TEC::CCM::EnableInterpolation(const size_t ROI_index)
{
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("Error! ROI index out of range!");
    }

    V.ROIs.at(ROI_index).interpolator.EnableInterpolation();
    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

void TEC::CCM::EnableInterpolation()
{
    for (auto &roi : V.ROIs) { roi.interpolator.EnableInterpolation(); }
    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

const TEC::FitCont TEC::CCM::GetCorrectionFit(const double time)
{
    if (fForceRebuildInterpolators) { this->BuildInterpolators(); }

    if (fCorrectionFits.find(time) != fCorrectionFits.end())
    {
        return fCorrectionFits.at(time);
    }

    FitCont fit_result;

    std::vector<double> x(V.ROIs.size());
    std::vector<double> y(V.ROIs.size());

    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        auto *interpolator = &V.ROIs[ROI_index].interpolator;
        if (interpolator->IsValueValid(time))
        {
            x.emplace_back(V.ROIs[ROI_index].desired_energy);
            y.emplace_back(V.ROIs[ROI_index].desired_energy - interpolator->Eval(time));
        }
    }

    std::pair<TF1 *, std::string> *fcn = this->FindCorrectionFunction((int)x.size());

    // this might be confusing, but sometimes we can have a time where no ROI is valid -
    // e.g. when we have a gap in the data (e.g on run change)
    if (strcmp(fcn->first->GetName(), EMPTY_FUNCTION_NAME.c_str()) == 0)
    {
        fit_result.functionUsed = EMPTY_FUNCTION_NAME;
        return fit_result;
    }
    TGraph gr((int)x.size(), x.data(), y.data());
    gr.Fit(fcn->first, fcn->second.c_str());
    fit_result.functionUsed = fcn->first->GetName();
    fit_result.coef.reserve((size_t)fcn->first->GetNpar());
    for (int i = 0; i < fcn->first->GetNpar(); i++)
    {
        fit_result.coef.emplace_back(fcn->first->GetParameter(i));
    }

    return fit_result;
}

// void TEC::CCM::CalculateCorrectionFits(int time_subdivision)
// {
//     fCorrectionFits.clear();
//     if (V.ROIs.size() == 0)
//     {
//         throw std::runtime_error("No ROIs added to the CCM object");
//     }
//     if (time_subdivision < 1)
//     {
//         throw std::runtime_error("Error! Time subdivision must be at least 1!");
//     }

//     const double step = V.TEMAT->GetXaxis()->GetBinWidth(1) /
//     (double)(time_subdivision);
//     // const double step = V.TEMAT->GetXaxis()->GetBinWidth(1) /
//     // (double)(time_subdivision+1);
//     double time;
//     for (size_t time_index = 0; time_index < V.time_bins; time_index++)
//     {
//         time = V.TEMAT->GetXaxis()->GetBinCenter((int)time_index + 1);
//         fCorrectionFits.insert_or_assign(time, this->GetCorrectionFit(time));

//         time = V.TEMAT->GetXaxis()->GetBinCenter((int)time_index + 1) - step;
//         while (time > V.TEMAT->GetXaxis()->GetBinLowEdge((int)time_index + 1))
//         {
//             fCorrectionFits.insert_or_assign(time, this->GetCorrectionFit(time));
//             time -= step;
//         }

//         time = V.TEMAT->GetXaxis()->GetBinCenter((int)time_index + 1) + step;

//         while (time < V.TEMAT->GetXaxis()->GetBinUpEdge((int)time_index + 1))
//         {
//             fCorrectionFits.insert_or_assign(time, this->GetCorrectionFit(time));
//             time += step;
//         }
//     }
//     fFitDone = true;
// }

void TEC::CCM::CalculateCorrectionFits(int time_subdivision)
{
    fCorrectionFits.clear();
    if (V.ROIs.size() == 0)
    {
        throw std::runtime_error("No ROIs added to the CCM object");
    }
    if (time_subdivision < 1)
    {
        throw std::runtime_error("Error! Time subdivision must be at least 1!");
    }

    const double step = V.TEMAT->GetXaxis()->GetBinWidth(1) / (double)(time_subdivision);
    double       center;

    for (size_t time_index = 0; time_index < V.time_bins; time_index++)
    {
        center = V.TEMAT->GetXaxis()->GetBinLowEdge((int)time_index + 1) + step / 2.0;
        while (center < V.TEMAT->GetXaxis()->GetBinUpEdge((int)time_index + 1))
        {
            fCorrectionFits.insert_or_assign(center, this->GetCorrectionFit(center));
            center += step;
        }
    }
    fFitDone = true;
}

void TEC::CCM::SaveFitTable(const std::string &data_filename,
                            const std::string &detector_name)
{
    if (!fFitDone)
    {
        // std::cerr << "Warning! PerformFits() function was not yet called, it "
        //              "will be called with default parameters "
        //              "(valid only, no spline)"
        //           << std::endl;
        this->CalculateCorrectionFits(1);
    }

    std::ofstream output;
    output.open(data_filename);
    output << "# detector_name\n" << detector_name << '\n';
    output << "# Timestamp range covered by this file\n";
    output << std::setprecision(20) << V.TEMAT->GetXaxis()->GetBinLowEdge(1) << " "
           << std::setprecision(20) << V.TEMAT->GetXaxis()->GetBinUpEdge((int)V.time_bins)
           << '\n';
    output << "# number of functions used \n";
    output << fCorrectionFunctions.size() << '\n';
    output << "# fcn_name number_of_parameters function_equation\n";
    for (const auto &f : fCorrectionFunctions)
    {
        output << f.first->GetName() << "\t" << f.first->GetNpar() << "\t"
               << f.first->GetExpFormula() << '\n';
    }
    output << "#############################################################" << '\n';
    output << "# TS_start TS_stop fcn_name par0 par1 ..." << '\n';

    for (size_t time_index = 0; time_index < V.time_bins; time_index++)
    {
        auto        time = this->GetMatrixTime(time_index - 1);
        const auto &fit  = fCorrectionFits[time];
        // this->GetCorrectionFit(V.TEMAT->GetXaxis()->GetBinCenter(time_index
        // + 1));
        output << V.TEMAT->GetXaxis()->GetBinLowEdge((int)time_index + 1) << '\t'
               << V.TEMAT->GetXaxis()->GetBinUpEdge((int)time_index + 1) << '\t';
        output << fit.functionUsed << '\t' << fit.coef.size() << '\t';
        for (const auto c : fit.coef) { output << c << '\t'; }
        output << '\n';
    }
}

std::shared_ptr<TH2> TEC::CCM::FixMatrix()
{

    if (!fFitDone) { this->CalculateCorrectionFits(1); }

    auto clone =
        dynamic_cast<TH2 *>(V.TEMAT->Clone(Form("%s_corrected", V.TEMAT->GetName())));
    if (!clone) { throw std::runtime_error("Error: Cloning TEMAT failed!"); }

    fFixedTEMAT = std::shared_ptr<TH2>(clone);
    fFixedTEMAT->Reset();

    Double_t bin_cont;
    Double_t total_ratio, en_low_edge, en_up_edge, new_en_low_edge, new_en_up_edge;

    Int_t          global_bin, bin_start, bin_end;
    const TAxis   *axis      = V.TEMAT->GetYaxis();
    const Double_t bin_width = axis->GetBinWidth(1);

    for (int time_index = 1; time_index <= V.TEMAT->GetXaxis()->GetNbins(); time_index++)
    {
        auto time = this->GetMatrixTime((size_t)(time_index - 1));
        if (fCorrectionFits.find(time) == fCorrectionFits.end())
        {
            throw std::runtime_error("Error! Fit for time " + std::to_string(time) +
                                     " not found!");
        }
        const auto &fit = fCorrectionFits[time];
        // load function with parameters
        TF1 *fcn{nullptr};
        for (size_t i = 0; i < fCorrectionFunctions.size(); i++)
        {
            if (fit.functionUsed == std::string(fCorrectionFunctions[i].first->GetName()))
            {
                fcn = fCorrectionFunctions[i].first;
                break;
            }
        }
        // this cannot happen, function is always at lest "EMPTY_FUNCTION"!
        if (fcn == nullptr)
        {
            throw std::runtime_error("No function found for time slice " +
                                     std::to_string(time_index) + " with name " +
                                     fit.functionUsed);
        }
        if (strcmp(fcn->GetName(), EMPTY_FUNCTION_NAME.c_str()) == 0) continue;
        for (unsigned int i = 0; i < fit.coef.size(); i++)
        {
            fcn->SetParameter((int)i, fit.coef[i]);
        }

        // set bins
        for (int en_bin = 1; en_bin <= axis->GetNbins(); en_bin++)
        {
            bin_cont = V.TEMAT->GetBinContent(V.TEMAT->GetBin(time_index, en_bin));
            if (bin_cont == 0.) continue;
            // we need to get the new corrected energies that are matching low
            // and up edges of the old bin
            en_low_edge     = axis->GetBinLowEdge(en_bin);
            en_up_edge      = axis->GetBinUpEdge(en_bin);
            new_en_low_edge = fcn->Eval(en_low_edge);
            new_en_up_edge  = fcn->Eval(en_up_edge);

            // the new_en will not cover range that might be multiple bins in
            // width, and there is also a difference in coverage - new low edge
            // will cover only partially the new bin, same for up edge so we
            // need to get ratios on how to distribute the content of the old
            // bin to the new bins

            bin_start = axis->FindBin(new_en_low_edge);
            bin_end   = axis->FindBin(new_en_up_edge);
            // in case we are covering just 1 bin, put the content in
            if (bin_start == bin_end)
            {
                global_bin = fFixedTEMAT->GetBin(time_index, bin_start);
                fFixedTEMAT->AddBinContent(global_bin, bin_cont);
                continue;
            }
            // get to total ratio of "bin equivalent widths" we are covering
            // with the new energy
            total_ratio = static_cast<double>(bin_end - bin_start - 1);
            total_ratio += (double)(axis->GetBinUpEdge(bin_start) - new_en_low_edge) /
                           (double)bin_width;
            total_ratio += (new_en_up_edge - axis->GetBinLowEdge(bin_end)) / bin_width;
            // set manually to new lower bin content
            {
                global_bin = fFixedTEMAT->GetBin(time_index, bin_start);
                fFixedTEMAT->AddBinContent(
                    global_bin,
                    (bin_cont * (axis->GetBinUpEdge(bin_start) - new_en_low_edge) /
                     bin_width) /
                        total_ratio);
            }
            // set manually to new upper bin content
            {
                global_bin = fFixedTEMAT->GetBin(time_index, bin_end);
                fFixedTEMAT->AddBinContent(
                    global_bin,
                    (bin_cont * (new_en_up_edge - axis->GetBinLowEdge(bin_end)) /
                     bin_width) /
                        total_ratio);
            }
            // in case there are more bins in between, fill them here (skipping
            // first already bin_start and bin_end)
            for (int ybin = bin_start + 1; ybin < bin_end; ybin++)
            {
                global_bin = fFixedTEMAT->GetBin(time_index, ybin);
                fFixedTEMAT->AddBinContent(global_bin, bin_cont / total_ratio);
            }
        }
    }

    fFixedTEMAT->SetEntries(V.TEMAT->GetEntries());
    fFixedTEMAT->ResetStats();
    return fFixedTEMAT;
}

std::shared_ptr<TH2> TEC::CCM::FixMatrix(const TH2 *input_mat)
{
    if (!fFitDone) { fCorrectionFits.clear(); }
    if (input_mat == nullptr)
    {
        throw std::runtime_error("TEC::CCM::FixMatrix: Error! Input matrix is nullptr!");
    }
    auto clone =
        dynamic_cast<TH2 *>(input_mat->Clone(Form("%s_corrected", V.TEMAT->GetName())));
    if (!clone) { throw std::runtime_error("Error: Cloning TEMAT failed!"); }

    std::shared_ptr<TH2> fixed_mat = std::shared_ptr<TH2>(clone);
    fixed_mat->SetName(Form("%s_fixed", input_mat->GetName()));
    fixed_mat->SetTitle(Form("%s_fixed", input_mat->GetTitle()));

    fixed_mat->SetDirectory(0);
    fixed_mat->Reset();

    Double_t bin_cont;
    Double_t total_ratio, en_low_edge, en_up_edge, new_en_low_edge, new_en_up_edge;

    Int_t global_bin, bin_start, bin_end;

    const TAxis   *axis               = fixed_mat->GetYaxis();
    const Double_t bin_width          = axis->GetBinWidth(1);
    const Double_t inverted_bin_width = 1. / (double)axis->GetBinWidth(1);

    for (int time_bin = 1; time_bin <= fixed_mat->GetXaxis()->GetNbins(); time_bin++)
    {
        auto time = fixed_mat->GetXaxis()->GetBinCenter(time_bin);

        // get fit
        const FitCont fit = this->GetCorrectionFit(time);

        // get fit function
        TF1 *fcn{nullptr};
        for (size_t i = 0; i < fCorrectionFunctions.size(); i++)
        {
            if (fit.functionUsed == std::string(fCorrectionFunctions[i].first->GetName()))
            {
                fcn = fCorrectionFunctions[i].first;
                break;
            }
        }
        // this cannot happen, function is always at lest "EMPTY_FUNCTION"!
        if (fcn == nullptr)
        {
            throw std::runtime_error("No function found for time slice " +
                                     std::to_string(time_bin) + " with name " +
                                     fit.functionUsed);
        }
        if (strcmp(fcn->GetName(), "EMPTY_FUNCTION") == 0)
        {
            std::cout << "Using empty function for time bin " << time_bin << std::endl;
            continue;
        }
        for (unsigned int i = 0; i < fit.coef.size(); i++)
        {
            fcn->SetParameter((int)i, fit.coef[i]);
        }

        // set bins
        for (int en_bin = 1; en_bin <= axis->GetNbins(); en_bin++)
        {
            bin_cont = input_mat->GetBinContent(input_mat->GetBin(time_bin, en_bin));
            if (bin_cont == 0) continue;
            // we need to get the new corrected energies that are matching low
            // and up edges of the old bin
            en_low_edge     = axis->GetBinLowEdge(en_bin);
            en_up_edge      = axis->GetBinUpEdge(en_bin);
            new_en_low_edge = fcn->Eval(en_low_edge);
            new_en_up_edge  = fcn->Eval(en_up_edge);

            // the new_en will not cover range that might be multiple bins in
            // width, and there is also a difference in coverage - new low edge
            // will cover only partially the new bin, same for up edge so we
            // need to get ratios on how to distribute the content of the old
            // bin to the new bins

            bin_start = axis->FindBin(new_en_low_edge);
            bin_end   = axis->FindBin(new_en_up_edge);
            // in case we are covering just 1 bin, put the content in
            if (bin_start == bin_end)
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_start);
                fixed_mat->AddBinContent(global_bin, bin_cont);
                continue;
            }
            // get to total ratio of "bin equivalent widths" we are covering
            // with the new energy
            // number of whole bins covered
            total_ratio = (new_en_up_edge - new_en_low_edge) * inverted_bin_width;
            if (total_ratio < 1.0) { total_ratio = 1.0; }

            // total_ratio = static_cast<double>(bin_end - bin_start - 1);
            // total_ratio += (double)(axis->GetBinUpEdge(bin_start) - new_en_low_edge) *
            //                inverted_bin_width;
            // total_ratio += (new_en_up_edge - axis->GetBinLowEdge(bin_end)) *
            //                inverted_bin_width;
            // set manually to new lower bin content
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_start);
                fixed_mat->AddBinContent(
                    global_bin,
                    (bin_cont * (axis->GetBinUpEdge(bin_start) - new_en_low_edge) *
                     inverted_bin_width) /
                        total_ratio);
            }
            // set manually to new upper bin content
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_end);
                fixed_mat->AddBinContent(
                    global_bin,
                    (bin_cont * (new_en_up_edge - axis->GetBinLowEdge(bin_end)) *
                     inverted_bin_width) /
                        total_ratio);
            }
            // in case there are more bins in between, fill them here (skipping
            // first already bin_start and bin_end)
            for (int ybin = bin_start + 1; ybin < bin_end; ybin++)
            {
                global_bin = fixed_mat->GetBin(time_bin, ybin);
                fixed_mat->AddBinContent(global_bin, bin_cont / total_ratio);
            }
        }
    }

    fixed_mat->SetEntries(input_mat->GetEntries());
    fixed_mat->ResetStats();
    return fixed_mat;
}

void TEC::CCM::FixTree(const std::string &tfilename,
                       const std::string &treename,
                       const std::string &e_branchname,
                       const std::string &ts_branchname,
                       const bool         valid_only,
                       const int          time_subdivision)
{
    std::cout << "Fix Tree not implemented yet!\n";
}

//     // input checks
//     TFile *old_file = TFile::Open(tfilename.c_str(), "READ");
//     if (!old_file || old_file->IsZombie())
//     {
//         std::cerr << "Error opening file\n";
//         return;
//     }
//     TTree *old_tree = (TTree *)old_file->Get(treename.c_str());
//     if (!old_tree)
//     {
//         std::cerr << "Error getting tree\n";
//         return;
//     }
//     TBranch *branch_e = old_tree->GetBranch(e_branchname.c_str());
//     TBranch *branch_t = old_tree->GetBranch(ts_branchname.c_str());
//     if (!branch_e || !branch_t)
//     {
//         std::cerr << "Branch does not exist\n";
//         return;
//     }

//     // prepare new root file

//     // copy all branches but the corrected ones
//     old_tree->SetBranchStatus("*", 1);
//     old_tree->SetBranchStatus(e_branchname.c_str(), 0);

//     // create new file and clone the tree
//     std::string new_file =
//         "corrected_" + tfilename.substr(tfilename.find_last_of('/') + 1);
//     TFile *new_tfile = new TFile(new_file.c_str(), "RECREATE");
//     TTree *new_tree  = old_tree->CloneTree();

//     // add missing branches in the new tree
//     double   new_e, old_e;
//     Long64_t timestamp;
//     TBranch *new_branch_e = new_tree->Branch(e_branchname.c_str(), &new_e);

//     // now we need to iterate through the old tree to get the original energy
//     // and timestamp, and store corrected values in a new tree
//     old_tree->SetBranchStatus("*", 0);
//     old_tree->SetBranchStatus(e_branchname.c_str(), 1);
//     old_tree->SetBranchStatus(ts_branchname.c_str(), 1);
//     old_tree->SetBranchAddress(e_branchname.c_str(), &old_e);
//     old_tree->SetBranchAddress(ts_branchname.c_str(), &timestamp);

//     std::vector<TGraph>   roi_graphs;
//     std::vector<TSpline3> roi_splines;

//     // replace missing poi
//     for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
//     {
//         roi_graphs.emplace_back();
//         roi_graphs.back().SetName(Form("graph_ROI_%i", ROI_index));
//         roi_graphs.back().SetTitle(Form(
//             "graph of ROI %i; time slice index; energy shift [energy units]",
//             ROI_index));
//         for (size_t time_index = 0; time_index < V.time_bins; time_index++)
//         {
//             if (valid_only && !ResVec[ROI_index][time_index].isValid) continue;
//             roi_graphs.back().AddPoint(this->GetMatrixTime(time_index),
//                                        V.ROIs[ROI_index].desired_energy -
//                                            ResVec[ROI_index][time_index].energy_shift);
//         }
//     }
//     for (auto &gr : roi_graphs)
//     {
//         roi_splines.emplace_back(Form("spline_%s", gr.GetName()), &gr);
//     }
//     // prepare some variables
//     const double   Xmax              = V.TEMAT->GetXaxis()->GetXmax();
//     const double   Xmin              = V.TEMAT->GetXaxis()->GetXmin();
//     const long     Nbins             = V.TEMAT->GetXaxis()->GetNbins();
//     const double   ts_half_bin_width = (Xmax - Xmin) / (double)(2 * Nbins);
//     const Long64_t ts_section_width =
//         (Long64_t)(Xmax - Xmin) / (Nbins * time_subdivision);
//     Long64_t next_section_ts_high = ts_section_width;
//     Long64_t next_section_ts_low  = 0;

//     std::pair<TF1 *, std::string> *correction_fcn{nullptr};

//     const Long64_t nentries = old_tree->GetEntries();
//     // std::cout << "***** Fixing " << treename << " of " << tfilename << "\n";
//     for (Long64_t entry = 0; entry < nentries; entry++)
//     {
//         // if (entry % 10000 == 0)
//         // {
//         //     std::cout << std::setprecision(3) << entry / (double)nentries * 100.
//         //               << " %    \r";
//         // }
//         old_tree->GetEntry(entry);
//         // check if function is set
//         if (timestamp > next_section_ts_high || timestamp < next_section_ts_low ||
//             !correction_fcn)
//         {
//             // find new timestamp section
//             while (timestamp > next_section_ts_high)
//             {
//                 next_section_ts_low = next_section_ts_high;
//                 next_section_ts_high += ts_section_width;
//             }
//             while (timestamp < next_section_ts_low)
//             {
//                 next_section_ts_high = next_section_ts_low;
//                 next_section_ts_low -= ts_section_width;
//             }

//             // get the parameters
//             // we need to take care not to extrapolate the spline by checking
//             // the boundaries we are playing with the ts_half_bin_width value,
//             // because as it is, splines uses knots that are in the center of
//             // the "time_bin", so we need to consider also events that are less
//             // than half of the bin width to still belong to the same bin and
//             // get the its parameters
//             TGraph gr;
//             uint   roi_index = 0;
//             for (const TSpline3 &spline : roi_splines)
//             {
//                 // std::cout << spline.GetXmin() << " " << spline.GetXmax() <<
//                 // std::endl;

//                 if ((double)timestamp < spline.GetXmin())
//                 {
//                     if ((double)timestamp + ts_half_bin_width > spline.GetXmin())
//                     {
//                         // add first point of the spline to compensate for the
//                         // different time of the left edge and center of the
//                         // time bin
//                         gr.AddPoint(V.ROIs[roi_index].desired_energy,
//                                     spline.Eval(spline.GetXmin()));
//                     }
//                     continue;
//                 }
//                 if ((double)timestamp > spline.GetXmax())
//                 {
//                     if ((double)timestamp < spline.GetXmax() + ts_half_bin_width)
//                     {
//                         gr.AddPoint(V.ROIs[roi_index].desired_energy,
//                                     spline.Eval(spline.GetXmax()));
//                     }
//                     continue;
//                 }

//                 gr.AddPoint(V.ROIs[roi_index].desired_energy,
//                             spline.Eval((double)timestamp));

//                 roi_index++;
//             }
//             if (gr.GetN() == 0) { correction_fcn = this->FindCorrectionFunction(0); }
//             else
//             {
//                 correction_fcn = this->FindCorrectionFunction(gr.GetN());
//                 gr.Fit(correction_fcn->first, correction_fcn->second.c_str());
//             }
//         }
//         // new_e = 40.;
//         new_e = correction_fcn->first->Eval(old_e);
//         // std::cout << correction_fcn->first->GetName() << std::setprecision(7)
//         // << " " << new_e << " " << old_e
//         //   << std::endl;
//         new_branch_e->Fill();
//     }
//     std::cout << "\n"
//               << "Done\n";

//     new_tfile->cd();
//     new_tree->Write();
//     new_tfile->Close();

//     // clean up
//     old_file->Close();
// }

std::pair<TF1 *, std::string> *TEC::CCM::FindCorrectionFunction(const int nrois)
{
    for (auto &f : fCorrectionFunctions)
    {
        if (nrois >= f.first->GetNpar()) { return &f; }
    }
    // we failed to find the function (should not happen since we have one with
    // zero parameters!)
    throw std::runtime_error(
        "Error! No suitable function found for given number of ROIs!");
    return nullptr;
}

const TEC::ResCont *TEC::CCM::GetResultContainer(const size_t ROI_no,
                                                 const size_t time_index) const noexcept
{
    if (ROI_no >= V.ROIs.size() || time_index >= V.time_bins) return nullptr;
    return &ResVec[ROI_no][time_index];
}
void TEC::CCM::SetResultStatus(const size_t ROI_no,
                               const size_t time_index,
                               const bool   valid)
{
    if (ROI_no >= V.ROIs.size() || time_index >= V.time_bins)
    {
        throw std::runtime_error("ROI_no or time_index are not correct");
    }
    ResVec[ROI_no][time_index].isValid = valid;

    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

void TEC::CCM::UseGaussianResult()
{
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        for (size_t time_index = 0; time_index < V.time_bins; time_index++)
        {
            ResVec[roi_index][time_index].bin_shift =
                ResVec[roi_index][time_index].gfit_mu +
                V.ROIs[roi_index].base_shift_value;
            ResVec[roi_index][time_index].energy_shift =
                V.TEMAT->GetYaxis()->GetBinWidth(1) *
                ResVec[roi_index][time_index].bin_shift;
            if (time_index == 1848)
            {
                std::cout << "gauss " << ResVec[roi_index][time_index].gfit_mu << " "
                          << ResVec[roi_index][time_index].bin_shift << " "
                          << ResVec[roi_index][time_index].energy_shift << std::endl;
            }
        }
    }
    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

void TEC::CCM::UsePolynomialResult()
{
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        for (size_t time_index = 0; time_index < V.time_bins; time_index++)
        {

            ResVec[roi_index][time_index].bin_shift =
                ResVec[roi_index][time_index].poly_shift +
                V.ROIs[roi_index].base_shift_value;
            ResVec[roi_index][time_index].energy_shift =
                V.TEMAT->GetYaxis()->GetBinWidth(1) *
                ResVec[roi_index][time_index].bin_shift;
            if (time_index == 1848)
            {
                std::cout << "poly2 " << ResVec[roi_index][time_index].poly_shift << " "
                          << ResVec[roi_index][time_index].bin_shift << " "
                          << ResVec[roi_index][time_index].energy_shift << std::endl;
            }
        }
    }
    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

void TEC::CCM::UseMaxDPResult()
{
    for (uint roi_index = 0; roi_index < V.ROIs.size(); roi_index++)
    {
        for (size_t time_index = 0; time_index < V.time_bins; time_index++)
        {
            auto maxdp_iterator =
                std::max_element(ResVec[roi_index][time_index].dp_vec.begin(),
                                 ResVec[roi_index][time_index].dp_vec.end());

            int maxdp_shift = std::distance(ResVec[roi_index][time_index].dp_vec.begin(),
                                            maxdp_iterator);
            // std::cout << "MaxDP shift: "
            //           << maxdp_shift + V.ROIs[roi_index].base_shift_value << " time "
            //           << time_index << std::endl;
            ResVec[roi_index][time_index].bin_shift =
                (maxdp_shift + V.ROIs[roi_index].base_shift_value);
            ResVec[roi_index][time_index].energy_shift =
                V.TEMAT->GetYaxis()->GetBinWidth(1) *
                ResVec[roi_index][time_index].bin_shift;
            if (time_index == 1848)
            {
                std::cout << "maxdp " << maxdp_shift << " "
                          << ResVec[roi_index][time_index].bin_shift << " "
                          << ResVec[roi_index][time_index].energy_shift << std::endl;
            }
        }
    }
    fForceRebuildInterpolators = true;
    fFitDone                   = false;
}

void TEC::CCM::SetReferenceVector(const unsigned int        ROI_index,
                                  const std::vector<float> &own_reference_vector)
{
    if (ROI_index > V.ROIs.size())
    {
        throw std::runtime_error("ROI index out of bounds");
    }

    int vsize = V.ROIs[ROI_index].bin_window_high - V.ROIs[ROI_index].bin_window_low;

    if ((int)own_reference_vector.size() != vsize)
    {
        throw std::runtime_error(
            "Reference vector size does not match expected sample vector size");
    }
    V.sample_vector[ROI_index].clear();
    V.sample_vector[ROI_index] = own_reference_vector;

    this->Normalize(V.sample_vector[ROI_index]);
}

void TEC::CCM::SetReferenceProjection(const TH1 *projection)
{
    V.sample_vector.clear();

    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        std::vector<float> vec{};
        for (int e_bin = V.ROIs[ROI_index].bin_window_low;
             e_bin < V.ROIs[ROI_index].bin_window_high; e_bin++)
        {
            vec.push_back((float)projection->GetBinContent(e_bin));
        }
        this->Normalize(vec);
        V.sample_vector.emplace_back(std::move(vec));
    }
}

void TEC::CCM::CheckReferenceVector(const size_t ROI_index)
{
    int vsize = V.ROIs[ROI_index].bin_window_high - V.ROIs[ROI_index].bin_window_low;

    if ((int)V.sample_vector[ROI_index].size() != vsize)
    {
        throw std::runtime_error("Reference vector not set for ROI " +
                                 std::to_string(ROI_index));
    }
}

void TEC::CCM::CheckReferenceVectors()
{
    for (uint ROI_index = 0; ROI_index < V.ROIs.size(); ROI_index++)
    {
        this->CheckReferenceVector(ROI_index);
    }
}

std::unique_ptr<TGraph> TEC::CCM::GetDotProductGraph(const size_t roi_index,
                                                     const int    time_bin)
{
    if (roi_index >= V.ROIs.size())
    {
        throw std::runtime_error("ROI index out of bounds");
    }
    std::unique_ptr<TGraph> gr = std::make_unique<TGraph>();

    const auto rc = this->GetResultContainer(roi_index, time_bin);
    if (rc == nullptr) { throw std::runtime_error("Result container is nullptr"); }
    double sum = 0;
    for (int i = 0; i < rc->dp_vec.size(); i++)
    {
        sum += rc->dp_vec.at(i);
        gr->AddPoint(i + V.ROIs[roi_index].base_shift_value, rc->dp_vec.at(i));
    }
    gr->SetName(Form("dot_product_graph_%i_%i", roi_index, time_bin));
    gr->SetTitle(Form("Dot product graph for ROI %i, time bin %i", roi_index, time_bin));
    gr->GetXaxis()->SetTitle(
        Form("Bin shift of ROI %i at time bin %i", roi_index, time_bin));

    std::cout << "Sum of dot product vector: " << sum << std::endl;
    return gr;
}

std::unique_ptr<TGraph> TEC::CCM::GetInterpolationGraph(const size_t ROI_index,
                                                        const int    subdivide,
                                                        const bool   valid_only)
{
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("ROI index out of bounds");
    }

    auto                   *interpolator = &V.ROIs[ROI_index].interpolator;
    std::unique_ptr<TGraph> gr           = std::make_unique<TGraph>();

    const double x_start = interpolator->GetFirstX();
    const double x_end   = interpolator->GetLastX();
    const double step =
        (x_end - x_start) / ((double)interpolator->GetNPoints() * subdivide);
    // double x_start = interpolator->GetFirstX() - step * subdivide / 2.;
    // double x_end   = interpolator->GetLastX() + step * subdivide / 2.;

    double y;
    for (double x = x_start; x < x_end; x += step)
    {
        if (valid_only && !interpolator->IsValueValid(x)) continue;
        y = interpolator->Eval(x);
        gr->AddPoint(x, y);
    }
    return gr;
}

void TEC::CCM::SmoothShifts(const SmootherType smoother,
                            const double       smoother_parameter,
                            const size_t       ROI_index)
{
    // check if we are in bounds
    if (ROI_index >= V.ROIs.size())
    {
        throw std::runtime_error("ROI index out of bounds");
    }

    if (smoother == SmootherType::NONE) return;

    // define shortcut function to get smoothed graph
    auto smoothing_function = [&](TGraph &data, CCMInterpolator *interpolator) -> void {
        TGraphSmooth gs;
        TGraph      *gr{nullptr};
        switch (smoother)
        {
        case SmootherType::LOWESS:
            gr = gs.SmoothLowess(&data, "", smoother_parameter);
            break;
        case SmootherType::KERNEL:
            gr = gs.SmoothKern(&data, "", smoother_parameter);
            break;
        case SmootherType::SUPER:
            gr = gs.SmoothSuper(&data, "", smoother_parameter);
            break;
        default:
            throw std::runtime_error("Impossible option");
            // case SmootherType::NONE: throw std::runtime_error("Impossible option");
        }
        for (int point = 0; point < gr->GetN(); point++)
        {
            interpolator->AddPoint(gr->GetPointX(point), gr->GetPointY(point), true);
        }
    };

    fFitDone = false;
    TGraph gr_data;

    auto *interpolator = &V.ROIs[ROI_index].interpolator;
    interpolator->ClearPoints();
    for (size_t time_index = 0; time_index < V.time_bins; time_index++)
    {
        if (ResVec[ROI_index][time_index].isValid)
        {
            gr_data.AddPoint(V.TEMAT->GetXaxis()->GetBinCenter((int)time_index + 1),
                             ResVec[ROI_index][time_index].energy_shift);
        }
        else
        {
            // do not smooth if we have too little points. Just copy points to
            // interpolator
            if (gr_data.GetN() < MINIMUM_SMOOTHING_POINTS)
            {
                for (int point = 0; point < gr_data.GetN(); point++)
                {
                    interpolator->AddPoint(gr_data.GetPointX(point),
                                           gr_data.GetPointY(point), true);
                }
            }
            else { smoothing_function(gr_data, interpolator); }
            interpolator->AddPoint(V.TEMAT->GetXaxis()->GetBinCenter((int)time_index + 1),
                                   ResVec[ROI_index][time_index].energy_shift, false);
            gr_data = TGraph();
        }
    }

    if (gr_data.GetN() < MINIMUM_SMOOTHING_POINTS)
    {
        for (int point = 0; point < gr_data.GetN(); point++)
        {
            interpolator->AddPoint(gr_data.GetPointX(point), gr_data.GetPointY(point),
                                   true);
        }
    }
    else { smoothing_function(gr_data, interpolator); }
}

void TEC::CCM::SmoothShifts(const SmootherType smoother, const double smoother_parameter)
{
    for (size_t i = 0; i < V.ROIs.size(); i++)
    {
        this->SmoothShifts(smoother, smoother_parameter, i);
    }
    return;
}