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

// Root
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"

#include "CCM.h"
// #include "CheckCCM.h"
#include "Cross_correlation.h"
#include "variables.h"

double useless_function(double *x, double *par)
{
    return x[0];
}

CCM::CCM(const TH2D &matrix, const std::vector<Region_of_interest> &_ROIs, const double reference_time_low,
         const double reference_time_high)
    : CCM(matrix, _ROIs)
{

    // create sample_vector
    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        this->CreateReferenceVector(ROI_index, reference_time_low, reference_time_high);
    }
}

CCM::CCM(const TH2D &matrix, const std::vector<Region_of_interest> &_ROIs)
{

    V.TEMAT = dynamic_cast<TH2D *>(matrix.Clone(Form("%s_clone", matrix.GetName())));
    fXbins = matrix.GetXaxis()->GetNbins();
    fYbins = matrix.GetYaxis()->GetNbins();

    for (const auto &roi : _ROIs)
    {
        V.ROIs.emplace_back(roi);
    }

    // V.ROIs = _ROIs;
    V.number_of_ROIs = _ROIs.size();
    V.total_tasks = fXbins * V.number_of_ROIs;
    V.time_bins = fXbins;

    // reserve space for result container structure
    ResVec = new ResCont *[V.number_of_ROIs];
    for (int ROI = 0; ROI < V.number_of_ROIs; ROI++)
    {
        ResVec[ROI] = new ResCont[V.time_bins];
    }
    // reserve space for fit container - fit results and other variables
    FitVec = new FitCont[V.time_bins];
    this->CopyMatrixContent(V.TEMAT);
}

CCM::~CCM()
{
    if (fFixedTEMAT)
        delete fFixedTEMAT;

    for (int ROI = 0; ROI < V.number_of_ROIs; ROI++)
    {
        delete[] ResVec[ROI];
    }
    delete[] ResVec;
    delete[] FitVec;

    for (auto i = 0; i < V.number_of_ROIs; i++)
    {
        for (auto j = 0; j < fXbins; j++)
        {
            delete[] V.TEMATarr[i][j];
        }
        delete[] V.TEMATarr[i];
    }
    delete[] V.TEMATarr;
    delete V.TEMAT;

    for (auto &f : fCorrectionFunctions)
    {
        delete f.first;
        f.first = nullptr;
    }
}

void CCM::SetCorrectionFunction(const TF1 &fcn, const std::string &fit_options)
{
    if (!fcn.GetFormula())
    {
        throw std::runtime_error("Error! Function must be created using TFormula, e.i. using the text string, "
                                 "otherwise it cannot be saved in a ROOT file and used in the CCM object!");
    }
    if (fCorrectionFunctions.size() > 0)
    {
        throw std::runtime_error("Error! Primary function needs to be set first, only one primary function can be set");
    }
    TF1 *fclone = (TF1 *)fcn.Clone();
    // TF1 *fclone = new TF1("gain_fcn", "[0]*x", 0, 4000);

    fclone->SetNpx(V.TEMAT->GetYaxis()->GetNbins() * 100);
    fclone->SetRange(V.TEMAT->GetYaxis()->GetBinLowEdge(1),
                     V.TEMAT->GetYaxis()->GetBinUpEdge(V.TEMAT->GetYaxis()->GetNbins()));
    fclone->Update();
    fCorrectionFunctions.emplace_back(std::make_pair(fclone, fit_options + "NQ"));
    // add empty function for cases when no ROI is valid - this eseential does nothing
    TF1 *empty_function = new TF1("empty_function", "x", V.TEMAT->GetYaxis()->GetBinLowEdge(1),
                                  V.TEMAT->GetYaxis()->GetBinUpEdge(V.TEMAT->GetYaxis()->GetNbins()));
    fCorrectionFunctions.emplace_back(std::make_pair(empty_function, "NQ"));
}

void CCM::SetFallbackCorrectionFunction(const TF1 &fcn, const std::string &fit_options)
{
    if (!fcn.GetFormula())
    {
        throw std::runtime_error("Error! Function must be created using TFormula, e.i. using the text string, "
                                 "otherwise it cannot be saved in a ROOT file and used in the CCM object!");
    }
    if (fCorrectionFunctions.size() == 0)
    {
        throw std::runtime_error("Error! Primary function not set, it needs to be set first!");
    }

    auto it = std::find_if(fCorrectionFunctions.begin(), fCorrectionFunctions.end(),
                           [&fcn](const auto &f) { return fcn.GetNpar() == f.first->GetNpar(); });
    if (it != fCorrectionFunctions.end())
    {
        throw std::runtime_error("Error! Function with the same number of parameters already exists!");
    }

    for (const auto &f : fCorrectionFunctions)
    {
        if (strcmp(f.first->GetName(), fcn.GetName()) == 0)
        {
            throw std::runtime_error("Error! You are trying to add function with same name of an existing function!");
        }
    }
    TF1 *fclone = (TF1 *)fcn.Clone();
    fclone->SetNpx(V.TEMAT->GetYaxis()->GetNbins() * 100);
    fclone->SetRange(V.TEMAT->GetYaxis()->GetBinLowEdge(1),
                     V.TEMAT->GetYaxis()->GetBinUpEdge(V.TEMAT->GetYaxis()->GetNbins()));
    fCorrectionFunctions.emplace_back(std::make_pair(fclone, fit_options + "NQ"));

    std::sort(fCorrectionFunctions.begin(), fCorrectionFunctions.end(),
              [](const auto &a, const auto &b) { return a.first->GetNpar() > b.first->GetNpar(); });
}

TGraph CCM::GetROIShifts(const int roi_index, const bool valid_only)
{
    TGraph gr;
    gr.SetBit(TGraph::kIsSortedX);
    gr.SetName(Form("shift_ROI_%i", roi_index));

    for (int time = 0; time < V.time_bins; time++)
    {
        if (valid_only && !ResVec[roi_index][time].isValid)
        {
            continue;
        }
        gr.AddPoint(GetMatrixTime(time), ResVec[roi_index][time].bin_shift);
    }
    return gr;
}

TGraph CCM::GetShiftProfile(const int time_bin, const bool valid_only)
{
    TGraph gr;
    gr.SetBit(TGraph::kIsSortedX);
    gr.SetTitle(Form("%i", time_bin));
    for (int roi_index = 0; roi_index < V.number_of_ROIs; roi_index++)
    {
        if (ResVec[roi_index][time_bin].isValid)
        {
            if (valid_only && !ResVec[roi_index][time_bin].isValid)
            {
                continue;
            }
            gr.AddPoint(V.ROIs[roi_index].desired_energy, ResVec[roi_index][time_bin].bin_shift);
        }
    }
    gr.SetMarkerColor(kRed);
    gr.SetMarkerStyle(20);
    return gr;
}

void CCM::CopyMatrixContent(TH2D *matrix)
{
    // V.TEMAT = matrix;
    // fXbins = V.TEMAT->GetXaxis()->GetNbins();
    // V.time_bins = fXbins;

    V.TEMATarr = new double **[V.number_of_ROIs];

    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        V.TEMATarr[ROI_index] = new double *[fXbins];
        for (uint j = 0; j < fXbins; j++)
        {
            V.TEMATarr[ROI_index][j] = new double[V.ROIs[ROI_index].displacement_range];
        }
    }

    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        for (int x = 0; x < fXbins; x++)
        {
            for (int y = 0; y < V.ROIs[ROI_index].displacement_range; y++)
            {
                V.TEMATarr[ROI_index][x][y] =
                    V.TEMAT->GetBinContent(x + 1, y + V.ROIs[ROI_index].bin_window_low +
                                                      V.ROIs[ROI_index].bin_displacement_low); // y+ 1?? - check
            }
        }
    }
}

void CCM::CreateReferenceVector(const int ROI_index, const double sample_time_low, const double sample_time_high)
{
    std::vector<double> vec(V.ROIs[ROI_index].vector_dimension);

    int vector_iterator{0};
    // loop over energy

    int tbin_start = V.TEMAT->GetXaxis()->FindBin(sample_time_low);
    int tbin_end = V.TEMAT->GetXaxis()->FindBin(sample_time_high);

    for (int e_bin = V.ROIs[ROI_index].bin_window_low; e_bin < V.ROIs[ROI_index].bin_window_high;
         e_bin++, vector_iterator++)
    {
        vec[vector_iterator] = 0;
        // loop over time, since sample vector should be formed of more than 1 time-bin width of TEMAT
        for (int t_bin = tbin_start; t_bin <= tbin_end; t_bin++)
        {
            vec[vector_iterator] += V.TEMAT->GetBinContent(t_bin, e_bin);
        }
    }

    this->Normalize(vec);

    V.sample_vector.emplace_back(std::move(vec));
}

void CCM::Normalize(std::vector<double> &v)
{
    double norm = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    if (norm <= 0)
    {
        throw std::runtime_error("Normalization ERROR, sample vector cannot consist of zeros!");
    }
    norm = sqrt(norm);
    std::transform(v.begin(), v.end(), v.begin(), [norm](double &x) { return x / norm; });
}

void CCM::CalculateEnergyShifts(const unsigned int threads)
{
    this->CheckReferenceVectors();
    fNthreads = threads;

    std::vector<std::thread> t;
    t.reserve(fNthreads);

    std::cout << "***** Starting calculations of Cross-correlation corrections  *****" << std::endl;

    fThreadTask = -1; // so that the first job start at zero

    CrossCorrel *starting_class[fNthreads];
    for (int i = 0; i < fNthreads; i++)
    {
        starting_class[i] = new CrossCorrel(&V, ResVec);
    }

    std::mutex mtx_task;
    std::mutex mtx_fit;
    // start threads
    if (fNthreads > V.total_tasks)
    {
        for (int i = 0; i < V.total_tasks; i++)
        {
            // pass thread task, VarManger variable, mutex and ResVec as reference
            t.push_back(std::thread(&CrossCorrel::Process, starting_class[i], i + 1, &fThreadTask, std::ref(mtx_task),
                                    std::ref(mtx_fit)));
        }
    }
    else
    {
        for (int i = 0; i < fNthreads; i++)
        {
            // pass thread task, VarManger variable, mutex and ResVec as reference
            t.push_back(std::thread(&CrossCorrel::Process, starting_class[i], i + 1, &fThreadTask, std::ref(mtx_task),
                                    std::ref(mtx_fit)));
        }
    }
    for (auto &th : t)
        th.join();

    std::cout << std::endl;
    std::cout << "Done" << std::endl;
}

void CCM::SaveToRootFile(const std::string &outroot_file)
{
    TFile file(outroot_file.c_str(), "RECREATE");
    {
        TH2D *mat = dynamic_cast<TH2D *>(V.TEMAT->Clone("CCM_matrix"));
        mat->Write();
        if (fFixedTEMAT)
            fFixedTEMAT->Write();
    }
    for (int roi = 0; roi < V.number_of_ROIs; roi++)
    {
        ResCont results;
        double time;

        // TTree *t = new TTree(Form("ROI_%i", roi), Form("CCM tree of ROI %i", roi));
        TTree *t = new TTree(Form("ROI_%i", roi), Form("CCM tree of ROI %i", roi));
        t->Branch("fit_valid", &results.isValid);
        t->Branch("bin_shift", &results.bin_shift);
        t->Branch("energy_shift", &results.energy_shift);
        t->Branch("dot_product", &results.dp);
        t->Branch("dp_vector", &results.dp_vec);
        t->Branch("gfit_chi2", &results.gfit_chi2);
        t->Branch("gfit_sigma", &results.gfit_sigma);
        t->Branch("gfit_mu", &results.gfit_sigma);
        t->Branch("time", &time);

        for (int time_bin = 0; time_bin < V.time_bins; time_bin++)
        {
            time = this->GetMatrixTime(time_bin);
            results = ResVec[roi][time_bin];
            t->Fill();
        }
        file.cd();
        t->Write();
    }

    {
        std::string fcn_name;
        std::vector<double> fcn_coefs;
        double time;
        TTree *t = new TTree("fit_coef", "Fit coefficients for given time, parameters = 0.");

        t->Branch("time", &time);
        t->Branch("fcn_name", &fcn_name);
        t->Branch("fcn_coef", &fcn_coefs);

        for (int time_bin = 0; time_bin < V.time_bins; time_bin++)
        {
            time = this->GetMatrixTime(time_bin);
            // first zero out params
            fcn_name = FitVec[time_bin].functionUsed;
            fcn_coefs = FitVec[time_bin].coef;
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

void CCM::SaveShiftTable(const std::string &table_filename)
{

    std::ofstream output;
    output.open(table_filename);
    // table head
    output << "# number of ROIs (regions of interest) in the matrix \n";
    output << V.number_of_ROIs << '\n';
    output << "# desired energies for each ROI \n";
    for (const auto &r : V.ROIs)
    {
        output << std::setprecision(15) << r.desired_energy << " ";
    }
    output << '\n';
    output << "################ \n";
    output << "# time_bin time_start time_end ";
    for (int roi = 0; roi < V.number_of_ROIs; roi++)
    {
        output << "ROI" << roi << "_ok(bool) "
               << "ROI" << roi << "_shift(double) ";
    }
    output << "\n";

    for (int time_slice_index = 0; time_slice_index < V.time_bins; time_slice_index++)
    {
        output << time_slice_index + 1 << std::setprecision(15) << '\t';
        output << std::setprecision(15) << V.TEMAT->GetXaxis()->GetBinLowEdge(time_slice_index + 1) << '\t';
        output << std::setprecision(15) << V.TEMAT->GetXaxis()->GetBinUpEdge(time_slice_index + 1) << '\t';
        for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
        {
            // check if shift is converted to energy units
            output << '\t' << ResVec[ROI_index][time_slice_index].isValid << '\t'
                   << ResVec[ROI_index][time_slice_index].energy_shift;
        }
        output << '\n';
    }
}

// void CCM::BuildInterpolators(const bool valid_only)
// {
//     for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
//     {
//         bool graph_defined = std::find(fInterpolator.begin(), fInterpolator.end(),
//                                        [&](const CCMInterpolator &x) { return x.GetROI() == ROI_index; }

//                                        ) == fInterpolator.end();
//         if (graph_defined)
//             continue;
//         fInterpolator.emplace_back(ROI_index);
//         for (int time_index = 0; time_index < V.time_bins; time_index++)
//         {
//             if (valid_only)
//             {
//                 fInterpolator.back().AddPoint(V.TEMAT->GetXaxis()->GetBinCenter(time_index + 1),
//                                               ResVec[ROI_index][time_index].energy_shift,
//                                               ResVec[ROI_index][time_index].isValid);
//             }
//             else
//             {
//                 fInterpolator.back().AddPoint(V.TEMAT->GetXaxis()->GetBinCenter(time_index + 1),
//                                               ResVec[ROI_index][time_index].energy_shift, true);
//             }
//         }
//     }
// }

// void CCM::PerformFits(const bool valid_only, const bool use_spline)
// {
//     if (V.number_of_ROIs == 0)
//     {
//         throw std::runtime_error("No ROIs added to the CCM object");
//     }

//     int nrois;
//     double x[V.number_of_ROIs];
//     double y[V.number_of_ROIs];

//     std::vector<TGraph> roi_graphs;
//     std::vector<TSpline3> roi_splines;

//     // prepare splines if needed
//     if (use_spline)
//     {
//         for (int time_index = 0; time_index < V.time_bins; time_index++)
//         {
//             roi_graphs.emplace_back();
//             roi_graphs.back().SetName(Form("ROI_%i", roi_graphs.size() - 1));
//             roi_graphs.back().SetTitle(
//                 Form("ROI_%i; time slice index; energy shift [energy units]", roi_graphs.size() - 1));
//             for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
//             {
//                 if (valid_only && !ResVec[ROI_index][time_index].isValid)
//                     continue;
//                 roi_graphs.back().AddPoint(time_index, V.ROIs[ROI_index].desired_energy +
//                                                            ResVec[ROI_index][time_index].energy_shift);
//             }
//         }
//         for (auto &gr : roi_graphs)
//         {
//             roi_splines.emplace_back(Form("spline_%s", gr.GetName()), &gr);
//         }
//     }

//     for (int time_index = 0; time_index < V.time_bins; time_index++)
//     {
//         nrois = 0;

//         for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
//         {
//             // ResVec[ROI_index][time_index].print();
//             if (valid_only && ResVec[ROI_index][time_index].isValid)
//             {
//                 x[nrois] = V.ROIs[ROI_index].desired_energy;
//                 y[nrois] = V.ROIs[ROI_index].desired_energy - ResVec[ROI_index][time_index].energy_shift;
//                 nrois++; // get number of good points
//             }
//             else if (use_spline)
//             {
//                 x[nrois] = V.ROIs[ROI_index].desired_energy;
//                 y[nrois] = V.ROIs[ROI_index].desired_energy - roi_splines[ROI_index].Eval(time_index);
//                 nrois++;
//             }
//         }

//         std::pair<TF1 *, std::string> *fcn = nullptr; // TF1 *fcn{nullptr};
//         for (auto &f : fCorrectionFunctions)
//         {
//             if (nrois >= f.first->GetNpar())
//             {
//                 fcn = &f;
//                 break;
//             }
//         }
//         // if (fcn->first == nullptr)
//         // {
//         //     auto min =
//         //         std::min_element(fCorrectionFunctions.begin(), fCorrectionFunctions.end(),
//         //                          [](const auto &a, const auto &b) { return a.first->GetNpar() <
//         b.first->GetNpar();
//         //                          });
//         //     std::cerr
//         //         << "No suitable function found for time slice " << time_index << ", this slice has in total " <<
//         //         nrois
//         //         << " valid ROIs available to fit, but smallest minimum degree of freedom from provided functions
//         is "
//         //         << (*min).first->GetNpar() << std::endl;
//         //     std::cerr << "skipping slice" << std::endl;
//         //     continue;
//         // }

//         if (strcmp(fcn->first->GetName(), "empty_function") == 0)
//         {
//             FitVec[time_index].functionUsed = "empty_function";
//             continue;
//         }
//         TGraph gr(nrois, x, y);
//         gr.Fit(fcn->first, fcn->second.c_str());
//         FitVec[time_index].functionUsed = fcn->first->GetName();
//         for (int i = 0; i < fcn->first->GetNpar(); i++)
//         {
//             // std::cout << std::setprecision(6) << "par[" << i << "] = " << fcn->first->GetParameter(i) <<
//             std::endl; FitVec[time_index].coef.emplace_back(fcn->first->GetParameter(i));
//         }
//     }

//     // SaveFitTable(&V, FitVec);
//     fFitDone = true;
// }

void CCM::PerformFits(const bool valid_only, const bool use_spline)
{
    if (V.number_of_ROIs == 0)
    {
        throw std::runtime_error("No ROIs added to the CCM object");
    }

    int nrois;
    double x[V.number_of_ROIs];
    double y[V.number_of_ROIs];

    std::vector<TGraph> roi_graphs;
    std::vector<TSpline3> roi_splines;

    // prepare splines if needed
    if (use_spline)
    {
        for (int time_index = 0; time_index < V.time_bins; time_index++)
        {
            roi_graphs.emplace_back();
            roi_graphs.back().SetName(Form("ROI_%i", roi_graphs.size() - 1));
            roi_graphs.back().SetTitle(
                Form("ROI_%i; time slice index; energy shift [energy units]", roi_graphs.size() - 1));
            for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
            {
                if (valid_only && !ResVec[ROI_index][time_index].isValid)
                    continue;
                roi_graphs.back().AddPoint(time_index, V.ROIs[ROI_index].desired_energy +
                                                           ResVec[ROI_index][time_index].energy_shift);
            }
        }
        for (auto &gr : roi_graphs)
        {
            roi_splines.emplace_back(Form("spline_%s", gr.GetName()), &gr);
        }
    }

    for (int time_index = 0; time_index < V.time_bins; time_index++)
    {
        nrois = 0;

        for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
        {
            // ResVec[ROI_index][time_index].print();
            if (valid_only && ResVec[ROI_index][time_index].isValid)
            {
                x[nrois] = V.ROIs[ROI_index].desired_energy;
                y[nrois] = V.ROIs[ROI_index].desired_energy - ResVec[ROI_index][time_index].energy_shift;
                nrois++; // get number of good points
            }
            else if (use_spline)
            {
                x[nrois] = V.ROIs[ROI_index].desired_energy;
                y[nrois] = V.ROIs[ROI_index].desired_energy - roi_splines[ROI_index].Eval(time_index);
                nrois++;
            }
        }

        std::pair<TF1 *, std::string> *fcn = nullptr; // TF1 *fcn{nullptr};
        for (auto &f : fCorrectionFunctions)
        {
            if (nrois >= f.first->GetNpar())
            {
                fcn = &f;
                break;
            }
        }
        // if (fcn->first == nullptr)
        // {
        //     auto min =
        //         std::min_element(fCorrectionFunctions.begin(), fCorrectionFunctions.end(),
        //                          [](const auto &a, const auto &b) { return a.first->GetNpar() < b.first->GetNpar();
        //                          });
        //     std::cerr
        //         << "No suitable function found for time slice " << time_index << ", this slice has in total " <<
        //         nrois
        //         << " valid ROIs available to fit, but smallest minimum degree of freedom from provided functions is "
        //         << (*min).first->GetNpar() << std::endl;
        //     std::cerr << "skipping slice" << std::endl;
        //     continue;
        // }

        if (strcmp(fcn->first->GetName(), "empty_function") == 0)
        {
            FitVec[time_index].functionUsed = "empty_function";
            continue;
        }
        TGraph gr(nrois, x, y);
        gr.Fit(fcn->first, fcn->second.c_str());
        FitVec[time_index].functionUsed = fcn->first->GetName();
        for (int i = 0; i < fcn->first->GetNpar(); i++)
        {
            // std::cout << std::setprecision(6) << "par[" << i << "] = " << fcn->first->GetParameter(i) << std::endl;
            FitVec[time_index].coef.emplace_back(fcn->first->GetParameter(i));
        }
    }

    // SaveFitTable(&V, FitVec);
    fFitDone = true;
}

void CCM::SaveFitTable(const std::string &data_filename, const std::string &detector_name)
{
    if (!fFitDone)
    {
        std::cerr << "Warning! PerformFits() function was not yet called, it will be called with default parameters "
                     "(valid only, no spline)"
                  << std::endl;
        this->PerformFits();
    }

    std::ofstream output;
    output.open(data_filename);
    output << "# detector_name\n" << detector_name << '\n';
    output << "# Timestamp range covered by this file\n";
    output << std::setprecision(20) << V.TEMAT->GetXaxis()->GetBinLowEdge(1) << " " << std::setprecision(20)
           << V.TEMAT->GetXaxis()->GetBinUpEdge(V.time_bins) << '\n';
    output << "# number of functions used \n";
    output << fCorrectionFunctions.size() << '\n';
    output << "# fcn_name number_of_parameters function_equation\n";
    for (const auto &f : fCorrectionFunctions)
    {
        output << f.first->GetName() << "\t" << f.first->GetNpar() << "\t" << f.first->GetExpFormula() << '\n';
    }
    output << "#############################################################" << '\n';
    output << "# TS_start TS_stop fcn_name par0 par1 ..." << '\n';

    for (int time_index = 0; time_index < V.time_bins; time_index++)
    {
        output << V.TEMAT->GetXaxis()->GetBinLowEdge(time_index + 1) << '\t'
               << V.TEMAT->GetXaxis()->GetBinUpEdge(time_index + 1) << '\t';
        output << FitVec[time_index].functionUsed << '\t' << FitVec[time_index].coef.size() << '\t';
        for (int i = 0; i < FitVec[time_index].coef.size(); i++)
        {
            output << FitVec[time_index].coef[i] << '\t';
        }
        output << '\n';
    }
}

TH2D *CCM::FixMatrix()
{
    if (!fFitDone)
    {
        std::cerr << "Warning! PerformFits() function was not yet called, it will be called with default parameters "
                     "(valid only, no spline)"
                  << std::endl;
        this->PerformFits();
    }

    fFixedTEMAT = dynamic_cast<TH2D *>(V.TEMAT->Clone(Form("%s_corrected", V.TEMAT->GetName())));
    fFixedTEMAT->Reset();

    Double_t bin_cont;
    Double_t total_ratio, en_low_edge, en_up_edge, new_en_low_edge, new_en_up_edge;

    Int_t global_bin, bin_start, bin_end;
    const Double_t bin_width = V.TEMAT->GetYaxis()->GetBinWidth(1);

    const TAxis *axis = V.TEMAT->GetYaxis();

    for (int time_bin = 1; time_bin <= V.TEMAT->GetXaxis()->GetNbins(); time_bin++)
    {
        // load function with parameters
        TF1 *fcn{nullptr};
        for (int i = 0; i < fCorrectionFunctions.size(); i++)
        {
            if (FitVec[time_bin - 1].functionUsed == std::string(fCorrectionFunctions[i].first->GetName()))
            {
                fcn = fCorrectionFunctions[i].first;
                break;
            }
        }
        if (fcn == nullptr)
        {
            throw std::runtime_error("No function found for time slice " + std::to_string(time_bin) + " with name " +
                                     FitVec[time_bin - 1].functionUsed);
        }
        for (unsigned int i = 0; i < FitVec[time_bin - 1].coef.size(); i++)
        {
            fcn->SetParameter(i, FitVec[time_bin - 1].coef[i]);
        }
        if (strcmp(fcn->GetName(), "useless_function") == 0)
            continue;

        // set bins
        for (int en_bin = 1; en_bin <= axis->GetNbins(); en_bin++)
        {
            bin_cont = V.TEMAT->GetBinContent(V.TEMAT->GetBin(time_bin, en_bin));
            if (bin_cont == 0)
                continue;
            // we need to get the new corrected energies that are matching low and up edges of the old bin
            en_low_edge = axis->GetBinLowEdge(en_bin);
            en_up_edge = axis->GetBinUpEdge(en_bin);
            new_en_low_edge = fcn->Eval(en_low_edge);
            new_en_up_edge = fcn->Eval(en_up_edge);

            // the new_en will not cover range that might be multiple bins in width, and there is also a difference in
            // coverage - new low edge will cover only partially the new bin, same for up edge so we need to get ratios
            // on how to distribute the content of the old bin to the new bins

            bin_start = axis->FindBin(new_en_low_edge);
            bin_end = axis->FindBin(new_en_up_edge);
            // in case we are covering just 1 bin, put the content in
            if (bin_start == bin_end)
            {
                global_bin = fFixedTEMAT->GetBin(time_bin, bin_start);
                fFixedTEMAT->AddBinContent(global_bin, bin_cont);
                continue;
            }
            // get to total ratio of "bin equivalent widths" we are covering with the new energy
            total_ratio = static_cast<double>(bin_end - bin_start - 1);
            // std::cout << std::setprecision(6) << "total ratio " << total_ratio << " ";
            total_ratio += (double)(axis->GetBinUpEdge(bin_start) - new_en_low_edge) / (double)bin_width;
            // std::cout << std::setprecision(6) << total_ratio << " ";
            total_ratio += (new_en_up_edge - axis->GetBinLowEdge(bin_end)) / bin_width;
            // std::cout << std::setprecision(6) << total_ratio << std::endl;
            // set manually to new lower bin content
            {
                global_bin = fFixedTEMAT->GetBin(time_bin, bin_start);
                fFixedTEMAT->AddBinContent(global_bin,
                                           (bin_cont * (axis->GetBinUpEdge(bin_start) - new_en_low_edge) / bin_width) /
                                               total_ratio);
            }
            // set manually to new upper bin content
            {
                global_bin = fFixedTEMAT->GetBin(time_bin, bin_end);
                fFixedTEMAT->AddBinContent(
                    global_bin, (bin_cont * (new_en_up_edge - axis->GetBinLowEdge(bin_end)) / bin_width) / total_ratio);
            }
            // in case there are more bins in between, fill them here (skipping first already bin_start and bin_end)
            for (int ybin = bin_start + 1; ybin < bin_end; ybin++)
            {
                global_bin = fFixedTEMAT->GetBin(time_bin, ybin);
                fFixedTEMAT->AddBinContent(global_bin, bin_cont / total_ratio);
            }
        }
    }

    fFixedTEMAT->SetEntries(V.TEMAT->GetEntries());
    fFixedTEMAT->ResetStats();
    return fFixedTEMAT;
}

TH2D *CCM::FixMatrix(const TH2D *input_mat, const bool valid_only)
{
    if (!fFitDone)
    {
        std::cerr << "Warning! PerformFits() function was not yet called, it will be called with default parameters "
                     "(valid only, no spline)"
                  << std::endl;
        this->PerformFits();
    }

    TH2D *fixed_mat = dynamic_cast<TH2D *>(input_mat->Clone(Form("%s_corrected", input_mat->GetName())));
    fixed_mat->Reset();
    fixed_mat->SetDirectory(0);

    Double_t bin_cont;
    Double_t total_ratio, en_low_edge, en_up_edge, new_en_low_edge, new_en_up_edge;

    Int_t global_bin, bin_start, bin_end;
    const Double_t bin_width = fixed_mat->GetYaxis()->GetBinWidth(1);

    const TAxis *axis = fixed_mat->GetYaxis();

    std::vector<TGraph> roi_graphs;
    std::vector<TSpline3> roi_splines;
    // replace missing poi
    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        roi_graphs.emplace_back();
        roi_graphs.back().SetName(Form("graph_ROI_%i", ROI_index));
        roi_graphs.back().SetTitle(Form("graph of ROI %i; time slice index; energy shift [energy units]", ROI_index));
        for (int time_index = 0; time_index < V.time_bins; time_index++)
        {
            if (valid_only && !ResVec[ROI_index][time_index].isValid)
                continue;
            roi_graphs.back().AddPoint(this->GetMatrixTime(time_index),
                                       V.ROIs[ROI_index].desired_energy - ResVec[ROI_index][time_index].energy_shift);
        }
    }
    for (auto &gr : roi_graphs)
    {
        roi_splines.emplace_back(Form("spline_%s", gr.GetName()), &gr);
    }

    double timestamp;
    TF1 *fcn{nullptr};
    std::pair<TF1 *, std::string> *correction_fcn{nullptr};

    for (int time_bin = 1; time_bin <= fixed_mat->GetXaxis()->GetNbins(); time_bin++)
    {
        timestamp = input_mat->GetXaxis()->GetBinCenter(time_bin);
        // load function with parameters
        {
            TGraph gr;
            int roi_index = 0;
            for (const TSpline3 &spline : roi_splines)
            {
                if (timestamp < spline.GetXmin())
                {
                    gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(spline.GetXmin()));
                    continue;
                }
                if (timestamp > spline.GetXmax())
                {
                    gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(spline.GetXmax()));
                    continue;
                }

                gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(timestamp));
                roi_index++;
            }
            if (gr.GetN() == 0)
            {
                correction_fcn = this->FindCorrectionFunction(0);
            }
            else
            {
                correction_fcn = this->FindCorrectionFunction(gr.GetN());
                gr.Fit(correction_fcn->first, correction_fcn->second.c_str());
            }
        }

        // set bins
        for (int en_bin = 1; en_bin <= axis->GetNbins(); en_bin++)
        {
            bin_cont = input_mat->GetBinContent(input_mat->GetBin(time_bin, en_bin));
            if (bin_cont == 0)
                continue;
            // we need to get the new corrected energies that are matching low and up edges of the old bin
            en_low_edge = axis->GetBinLowEdge(en_bin);
            en_up_edge = axis->GetBinUpEdge(en_bin);
            new_en_low_edge = correction_fcn->first->Eval(en_low_edge);
            new_en_up_edge = correction_fcn->first->Eval(en_up_edge);

            // the new_en will not cover range that might be multiple bins in width, and there is also a difference in
            // coverage - new low edge will cover only partially the new bin, same for up edge so we need to get ratios
            // on how to distribute the content of the old bin to the new bins

            bin_start = axis->FindBin(new_en_low_edge);
            bin_end = axis->FindBin(new_en_up_edge);
            // in case we are covering just 1 bin, put the content in
            if (bin_start == bin_end)
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_start);
                fixed_mat->AddBinContent(global_bin, bin_cont);
                continue;
            }
            // get to total ratio of "bin equivalent widths" we are covering with the new energy
            total_ratio = static_cast<double>(bin_end - bin_start - 1);
            // std::cout << std::setprecision(6) << "total ratio " << total_ratio << " ";
            total_ratio += (double)(axis->GetBinUpEdge(bin_start) - new_en_low_edge) / (double)bin_width;
            // std::cout << std::setprecision(6) << total_ratio << " ";
            total_ratio += (new_en_up_edge - axis->GetBinLowEdge(bin_end)) / bin_width;
            // std::cout << std::setprecision(6) << total_ratio << std::endl;
            // set manually to new lower bin content
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_start);
                fixed_mat->AddBinContent(global_bin,
                                         (bin_cont * (axis->GetBinUpEdge(bin_start) - new_en_low_edge) / bin_width) /
                                             total_ratio);
            }
            // set manually to new upper bin content
            {
                global_bin = fixed_mat->GetBin(time_bin, bin_end);
                fixed_mat->AddBinContent(
                    global_bin, (bin_cont * (new_en_up_edge - axis->GetBinLowEdge(bin_end)) / bin_width) / total_ratio);
            }
            // in case there are more bins in between, fill them here (skipping first already bin_start and bin_end)
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

void CCM::FixTree(const std::string &tfilename, const std::string &treename, const std::string &e_branchname,
                  const std::string &ts_branchname, const bool valid_only, const int time_subdivision)
{

    // input checks
    TFile *old_file = TFile::Open(tfilename.c_str(), "READ");
    if (!old_file || old_file->IsZombie())
    {
        std::cerr << "Error opening file\n";
        return;
    }
    TTree *old_tree = (TTree *)old_file->Get(treename.c_str());
    if (!old_tree)
    {
        std::cerr << "Error getting tree\n";
        return;
    }
    TBranch *branch_e = old_tree->GetBranch(e_branchname.c_str());
    TBranch *branch_t = old_tree->GetBranch(ts_branchname.c_str());
    if (!branch_e || !branch_t)
    {
        std::cerr << "Branch does not exist\n";
        return;
    }

    // if (!fFitDone)
    // {
    //     std::cerr << "Warning! PerformFits() function was not yet called, it will be called with default parameters "
    //                  "(valid only, no spline)"
    //               << std::endl;
    //     this->PerformFits();
    // }

    // prepare new root file

    // copy all branches but the corrected ones
    old_tree->SetBranchStatus("*", 1);
    old_tree->SetBranchStatus(e_branchname.c_str(), 0);

    // create new file and clone the tree
    std::string new_file = "corrected_" + tfilename.substr(tfilename.find_last_of('/') + 1);
    TFile *new_tfile = new TFile(new_file.c_str(), "RECREATE");
    TTree *new_tree = old_tree->CloneTree();

    // add missing branches in the new tree
    double new_e, old_e;
    Long64_t timestamp;
    TBranch *new_branch_e = new_tree->Branch(e_branchname.c_str(), &new_e);

    // now we need to iterate through the old tree to get the original energy and timestamp, and store corrected values
    // in a new tree
    old_tree->SetBranchStatus("*", 0);
    old_tree->SetBranchStatus(e_branchname.c_str(), 1);
    old_tree->SetBranchStatus(ts_branchname.c_str(), 1);
    old_tree->SetBranchAddress(e_branchname.c_str(), &old_e);
    old_tree->SetBranchAddress(ts_branchname.c_str(), &timestamp);

    std::vector<TGraph> roi_graphs;
    std::vector<TSpline3> roi_splines;

    // replace missing poi
    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        roi_graphs.emplace_back();
        roi_graphs.back().SetName(Form("graph_ROI_%i", ROI_index));
        roi_graphs.back().SetTitle(Form("graph of ROI %i; time slice index; energy shift [energy units]", ROI_index));
        for (int time_index = 0; time_index < V.time_bins; time_index++)
        {
            if (valid_only && !ResVec[ROI_index][time_index].isValid)
                continue;
            roi_graphs.back().AddPoint(this->GetMatrixTime(time_index),
                                       V.ROIs[ROI_index].desired_energy - ResVec[ROI_index][time_index].energy_shift);
        }
    }
    for (auto &gr : roi_graphs)
    {
        roi_splines.emplace_back(Form("spline_%s", gr.GetName()), &gr);
    }
    // prepare some variables
    const double Xmax = V.TEMAT->GetXaxis()->GetXmax();
    const double Xmin = V.TEMAT->GetXaxis()->GetXmin();
    const long Nbins = V.TEMAT->GetXaxis()->GetNbins();
    const double ts_half_bin_width = (Xmax - Xmin) / (double)(2. * Nbins);
    const Long64_t ts_section_width = (Xmax - Xmin) / (Nbins * time_subdivision);
    Long64_t next_section_ts_high = ts_section_width;
    Long64_t next_section_ts_low = 0;

    std::pair<TF1 *, std::string> *correction_fcn{nullptr};

    const Long64_t nentries = old_tree->GetEntries();
    std::cout << "***** Fixing " << treename << " of " << tfilename << "\n";
    for (Long64_t entry = 0; entry < nentries; entry++)
    {
        if (entry % 10000 == 0)
        {
            std::cout << std::setprecision(3) << entry / (double)nentries * 100. << " %    \r";
        }
        old_tree->GetEntry(entry);
        // check if function is set
        if (timestamp > next_section_ts_high || timestamp < next_section_ts_low || !correction_fcn)
        {
            // find new timestamp section
            while (timestamp > next_section_ts_high)
            {
                next_section_ts_low = next_section_ts_high;
                next_section_ts_high += ts_section_width;
            }
            while (timestamp < next_section_ts_low)
            {
                next_section_ts_high = next_section_ts_low;
                next_section_ts_low -= ts_section_width;
            }

            // get the parameters
            // we need to take care not to extrapolate the spline by checking the boundaries
            // we are playing with the ts_half_bin_width value, because as it is, splines uses knots that are in the
            // center of the "time_bin", so we need to consider also events that are less than half of the bin width to
            // still belong to the same bin and get the its parameters
            TGraph gr;
            int roi_index = 0;
            for (const TSpline3 &spline : roi_splines)
            {
                // std::cout << spline.GetXmin() << " " << spline.GetXmax() << std::endl;

                if (timestamp < spline.GetXmin())
                {
                    if (timestamp + ts_half_bin_width > spline.GetXmin())
                    {
                        // add first point of the spline to compensate for the different time of the left edge and
                        // center of the time bin
                        gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(spline.GetXmin()));
                    }
                    continue;
                }
                if (timestamp > spline.GetXmax())
                {
                    if (timestamp < spline.GetXmax() + ts_half_bin_width)
                    {
                        gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(spline.GetXmax()));
                    }
                    continue;
                }

                gr.AddPoint(V.ROIs[roi_index].desired_energy, spline.Eval(timestamp));

                roi_index++;
            }
            if (gr.GetN() == 0)
            {
                correction_fcn = this->FindCorrectionFunction(0);
            }
            else
            {
                correction_fcn = this->FindCorrectionFunction(gr.GetN());
                gr.Fit(correction_fcn->first, correction_fcn->second.c_str());
            }
        }
        // new_e = 40.;
        new_e = correction_fcn->first->Eval(old_e);
        // std::cout << correction_fcn->first->GetName() << std::setprecision(7) << " " << new_e << " " << old_e
        //   << std::endl;
        new_branch_e->Fill();
    }
    std::cout << "\n"
              << "Done\n";

    new_tfile->cd();
    new_tree->Write();
    new_tfile->Close();

    // clean up
    old_file->Close();
}

std::pair<TF1 *, std::string> *CCM::FindCorrectionFunction(const int nrois)
{
    for (auto &f : fCorrectionFunctions)
    {
        if (nrois >= f.first->GetNpar())
        {
            return &f;
        }
    }
    // we failed to find the function (should not happen since we have one with zero parameters!)
    throw std::runtime_error("Error! No suitable function found for given number of ROIs!");
    return nullptr;
}

const ResCont *CCM::GetResultContainer(const int ROI_no, const int time_index) const noexcept
{
    if (ROI_no >= V.number_of_ROIs || time_index >= V.time_bins)
        return nullptr;
    return &ResVec[ROI_no][time_index];
}
void CCM::SetInvalidResult(const int ROI_no, const int time_index)
{
    if (ROI_no >= V.number_of_ROIs || time_index >= V.time_bins)
        throw std::runtime_error("ROI_no or time_index are not correct");
    ResVec[ROI_no][time_index].isValid = false;
}

void CCM::UseGaussianResult()
{
    for (int roi_index = 0; roi_index < V.number_of_ROIs; roi_index++)
    {
        for (int time_index = 0; time_index < V.time_bins; time_index++)
        {
            ResVec[roi_index][time_index].bin_shift =
                ResVec[roi_index][time_index].gfit_mu + V.ROIs[roi_index].base_shift_value;
            ResVec[roi_index][time_index].energy_shift =
                V.TEMAT->GetYaxis()->GetBinWidth(1) * ResVec[roi_index][time_index].bin_shift;
        }
    }
}

void CCM::SetReferenceVector(const unsigned int ROI_index, const std::vector<double> &own_reference_vector)
{
    if (ROI_index > V.number_of_ROIs)
    {
        throw std::runtime_error("ROI index out of bounds");
    }

    int vsize = V.ROIs[ROI_index].bin_window_high - V.ROIs[ROI_index].bin_window_low;

    if (own_reference_vector.size() != vsize)
    {
        throw std::runtime_error("Reference vector size does not match expected sample vector size");
    }
    V.sample_vector[ROI_index].clear();
    V.sample_vector[ROI_index] = own_reference_vector;

    this->Normalize(V.sample_vector[ROI_index]);
}

void CCM::SetReferenceProjection(const TH1 *projection)
{
    V.sample_vector.clear();

    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        std::vector<double> vec{};
        for (int e_bin = V.ROIs[ROI_index].bin_window_low; e_bin < V.ROIs[ROI_index].bin_window_high; e_bin++)
        {
            vec.push_back(projection->GetBinContent(e_bin));
        }
        this->Normalize(vec);
        V.sample_vector.emplace_back(std::move(vec));
    }
}

void CCM::CheckReferenceVector(int ROI_index)
{
    int vsize = V.ROIs[ROI_index].bin_window_high - V.ROIs[ROI_index].bin_window_low;

    if (V.sample_vector[ROI_index].size() != vsize)
    {
        throw std::runtime_error("Reference vector not set for ROI " + std::to_string(ROI_index));
    }
}

void CCM::CheckReferenceVectors()
{
    for (int ROI_index = 0; ROI_index < V.number_of_ROIs; ROI_index++)
    {
        this->CheckReferenceVector(ROI_index);
    }
}