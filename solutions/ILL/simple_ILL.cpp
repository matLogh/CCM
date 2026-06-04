#include <algorithm>
#include <array>
#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Root
#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TTree.h"

#include "CCM.h"
#include "Cross_correlation.h"
#include "variables.h"

using namespace std::chrono;
using namespace TEC;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

const double reference_time_bgn = 247800;
const double reference_time_end = 247804;

// const std::array<double, 5>              ROIarr{2746., 2754., -90., 10., 2749.};
std::vector<TEC::RegionOfInterest> gROIarrs;

std::shared_ptr<TH2> get_matrix()
{
    TFile               *f = TFile::Open(std::string("~/data/frank.root").c_str());
    std::shared_ptr<TH2> mat((TH2F *)f->Get("RunVsEnergy_247253_248623/RunVsEnergy_det7"));

    return mat;
}

void save_shift_profile_views(CCM &fix, const std::string &filename = "simpleILL_shiftProfiles.root")
{
    TFile output(filename.c_str(), "RECREATE");
    if (!output.IsOpen()) { throw std::runtime_error("Cannot create " + filename); }

    TDirectory *profile_dir = output.mkdir("shift_profiles");
    if (!profile_dir) { throw std::runtime_error("Cannot create shift_profiles directory in " + filename); }
    profile_dir->cd();

    const int n_time_bins = static_cast<int>(fix.GetNumberOfTimeIndices());
    for (int xbin = 1; xbin <= n_time_bins; ++xbin)
    {
        auto profile = fix.GetShiftProfile(xbin - 1);
        profile->SetName(Form("shift_profile_xbin_%04d", xbin));
        profile->SetTitle(Form("Shift profile for TEMAT X bin %d;ROI energy;bin shift", xbin));
        profile->Write();
    }

    output.cd();
    output.Close();
}

void save_shift_profile_animation(CCM &fix, const std::string &filename = "simpleILL_shiftProfiles.gif")
{
    TCanvas canvas("c_shift_profile_animation", "Shift profile animation", 900, 650);
    canvas.SetCrosshair(1);

    const auto matrix      = fix.GetInputMatrix();
    const int  n_time_bins = static_cast<int>(fix.GetNumberOfTimeIndices());
    for (int xbin = 1; xbin <= n_time_bins; ++xbin)
    {
        auto profile = fix.GetShiftProfile(xbin - 1);
        profile->SetTitle(Form("Shift profile for TEMAT X bin %d, x = %.0f;ROI energy;bin shift", xbin, matrix->GetXaxis()->GetBinCenter(xbin)));
        profile->SetMarkerColor(kBlue);
        profile->SetMarkerStyle(20);
        profile->SetMarkerSize(1.0);
        profile->SetLineColor(kBlue);

        canvas.cd();
        profile->Draw("ALP");
        canvas.Modified();
        canvas.Update();

        const std::string print_name = filename + (xbin == n_time_bins ? "++" : "+");
        canvas.Print(print_name.c_str());
    }
}

int main(int argc, char **argv)
{
    TApplication app("app", &argc, argv);

    std::array<double, 5> ROIarr{19600, 19300, 19800, -400, 400};

    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::shared_ptr<TH2> TEMAT;

    {
        std::cout << "Fetching matrix...                                   " << std::flush;
        TEMAT         = get_matrix();
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    gROIarrs.emplace_back(TEMAT, 19300, 19800, -400, 400, 19600);
    gROIarrs.emplace_back(TEMAT, 3800, 3900, -150, 150, 3845);
    gROIarrs.emplace_back(TEMAT, 16800, 17000, -300, 300, 16900);
    gROIarrs.emplace_back(TEMAT, 10250, 10450, -200, 200, 10350);
    gROIarrs.emplace_back(TEMAT, 6650, 6750, -100, 100, 6708);
    gROIarrs.emplace_back(TEMAT, 1750, 1850, -100, 100, 1800);

    std::cout << "Constructing CCM object...                           " << std::flush;
    TF1 fcn("lin_fcn", "[0] + [1]*x", 0, 32000);
    // TF1 fcn("gain_fcn", "[0]*x", 0, 32000);
    // create CCM object
    CCM fix(TEMAT, gROIarrs, reference_time_bgn, reference_time_end);
    fix.SetCorrectionFunction(fcn, "");
    {
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Calculating energy shifts...                         " << std::flush;
        fix.CalculateEnergyShifts(8);
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Switching precise shift calculation to polynomial... " << std::flush;
        // fix.UseMaxDPResult();
        fix.UsePolynomialResult();
        // fix.UseGaussianResult();
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    std::shared_ptr<TH1F> profile_rchi2  = std::make_shared<TH1F>("shift_profile_rchi2", "rchi2 of shift profiles", 10000, 0, 100);
    std::shared_ptr<TH1F> profile_maxres = std::make_shared<TH1F>("shift_profile_maxres", "max residual of shift profiles", 10000, 0, 100);

    // evaluate "linarity" of the ROI shifts
    std::vector<size_t>    profiles_to_remove;
    const float            max_residual_threshold = 3.f;
    const float            max_dp_threshold       = 0.6f;
    const std::vector<int> good_roi_indices       = {0, 1}; // index of ROI that is expected to be linear and will be used for correction

    {
        std::cout << "Evaluating linearity of shifts...                   " << std::flush;
        const size_t n_time_bins = fix.GetNumberOfTimeIndices();
        const size_t n_rois      = fix.GetNumberOfROIs();
        for (size_t t_index = 0; t_index < n_time_bins; t_index++)
        {
            auto          profile = fix.GetShiftProfile(t_index);
            TFitResultPtr fit     = profile->Fit(&fcn, "Q0S");
            float         rchi2   = fit->Chi2() / fit->Ndf();
            profile_rchi2->Fill(rchi2);
            float max_res = 0;
            for (int i = 1; i <= profile->GetN(); i++)
            {
                float res = std::abs(profile->GetY()[i - 1] - fcn.Eval(profile->GetX()[i - 1]));
                if (res > max_res) { max_res = res; }
            }
            profile_maxres->Fill(max_res);
            if (max_res > max_residual_threshold) { profiles_to_remove.push_back(t_index); }

            // check max dp and decide if the profile should be removed or not
            bool is_valid = true;
            for (auto roi_index : good_roi_indices)
            {
                auto *res = fix.GetResultContainer(roi_index, t_index);
                if (res->dp > max_dp_threshold) { is_valid = false; }
            }
            if (!is_valid) { profiles_to_remove.push_back(t_index); }
        }
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
    }

    // remove extra ROIs
    {
        const size_t n_time_bins = fix.GetNumberOfTimeIndices();
        const size_t n_rois      = fix.GetNumberOfROIs();
        for (size_t t_index = 0; t_index < n_time_bins; t_index++)
        {
            for (size_t roi_index = 0; roi_index < n_rois; roi_index++)
            {
                if (std::find(good_roi_indices.begin(), good_roi_indices.end(), roi_index) == good_roi_indices.end())
                {
                    fix.SetResultStatus(roi_index, t_index, false);
                }
            }
            if (std::find(profiles_to_remove.begin(), profiles_to_remove.end(), t_index) != profiles_to_remove.end())
            {
                for (int good_roi_index : good_roi_indices) { fix.SetResultStatus(good_roi_index, t_index, false); }
            }
        }
    }

    {
        std::cout << "Calculating correction fits...                       " << std::flush;
        fix.CalculateCorrectionFits();
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    // remove non linear profiles - give them -1 correction gain
    {
        fix.RemoveInvalidProjections();
    }

    std::shared_ptr<TH2> TEMAT_fixed;
    {
        std::cout << "Fixing input matrix...                               " << std::flush;
        TEMAT_fixed   = fix.FixMatrix();
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Creating output files (ROOT file, shift&fit table)..." << std::flush;
        fix.SaveToRootFile();
        fix.SaveShiftTable("simpleILL_shiftTable.txt");
        fix.SaveFitTable("simpleILL_fitTable.txt");
        save_shift_profile_views(fix);
        // save_shift_profile_animation(fix);
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    std::cout << "Creating and drawing nice plots, histograms etc...   " << std::flush;
    auto *proj_old = TEMAT->ProjectionY();
    proj_old->SetLineColor(kRed);

    auto *proj_fixed = TEMAT_fixed->ProjectionY();
    proj_fixed->SetTitle("TEMAT energy projection");
    proj_fixed->SetLineColor(kBlue);

    TCanvas c0("c_mat", "Matrices - simple_example", 1200, 600);
    c0.SetCrosshair(1);
    c0.Divide(1, 2);
    c0.cd(1);
    TEMAT->GetYaxis()->SetRangeUser(ROIarr[0] + ROIarr[2], ROIarr[1] + ROIarr[3]);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    c0.cd(2);
    TEMAT_fixed->GetYaxis()->SetRangeUser(ROIarr[0] + ROIarr[2], ROIarr[1] + ROIarr[3]);
    TEMAT_fixed->SetTitle("Corrected matrix");
    TEMAT_fixed->Draw("COLZ");

    TCanvas c1("c_hist", "Histograms - simple_example", 800, 600);
    c1.SetLogy();
    c1.SetCrosshair(1);
    proj_fixed->Draw();
    proj_old->Draw("SAME");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(proj_fixed, "After correction", "l");
    legend->AddEntry(proj_old, "Before correction", "l");
    legend->Draw("SAME");

    TCanvas c2("c_shifts", "Shifts - simple_example", 800, 600);
    auto    shifts = fix.GetROIShifts(0);
    shifts->SetTitle("ROI 0: calculated shifts");
    shifts->SetMarkerColor(kBlue);
    shifts->SetMarkerStyle(20);
    shifts->SetMarkerSize(0.5);
    shifts->Draw("ALP");

    TCanvas c3("c_interpolation", "Interpolation - simple_example", 800, 600);
    auto    interpol = fix.GetInterpolationGraph(0, 10);
    interpol->SetTitle("ROI 0: 10x interpolation");
    interpol->SetMarkerColor(kRed);
    interpol->SetMarkerStyle(20);
    interpol->SetMarkerSize(0.5);
    interpol->Draw("ALP");

    TCanvas c4("c_shifts", "Shifts - simple_example", 800, 600);
    c4.Divide(1, 2);
    c4.SetCrosshair(1);
    c4.cd(1);
    TEMAT->GetYaxis()->SetRangeUser(ROIarr[0] + ROIarr[2], ROIarr[1] + ROIarr[3]);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    c4.cd(2);
    shifts->Draw("ALP");

    TCanvas c5("c_dot_product", "Dot product - simple_example", 800, 600);
    c5.SetCrosshair(1);
    auto dot_product = fix.GetDotProductGraph(0, 0);
    dot_product->SetMarkerColor(kBlue);
    dot_product->SetMarkerStyle(20);
    dot_product->SetMarkerSize(0.5);
    dot_product->Draw("ALP");

    TCanvas c6("c_shift_profile", "Shift profile - simple_example", 800, 600);
    c6.SetCrosshair(1);
    auto shift_profile = fix.GetShiftProfile(1);
    shift_profile->SetMarkerColor(kBlue);
    shift_profile->SetMarkerStyle(20);
    shift_profile->SetMarkerSize(0.5);
    shift_profile->Draw("ALP");

    TCanvas c7("c_shift_profile_rchi2", "Rchi2 of shift profiles - simple_example", 800, 600);
    c7.SetCrosshair(1);
    profile_rchi2->SetLineColor(kRed);
    profile_rchi2->Draw("hist");

    TCanvas c8("c_shift_profile_max_res", "Max residual of shift profiles - simple_example", 800, 600);
    c8.SetCrosshair(1);
    profile_maxres->SetLineColor(kRed);
    profile_maxres->Draw("hist");

    {
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();

        duration = duration_cast<seconds>(high_resolution_clock::now() - t0).count();
        std::cout << "\nTOTAL DURATION OF " << duration << " seconds" << std::endl;
    }

    app.Run();
    return 0;
}
