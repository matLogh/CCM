// #include <array>
#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Root
#include "TApplication.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"

#include "CCM.h"
// #include "CheckCCM.h"
#include "Cross_correlation.h"
#include "TheuerkaufPeak.hpp"
#include "variables.h"

using namespace std::chrono;
using namespace TEC;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

struct ccm_settings
{
    double cost{std::numeric_limits<double>::max()};
    uint temat_rebin_x{1};
    uint temat_rebin_y{1};
    // std::vector<Region_of_interest> ROIs{};
    bool use_gaussian{true}; // gaussian or polynomial dot product fit
    bool valid_only{true};
    std::string interpolator_type{"AKIMA"};
    bool interpolator_smoothing{false};
    std::string smoother_type{""};
    double smoother_par0{0};
    double smoother_par1{0};
    double smoother_par2{0};

    void print_header(std::ostream &os)
    {
        os << "cost, temat_rebin_x, temat_rebin_y, ROIs, use_gaussian, valid_only, interpolator_type, "
              "interpolator_smoothing, smoother_type, smoother_par0, smoother_par1, smoother_par2 \n ";
    }

    void print_values(std::ostream &os)
    {
        os << this->cost << " " << this->temat_rebin_x << " " << this->temat_rebin_y << " " << std::boolalpha
           << this->use_gaussian << " " << this->valid_only << " " << this->interpolator_type << " "
           << this->interpolator_smoothing << " " << this->smoother_type << " " << this->smoother_par0 << " "
           << this->smoother_par1 << " " << this->smoother_par2 << "\n";
    }
};

CCM *GLOBAL_MINIMIZATION_CCM{nullptr};

double get_fwfm(TH1 *histo, const double center, const double min, const double max)
{
    TheuerkaufFitter fitter(min, max);
    fitter.AddPeak(center, true, false, false);
    fitter.Fit(histo, "OUTPUT_NONE");
    return fitter.GetPeak(0)->GetFWxM(5);
}

double ccm_optimizer(CCM input_ccm)
{
    return 0;
}

std::vector<ccm_settings> set_minimizer_grid()
{
    std::vector<ccm_settings> grid;
    const std::vector<int> rebin_values = {5};
    // const std::vector<int> rebin_values = {1, 2, 5, 10};
    for (auto rebin_x : rebin_values)
    {
        for (auto rebin_y : rebin_values)
        {
            ccm_settings s;
            s.temat_rebin_x = rebin_x;
            s.temat_rebin_y = rebin_y;
            grid.push_back(s);
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.use_gaussian = false;
            //     grid.push_back(st);
            // }
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.interpolator_type = "CSPLINE";
            //     grid.push_back(st);
            // }
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.interpolator_type = "CSPLINE_PERIODIC";
            //     grid.push_back(st);
            // }
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.interpolator_type = "POLYNOMIAL";
            //     grid.push_back(st);
            // }
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.interpolator_type = "AKIMA_PERIODIC";
            //     grid.push_back(st);
            // }
            {
                ccm_settings st;
                st = s;
                st.interpolator_smoothing = true;
                st.smoother_type = "KernelSmoother";
                // st.smoother_par0 = 1.;
                // grid.push_back(st);

                // st.smoother_par0 = 2.;
                // grid.push_back(st);

                // st.smoother_par0 = 5.;
                // grid.push_back(st);

                st.smoother_par0 = 20.;
                grid.push_back(st);
            }
            // {
            //     ccm_settings st;
            //     st = s;
            //     st.interpolator_smoothing = true;
            //     st.smoother_type = "SmoothShifts_Lowess";
            //     st.smoother_par0 = 0.1;
            //     grid.push_back(st);

            //     st.smoother_par0 = 0.2;
            //     grid.push_back(st);

            //     st.smoother_par0 = 0.4;
            //     grid.push_back(st);

            //     st.smoother_par0 = 0.6;
            //     grid.push_back(st);

            //     st.smoother_par0 = 0.8;
            //     grid.push_back(st);

            //     st.smoother_par0 = 0.9;
            //     grid.push_back(st);
            // }
        }
    }
    // grid.back()
    return grid;
}

int main(int argc, char **argv)
{
    TFile savefile("savefile.root", "RECREATE");
    // TApplication theApp("App", 0, 0);
    auto settings = set_minimizer_grid();

    TFile *matfile = TFile::Open("/home/mbalogh/to_win/out_huge_00A_1010.root", "READ");
    TH2D *TEMAT_original = (TH2D *)matfile->Get("hE0_TS_00A");
    TF1 fcn("gain_fcn", "[0]*x", 0, 4000);
    double reference_time_bgn = 3100;
    double reference_time_end = 3200;

    int counter = 0;
    for (auto &s : settings)
    {
        std::cout << "Working on settings: " << ++counter << "/" << settings.size() << std::endl;
        TH2D *TEMAT = (TH2D *)TEMAT_original->Clone(Form("hE0_TS_00A_rebinned_%i", counter));
        TEMAT->RebinX(s.temat_rebin_x);
        TEMAT->RebinY(s.temat_rebin_y);

        std::vector<Region_of_interest> ROIs;
        ROIs.emplace_back(Region_of_interest(*TEMAT, 2200., 2250., -20., 20., 2223.)); // ROI1 is the region of interest

        CCM ccm_fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
        ccm_fix.SetCorrectionFunction(fcn, "");
        ccm_fix.CalculateEnergyShifts(8);
        if (s.use_gaussian)
            ccm_fix.UseGaussianResult();
        else
            ccm_fix.UsePolynomialResult();

        std::cout << std::endl << "smoothing " << s.interpolator_smoothing << std::endl;
        if (!s.interpolator_type.empty())
        {
            for (int roi_index = 0; roi_index < ROIs.size(); roi_index++)
            {
                ccm_fix.ConfigureShiftInterpolator(roi_index, s.interpolator_type, s.valid_only);
                ccm_fix.EnableInterpolation(roi_index);
            }
        }
        if (s.interpolator_smoothing)
        {
            for (int roi_index = 0; roi_index < ROIs.size(); roi_index++)
            {
                if (s.smoother_type == "KernelSmoother")
                {
                    std::cout << "KernelSmoother " << s.smoother_par0 << std::endl;
                    ccm_fix.SmoothShifts_KernelSmoother(roi_index, s.smoother_par0);
                }
                else if (s.smoother_type == "SmoothShifts_Lowess")
                {
                    std::cout << "SmoothShifts_Lowess " << s.smoother_par0 << std::endl;

                    ccm_fix.SmoothShifts_Lowess(roi_index, s.smoother_par0);
                }
            }
            ccm_fix.CalculateCorrectionFits();

            TApplication theApp("App", 0, 0);
            ccm_fix.GetInterpolationGraph(0)->Draw("ALP");
            auto *TEMAT_large_new = ccm_fix.FixMatrix(TEMAT_original);
            auto *proj = TEMAT_large_new->ProjectionY();

            s.cost = get_fwfm(proj, 7917, 7880, 7940);
            std::cout << "Cost: " << std::setprecision(8) << s.cost << std::endl;
            theApp.Run();
        }

        auto *TEMAT_large_new = ccm_fix.FixMatrix(TEMAT_original);
        auto *proj = TEMAT_large_new->ProjectionY();

        s.cost = get_fwfm(proj, 7917, 7880, 7940);
        std::cout << "Cost: " << std::setprecision(8) << s.cost << std::endl;
        savefile.cd();
        proj->Write(Form("proj_%i", counter));
        ccm_fix.GetInterpolationGraph(0)->Write(Form("interpol_%i", counter));
        // TEMAT_large_new->Write(Form("TEMAT_large_new_%i", counter));
        delete TEMAT;
        delete TEMAT_large_new;
        // delete proj;
    }

    auto min_it = std::min_element(settings.begin(), settings.end(),
                                   [](const ccm_settings &a, const ccm_settings &b) { return a.cost < b.cost; });

    std::cout << "Best settings: " << std::endl;
    min_it->print_values(std::cout);

    std::fstream outfile("minimizer_grid.dat", std::ios::out);
    settings[0].print_header(outfile);

    for (auto &s : settings)
    {
        s.print_values(outfile);
    }

    TApplication theApp("App", 0, 0);

    // return 1;

    TH2D *TEMAT = (TH2D *)TEMAT_original->Clone(Form("hE0_TS_00A_rebinned_%i", counter));
    TEMAT->RebinX(2);
    TEMAT->RebinY(2);

    std::vector<Region_of_interest> ROIs;
    ROIs.emplace_back(Region_of_interest(*TEMAT, 2200., 2250., -20., 20., 2223.)); // ROI1 is the region of interest

    CCM ccm_fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
    ccm_fix.SetCorrectionFunction(fcn, "");
    ccm_fix.CalculateEnergyShifts(8);
    ccm_fix.UseGaussianResult();

    ccm_fix.CalculateCorrectionFits(5);

    auto *TEMAT_new = ccm_fix.FixMatrix();
    auto *proj = TEMAT_new->ProjectionY();
    proj->SetName("proj1");
    proj->SetLineColor(kRed);

    auto *TEMAT_large_new = ccm_fix.FixMatrix(TEMAT_original);
    // TEMAT_large_new->Draw("COLZ");
    // new TCanvas();
    // TEMAT->Draw("colz");
    // theApp.Run();
    // auto *TEMAT_large_new = ccm_fix.FixMatrix(TEMAT);
    TEMAT_large_new->RebinX(5);
    TEMAT_large_new->RebinY(2);
    auto proj_large = TEMAT_large_new->ProjectionY();
    proj_large->SetName("proj2");
    proj_large->SetLineColor(kGreen);

    new TCanvas();
    proj->Draw("");
    proj_large->Draw("SAME");
    // TEMAT_original->ProjectionY()->Draw("SAME");
    TEMAT->ProjectionY()->Draw("SAME");
    // new TCanvas();
    // projections.emplace_back(proj);

    double fwfm_proj = get_fwfm(proj, 7917, 7880, 7940);
    double fwfm_proj_large = get_fwfm(proj_large, 7917, 7880, 7940);

    std::cout << std::setprecision(8) << "proj " << fwfm_proj << " proj_large " << fwfm_proj_large << std::endl;

    auto *gr_interpol = ccm_fix.GetInterpolationGraph(0);
    // gr_interpol->Draw("ALP");

    ccm_fix.SmoothShifts_KernelSmoother(0, 20.);
    auto *TEMAT_large_smooth = ccm_fix.FixMatrix(TEMAT_original);
    TEMAT_large_smooth->RebinX(5);
    TEMAT_large_smooth->RebinY(2);
    auto proj_fixed_smooth = TEMAT_large_smooth->ProjectionY();
    proj_fixed_smooth->SetName("proj_fixed_smooth");
    proj_fixed_smooth->SetLineColor(kOrange);
    proj_fixed_smooth->Draw("SAME");

    double fwfm_proj_smooth = get_fwfm(proj_fixed_smooth, 7917, 7880, 7940);
    std::cout << "fwfm_proj_smooth " << fwfm_proj_smooth << std::endl;

    auto *gr_smooth = ccm_fix.GetInterpolationGraph(0);
    // new TCanvas();
    // gr_smooth->Draw("ALP");
    // theApp.Run();

    gr_smooth->SetMarkerStyle(23);
    gr_smooth->SetLineColor(4);
    gr_smooth->SetLineWidth(1);
    gr_smooth->SetFillStyle(0);
    gr_smooth->SetDrawOption("ACP");

    // for (int i = 20; i < ccm_fix.GetNumberOfTimeIndices() / 2.; i++)
    // {
    //     ccm_fix.SetInvalidResult(0, i);
    // }

    // ccm_fix.SaveToRootFile("elia_roi1.root");

    auto *gr_points = ccm_fix.GetROIShifts(0);
    gr_points->SetMarkerStyle(22);
    gr_points->SetLineColor(4);
    gr_points->SetLineWidth(1);
    gr_points->SetFillStyle(0);
    gr_points->SetDrawOption("P");

    gr_interpol->SetMarkerStyle(22);
    gr_interpol->SetMarkerColor(2);
    gr_interpol->SetDrawOption("AP");
    gr_interpol->SetLineColor(3);
    gr_interpol->SetLineWidth(2);
    gr_interpol->SetFillStyle(0);

    auto mg = new TMultiGraph;

    mg->Add(gr_points);
    mg->Add(gr_smooth);
    // mg->Add(gr_interpol);
    new TCanvas();

    mg->Draw("APL");

    theApp.Run();

    return 0;
}