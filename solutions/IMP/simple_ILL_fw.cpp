#include <chrono> //measure time
#include <algorithm>
#include <array>
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
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSystem.h"
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
    TFile               *f = TFile::Open(std::string("HistAll.root").c_str());
    std::shared_ptr<TH2> mat((TH2D *)f->Get("hetl10"));

    return mat;
}

std::unique_ptr<TH2D> make_shift_profile_heatmap(CCM &fix)
{
    std::unique_ptr<TGraph> first_profile;
    const int               n_time_bins = static_cast<int>(fix.GetNumberOfTimeIndices());
    for (int xbin = 1; xbin <= n_time_bins; ++xbin)
    {
        first_profile = fix.GetShiftProfile(xbin - 1);
        if (first_profile->GetN() > 0) { break; }
    }

    const int n_points = first_profile ? first_profile->GetN() : 0;
    if (n_points == 0) { throw std::runtime_error("Cannot build shift-profile heatmap from empty shift profiles"); }

    std::vector<std::pair<double, double>> profile_points;
    profile_points.reserve(static_cast<size_t>(n_points));
    for (int point = 0; point < n_points; ++point)
    {
        double energy;
        double shift;
        first_profile->GetPoint(point, energy, shift);
        profile_points.emplace_back(energy, shift);
    }
    std::sort(profile_points.begin(), profile_points.end(), [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

    std::vector<double> energies(static_cast<size_t>(n_points));
    for (int point = 0; point < n_points; ++point) { energies.at(static_cast<size_t>(point)) = profile_points.at(static_cast<size_t>(point)).first; }

    std::vector<double> energy_edges(static_cast<size_t>(n_points) + 1);
    if (n_points == 1)
    {
        energy_edges.at(0) = energies.front() - 0.5;
        energy_edges.at(1) = energies.front() + 0.5;
    }
    else
    {
        energy_edges.front() = energies.front() - 0.5 * (energies.at(1) - energies.front());
        for (int point = 1; point < n_points; ++point)
        {
            energy_edges.at(static_cast<size_t>(point)) = 0.5 * (energies.at(static_cast<size_t>(point - 1)) + energies.at(static_cast<size_t>(point)));
        }
        energy_edges.back() = energies.back() + 0.5 * (energies.back() - energies.at(static_cast<size_t>(n_points - 2)));
    }

    auto heatmap = std::make_unique<TH2D>("shift_profile_heatmap",
                                               "Shift profiles;ROI energy;TEMAT X bin;bin shift",
                                               n_points,
                                               energy_edges.data(),
                                               n_time_bins,
                                               0.5,
                                               static_cast<double>(n_time_bins) + 0.5);
    heatmap->SetDirectory(nullptr);

    for (int xbin = 1; xbin <= n_time_bins; ++xbin)
    {
        auto profile = fix.GetShiftProfile(xbin - 1);
        for (int point = 0; point < profile->GetN(); ++point)
        {
            double energy;
            double shift;
            profile->GetPoint(point, energy, shift);
            heatmap->SetBinContent(heatmap->GetXaxis()->FindBin(energy), xbin, shift);
        }
    }

    return heatmap;
}

void save_shift_profile_views(CCM &fix, const std::string &filename = "simpleIMP_shiftProfiles.root")
{
    TFile output(filename.c_str(), "RECREATE");
    if (!output.IsOpen()) { throw std::runtime_error("Cannot create " + filename); }

    auto heatmap = make_shift_profile_heatmap(fix);
    heatmap->Write();

    TCanvas heatmap_canvas("c_shift_profile_heatmap_lego", "Shift profile heatmap - lego", 1000, 700);
    heatmap_canvas.SetTheta(35);
    heatmap_canvas.SetPhi(35);
    heatmap->Draw("LEGO2Z");
    heatmap_canvas.Write();

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

void save_shift_profile_animation(CCM &fix, const std::string &filename = "simpleIMP_shiftProfiles.gif")
{
    gSystem->Unlink(filename.c_str());

    TCanvas canvas("c_shift_profile_animation", "Shift profile animation", 900, 650);
    canvas.SetCrosshair(1);

    const auto matrix      = fix.GetInputMatrix();
    const int  n_time_bins = static_cast<int>(fix.GetNumberOfTimeIndices());
    for (int xbin = 1; xbin <= n_time_bins; ++xbin)
    {
        auto profile = fix.GetShiftProfile(xbin - 1);
        profile->SetTitle(Form("Shift profile for TEMAT X bin %d, x = %.0f;ROI energy;bin shift",
                               xbin,
                               matrix->GetXaxis()->GetBinCenter(xbin)));
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

    //gROIarrs.emplace_back(TEMAT, 19300, 19800, -400, 400, 19600);
    gROIarrs.emplace_back(TEMAT, 18400, 18600, -400, 400, 18510);
    //gROIarrs.emplace_back(TEMAT, 16800, 17000, -300, 300, 16900);
    //gROIarrs.emplace_back(TEMAT, 10250, 10450, -200, 200, 10350);
    //gROIarrs.emplace_back(TEMAT, 6650, 6750, -100, 100, 6708);
    //gROIarrs.emplace_back(TEMAT, 3800, 3900, -100, 100, 3845);
    //gROIarrs.emplace_back(TEMAT, 1750, 1850, -100, 100, 1800);

    std::cout << "Constructing CCM object...                           " << std::flush;
    TF1 fcn("gain_fcn", "[0]*x", 0, 32000);
    // create CCM object
    CCM fix(TEMAT, gROIarrs, reference_time_bgn, reference_time_end);
    fix.SetCorrectionFunction(fcn, "");
    std::cout << "Number of TEMAT X bins: " << fix.GetNumberOfTimeIndices() << std::endl;
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

    {
        std::cout << "Calculating correction fits...                       " << std::flush;
        fix.CalculateCorrectionFits();
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6 << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
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
        fix.SaveToRootFile("simpleIMP_fw_output.root");
        fix.SaveShiftTable("simpleIMP_fw_shiftTable.txt");
        fix.SaveFitTable("simpleIMP_fw_fitTable.txt");
        save_shift_profile_views(fix, "simpleIMP_fw_shiftProfiles.root");
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

    TCanvas c7("c_shift_profile_heatmap", "Shift profile heatmap - simple_example", 1000, 700);
    c7.SetCrosshair(1);
    c7.SetTheta(35);
    c7.SetPhi(35);
    auto shift_profile_heatmap = make_shift_profile_heatmap(fix);
    shift_profile_heatmap->Draw("LEGO2Z");

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
