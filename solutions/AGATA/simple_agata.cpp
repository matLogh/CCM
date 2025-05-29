#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Root
#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
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

const double reference_time_bgn = 3700;
const double reference_time_end = 3710;

const std::array<double, 5> ROIarr{2746., 2754., -90., 10., 2749.};

std::shared_ptr<TH2> get_matrix()
{
    TFile *f = TFile::Open(
        (std::string(DATA_PATH) + "/example_data/agata_LNLexp035_0006_07B.root").c_str());
    std::shared_ptr<TH2> mat((TH2F *)f->Get("hE0_TS_07B"));

    return mat;
}

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);

    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::shared_ptr<TH2> TEMAT;

    {
        std::cout << "Fetching matrix...                                   "
                  << std::flush;
        TEMAT                                = get_matrix();
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    std::cout << "Constructing CCM object...                           " << std::flush;
    TF1                           fcn("gain_fcn", "[0]*x", 0, 32000);
    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(
        RegionOfInterest(TEMAT, ROIarr[0], ROIarr[1], ROIarr[2], ROIarr[3], ROIarr[4]));
    // create CCM object
    CCM fix(TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix.SetCorrectionFunction(fcn, "");

    auto duration =
        duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
    std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
              << " seconds" << std::endl;

    {
        std::cout << "Calculating energy shifts...                         "
                  << std::flush;
        fix.CalculateEnergyShifts(8);
        auto duration =
            duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Switching precise shift calculation to polynomial... "
                  << std::flush;
        // fix.UseMaxDPResult();
        fix.UsePolynomialResult();
        // fix.UseGaussianResult();
        auto duration =
            duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Calculating correction fits...                       "
                  << std::flush;
        fix.CalculateCorrectionFits();
        auto duration =
            duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    std::shared_ptr<TH2> TEMAT_fixed;
    {
        std::cout << "Fixing input matrix...                               "
                  << std::flush;
        TEMAT_fixed = fix.FixMatrix();
        auto duration =
            duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
        t1 = high_resolution_clock::now();
    }

    {
        std::cout << "Creating output files (ROOT file, shift&fit table)..."
                  << std::flush;
        fix.SaveToRootFile();
        fix.SaveShiftTable("simpleAgata_shiftTable_0006_07B.txt");
        fix.SaveFitTable("simpleAgata_fitTable_0006_07B.txt");
        auto duration =
            duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
        std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
                  << " seconds" << std::endl;
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
    auto dot_product = fix.GetDotProductGraph(0, 1849);
    dot_product->SetMarkerColor(kBlue);
    dot_product->SetMarkerStyle(20);
    dot_product->SetMarkerSize(0.5);
    dot_product->Draw("ALP");

    duration = duration_cast<microseconds>(high_resolution_clock::now() - t1).count();
    std::cout << "done in " << std::setprecision(2) << (double)duration / 1e6
              << " seconds" << std::endl;
    t1 = high_resolution_clock::now();

    duration = duration_cast<seconds>(high_resolution_clock::now() - t0).count();
    std::cout << "\nTOTAL DURATION OF " << duration << " seconds" << std::endl;

    app.Run();
    return 0;
}