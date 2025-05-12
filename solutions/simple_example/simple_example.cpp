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

int    bins       = 1000;
double time_start = 22E12;
double time_end   = 50E12;

double reference_time_bgn = 1300;
double reference_time_end = 1320;

std::shared_ptr<TH2F> get_matrix()
{
    TFile *f = TFile::Open(
        (std::string(DATA_PATH) + "/example_data/simple_example.root").c_str());
    std::shared_ptr<TH2F> mat((TH2F *)f->Get("temat"));

    return mat;
}

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);
    auto         TEMAT = get_matrix();

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(TEMAT, 1445., 1475., -35., 25., 1462.));

    // create CCM object
    CCM fix(TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix.SetCorrectionFunction(fcn, "");
    fix.CalculateEnergyShifts(8);
    fix.CalculateCorrectionFits();
    auto TEMAT_fixed = fix.FixMatrix();
    fix.SaveToRootFile();

    auto *proj_old = TEMAT->ProjectionY();
    proj_old->SetLineColor(kRed);

    auto *proj_fixed = TEMAT_fixed->ProjectionY();
    proj_fixed->SetTitle("TEMAT energy projection");
    proj_fixed->SetLineColor(kBlue);

    TCanvas c0("c_mat", "Matrices - simple_example", 1200, 600);
    c0.Divide(1, 2);
    c0.cd(1);
    TEMAT->GetYaxis()->SetRangeUser(1400, 1520);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    c0.cd(2);
    TEMAT_fixed->GetYaxis()->SetRangeUser(1400, 1520);
    TEMAT_fixed->SetTitle("Corrected matrix");
    TEMAT_fixed->Draw("COLZ");

    TCanvas c1("c_hist", "Histograms - simple_example", 800, 600);
    c1.SetLogy();
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

    TCanvas cx("c_shifts", "Shifts - simple_example", 800, 600);
    cx.Divide(1, 2);
    cx.cd(1);
    TEMAT->GetYaxis()->SetRangeUser(1400, 1520);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    cx.cd(2);
    shifts->Draw("ALP");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration                        = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds"
              << std::endl;

    app.Run();
    return 0;
}