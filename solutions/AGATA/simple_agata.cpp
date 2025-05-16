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

const double reference_time_bgn = 3710;
const double reference_time_end = 3711;

std::shared_ptr<TH2F> get_matrix()
{
    TFile *f = TFile::Open(
        (std::string(DATA_PATH) + "/example_data/agata_LNLexp035_0006_07B.root").c_str());
    std::shared_ptr<TH2F> mat((TH2F *)f->Get("hE0_TS_07B"));

    return mat;
}

int main(int argc, char **argv)
{
    TApplication         app("app", 0, 0);
    auto                 temat = get_matrix();
    std::shared_ptr<TH2> TEMAT(temat->Rebin2D(1, 1));

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    TF1 fcn("gain_fcn", "[0]*x", 0, 32000);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(TEMAT, 2740, 2760, -100., 100., 2749. * 0.98));

    // create CCM object
    CCM fix(TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix.SetCorrectionFunction(fcn, "");
    fix.CalculateEnergyShifts(8);
    // fix.UseMaxDPResult();
    fix.UsePolynomialResult();
    // fix.UseGaussianResult();
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
    TEMAT->GetYaxis()->SetRangeUser(2600, 2800);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    c0.cd(2);
    TEMAT_fixed->GetYaxis()->SetRangeUser(2600, 2800);
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
    TEMAT->GetYaxis()->SetRangeUser(2600, 2800);
    TEMAT->SetTitle("Original matrix");
    TEMAT->Draw("COLZ");
    cx.cd(2);
    shifts->Draw("ALP");

    TCanvas cxx("c_dot_product", "Dot product - simple_example", 800, 600);
    auto    dot_product = fix.GetDotProductGraph(0, 1848);
    dot_product->Draw("ALP");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration                        = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds"
              << std::endl;

    app.Run();
    return 0;
}