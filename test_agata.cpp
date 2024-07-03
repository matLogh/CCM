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
#include "variables.h"

using namespace std::chrono;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

int main(int argc, char **argv)
{
    TApplication theApp("App", 0, 0);

    // high_resolution_clock::time_point t1 = high_resolution_clock::now();
    // high_resolution_clock::time_point t2 = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(t2 - t1).count();
    // std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds" << std::endl;

    TFile *matfile = TFile::Open("/home/mbalogh/to_win/out_huge_00A_1010.root", "READ");
    TH2D *TEMAT = (TH2D *)matfile->Get("hE0_TS_00A");
    // TEMAT->RebinX(5);
    // TEMAT->RebinY(2);

    TF1 fcn("gain_fcn", "[0]*x", 0, 4000);

    std::vector<Region_of_interest> ROIs;
    ROIs.emplace_back(Region_of_interest(*TEMAT, 2200., 2250., -20., 20., 2223.)); // ROI1 is the region of interest
    // Set reference vector
    // double scale_factor = (time_end - time_start) / (double)matrix_time_bins;
    double reference_time_bgn = 3100;
    double reference_time_end = 3200;

    // create CCM object

    std::vector<TH1 *> projections;
    // t->Draw();

    CCM fix1(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix1.SetCorrectionFunction(fcn, "");
    fix1.CalculateEnergyShifts(8);
    fix1.UseGaussianResult();

    auto *gr_interpol = fix1.GetInterpolationGraph(0);
    gr_interpol->Draw("ALP");

    fix1.SmoothShifts_KernelSmoother(0, 8.);

    auto *gr_smooth_lowess = fix1.GetInterpolationGraph(0);
    new TCanvas();
    gr_smooth_lowess->Draw("ALP");
    // theApp.Run();

    gr_smooth_lowess->SetMarkerStyle(23);
    gr_smooth_lowess->SetLineColor(4);
    gr_smooth_lowess->SetLineWidth(1);
    gr_smooth_lowess->SetFillStyle(0);
    gr_smooth_lowess->SetDrawOption("ACP");

    // for (int i = 20; i < fix1.GetNumberOfTimeIndices() / 2.; i++)
    // {
    //     fix1.SetInvalidResult(0, i);
    // }

    fix1.CalculateCorrectionFits();
    // auto *TEMAT_new = fix1.FixMatrix();
    // auto *TEMAT_new = fix1.FixMatrix();
    // auto *TEMAT_new = fix1.FixMatrix(TEMAT_full);
    // auto *proj = TEMAT_new->ProjectionY();
    // proj->SetName("proj1");
    // proj->SetLineColor(kRed);
    // TEMAT->ProjectionY()->Draw();
    // proj->Draw("SAME");
    // projections.emplace_back(proj);
    // fix1.SaveToRootFile("elia_roi1.root");
    // new TCanvas();
    auto *gr_points = fix1.GetROIShifts(0);
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
    mg->Add(gr_smooth_lowess);
    // mg->Add(gr_interpol);
    new TCanvas();

    mg->Draw("APL");

    theApp.Run();

    return 0;
}