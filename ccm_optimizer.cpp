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

    TFile *matfile = TFile::Open("/home/mbalogh/to_win/out_ext_1010.root", "READ");
    TH2D *TEMAT = (TH2D *)matfile->Get("coreE0_TS/hE0_TS_0");
    TEMAT->RebinX(20);

    TF1 fcn("gain_fcn", "[0]*x", 0, 4000);

    std::vector<Region_of_interest> ROIs;
    ROIs.emplace_back(Region_of_interest(*TEMAT, 2200., 2250., -20., 20., 2223.)); // ROI1 is the region of interest
    // Set reference vector
    // double scale_factor = (time_end - time_start) / (double)matrix_time_bins;
    double reference_time_bgn = 3100;
    double reference_time_end = 3200;

    // create CCM object

    std::vector<TH1 *> projections;
    auto *t = new TBrowser;
    // t->Draw();

    CCM fix1(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix1.SetCorrectionFunction(fcn, "");
    fix1.CalculateEnergyShifts(1);
    fix1.PerformFits();
    auto *TEMAT_new = fix1.FixMatrix();
    auto *proj = TEMAT_new->ProjectionY();
    proj->SetName("proj1");
    proj->SetLineColor(kRed);
    proj->Draw("SAME");
    projections.emplace_back(proj);
    fix1.SaveToRootFile("elia_roi1.root");

    ROIs.emplace_back(Region_of_interest(*TEMAT, 7600., 7700., -20., 20., 7639.)); // ROI1 is the region of interest

    CCM fix2(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix2.SetCorrectionFunction(fcn, "");
    fix2.CalculateEnergyShifts(1);
    fix2.PerformFits();
    TEMAT_new = fix2.FixMatrix();
    proj = TEMAT_new->ProjectionY();
    proj->SetName("proj2");
    proj->SetLineColor(kGreen);
    proj->Draw("SAME");
    projections.emplace_back(proj);
    fix2.SaveToRootFile("elia_roi2.root");

    // TEMAT->GetYaxis()->SetRangeUser(7600, 7700);
    // TEMAT_new->GetYaxis()->SetRangeUser(7600, 7700);
    projections.emplace_back(TEMAT->ProjectionY());
    projections.back()->Draw("SAME");

    TFile *outfile = TFile::Open("out.root", "RECREATE");
    for (auto &p : projections)
    {
        p->Write();
    }

    theApp.Run();

    return 0;
}