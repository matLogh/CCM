#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Root
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"

#include "CCM.h"
#include "Cross_correlation.h"
#include "variables.h"

using namespace std::chrono;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

int    bins       = 1000;
double time_start = 22E12;
double time_end   = 50E12;

double reference_time_bgn = 3.919e13;
double reference_time_end = 3.9468e13;

TH2D *get_matrix(const int t_bins, const int e_bins)
{
    TFile      *file = TFile::Open("ccm_matrix.root", "UPDATE");
    std::string mat_name =
        "temat_t_" + std::to_string(t_bins) + "_e_" + std::to_string(e_bins);
    TH2D *mat = (TH2D *)file->Get(mat_name.c_str());
    if (mat)
    {
        mat->SetDirectory(0);
        file->Close();
        return mat;
    }
    // Get tree from rootfiles
    TChain *tree = new TChain("data_tree");
    // char str[100];
    for (int i = 3; i < 14; i++)
    {
        // this is a bit chaotic due to size limitation of the github repository
        std::string path = DATA_PATH;
        path += "/example_data/DecayGammaSpectroscopy_timeUnstable_" + std::to_string(i) +
                ".root";
        tree->Add(path.c_str());
    }

    Double_t e;
    Long64_t t;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("e", 1);
    tree->SetBranchStatus("t", 1);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("t", &t);

    // Fill time-energy matrix (TEMAT) with tree events
    // Binning of the TEMAT (time vs energy matrix) is crucial!

    std::cout << "creating matrix " << mat_name << std::endl;
    const Long64_t nentries = tree->GetEntries();
    mat = new TH2D(mat_name.c_str(), mat_name.c_str(), t_bins, time_start, time_end,
                   e_bins, 0, 1600);
    for (Long64_t entry = 0; entry < nentries; entry++)
    {
        tree->GetEntry(entry);
        mat->Fill(t, e);
    }
    file->cd();
    mat->Write();
    mat->SetDirectory(0);
    file->Close();
    return mat;
}

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);
    TH2D        *TEMAT = get_matrix(1000, 3200);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(*TEMAT, 1445., 1475., -35., 25., 1480.));

    // create CCM object
    CCM fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
    fix.SeTECectionFunction(fcn, "");
    fix.CalculateEnergyShifts(8);
    fix.PerformFits();
    auto *TEMAT_fixed = fix.FixMatrix();
    fix.SaveToRootFile();

    auto *proj_old = TEMAT->ProjectionY();
    proj_old->SetLineColor(kRed);

    auto *proj_fixed = TEMAT_fixed->ProjectionY();
    proj_fixed->SetLineColor(kBlue);

    proj_fixed->Draw();
    proj_old->Draw("SAME");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration                        = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds"
              << std::endl;

    app.Run();
    return 0;
}