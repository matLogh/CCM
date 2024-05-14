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
// #include "CheckCCM.h"
#include "Cross_correlation.h"
#include "variables.h"

#include "TheuerkaufPeak.hpp"

using namespace std::chrono;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

int bins = 1000;
double time_start = 22E12;
double time_end = 50E12;

double reference_time_bgn = 3.919e13;
double reference_time_end = 3.9468e13;

TH2D *get_matrix(const int t_bins, const int e_bins)
{
    TFile *file = TFile::Open("ccm_matrix.root", "UPDATE");
    std::string mat_name = "temat_t_" + std::to_string(t_bins) + "_e_" + std::to_string(e_bins);
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
        path += "/data/DecayGammaSpectroscopy_timeUnstable_" + std::to_string(i) + ".root";
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
    mat = new TH2D(mat_name.c_str(), mat_name.c_str(), t_bins, time_start, time_end, e_bins, 0, 1600);
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

void build_fix_tree()
{
    TFile *file = TFile::Open("tree_fix_matrix.root", "RECREATE");
    // Get tree from rootfiles
    TChain *tree = new TChain("data_tree");
    // char str[100];
    for (int i = 3; i < 14; i++)
    {
        // this is a bit chaotic due to size limitation of the github repository
        //  sprintf(str,"data/DecayGammaSpectroscopy_timeUnstable_%i.root",i);
        std::string path = DATA_PATH;
        path = "corrected_DecayGammaSpectroscopy_timeUnstable_" + std::to_string(i) + ".root";
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

    std::cout << "creating matrix " << std::endl;
    Long64_t nentries = tree->GetEntries();
    TH2D *TEMAT = new TH2D("dte", "distorted T vs E matrix", bins, time_start, time_end, 3200, 0, 1600);
    for (Long64_t entry = 0; entry < nentries; entry++)
    {
        tree->GetEntry(entry);
        TEMAT->Fill(t, e);
    }
    file->cd();
    TEMAT->Write();
    file->Close();
}

void example_fix_tree()
{
    auto *TEMAT = get_matrix(1000, 3200);
    std::vector<Region_of_interest> ROIs;
    ROIs.emplace_back(Region_of_interest(*TEMAT, 1445., 1475., -35., 25., 1480.)); // ROI1 is the region of interest
    CCM fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);

    TF1 fcn("gain_fcn", "[0]*x", 0, 4000);

    // fix.SelectFitFunction(1, 1);
    // START CCM
    fix.SetCorrectionFunction(fcn, "");
    fix.CalculateEnergyShifts(8);
    fix.SaveShiftTable();
    fix.PerformFits();
    fix.FixMatrix();
    fix.SaveToRootFile();

    for (int i = 3; i < 14; i++)
    {
        std::string filename =
            "/home/mbalogh/CCM/data/DecayGammaSpectroscopy_timeUnstable_" + std::to_string(i) + ".root";
        fix.FixTree(filename.c_str(), "data_tree", "e", "t", true, 10);
    }
    build_fix_tree();
}
#include "fstream"

int main(int argc, char **argv)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    const std::vector<double> tbins{10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000};
    const std::vector<double> ebins{1600, 3200, 4800, 6400};

    std::ofstream output;
    output.open("fwhm.dat");

    TheuerkaufFitter fitter(1400, 1500);
    fitter.SetBacgroundPoly(3);
    fitter.AddPeak(1456, true, false, false);

    int best_time_hinning;
    int best_energy_binning;
    double best_fwhm = std::numeric_limits<double>::max();

    for (double t : tbins)
    {
        for (double e : ebins)
        {
            std::cout << "T " << t << " E " << e << std::endl;
            TH2D *TEMAT = get_matrix(t, e);

            std::vector<Region_of_interest> ROIs;
            ROIs.emplace_back(
                Region_of_interest(*TEMAT, 1445., 1475., -35., 25., 1480.)); // ROI1 is the region of interest

            CCM fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);
            TF1 fcn("gain_fcn", "[0]*x", 0, 4000);

            fix.SetCorrectionFunction(fcn, "");
            fix.CalculateEnergyShifts(8);

            // fix.SaveShiftTable();
            fix.PerformFits();
            auto *fixed_matrix = fix.FixMatrix();
            auto *proj = fixed_matrix->ProjectionY();
            // proj->Draw();
            fix.SaveFitTable("fit_table.dat", "detector_name");
            fitter.Fit(proj, "OUTPUT_NONE");
            output << t << " " << e << " " << fitter.GetPeak(0)->GetFWHM() << "\n";
            // fitter.Analyze(proj);
            if (fitter.GetPeak(0)->GetFWHM() < best_fwhm)
            {
                best_fwhm = fitter.GetPeak(0)->GetFWHM();
                best_time_hinning = t;
                best_energy_binning = e;
            }
            delete TEMAT;
        }
    }

    std::cout << "Minimum FWHM " << best_fwhm << " is reached for " << best_time_hinning << " time bins and "
              << best_energy_binning << " energy bins" << std::endl;

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds" << std::endl;
    return 0;
}