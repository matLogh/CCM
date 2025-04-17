#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Root
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>

#include "CCM.h"

using namespace std::chrono;
using namespace TEC;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

int    bins       = 1000;
double time_start = 22E12;
double time_end   = 50E12;

const std::string INPUT_ROOT_FILE =
    // "/home/mbalogh/data/56ni/temat_merged_compressed.root";
    "/home/mbalogh/data/56ni/temat_merged.root";

// const std::string MATRIX_NAME = "strips/summed_compressed_temat_strip_0";
const std::string MATRIX_NAME = "strips/summed_temat_strip_0";

TH2D *get_matrix()
{
    TFile *infile = TFile::Open(INPUT_ROOT_FILE.c_str(), "READ");

    if (!infile || infile->IsZombie())
    {
        throw std::runtime_error(
            Form("Error: Could not open input ROOT file: %s", INPUT_ROOT_FILE.c_str()));
    }

    TH2D *mat = (TH2D *)infile->Get(MATRIX_NAME.c_str());
    // std::cout << mat->GetName() << std::endl;
    if (mat)
    {
        mat->SetDirectory(0);
        infile->Close();
        return mat;
    }
    else
    {
        throw std::runtime_error(
            Form("Error: Could not open matrix %s", MATRIX_NAME.c_str()));
    }
}

TH1D *get_model_histogram(TH2D *matrix, double time_start, double time_end)
{
    int bin_start = matrix->GetXaxis()->FindBin(time_start);
    int bin_end   = matrix->GetXaxis()->FindBin(time_end);

    TH1D *proj = matrix->ProjectionY("model_projection", bin_start, bin_end);
    proj->SetName("model");
    proj->SetTitle("model");
    proj->SetLineColor(kBlack);
    return proj;
}

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);
    TH2D        *TEMAT = get_matrix();

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(*TEMAT, 5500, 8000, -100, 1500, 6166));
    // ROIs.emplace_back(RegionOfInterest(*TEMAT, 3000, 4300, -500, 1500, 3900));

    const double reference_time_bgn = 4.9e12;
    const double reference_time_end = 5.35e12;

    // create CCM object

    for (int i = 0; i < 1000; i++)
    {
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        CCM fix(*TEMAT, ROIs, reference_time_bgn, reference_time_end);

        fix.SetCorrectionFunction(fcn, "");
        fix.CalculateEnergyShifts(12);

        fix.DisableInterpolation(0);
        fix.CalculateCorrectionFits();
        auto *TEMAT_fixed = fix.FixMatrix();

        auto *proj_old =
            TEMAT->ProjectionY("proj_old", 1030, TEMAT->GetXaxis()->GetNbins());
        proj_old->SetLineColor(kRed);

        auto *proj_fixed =
            TEMAT_fixed->ProjectionY("proj_fixed", 1030, TEMAT->GetXaxis()->GetNbins());
        proj_fixed->SetLineColor(kBlue);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds"
                  << std::endl;
    }
    return 0;
}