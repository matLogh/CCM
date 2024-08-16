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

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);

    TFile *f = TFile::Open((std::string(DATA_PATH) + "/example_data/aoq_vs_xfp.root").c_str());
    TH2D *mat = (TH2D *)f->Get("Aoqnew_xfp_38");

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // define your correction functions here, but remember to define them using TFormula, otherwise you will get an
    // error

    // TF1 fcn1("offset_fcn", "[0] + x", 0, 1);
    // TF1 analytical_fcn("analytical_fcn", "[0]*sqrt([1]*[1] - 1)", 0, 1);
    // TF1 fcn2("lin_fcn", "[0] + [1]*x", 0, 1);
    // TF1 fcn3("quad_rcn", "sqrt(x)*[0] + [1] + [2]*x", 0, 1);

    // define ROIs

    std::vector<double> peaks{

        267.77102, 270.61866, 273.31970, 276.20607, 278.90865, 281.79899, 284.74140, 287.62261, 290.64785,
        293.73331, 296.92020, 300.04220, 303.19280, 306.63186, 310.07831, 313.41412, 317.00486, 320.54116};

    std::vector<Region_of_interest> ROIs;

    for (int i = 0; i < peaks.size(); i++)
    {
        ROIs.emplace_back(Region_of_interest(*mat, peaks[i] - 2, peaks[i] + 2, -1.5, 1.5, peaks[i]));
    }

    // create CCM object
    CCM fix(*mat, ROIs, 727., 731.);

    // Set correction functions, you can have multiple functions as a fallback scenarios

    // fix.SeTECectionFunction(fcn1, "");
    // fix.SeTECectionFunction(analytical_fcn, "");
    fix.SetCorrectionFunction(TF1("pol1", "pol1", 0, 1), "");
    // fix.SetFallbackCorrectionFunction(TF1("pol2", "pol2", 0, 1), "");
    // fix.SetFallbackCorrectionFunction(TF1("pol3", "pol3", 0, 1), "");
    // fix.SetFallbackCorrectionFunction(TF1("pol4", "pol4", 0, 1), "");
    // fix.SetFallbackCorrectionFunction(fcn2, "");
    // fix.SetFallbackCorrectionFunction(fcn1, "");

    // calculate the offsets
    fix.CalculateEnergyShifts(8);

    // just a way how the final dot product is obtained
    fix.UseGaussianResult();
    // shift table
    fix.SaveShiftTable();

    // filter for valid results
    for (int t = 0; t < fix.GetNumberOfTimeIndices(); t++)
    {
        int goodroi = fix.GetNumberOfROIs();

        for (int roi = 0; roi < fix.GetNumberOfROIs(); roi++)
        {
            const auto *result = fix.GetResultContainer(roi, t);
            if (result->gfit_sigma < 1 || result->gfit_sigma > 5 || result->dp < 0.6)
            {
                fix.SetInvalidResult(roi, t);
                goodroi--;
            }
        }
    }
    // calculate the calibration functions
    fix.CalculateCorrectionFits();
    fix.SaveFitTable();

    // get corrected matrix
    auto *mat_fixed = fix.FixMatrix();
    fix.SaveToRootFile("prisma_abberations.root");

    // just to draw some shift profiles: for given "time" (x-axis slice) draw the ROI_center vs the calculated shift
    // this is useful to understand what correction function you should use
    fix.GetShiftProfile(497)->Draw("ALP");
    new TCanvas();
    fix.GetShiftProfile(440)->Draw("ALP");
    new TCanvas();
    fix.GetShiftProfile(540)->Draw("ALP");
    new TCanvas();

    mat_fixed->Draw("COLZ");
    new TCanvas();
    mat->Draw("COLZ");
    new TCanvas();
    auto *proj_old = mat->ProjectionY();
    proj_old->SetLineColor(kRed);

    auto *proj_fixed = mat_fixed->ProjectionY();
    proj_fixed->SetLineColor(kBlue);

    proj_fixed->Draw();
    proj_old->Draw("SAME");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds" << std::endl;

    app.Run();
    return 0;
}