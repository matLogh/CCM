#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include "CCM.h"

#include "common.cpp"
#include <thread>

using namespace TEC;

std::vector<int>         gRUNLIST;
std::vector<std::string> gCRYSTALLIST;
std::string              gDIR            = "timeEvo";
float                    gSHIFTTHRESHOLD = 0.5; // keV
std::vector<float>       gROIarr;
bool                     gDrawCanvases = false;

std::vector<std::shared_ptr<TObject>> guiObjects;

std::array<double, 6> calculate_statistics(const std::vector<float> &values)
{
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Standard deviation
    double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
    double stdev  = std::sqrt(sq_sum / values.size() - mean * mean);

    // RMS
    double rms = std::sqrt(sq_sum / values.size());

    // Max absolute value
    double max_abs =
        *std::max_element(values.begin(), values.end(),
                          [](double a, double b) { return std::abs(a) < std::abs(b); });

    // Median
    std::vector<float> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    double median = sorted[sorted.size() / 2];

    double above_threshold =
        std::count_if(values.begin(), values.end(),
                      [](double value) { return std::abs(value) > gSHIFTTHRESHOLD; });

    return {above_threshold, mean, stdev, rms, max_abs, median};
}

std::array<float, 2> get_ref_time(std::shared_ptr<TH2> TEMAT)
{
    static const int     bin_window_width = 5;
    std::array<float, 2> ref_time{-1., -1.};

    int nbinsX = TEMAT->GetNbinsX();
    int nbinsY = TEMAT->GetNbinsY();
    for (int start_bin = nbinsX / 4; start_bin < nbinsX - bin_window_width; start_bin++)
    {
        auto one_bin_integral = TEMAT->Integral(start_bin, start_bin + 1, 1, nbinsY);
        auto width_integral =
            TEMAT->Integral(start_bin, start_bin + bin_window_width, 1, nbinsY);

        // lets make sure that there are some events in the time bin and if the width
        // contain at least the same amount of events
        if (one_bin_integral > 0 &&
            width_integral / bin_window_width * 1.2 > one_bin_integral)
        {
            ref_time.at(0) = TEMAT->GetXaxis()->GetBinLowEdge(start_bin);
            ref_time.at(1) =
                TEMAT->GetXaxis()->GetBinUpEdge(start_bin + bin_window_width);
            break;
        }
        {
            ref_time.at(0) = TEMAT->GetXaxis()->GetBinLowEdge(start_bin);
            ref_time.at(1) =
                TEMAT->GetXaxis()->GetBinUpEdge(start_bin + bin_window_width);
            break;
        }
    }

    return ref_time;
}

bool detect_time_evolution(std::shared_ptr<TH2> temat,
                           std::vector<float>  &over_threshold_values)
{

    std::shared_ptr<TH2> TEMAT(temat->Rebin2D(2, 2, Form("%s_rebin", temat->GetName())));
    TEMAT->SetDirectory(0);
    auto ref_time = get_ref_time(TEMAT);
    if (ref_time.at(0) < 0 || ref_time.at(1) < 0)
    {
        std::cerr << "Error! Could not find suitable reference time in the matrix "
                  << TEMAT->GetName() << " in " << gDirectory->GetPath() << std::endl;
        return false;
    }

    // Create a vector of RegionOfInterest objects
    std::vector<RegionOfInterest> rois;
    rois.emplace_back(RegionOfInterest(TEMAT, gROIarr.at(1), gROIarr.at(2), gROIarr.at(3),
                                       gROIarr.at(4), gROIarr.at(0)));
    CCM ccm(TEMAT, rois, ref_time.at(0), ref_time.at(1));
    ccm.CalculateEnergyShifts(8);

    over_threshold_values.clear();
    for (int i = 0; i < ccm.GetNumberOfTimeIndices(); i++)
    {
        auto res = ccm.GetResultContainer(0, i);
        if (res->energy_shift > gSHIFTTHRESHOLD)
        {
            over_threshold_values.push_back(res->energy_shift);
        }
    }

    if (over_threshold_values.size() > 0)
    {
        if (gDrawCanvases)
        {
            std::string canvas_name = "canvas_" + std::string(TEMAT->GetName());
            std::string canvas_title =
                "TimeEvo DETECTED in " + std::string(TEMAT->GetName());
            std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>(
                canvas_name.c_str(), canvas_title.c_str(), 1600, 900);
            c->SetCrosshair(1);
            c->Divide(1, 2);

            c->cd(2);
            auto _m = ccm.GetInputMatrix();
            _m->GetYaxis()->SetRangeUser(gROIarr.at(1), gROIarr.at(2));
            _m->Draw("COLZ");
            guiObjects.emplace_back(_m);
            c->cd(1);
            auto gr = ccm.GetROIShifts(0, false);
            gr->GetXaxis()->SetRangeUser(_m->GetXaxis()->GetXmin(),
                                         _m->GetXaxis()->GetXmax());
            gr->Draw("ALP");

            guiObjects.emplace_back(gr.release());
            guiObjects.emplace_back(c);
        }
        return true;
    }
    return false;
}

void print_help()
{
    std::cout << "This code uses some assumptions and hardcoded values to run CCM - "
                 "namely it uses a reference time in around 1/3 of the matrix that is 2x "
                 "in energy and 2x in time\n";
    std::cout << "If a shift larger than set threshold (default 0.5) is detected, a "
                 "TimeEvo is reported.\n";
    std::cout << "Two reports are made: one while code is running showing also some "
                 "statistics, and one at the end sorted by run number.\n";
    std::cout << "You can use --draw option to draw matrices with calculated shifts "
                 "while the code is running.\n\n";

    std::cout << "Usage: detectTimeEvo_AGATA [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h                 Show this help message\n";
    std::cout << "  --crystal <1> <...>        Specify the crystal(s) name\n";
    std::cout << "  --allcrys                  Run for all crystals of EXP_035\n";
    std::cout << "  --run <1> <...>            Specify the run number\n";
    std::cout << "  --shift_threshold <1>      Energy threshold, if energy shift value "
                 "threshold is found the timeEvo is reported (default 0.5)\n";
    std::cout << "  --draw                     Use this flag to enable drawing the "
                 "matrices that are over the set threshold\n";

    std::cout << "  --ROI <1> <2> <3> <4> <5>  Specify the Region of "
                 "Interest (ROI) as:\n"
              << "                                <1> - desired energy of "
                 "the ROI\n"
              << "                                <2> - left edge of ROI n\n"
              << "                                <3> - right edge of ROI\n"
              << "                                <4> - shift ROI by "
                 "maximum of <4> to "
                 "the LEFT (neg value!)\n"
              << "                                <5> - shift ROI by "
                 "maximum of <5> to "
                 "the RIGHT\n";

    std::cout
        << "  --ROIsource <1>            Define ROI for calibration sources. Currently "
           "recognized are: 60Co, 66Ga, 133Ba, 226Ra \n";

    std::cout << "  --dir <1>                  Set directory in which to search for "
                 "matrices\n";

    std::cout << std::endl << std::endl;
}

void parse_args(int argc, char **argv)
{
    if (argc < 2)
    {
        print_help();
        throw std::invalid_argument("No arguments provided");
    }
    for (int i = 1; i < argc; ++i)
    {
        // gCOST_PEAK
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h")
        {
            print_help();
            exit(0);
        }
        else if (arg == "--run")
        {
            if (i + 1 < argc)
            {
                while (i + 1 < argc && std::isdigit(argv[i + 1][0]))
                {
                    gRUNLIST.emplace_back(std::stoi(argv[++i]));
                }
                if (gRUNLIST.empty())
                {
                    throw std::invalid_argument(
                        "--run must be followed by at least one integer value");
                }
            }
            else { throw std::invalid_argument("Missing value for --run"); }
        }
        else if (arg == "--crys" || arg == "--crystal" || arg == "--crystals")
        {
            auto _crys = parse_space_separated_crystals(i, argc, argv);
            // avoid duplicates
            for (const auto &cry : _crys)
            {
                if (std::find(gCRYSTALLIST.begin(), gCRYSTALLIST.end(), cry) ==
                    gCRYSTALLIST.end())
                {
                    gCRYSTALLIST.emplace_back(cry);
                }
            }
        }
        else if (arg == "--allcrys")
        {
            std::vector<std::string> _c = {
                "00A", "00B", "00C", "01A", "01C", "02A", "02B", "02C",
                "04A", "04B", "04C", "05B", "05C", "06A", "06B", "06C",
                "07A", "07B", "08A", "08B", "09A", "09B", "09C", "10A",
                "10B", "10C", "11A", "11B", "11C", "14A", "14B", "14C"};
            for (const auto &cry : _c)
            {
                if (std::find(gCRYSTALLIST.begin(), gCRYSTALLIST.end(), cry) ==
                    gCRYSTALLIST.end())
                {
                    gCRYSTALLIST.emplace_back(cry);
                }
            }
        }
        else if (arg == "--dir")
        {
            if (i + 1 < argc) { gDIR = argv[++i]; }
            else { throw std::invalid_argument("Missing value for --dir"); }
        }
        else if (arg == "--shift_threshold")
        {
            if (i + 1 < argc)
            {
                try
                {
                    gSHIFTTHRESHOLD = std::stof(argv[++i]);
                }
                catch (const std::invalid_argument &)
                {
                    throw std::invalid_argument(
                        "Invalid float value for --shift_threshold");
                }
            }
            else { throw std::invalid_argument("Missing value for --shift_threshold"); }
        }
        else if (arg == "--ROI")
        {
            gROIarr = parse_space_separated_floats(i, argc, argv, 5);
            if (gROIarr.size() != 5)
            {
                throw std::invalid_argument(
                    "Invalid number of arguments for --ROI. Expected 5 values.");
            }
        }
        else if (arg == "--ROIsource")
        {
            std::vector<float> peak;
            if (i + 1 < argc) { parse_ROI_source(argv[++i], gROIarr, peak); }
            else { throw std::invalid_argument("Missing value for --ROIsource"); }
        }
        else if (arg == "--draw") { gDrawCanvases = true; }

        else
        {
            print_help();
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }
    if (gROIarr.size() != 5)
    {
        throw std::invalid_argument(
            "Invalid number of arguments for --ROI. Expected 5 values.");
        exit(10);
    }

    std::cout << "Parameters used are:" << std::endl;
    std::cout << "Run number(s):    ";
    for (const auto &run : gRUNLIST) { std::cout << run << " "; }
    std::cout << std::endl;
    std::cout << "Crystals:         ";
    for (const auto &cry : gCRYSTALLIST) { std::cout << cry << " "; }
    std::cout << std::endl;
    std::cout << "Data directory:        " << gDIR << std::endl;
    std::cout << "Energy shift threshold: " << gSHIFTTHRESHOLD << std::endl;
    std::cout << "ROI:              ";
    for (const auto &roi : gROIarr) { std::cout << roi << " "; }
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
}

int main(int argc, char **argv)
{
    TApplication app("app", 0, 0);
    parse_args(argc, argv);

    std::vector<std::pair<int, std::string>> detected_timeEvo;

    for (const auto &run : gRUNLIST)
    {
        for (const auto &crystal : gCRYSTALLIST)
        {
            gSystem->ProcessEvents();

            auto  rfname = get_rootfilename(gDIR, run, crystal);
            TFile matfile(rfname.c_str(), "READ");
            if (matfile.IsZombie())
            {
                throw std::runtime_error("Error! could not open/find the ROOT file " +
                                         rfname + " file");
            }
            std::string matrix_name = "hE0_TS_" + crystal;
            TH2        *raw         = (TH2 *)matfile.Get(matrix_name.c_str());
            if (!raw) throw std::runtime_error("Matrix not found: " + matrix_name);
            std::shared_ptr<TH2> TEMAT(raw);
            TEMAT->SetDirectory(0);
            std::vector<float> over_threshold_value;

            if (detect_time_evolution(TEMAT, over_threshold_value))
            {
                if (detected_timeEvo.size() == 0)
                {
                    std::cout << "TimeEvolution detected!" << std::endl;
                }
                auto stats = calculate_statistics(over_threshold_value);
                std::cout << std::fixed << std::setprecision(2) << "  Run " << run
                          << " crys " << crystal << ": "
                          << "above thr: " << (int)stats.at(0) << " mean: " << stats.at(1)
                          << " stdev: " << stats.at(2) << " rms: " << stats.at(3)
                          << " max_abs: " << stats.at(4) << " median: " << stats.at(5)
                          << std::endl;
                detected_timeEvo.emplace_back(std::make_pair(run, crystal));
            }
            matfile.Close();
        }
    }

    if (detected_timeEvo.size() != 0)
    {
        std::cout << "\n\n    *****Summary*****\n\nTimeEvolution shift above "
                  << gSHIFTTHRESHOLD << " threshold was detected for:";

        std::sort(detected_timeEvo.begin(), detected_timeEvo.end(),
                  [](const auto &a, const auto &b) {
                      if (a.first != b.first) return a.first < b.first;
                      else
                      {
                          int num_a = std::stoi(a.second.substr(0, 2));
                          int num_b = std::stoi(b.second.substr(0, 2));
                          if (num_a != num_b) return num_a < num_b;
                          else { return a.second[2] < b.second[2]; }
                      }
                  });
        int _this_run = -1;
        detected_timeEvo.front().first;
        for (const auto &suspects : detected_timeEvo)
        {
            if (_this_run != suspects.first)
            {
                _this_run = suspects.first;
                std::cout << "\n  Run " << _this_run << ":";
            }
            std::cout << " " << suspects.second;
        }
        std::cout << std::endl;
    }
    if (guiObjects.size() != 0) { app.Run(); }

    return 0;
}