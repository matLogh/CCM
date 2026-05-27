#include <algorithm>
#include <array>
#include <atomic>
#include <csignal>
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

std::string        gROOTFILE      = "";
int                gSTART_RUN     = -1;
int                gEND_RUN       = -1;
int                gDETECTOR      = -1;
bool               gALLDET        = false;
std::vector<int>   gDETECTORLIST;
std::string        gROOT_DIR_NAME = "";
double             gSHIFTTHRESHOLD = 0.5; // keV
std::vector<float> gROIarr;
bool               gDrawCanvases = false;

std::vector<std::shared_ptr<TObject>> guiObjects;

std::atomic<bool> gStopExecution{false};
std::atomic<bool> gTerminate_program{false};

std::string get_ill_directory_name()
{
    return "RunVsEnergy_" + std::to_string(gSTART_RUN) + "_" + std::to_string(gEND_RUN);
}

std::string get_ill_matrix_name(const int detector)
{
    return "RunVsEnergy_det" + std::to_string(detector);
}

void signal_handler(int signal)
{
    if (signal == SIGINT)
    {
        // if we are drawing canvases we will stop execution loop and jump to running
        // TApplication
        if (!gStopExecution.load() && gDrawCanvases)
        {
            gStopExecution = true;
            std::cout << "\nCtrl+C detected. Stopping execution loop and no new "
                         "histograms will be drawn!"
                      << std::endl;
        }
        // if we are not drawing canvases or ctrl+c was pressed twice code terminates
        else
        {
            std::cout << "\nCtrl+C detected again. Terminating program..." << std::endl;
            exit(0);
        }
    }
}

void setup_signal_handling() { std::signal(SIGINT, signal_handler); }

std::array<double, 6> calculate_statistics(const std::vector<float> &values)
{
    double mean = std::accumulate(values.begin(), values.end(), 0.0) /
                  static_cast<double>(values.size());

    // Standard deviation
    double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
    double stdev  = std::sqrt(sq_sum / static_cast<double>(values.size()) - mean * mean);

    // RMS
    double rms = std::sqrt(sq_sum / static_cast<double>(values.size()));

    // Max absolute value
    double max_abs =
        *std::max_element(values.begin(), values.end(),
                          [](double a, double b) { return std::abs(a) < std::abs(b); });

    // Median
    std::vector<float> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    double median = sorted[sorted.size() / 2];

    double above_threshold =
        static_cast<double>(std::count_if(values.begin(), values.end(), [](double value) {
            return std::abs(value) > gSHIFTTHRESHOLD;
        }));

    return {above_threshold, mean, stdev, rms, max_abs, median};
}

std::array<float, 2> get_ref_time(std::shared_ptr<TH2> TEMAT)
{
    static const int     bin_window_width = 5;
    std::array<float, 2> ref_time{-1., -1.};

    int nbinsX = TEMAT->GetNbinsX();
    // int nbinsY = TEMAT->GetNbinsY();
    for (int start_bin = nbinsX / 3; start_bin < nbinsX - bin_window_width; start_bin++)
    {
        auto one_bin_integral = TEMAT->Integral(
            start_bin, start_bin, TEMAT->GetYaxis()->FindBin(gROIarr.at(1)),
            TEMAT->GetYaxis()->FindBin(gROIarr.at(2)));
        auto width_integral = TEMAT->Integral(start_bin, start_bin + bin_window_width,
                                              TEMAT->GetYaxis()->FindBin(gROIarr.at(1)),
                                              TEMAT->GetYaxis()->FindBin(gROIarr.at(2)));

        // lets make sure that there are some events in the time bin and if the width
        // contain at least the same amount of events
        if (one_bin_integral > 0 &&
            width_integral / (double)bin_window_width * 0.9 > one_bin_integral)
        {
            ref_time.at(0) =
                static_cast<float>(TEMAT->GetXaxis()->GetBinLowEdge(start_bin));
            ref_time.at(1) = static_cast<float>(
                TEMAT->GetXaxis()->GetBinUpEdge(start_bin + bin_window_width));
            break;
        }
    }

    return ref_time;
}

bool detect_time_evolution(std::shared_ptr<TH2> temat,
                           std::vector<float>  &over_threshold_values)
{

    std::shared_ptr<TH2> TEMAT(temat->Rebin2D(1, 1, Form("%s_rebin", temat->GetName())));
    TEMAT->SetDirectory(0);
    auto ref_time = get_ref_time(TEMAT);
    if (ref_time.at(0) < 0 || ref_time.at(1) < 0)
    {
        if (!gALLDET)
        {
            std::cerr << "Error! Could not find suitable reference time in the matrix "
                      << TEMAT->GetName() << " in " << gDirectory->GetPath() << std::endl;
        }
        return false;
    }

    // Create a vector of RegionOfInterest objects
    std::vector<RegionOfInterest> rois;
    rois.emplace_back(RegionOfInterest(TEMAT, gROIarr.at(1), gROIarr.at(2), gROIarr.at(3),
                                       gROIarr.at(4), gROIarr.at(0)));
    CCM ccm(TEMAT, rois, ref_time.at(0), ref_time.at(1));
    ccm.CalculateEnergyShifts(8);

    over_threshold_values.clear();
    for (std::size_t i = 0; i < ccm.GetNumberOfTimeIndices(); i++)
    {
        auto res = ccm.GetResultContainer(0, i);
        if (res->energy_shift > gSHIFTTHRESHOLD)
        {
            over_threshold_values.push_back(static_cast<float>(res->energy_shift));
        }
    }

    if (gStopExecution.load()) return false;
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
            gr->SetMarkerColor(kRed);
            gr->SetMarkerStyle(20);
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
    std::cout << "A report is made showing whether the detector has shifts above the "
                 "configured threshold.\n";
    std::cout << "You can use --draw option to draw matrices with calculated shifts "
                 "while the code is running.\n\n";

    std::cout << "Usage: detectTimeEvo_ILL [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h                 Show this help message\n";
    std::cout << "  --rootfile <1>             Direct path to the ROOT file including its name\n";
    std::cout << "  --start_run <1>            START run number; used in RunVsEnergy_START_END\n";
    std::cout << "  --end_run <1>              END run number; used in RunVsEnergy_START_END\n";
    std::cout << "  --detector <1>             Detector number; reads RunVsEnergy_detNUMBER\n";
    std::cout << "  --alldet                   Run all detectors from 0 to 40; missing spectra are skipped\n";
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
        else if (arg == "--rootfile" || arg == "--rfile")
        {
            if (i + 1 < argc) { gROOTFILE = argv[++i]; }
            else
            {
                throw std::invalid_argument("Missing value for --rootfile");
            }
        }
        else if (arg == "--start_run" || arg == "--start")
        {
            if (i + 1 < argc) { gSTART_RUN = std::stoi(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --start_run");
            }
        }
        else if (arg == "--end_run" || arg == "--end")
        {
            if (i + 1 < argc) { gEND_RUN = std::stoi(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --end_run");
            }
        }
        else if (arg == "--detector" || arg == "--det")
        {
            if (i + 1 < argc) { gDETECTOR = std::stoi(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --detector");
            }
        }
        else if (arg == "--alldet") { gALLDET = true; }
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
            else
            {
                throw std::invalid_argument("Missing value for --shift_threshold");
            }
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
            else
            {
                throw std::invalid_argument("Missing value for --ROIsource");
            }
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
    if (gROOTFILE.empty()) { throw std::invalid_argument("--rootfile is required"); }
    if (gSTART_RUN == -1) { throw std::invalid_argument("--start_run is required"); }
    if (gEND_RUN == -1) { throw std::invalid_argument("--end_run is required"); }
    if (!gALLDET && gDETECTOR == -1) { throw std::invalid_argument("--detector or --alldet is required"); }

    gROOT_DIR_NAME = get_ill_directory_name();
    if (gALLDET)
    {
        gDETECTORLIST.clear();
        for (int detector = 0; detector <= 40; detector++) { gDETECTORLIST.emplace_back(detector); }
    }
    else
    {
        gDETECTORLIST = {gDETECTOR};
    }

    std::cout << "Parameters used are:" << std::endl;
    std::cout << "Root file:              " << gROOTFILE << std::endl;
    std::cout << "ROOT directory:         " << gROOT_DIR_NAME << std::endl;
    std::cout << "Detector(s):            ";
    for (const auto detector : gDETECTORLIST) { std::cout << detector << " "; }
    std::cout << std::endl;
    std::cout << "Energy shift threshold: " << gSHIFTTHRESHOLD << std::endl;
    std::cout << "ROI:              ";
    for (const auto &roi : gROIarr) { std::cout << roi << " "; }
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
}

int main(int argc, char **argv)
{
    setup_signal_handling();

    TApplication app("app", 0, 0);
    parse_args(argc, argv);

    TFile matfile(gROOTFILE.c_str(), "READ");
    if (matfile.IsZombie())
    {
        throw std::runtime_error("Error! could not open/find the ROOT file " + gROOTFILE);
    }

    std::vector<std::pair<int, std::array<double, 6>>> detected_timeEvo;
    int                                                checked_detectors = 0;

    for (const auto detector : gDETECTORLIST)
    {
        if (gStopExecution.load()) break;
        gSystem->ProcessEvents();

        gDETECTOR                    = detector;
        const std::string matrix_name = get_ill_matrix_name(detector);
        const std::string matrix_path = gROOT_DIR_NAME + "/" + matrix_name;
        TH2              *raw         = static_cast<TH2 *>(matfile.Get(matrix_path.c_str()));
        if (!raw)
        {
            if (gALLDET) { continue; }
            throw std::runtime_error("Matrix not found: " + matrix_path);
        }

        std::shared_ptr<TH2> TEMAT(raw);
        TEMAT->SetDirectory(0);
        checked_detectors++;

        std::vector<float> over_threshold_value;
        const bool         detected = detect_time_evolution(TEMAT, over_threshold_value);

        if (detected)
        {
            const auto stats = calculate_statistics(over_threshold_value);
            std::cout << std::fixed << std::setprecision(2) << "  Detector " << detector << ": "
                      << "above thr: " << static_cast<int>(stats.at(0)) << " mean: " << stats.at(1) << " stdev: " << stats.at(2)
                      << " rms: " << stats.at(3) << " max_abs: " << stats.at(4) << " median: " << stats.at(5) << std::endl;
            detected_timeEvo.emplace_back(detector, stats);
        }
    }

    matfile.Close();

    if (!detected_timeEvo.empty())
    {
        std::cout << "\nTimeEvolution shift above " << gSHIFTTHRESHOLD << " threshold was detected for detector(s):";
        for (const auto &suspect : detected_timeEvo) { std::cout << " " << suspect.first; }
        std::cout << std::endl;
    }
    else if (checked_detectors == 0)
    {
        std::cout << "No detector spectra found in " << gROOT_DIR_NAME << std::endl;
    }
    else
    {
        std::cout << "No TimeEvolution detected above threshold " << gSHIFTTHRESHOLD << std::endl;
    }
    if (guiObjects.size() != 0) { app.Run(); }

    return 0;
}
