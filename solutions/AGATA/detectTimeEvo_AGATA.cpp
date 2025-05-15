#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH2.h>

#include "CCM.h"

#include "common.cpp"

using namespace TEC;

std::vector<int>         gRUNLIST;
std::vector<std::string> gCRYSTALLIST;
std::string              gDIR            = "timeEvo";
float                    gSHIFTTHRESHOLD = 0.5; // keV
std::vector<float>       gROIarr;
std::string              gMATRIX_NAME;

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

bool detect_time_evolution(const int          run,
                           const std::string &crystal,
                           int               &over_threshold_counter,
                           float             &average_over_threshold_value)
{
    auto  rfname = get_rootfilename(gDIR, run, crystal);
    TFile matfile(rfname.c_str(), "READ");
    if (matfile.IsZombie())
    {
        throw std::runtime_error("Error! could not open/find the ROOT file " + rfname +
                                 " file");
    }
    std::shared_ptr<TH2> TEMAT((TH2 *)matfile.Get(gMATRIX_NAME.c_str()));

    auto ref_time = get_ref_time(TEMAT);
    if (ref_time.at(0) < 0 || ref_time.at(1) < 0)
    {
        std::cerr << "Error! Could not find suitable reference time in the matrix "
                  << rfname << std::endl;
        return false;
    }

    // Create a vector of RegionOfInterest objects
    std::vector<RegionOfInterest> rois;
    rois.emplace_back(RegionOfInterest(TEMAT, gROIarr.at(1), gROIarr.at(2), gROIarr.at(3),
                                       gROIarr.at(4), gROIarr.at(0)));
    // Create a CCM object
    CCM ccm(TEMAT, rois, ref_time.at(0), ref_time.at(1));
    ccm.CalculateEnergyShifts(8);

    over_threshold_counter       = 0;
    average_over_threshold_value = 0.;
    for (int i = 0; i < ccm.GetNumberOfTimeIndices(); i++)
    {
        auto res = ccm.GetResultContainer(0, i);
        if (res->energy_shift > gSHIFTTHRESHOLD)
        {
            over_threshold_counter++;
            average_over_threshold_value += res->energy_shift;
        }
    }
    matfile.Close();

    if (over_threshold_counter > 0) return true;
    return false;
}

void print_help()
{
    std::cout << "Usage: detectTimeEvo_AGATA [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h                 Show this help message\n";
    std::cout << "  --crystal <1> <...>        Specify the crystal(s) name\n";
    std::cout << "  --allcrys                  Run for all crystals of EXP_035\n";
    std::cout << "  --run <1> <...>            Specify the run number\n";
    std::cout << "  --shift_threshold <1>      Energy threshold, if energy shift value "
                 "threshold is found the timeEvo is reported\n";
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
    std::cout << "  --matrix <1>               Specify the matrix name \n";
    std::cout << std::endl << std::endl;
}

void parse_args(int argc, char **argv)
{
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
        else if (arg == "--matrix" || arg == "--mat")
        {
            if (i + 1 < argc) { gMATRIX_NAME = argv[++i]; }
            else { throw std::invalid_argument("Missing value for --matrix"); }
        }
        else
        {
            print_help();
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    for (const auto &run : gRUNLIST)
    {
        for (const auto &crystal : gCRYSTALLIST)
        {
            // Add your processing logic here
        }
    }
}