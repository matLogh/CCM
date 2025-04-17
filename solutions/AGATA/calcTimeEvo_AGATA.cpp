// #include <array>
#include <chrono> //measure time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Root
#include <TApplication.h>
#include <TBrowser.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TTree.h>

#include "CCM.h"
// #include "CheckCCM.h"
#include "Cross_correlation.h"
#include "RegionOfInterest.h"
#include "TheuerkaufPeak.hpp"
#include "variables.h"
#include <thread>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include <TGraphSmooth.h>
#pragma GCC diagnostic pop

using namespace std::chrono;
using namespace TEC;

#include <functional>
#include <future>
#include <map>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

// Global variables for input parameters
std::string        gCRYSTAL = "";
int                gRUN     = -1;
std::vector<float> gROIarr;
std::vector<float> gREFERENCE_TIME;
std::vector<float> gFIT_PEAK;
int                gREGIME      = 0; // 0 - default grid search, 1 - empirical settings
std::string        gROOTFILE    = "";
std::string        gMATRIX_NAME = "";

const double MINUTES_TO_TIMESTAMPS = 6.0e9;

std::string fourCharInt(int I)
{
    std::stringstream ID;
    ID << std::setfill('0') << std::setw(4) << I;
    return ID.str();
}

int get_crystal_id(const std::string &input)
{
    if (input.size() != 3 || !isdigit(input[0]) || !isdigit(input[1]) ||
        !isalpha(input[2]))
    {
        throw std::invalid_argument(
            "Input must be a 3-character string with 2 digits followed by a letter.");
    }

    int number =
        (input[0] - '0') * 10 +
        (input[1] - '0'); // Combine the first two characters into a single integer
    char letter = std::toupper(input[2]); // Extract the third character as the letter

    int retval = number * 3;
    switch (letter)
    {
    case 'A': retval += 0; break;
    case 'B': retval += 1; break;
    case 'C': retval += 2; break;
    default: throw std::invalid_argument("Invalid letter. Only A, B, or C are allowed.");
    }

    return retval;
}

void write_timeevo_agata_file(CCM              &corrections,
                              const std::string fname             = "TimeEvoCC.conf",
                              int               seconds_per_point = 30)
{
    assert(!fname.empty());
    assert(seconds_per_point > 0);
    std::ofstream file(fname);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << fname << std::endl;
        return;
    }

    std::cout << "Corrections for postPSAfilter are being written to: " << fname
              << std::endl;

    auto         matrix        = corrections.GetInputMatrix();
    const double time_low_edge = matrix->GetXaxis()->GetBinLowEdge(1);
    const double time_up_edge =
        matrix->GetXaxis()->GetBinUpEdge(matrix->GetXaxis()->GetNbins());
    const double step = (double)seconds_per_point / 60.;

    double TS_start, TS_end;
    double gain;
    double time;

    // Write the header
    file << "# TS_start        TS_end          gain \n";

    TS_start = time_low_edge;
    TS_end   = TS_start + step;

    // std::cout << "TS_start: " << TS_start << " TS_end: " << TS_end << std::endl;
    // std::cout << "Time low edge: " << time_low_edge << std::endl;
    // std::cout << "Time up edge: " << time_up_edge << std::endl;

    while (TS_end < time_up_edge)
    {
        time           = (double)TS_start + step / 2.0;
        const auto fit = corrections.GetCorrectionFit(time);

        if (fit.coef.size() != 1)
        {
            file << std::setw(20) << (Long64_t)TS_start * MINUTES_TO_TIMESTAMPS
                 << std::setw(20) << (Long64_t)TS_end * MINUTES_TO_TIMESTAMPS
                 << std::setw(15) << 0. << "\n";
        }
        else
        {
            file << std::setw(20) << (Long64_t)TS_start * MINUTES_TO_TIMESTAMPS
                 << std::setw(15) << (Long64_t)TS_end * MINUTES_TO_TIMESTAMPS
                 << std::setw(20) << fit.coef.front() << "\n";
        }

        TS_start += step;
        TS_end += step;
    }

    file.close();
}

struct ccm_settings
{
    double cost{std::numeric_limits<double>::max()};
    uint   temat_rebin_x{1};
    uint   temat_rebin_y{1};
    // std::vector<RegionOfInterest> ROIs{};
    bool              use_gaussian{true}; // gaussian or polynomial dot product fit
    bool              valid_only{true};
    std::string       interpolator_type{""};
    bool              interpolator_smoothing{false};
    TEC::SmootherType smoother_type{TEC::SmootherType::NONE};
    double            smoother_par{0};

    static void print_header(std::ostream &os = std::cout)
    {
        os << std::setw(10) << "Cost" << std::setw(8) << "RebX" << std::setw(8) << "RebY"
           << std::setw(12) << "Gaussian" << std::setw(10) << "Valid" << std::setw(20)
           << "InterpType" << std::setw(12) << "Smoothing" << std::setw(15)
           << "SmootherType" << std::setw(12) << "SmoothPar" << std::setw(12) << "\n";
    }

    void print_values(std::ostream &os = std::cout) const
    {
        std::string smoother_type_str;
        switch (this->smoother_type)
        {
        case TEC::SmootherType::NONE: smoother_type_str = "NONE"; break;
        case TEC::SmootherType::KERNEL: smoother_type_str = "KERNEL"; break;
        case TEC::SmootherType::LOWESS: smoother_type_str = "LOWESS"; break;
        case TEC::SmootherType::SUPER: smoother_type_str = "SUPER"; break;
        }

        os << std::setw(10) << std::setprecision(6) << this->cost << std::setw(8)
           << this->temat_rebin_x << std::setw(8) << this->temat_rebin_y << std::setw(12)
           << std::boolalpha << this->use_gaussian << std::setw(10) << this->valid_only
           << std::setw(20) << this->interpolator_type << std::setw(12)
           << this->interpolator_smoothing << std::setw(15) << smoother_type_str
           << std::setw(12) << std::setprecision(5) << this->smoother_par << "\n";
    }
};

ccm_settings gSUPER_SETTINGS;

double get_fwfm(TH1 *histo, const double center, const double min, const double max)
{
    TheuerkaufFitter fitter(min, max);
    fitter.AddPeak(center, true, false, false);
    fitter.Fit(histo, "OUTPUT_NONE");
    return fitter.GetPeak(0)->GetFWxM(5);
}

std::string get_pointer_string(void *address)
{
    std::ostringstream oss;
    oss << address;
    return oss.str();
}

double fitfcn(double *x, double *par)
{
    // par[0] = gain
    // par[1] = offset
    return par[0] * x[0];
}

void run_ccm_super_settings(std::shared_ptr<TH2> TEMAT, const ccm_settings &settings)
{
    std::cout << "Running final corrections with super settings..." << std::endl;
    settings.print_header(std::cout);
    settings.print_values(std::cout);

    auto mr = TEMAT->Rebin2D(settings.temat_rebin_x, settings.temat_rebin_y,
                             Form("%%s_rebin_%ix_%iy", TEMAT->GetName(),
                                  settings.temat_rebin_x, settings.temat_rebin_y));

    if (!mr) { throw std::runtime_error("Error: Rebinning TEMAT failed!"); }
    std::shared_ptr<TH2> rTEMAT = std::shared_ptr<TH2>(mr);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(rTEMAT, gROIarr.at(1), gROIarr.at(2),
                                       gROIarr.at(3), gROIarr.at(4), gROIarr.at(0)));
    CCM ccm_fix(rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));

    std::string addressStr = "gain_fcn_" + get_pointer_string(&ccm_fix);
    TF1         fcn("gain_fcn", "[0]*x", 0, 32000);

    ccm_fix.SetCorrectionFunction(fcn, "");
    ccm_fix.CalculateEnergyShifts(4);

    if (settings.use_gaussian) { ccm_fix.UseGaussianResult(); }
    else { ccm_fix.UsePolynomialResult(); }

    if (settings.interpolator_type.empty() && !settings.interpolator_smoothing)
    {
        ccm_fix.DisableInterpolation();
    }
    if (!settings.interpolator_type.empty() && !settings.interpolator_smoothing)
    {
        ccm_fix.EnableInterpolation();
        ccm_fix.ConfigureShiftInterpolator(settings.interpolator_type,
                                           settings.valid_only);
    }

    if (settings.interpolator_smoothing)
    {
        ccm_fix.SmoothShifts(settings.smoother_type, settings.smoother_par);
    }

    write_timeevo_agata_file(ccm_fix, "TimeEvoCC.conf", 30);

    {
        std::string diagnostic_file_name =
            "corrected_timeEvo_r" + std::to_string(gRUN) + "_" + gCRYSTAL + ".root";
        TFile diagnostic_file(diagnostic_file_name.c_str(), "recreate");
        if (!diagnostic_file.IsOpen())
        {
            std::cerr << "Error opening file: " << diagnostic_file_name << std::endl;
            return;
        }
        diagnostic_file.cd();

        auto TEMAT_fixed = ccm_fix.FixMatrix(TEMAT.get());

        std::string proj_name = "projY_" + get_pointer_string(TEMAT_fixed.get());
        TH1        *proj      = TEMAT_fixed->ProjectionY(proj_name.c_str());
        auto        shifts    = ccm_fix.GetROIShifts(0);
        auto profile = ccm_fix.GetInterpolationGraph(0, settings.temat_rebin_x, true);

        TEMAT_fixed->GetYaxis()->SetRangeUser(gROIarr.at(1) + gROIarr.at(3),
                                              gROIarr.at(2) + gROIarr.at(4));
        TEMAT->GetYaxis()->SetRangeUser(gROIarr.at(1) + gROIarr.at(3),
                                        gROIarr.at(2) + gROIarr.at(4));

        TEMAT->Write();
        proj->Write();
        TEMAT_fixed->Write();
        shifts->Write();
        profile->Write();
    }
}

std::vector<ccm_settings> ccm_local_optimizer(const std::shared_ptr<TH2> original_TEMAT,
                                              std::shared_ptr<CCM>       ccm_fix,
                                              ccm_settings               settings,
                                              const std::function<double(TH1 *)> &costFcn)
{
    std::vector<ccm_settings> results;

    // get fcn unique name by setting it to the address of the CCM
    // ojbect
    std::string addressStr = "gain_fcn_" + get_pointer_string(ccm_fix.get());
    TF1         fcn("gain_fcn", "[0]*x", 0, 32000);
    // TF1 fcn(addressStr.c_str(), fitfcn, 0, 32000, 1);

    ccm_fix->SetCorrectionFunction(fcn, "");
    ccm_fix->CalculateEnergyShifts(1);

    // now, loop over all the settings for given CCM
    for (const bool use_gaussian : {true, false})
    {
        settings.use_gaussian = use_gaussian;
        if (use_gaussian) { ccm_fix->UseGaussianResult(); }
        else { ccm_fix->UsePolynomialResult(); }

        {
            settings.interpolator_smoothing = false;
            settings.interpolator_type      = "";
            settings.smoother_type          = TEC::SmootherType::NONE;
            settings.smoother_par           = -1.;
            ccm_fix->DisableInterpolation();
            auto        TEMAT_fixed = ccm_fix->FixMatrix(original_TEMAT.get());
            std::string proj_name   = "projY_" + get_pointer_string(TEMAT_fixed.get());
            TH1        *proj        = TEMAT_fixed->ProjectionY(proj_name.c_str());
            proj->SetDirectory(0);
            double cost   = costFcn(proj);
            settings.cost = cost;
            results.emplace_back(settings);
            delete proj;
        }

        {
            settings.interpolator_smoothing = false;
            settings.interpolator_type      = "akima";
            settings.smoother_type          = TEC::SmootherType::NONE;
            settings.smoother_par           = -1.;

            ccm_fix->EnableInterpolation();
            auto        TEMAT_fixed = ccm_fix->FixMatrix(original_TEMAT.get());
            std::string proj_name   = "projY_" + get_pointer_string(TEMAT_fixed.get());
            TH1        *proj        = TEMAT_fixed->ProjectionY(proj_name.c_str());
            proj->SetDirectory(0);
            double cost   = costFcn(proj);
            settings.cost = cost;
            results.emplace_back(settings);
            delete proj;
        }

        // for (const std::string interpolator :
        //      {"linear", "", "cspline", "cspline_periodic", "akima",
        //      "akima_periodic"})
        // {
        //     settings.interpolator_type = interpolator;
        //     if (interpolator.empty()) {
        //     ccm_fix->DisableInterpolation(); } else
        //     {
        //         ccm_fix->EnableInterpolation();
        //         ccm_fix->ConfigureShiftInterpolator(interpolator,
        //         true);
        //     }

        //     auto        TEMAT_fixed =
        //     ccm_fix->FixMatrix(original_TEMAT.get()); std::string
        //     proj_name   = "projY_" +
        //     get_pointer_string(TEMAT_fixed.get()); TH1        *proj
        //     = TEMAT_fixed->ProjectionY(proj_name.c_str());
        //     proj->SetDirectory(0);
        //     double cost   = costFcn(proj);
        //     settings.cost = cost;
        //     results.emplace_back(settings);
        //     delete proj;
        // }

        for (const auto smoother : {TEC::SmootherType::KERNEL, TEC::SmootherType::LOWESS,
                                    TEC::SmootherType::SUPER})
        {
            std::vector<double> smooth_param;
            if (smoother == TEC::SmootherType::LOWESS)
            {
                smooth_param = {.2, .4, .6, .8, 1.0};
                // smooth_param = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                // 0.8, 0.9, 1.0};
            }
            else { smooth_param = {1, 5, 10, 20, 50, 100, 200}; }

            for (const auto par : smooth_param)
            {
                settings.interpolator_smoothing = true;
                settings.smoother_type          = smoother;
                settings.smoother_par           = par;
                ccm_fix->SmoothShifts(smoother, par);

                auto        TEMAT_fixed = ccm_fix->FixMatrix(original_TEMAT.get());
                std::string proj_name = "projY_" + get_pointer_string(TEMAT_fixed.get());
                TH1        *proj      = TEMAT_fixed->ProjectionY(proj_name.c_str());
                proj->SetDirectory(0);
                double cost   = costFcn(proj);
                settings.cost = cost;
                results.emplace_back(settings);
                delete proj;
            }
        }
    }
    return results;
}

std::vector<ccm_settings> ccm_optimizer_global(
    const std::shared_ptr<TH2>          TEMAT,
    const std::vector<float>            reference_time,
    const std::function<double(TH1 *)> &costFcn)
{
    std::vector<ccm_settings> global_results;

    std::list<std::future<std::vector<ccm_settings>>> futures;

    ccm_settings s;
    s.valid_only = true;
    // rebinX loop
    for (const int rebinX : {1, 2, 5, 10})
    // for (const int rebinX : {1,10})
    {
        s.temat_rebin_x = rebinX;
        // rebinY loop
        for (const int rebinY : {1, 2, 5, 10})
        // for (const int rebinY : {10})
        {
            s.temat_rebin_y = rebinY;

            auto mr = TEMAT->Rebin2D(s.temat_rebin_x, s.temat_rebin_y,
                                     Form("%%s_rebin_%ix_%iy", TEMAT->GetName(),
                                          s.temat_rebin_x, s.temat_rebin_y));

            if (!mr) { throw std::runtime_error("Error: Rebinning TEMAT failed!"); }
            std::shared_ptr<TH2> rTEMAT = std::shared_ptr<TH2>(mr);

            std::vector<RegionOfInterest> ROIs;
            ROIs.emplace_back(RegionOfInterest(rTEMAT, 2200., 2250., -30., 30.,
                                               2223.)); // ROI1 is the region of interest

            std::shared_ptr<CCM> ccm_fix = std::make_shared<CCM>(
                rTEMAT, ROIs, reference_time.at(0), reference_time.at(1));

            auto res = ccm_local_optimizer(rTEMAT, ccm_fix, s, costFcn);

            global_results.insert(global_results.end(), res.begin(), res.end());
        }
    }

    return global_results;
}

// Function to print help message
void print_help()
{
    std::cout << "Usage: program [OPTIONS]\n\n"
              << "Options:\n"
              << "  --crystal [1]              Specify the crystal name "
                 "(e.g. 00A).\n"
              << "  --run [1]                  Specify the run number\n"
              << "  --ROI [1] [2] [3] [4] [5]  Specify the Region of "
                 "Interest (ROI) as:\n"
              << "                                [1] - desired energy of "
                 "the ROI\n"
              << "                                [2] - left edge of ROI n\n"
              << "                                [3] - right edge of ROI\n"
              << "                                [4] - shift ROI by "
                 "maximum of [4] to "
                 "the LEFT (neg value!)\n"
              << "                                [5] - shift ROI by "
                 "maximum of [5] to "
                 "the RIGHT\n"
              << "  --ref_time [1] [2]         Specify the reference time "
                 "interval \n"
              << "  --fit_peak [1] [2] [3]     If running in minimization "
                 "mode, specify "
                 "peak used \n"
              << "                                which FWFM is used to "
                 "find the optimal "
                 "parameters\n"
              << "                                [1] peak center\n"
              << "                                [2] left fit region \n"
              << "                                [3] right fit region \n"
              << "                             Specify the peak used to "
                 "find optimal "
                 "parameters \n"
              << "  --help                     Display this help message "
                 "and exit.\n"
              << "  --rootfile [1]             Specify the root file "
                 "name\n"
              << "  --matrix [1]               Specify the matrix name " << std::endl
              << std::endl;
}

// Helper function to parse space-separated floats
std::vector<float> parse_space_separated_floats(int &i, int argc, char **argv, int count)
{
    std::vector<float> result;
    for (int j = 0; j < count; ++j)
    {
        if (i + 1 < argc)
        {
            try
            {
                result.push_back(std::stof(argv[++i]));
            }
            catch (const std::invalid_argument &)
            {
                throw std::runtime_error("Invalid float value: " + std::string(argv[i]));
            }
        }
        else { throw std::runtime_error("Missing float value for parameter"); }
    }
    return result;
}

// Function to parse command-line arguments
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
        else if (arg == "--fit_peak")
        {
            gFIT_PEAK = parse_space_separated_floats(i, argc, argv, 3);
        }
        else if (arg == "--crystal")
        {
            if (i + 1 < argc) { gCRYSTAL = argv[++i]; }
            else { throw std::runtime_error("Missing value for --crystal"); }
        }
        else if (arg == "--rootfile")
        {
            if (i + 1 < argc) { gROOTFILE = argv[++i]; }
            else { throw std::runtime_error("Missing value for --rootfile"); }
        }
        else if (arg == "--matrix")
        {
            if (i + 1 < argc) { gMATRIX_NAME = argv[++i]; }
            else { throw std::runtime_error("Missing value for --matrix"); }
        }
        else if (arg == "--run")
        {
            if (i + 1 < argc)
            {
                try
                {
                    gRUN = std::stoi(argv[++i]);
                }
                catch (const std::invalid_argument &)
                {
                    throw std::runtime_error("Invalid integer value for --run");
                }
            }
            else { throw std::runtime_error("Missing value for --run"); }
        }
        else if (arg == "--ROI")
        {
            gROIarr =
                parse_space_separated_floats(i, argc, argv, 5); // Expecting 5 floats
        }
        else if (arg == "--ref_time")
        {
            gREFERENCE_TIME =
                parse_space_separated_floats(i, argc, argv, 2); // Expecting 2 floats
        }
        else
        {
            print_help();
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (gCRYSTAL.empty() && gMATRIX_NAME.empty())
    {
        print_help();
        throw std::runtime_error("Either --crystal or --matrix must be provided");
    }

    if (gCRYSTAL.empty() && gROOTFILE.empty())
    {
        print_help();
        throw std::runtime_error("Either --crystal or --rootfile must be provided");
    }

    // Validate required arguments
    // if (gCRYSTAL.empty()) { throw std::runtime_error("--crystal is
    // required"); }
    if (gRUN == -1 && (gMATRIX_NAME.empty() || gROOTFILE.empty()))
    {
        print_help();
        throw std::runtime_error("--run is required when matrix name "
                                 "and rootfile are not defined");
    }
    if (gROIarr.size() != 5)
    {
        print_help();
        throw std::runtime_error("--ROI must have exactly 5 float values");
    }
    if (gREFERENCE_TIME.size() != 2)
    {
        print_help();
        throw std::runtime_error("--ref_time must have exactly 2 float values");
    }

    // Set default root file and matrix name as produced by
    // matTimeEvo_AGATA.cpp
    if (!gCRYSTAL.empty()) { auto crysId = get_crystal_id(gCRYSTAL); }
    if (gMATRIX_NAME.empty()) { gMATRIX_NAME = "hE0_TS_" + gCRYSTAL; }
    if (gROOTFILE.empty())
    {
        gROOTFILE = "Out/run_" + fourCharInt(gRUN) + "/out_" + fourCharInt(gRUN) + "_" +
                    gCRYSTAL + ".root";
    }
}

int main(int argc, char **argv)
{

    gSUPER_SETTINGS.temat_rebin_x          = 1;
    gSUPER_SETTINGS.temat_rebin_y          = 1;
    gSUPER_SETTINGS.use_gaussian           = false;
    gSUPER_SETTINGS.valid_only             = true;
    gSUPER_SETTINGS.interpolator_type      = "akima";
    gSUPER_SETTINGS.interpolator_smoothing = true;
    gSUPER_SETTINGS.smoother_type          = TEC::SmootherType::KERNEL;
    gSUPER_SETTINGS.smoother_par           = 20;
    gSUPER_SETTINGS.cost                   = std::numeric_limits<double>::quiet_NaN();

    // this makes it slower!!!
    // ROOT::EnableImplicitMT();
    // ROOT::EnableThreadSafety();

    parse_args(argc, argv);

    // test();

    // std::string rfname  =
    // "/home/mbalogh/data/ccm_agata/test/out_huge_00A_1010.root";
    TFile *matfile = TFile::Open(gROOTFILE.c_str(), "READ");
    if (!matfile || matfile->IsZombie())
    {
        throw std::runtime_error("Error! could not open/find the " + gROOTFILE + " file");
    }
    std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(gMATRIX_NAME.c_str()));
    // TH2D *TEMAT_original = (TH2D *)matfile->Get(gMATRIX_NAME.c_str());
    if (!TEMAT_original)
    {
        throw std::runtime_error("Error! could not open/find the " + gMATRIX_NAME +
                                 " matrix");
    }

    run_ccm_super_settings(TEMAT_original, gSUPER_SETTINGS);
    return 0;

    std::cout << "Starting jobs..." << std::endl;
    auto start = high_resolution_clock::now();

    auto result = ccm_optimizer_global(
        std::shared_ptr<TH2>(TEMAT_original), gREFERENCE_TIME, [](TH1 *histo) {
            return get_fwfm(histo, gFIT_PEAK.at(0), gFIT_PEAK.at(1), gFIT_PEAK.at(2));
        });

    std::sort(result.begin(), result.end(),
              [](const auto &a, const auto &b) { return a.cost < b.cost; });

    auto stop     = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration: " << duration.count() << " ms" << std::endl;

    ccm_settings::print_header(std::cout);
    for (const auto &r : result) { r.print_values(std::cout); }

    std::string minimization_file =
        "CCMconf_r" + std::to_string(gRUN) + "_" + gCRYSTAL + ".txt";
    std::fstream out_file(minimization_file.c_str(), std::ios::out);
    if (!out_file)
    {
        std::cerr << "Error opening CCM conf file for writing." << std::endl;
    }
    ccm_settings::print_header(out_file);
    for (const auto &r : result) { r.print_values(out_file); }

    // write_timeevo_agata_file(ccm_fix, "TimeEvoCC.conf", 3600);
    return 0;
}

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int test()
{

    std::string rfname  = "/home/mbalogh/data/ccm_agata/test/out_huge_00A_1010.root";
    TFile      *matfile = TFile::Open(rfname.c_str(), "READ");
    if (!matfile || matfile->IsZombie())
    {
        throw std::runtime_error("Error! could not open/find the " + rfname + " file");
    }
    TH2D *TEMAT_original = (TH2D *)matfile->Get("hE0_TS_00A");
    if (!TEMAT_original)
    {
        throw std::runtime_error("Error! could not open/find the hE0_TS_00A matrix");
    }

    double reference_time[]{3103, 3119};

    TF1 fcn("gain_fcn", "[0]*x", 0, 4000);
    // // double reference_time_bgn = 3100;
    // // double reference_time_end = 3200;
    std::shared_ptr<TH2D> TEMAT(dynamic_cast<TH2D *>(
        TEMAT_original->Clone(Form("%s_rebinned", TEMAT_original->GetName()))));

    TEMAT->RebinX(2);
    TEMAT->RebinY(2);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(TEMAT, 2200., 2250., -30., 30., 2223.));

    CCM ccm_fix(TEMAT, ROIs, reference_time[0], reference_time[1]);
    ccm_fix.SetCorrectionFunction(fcn, "");
    ccm_fix.CalculateEnergyShifts(8);
    ccm_fix.UseGaussianResult();

    // ccm_fix.CalculateCorrectionFits();
    // for (int roi_index = 0; roi_index < ROIs.size(); roi_index++)
    // {
    //     // ccm_fix.ConfigureShiftInterpolator(roi_index,
    //     // ROOT::Math::Interpolation::Type::kAKIMA,
    //     //                                    true);
    //     ccm_fix.EnableInterpolation(roi_index);
    // }
    TApplication app("App", 0, 0);
    app.SetReturnFromRun(true);

    ccm_fix.ConfigureShiftInterpolator(ROOT::Math::Interpolation::Type::kAKIMA, true);
    ccm_fix.EnableInterpolation();

    auto  TEMAT_new = ccm_fix.FixMatrix();
    auto *proj      = TEMAT_new->ProjectionY();
    proj->RebinX(8);
    proj->SetName("proj1");
    proj->SetLineColor(kRed);

    auto TEMAT_large_new = ccm_fix.FixMatrix(TEMAT_original);

    // new TCanvas();
    // TEMAT_new->SetTitle("TEMAT_new_corrected but whatever");
    // TEMAT_new->GetYaxis()->SetRangeUser(2200, 2240);
    // TEMAT_new->Draw("COLZ");

    // new TCanvas();
    // TEMAT_large_new->RebinY(2);
    // TEMAT_large_new->RebinX(2);
    // TEMAT_large_new->GetYaxis()->SetRangeUser(2200, 2240);
    // TEMAT_large_new->Draw("COLZ");

    // new TCanvas();
    // TEMAT_original->SetTitle("TEMAT_originak");
    // TEMAT_original->RebinY(2);
    // TEMAT_original->RebinX(2);
    // TEMAT_original->GetYaxis()->SetRangeUser(2200, 2240);
    // TEMAT_original->Draw("COLZ");
    // app.Run();

    // TEMAT_large_new->RebinX(5);
    // TEMAT_large_new->RebinY(2);
    auto proj_large = TEMAT_large_new->ProjectionY();
    proj_large->RebinX(8);
    proj_large->SetName("proj2");
    proj_large->SetLineColor(kGreen);

    new TCanvas();
    proj->Draw("");
    proj_large->Draw("SAME");
    // TEMAT_original->ProjectionY()->Draw("SAME");
    auto *_p = TEMAT->ProjectionY();
    _p->RebinX(8);
    _p->Draw("SAME");
    // new TCanvas();
    // projections.emplace_back(proj);

    double fwfm_proj       = get_fwfm(proj, 7917, 7880, 7940);
    double fwfm_proj_large = get_fwfm(proj_large, 7917, 7880, 7940);

    std::cout << std::setprecision(8) << "proj " << fwfm_proj << " proj_large "
              << fwfm_proj_large << std::endl;

    auto gr_interpol = ccm_fix.GetInterpolationGraph(0);
    gr_interpol->SetName("gr_interpol");
    // gr_interpol->Draw("ALP");
    gr_interpol->SetMarkerStyle(22);
    gr_interpol->SetMarkerColor(2);
    gr_interpol->SetDrawOption("L");
    gr_interpol->SetLineColor(kBlue);
    gr_interpol->SetLineWidth(2);
    gr_interpol->SetFillStyle(0);

    // ccm_fix.SmoothShifts_KernelSmoother(0, 20.);
    ccm_fix.SmoothShifts(TEC::SmootherType::SUPER, 20.);
    auto TEMAT_large_smooth = ccm_fix.FixMatrix(TEMAT_original);
    TEMAT_large_smooth->RebinY(2);
    auto proj_fixed_smooth = TEMAT_large_smooth->ProjectionY();
    proj_fixed_smooth->SetName("proj_fixed_smooth");
    proj_fixed_smooth->SetLineColor(kOrange);
    proj_fixed_smooth->RebinX(8);
    proj_fixed_smooth->Draw("SAME");

    double fwfm_proj_smooth = get_fwfm(proj_fixed_smooth, 7917, 7880, 7940);
    std::cout << "fwfm_proj_smooth " << fwfm_proj_smooth << std::endl;

    ccm_fix.SmoothShifts(TEC::SmootherType::SUPER, 1.);

    auto gr_smooth = ccm_fix.GetInterpolationGraph(0);

    gr_smooth->SetMarkerStyle(23);
    gr_smooth->SetLineColor(kRed);
    gr_smooth->SetLineWidth(3);
    gr_smooth->SetFillStyle(0);
    gr_smooth->SetDrawOption("ACP");

    // for (int i = 20; i < ccm_fix.GetNumberOfTimeIndices() / 2.;
    // i++)
    // {
    //     ccm_fix.SetInvalidResult(0, i);
    // }

    // ccm_fix.SaveToRootFile("elia_roi1.root");

    auto gr_points = ccm_fix.GetROIShifts(0);
    gr_points->SetMarkerStyle(22);
    gr_points->SetLineColor(kBlack);
    gr_points->SetLineWidth(1);
    gr_points->SetFillStyle(0);
    gr_points->SetDrawOption("P");

    TGraphSmooth gs("normal");
    auto        *dio_cane = gs.SmoothKern(gr_points.get(), "normal", 10);
    dio_cane->SetLineColor(kViolet);
    dio_cane->SetLineWidth(2);
    dio_cane->SetDrawOption("C");

    auto mg = new TMultiGraph;

    mg->Add(gr_points.get(), "P");
    mg->Add(gr_smooth.get(), "C");
    mg->Add(gr_interpol.get(), "L");
    mg->Add(dio_cane, "C");
    new TCanvas("Graph_profiles", "Graph_profiles", 800, 600);

    mg->Draw("A");

    write_timeevo_agata_file(ccm_fix, "TimeEvoCC.conf", 3600);

    app.Run();

    return 0;
}