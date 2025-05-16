// #include <array>
#include <chrono> //measure time
#include <filesystem>
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

#include "common.cpp"

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

// global parameters used in grid search
const std::vector<int>    gRebinX{1, 2, 4};
const std::vector<int>    gRebinY{1, 2, 4, 8};
const std::vector<double> gSmooth_param_lowess{.2, .4, .6, .8, 1.0};
const std::vector<double> gSmooth_param_others{1, 5, 10, 20, 50, 100, 200};

// Global variables for input parameters
std::string        gCRYSTAL = "";
int                gRUN     = -1;
std::vector<float> gROIarr;
std::vector<float> gREFERENCE_TIME;
std::vector<float> gFIT_PEAK;
bool               gUSE_SUPER_SETTINGS{false};
std::string        gROOTFILE      = "";
std::string        gDIR           = "timeEvo";
std::string        gMATRIX_NAME   = "";
int                gREFERENCE_RUN = -1;
std::vector<float> gREFERENCE_VECTOR;
std::vector<int>   gCHAIN_RUNS;

const double MINUTES_TO_TIMESTAMPS = 6.0e9;

void write_timeevo_agata_file(std::shared_ptr<CCM> corrections,
                              const std::string    fname             = "TimeEvoCC.conf",
                              int                  seconds_per_point = 30)
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

    auto         matrix        = corrections->GetInputMatrix();
    const double time_low_edge = matrix->GetXaxis()->GetBinLowEdge(1);
    const double time_up_edge =
        matrix->GetXaxis()->GetBinUpEdge(matrix->GetXaxis()->GetNbins());
    const double step = (double)seconds_per_point / 60.;

    double TS_start, TS_end;
    double gain;
    double time;

    // Write the header
    file << "#" << std::setw(21) << "TS_start" << std::setw(22) << "TS_end"
         << std::setw(22) << "gain"
         << "\n";

    TS_start = time_low_edge;
    TS_end   = TS_start + step;

    // std::cout << "TS_start: " << TS_start << " TS_end: " << TS_end << std::endl;
    // std::cout << "Time low edge: " << time_low_edge << std::endl;
    // std::cout << "Time up edge: " << time_up_edge << std::endl;

    while (TS_end < time_up_edge)
    {
        time           = (double)TS_start + step / 2.0;
        const auto fit = corrections->GetCorrectionFit(time);

        if (fit.coef.size() != 1)
        {
            file << std::fixed << std::setprecision(0) << std::setw(22)
                 << (Long64_t)(TS_start * MINUTES_TO_TIMESTAMPS) << std::setw(22)
                 << (Long64_t)(TS_end * MINUTES_TO_TIMESTAMPS) << std::setw(15) << 0.0
                 << "\n";
        }
        else
        {
            file << std::fixed << std::setprecision(0) << std::setw(22)
                 << (Long64_t)(TS_start * MINUTES_TO_TIMESTAMPS) << std::setw(22)
                 << (Long64_t)(TS_end * MINUTES_TO_TIMESTAMPS) << std::fixed
                 << std::setprecision(6) << std::setw(22) << fit.coef.front() << "\n";
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

void run_ccm_super_settings(std::shared_ptr<TH2> TEMAT,
                            const ccm_settings  &settings,
                            std::string          output_conffilename)
{
    std::cout << "Running final corrections for run " << gRUN << " with super settings..."
              << std::endl;
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

    std::shared_ptr<CCM> ccm_fix = nullptr;

    if (gREFERENCE_VECTOR.size() != 0)
    {
        // we need to "invent" the reference time or it may throw error if the time is
        // outside of this matrix range
        float stupid_ref_start = rTEMAT->GetXaxis()->GetBinLowEdge(1);
        float stupid_ref_end =
            rTEMAT->GetXaxis()->GetBinUpEdge(rTEMAT->GetXaxis()->GetNbins());
        ccm_fix = std::make_shared<CCM>(rTEMAT, ROIs, stupid_ref_start, stupid_ref_end);
        ccm_fix->SetReferenceVector(0, gREFERENCE_VECTOR);
    }
    else
    {
        ccm_fix = std::make_shared<CCM>(rTEMAT, ROIs, gREFERENCE_TIME.at(0),
                                        gREFERENCE_TIME.at(1));
    }

    // CCM ccm_fix(rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));

    std::string addressStr = "gain_fcn_" + get_pointer_string(&ccm_fix);
    TF1         fcn(addressStr.c_str(), "[0]*x", 0, 32000);

    ccm_fix->SetCorrectionFunction(fcn, "");
    ccm_fix->CalculateEnergyShifts(8);

    if (settings.use_gaussian) { ccm_fix->UseGaussianResult(); }
    else { ccm_fix->UsePolynomialResult(); }

    if (settings.interpolator_type.empty() && !settings.interpolator_smoothing)
    {
        ccm_fix->DisableInterpolation();
    }
    if (!settings.interpolator_type.empty() && !settings.interpolator_smoothing)
    {
        ccm_fix->EnableInterpolation();
        ccm_fix->ConfigureShiftInterpolator(settings.interpolator_type,
                                            settings.valid_only);
    }

    if (settings.interpolator_smoothing)
    {
        ccm_fix->SmoothShifts(settings.smoother_type, settings.smoother_par);
    }

    write_timeevo_agata_file(ccm_fix, output_conffilename, 30);

    {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << gRUN;
        std::string diagnostic_file_name = gDIR +
                                           "/diagnostic/"
                                           "correctedTimeEvo_run_" +
                                           oss.str() + "_" + gCRYSTAL + ".root";
        if (can_create_file(diagnostic_file_name))
        {

            TFile diagnostic_file(diagnostic_file_name.c_str(), "recreate");
            if (!diagnostic_file.IsOpen())
            {
                std::cerr << "Error opening file: " << diagnostic_file_name << std::endl;
                return;
            }
            diagnostic_file.cd();

            // settings.print_values(std::cout);
            auto TEMAT_fixed = ccm_fix->FixMatrix(TEMAT.get());

            std::string proj_name = "projY_" + get_pointer_string(TEMAT_fixed.get());
            TH1        *proj      = TEMAT_fixed->ProjectionY(proj_name.c_str());
            auto        shifts    = ccm_fix->GetROIShifts(0);
            auto        profile =
                ccm_fix->GetInterpolationGraph(0, settings.temat_rebin_x, true);

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
            settings.print_values(std::cout);
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
            settings.print_values(std::cout);
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
                smooth_param = gSmooth_param_lowess;
            }
            else { smooth_param = gSmooth_param_others; }

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
                settings.print_values(std::cout);
                results.emplace_back(settings);
                delete proj;
            }
        }
    }
    return results;
}

std::vector<ccm_settings> ccm_optimizer_global(
    const std::shared_ptr<TH2> TEMAT, const std::function<double(TH1 *)> &costFcn)
{
    std::vector<ccm_settings> global_results;

    std::list<std::future<std::vector<ccm_settings>>> futures;

    ccm_settings s;
    s.valid_only = true;
    // rebinX loop
    for (const int rebinX : gRebinX)
    {
        s.temat_rebin_x = rebinX;
        // rebinY loop
        for (const int rebinY : gRebinY)
        {
            s.temat_rebin_y = rebinY;

            auto mr = TEMAT->Rebin2D(s.temat_rebin_x, s.temat_rebin_y,
                                     Form("%%s_rebin_%ix_%iy", TEMAT->GetName(),
                                          s.temat_rebin_x, s.temat_rebin_y));

            if (!mr) { throw std::runtime_error("Error: Rebinning TEMAT failed!"); }
            std::shared_ptr<TH2> rTEMAT = std::shared_ptr<TH2>(mr);

            std::vector<RegionOfInterest> ROIs;
            ROIs.emplace_back(RegionOfInterest(rTEMAT, gROIarr.at(1), gROIarr.at(2),
                                               gROIarr.at(3), gROIarr.at(4),
                                               gROIarr.at(0)));

            std::shared_ptr<CCM> ccm_fix = std::make_shared<CCM>(
                rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));

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
              << "Options:\n";
    std::cout << "  --help                     Display this help message "
                 "and exit.\n";
    std::cout << "  --crystal <1>              Specify the crystal name "
                 "(e.g. 00A).\n";
    std::cout << "  --run <1>                  Specify the run number\n";
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
    std::cout << "  --ref_time <1> <2>         Specify the reference time "
                 "interval \n";
    std::cout << "  --fit_peak <1> <2> <3>     If running in minimization "
                 "mode, specify "
                 "peak used \n"
              << "                                which FWFM is used to "
                 "find the optimal "
                 "parameters\n"
              << "                                <1> peak center\n"
              << "                                <2> left fit region \n"
              << "                                <3> right fit region \n"
              << "                             Specify the peak used to "
                 "find optimal "
                 "parameters \n"
              << "                             Note that this should be different peak "
                 "than one contained in ROI, otherwise you are risking overfitting\n";
    std::cout << "  --dir <1>                  Set directory in which to search for "
                 "matrices and where TimeEvoCC.conf files will be saved\n";
    std::cout << "  --rootfile <1>             Specify the root file "
                 "name\n";
    std::cout << "  --matrix <1>               Specify the matrix name \n";
    std::cout
        << "  --super_settings           Run corrections with hardcoded parameters \n";
    // std::cout << "  --reference_other_run <1>  If you want to use the reference/sample
    // "
    //              "vector from a previous run.\n";
    std::cout
        << "  --chain_runs <1> [...]     Specify the runs that are going to use the same "
           "reference vector as defined for --run\n";
    std::cout << std::endl << std::endl;
}

void set_reference_vector(const int ref_run, const int run)
{
    // std::string reference_root_file = "Out/run_" + fourCharInt(ref_run) + "/out_" +
    //                                   fourCharInt(ref_run) + "_" + gCRYSTAL + ".root";
    std::string reference_root_file = get_rootfilename(gDIR, ref_run, gCRYSTAL);

    TFile *matfile = TFile::Open(reference_root_file.c_str(), "READ");
    if (!matfile || matfile->IsZombie())
    {
        throw std::runtime_error(
            "Error! could not open/find the old/REFERENCE ROOT file " +
            reference_root_file + " file");
    }
    std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(gMATRIX_NAME.c_str()));

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(TEMAT_original, gROIarr.at(1), gROIarr.at(2),
                                       gROIarr.at(3), gROIarr.at(4), gROIarr.at(0)));
    CCM ccm_fix(TEMAT_original, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));
    gREFERENCE_VECTOR = ccm_fix.GetReferenceVector(0);
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
        else if (arg == "--crystal" || arg == "--crys")
        {
            if (i + 1 < argc) { gCRYSTAL = argv[++i]; }
            else { throw std::runtime_error("Missing value for --crystal"); }
        }
        else if (arg == "--chain_runs")
        {
            if (i + 1 < argc)
            {
                gCHAIN_RUNS.clear(); // Clear any previous values
                while (i + 1 < argc && std::isdigit(argv[i + 1][0]))
                {
                    try
                    {
                        gCHAIN_RUNS.push_back(std::stoi(argv[++i]));
                    }
                    catch (const std::invalid_argument &)
                    {
                        throw std::runtime_error(
                            "Invalid integer value for --chain_runs");
                    }
                }
                if (gCHAIN_RUNS.empty())
                {
                    throw std::runtime_error(
                        "--chain_runs must be followed by at least one integer");
                }
            }
            else { throw std::runtime_error("Missing value for --chain_runs"); }
        }

        else if (arg == "--dir")
        {
            if (i + 1 < argc) { gDIR = argv[++i]; }
            else { throw std::runtime_error("Missing value for --dir"); }
            if (gDIR.back() == '/' && gDIR.size() > 1) gDIR.pop_back();
        }
        else if (arg == "--rootfile" || arg == "--rfile")
        {
            if (i + 1 < argc) { gROOTFILE = argv[++i]; }
            else { throw std::runtime_error("Missing value for --rootfile"); }
        }
        else if (arg == "--matrix" || arg == "--mat")
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
        else if (arg == "--ROIsource")
        {
            std::vector<float> peak;
            if (i + 1 < argc)
            {
                parse_ROI_source(argv[++i], gROIarr, peak);
                if (gFIT_PEAK.empty()) { gFIT_PEAK = peak; }
                else
                {
                    std::cerr << "Warning: --fit_peak is overwriting peak defined by "
                                 "--ROIsource\n";
                }
            }
            else { throw std::runtime_error("Missing value for --ROIsource"); }
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
        else if (arg == "--super_settings") { gUSE_SUPER_SETTINGS = true; }
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
        throw std::runtime_error("--ROI must have exactly 5 float values, but it has " +
                                 std::to_string(gROIarr.size()));
    }
    if (gREFERENCE_TIME.size() != 2)
    {
        print_help();
        throw std::runtime_error("--ref_time must have exactly 2 float values");
    }

    if (gFIT_PEAK.size() != 3 && !gUSE_SUPER_SETTINGS)
    {
        print_help();
        throw std::runtime_error(
            "if running without --super_settings flag, grid search for optimal "
            "parameters is used and thus --fit_peak must be defined");
    }

    // Set default root file and matrix name as produced by
    // matTimeEvo_AGATA.cpp
    if (!gCRYSTAL.empty()) { auto crysId = get_crystal_id(gCRYSTAL); }
    if (gMATRIX_NAME.empty()) { gMATRIX_NAME = "hE0_TS_" + gCRYSTAL; }

    // set
    if (gROOTFILE.empty()) { gROOTFILE = get_rootfilename(gDIR, gRUN, gCRYSTAL); }

    // Print all parsed input parameters
    std::cout << "Parsed Input Parameters:" << std::endl;
    std::cout << "  Crystal: " << gCRYSTAL << std::endl;
    std::cout << "  Run: " << gRUN << std::endl;
    std::cout << "  Root File: " << gROOTFILE << std::endl;
    std::cout << "  Matrix Name: " << gMATRIX_NAME << std::endl;
    std::cout << "  ROI: ";
    for (const auto &val : gROIarr) { std::cout << val << " "; }
    std::cout << std::endl;
    std::cout << "  Reference Run: " << gREFERENCE_RUN << std::endl;
    std::cout << "  Reference Time: ";
    for (const auto &val : gREFERENCE_TIME) { std::cout << val << " "; }
    std::cout << std::endl;
    std::cout << "  Fit Peak: ";
    for (const auto &val : gFIT_PEAK) { std::cout << val << " "; }
    std::cout << std::endl;
    std::cout << "  Use Super Settings: " << std::boolalpha << gUSE_SUPER_SETTINGS
              << std::endl;
    std::cout << "Output configuration file will be saved to: "
              << get_conffilename(gDIR, gRUN, gCRYSTAL) << std::endl;
}

void run_chained_runs(const ccm_settings &optimal_settings)
{
    if (gCHAIN_RUNS.empty()) { return; }

    set_reference_vector(gRUN, gCHAIN_RUNS.at(0));
    for (const auto c_run : gCHAIN_RUNS)
    {

        // std::string c_rootfile = "run_" + fourCharInt(c_run) + "/Out/TimeEvo/out_" +
        //                          fourCharInt(c_run) + "_" + gCRYSTAL + ".root";

        std::string c_rootfile = get_rootfilename(gDIR, c_run, gCRYSTAL);

        TFile *matfile = TFile::Open(c_rootfile.c_str(), "READ");
        if (!matfile || matfile->IsZombie())
        {
            throw std::runtime_error("Error! could not open/find the " + c_rootfile +
                                     " file");
        }
        // matrix name is the same, we can use gMATRIX_NAME
        std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(gMATRIX_NAME.c_str()));
        if (!TEMAT_original)
        {
            throw std::runtime_error("Error! could not open/find the " + gMATRIX_NAME +
                                     " matrix");
        }

        // set conf path for output file
        std::string conf_filename = get_conffilename(gDIR, c_run, gCRYSTAL);
        gRUN                      = c_run;
        std::cout << std::endl << "Running chained run: " << c_run << std::endl;
        run_ccm_super_settings(TEMAT_original, optimal_settings, conf_filename);
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    // this makes it slower!!!
    // ROOT::EnableImplicitMT();
    // ROOT::EnableThreadSafety();

    TFile *matfile = TFile::Open(gROOTFILE.c_str(), "READ");
    if (!matfile || matfile->IsZombie())
    {
        throw std::runtime_error("Error! could not open/find the " + gROOTFILE + " file");
    }
    std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(gMATRIX_NAME.c_str()));
    if (!TEMAT_original)
    {
        throw std::runtime_error("Error! could not open/find the " + gMATRIX_NAME +
                                 " matrix");
    }

    if (gUSE_SUPER_SETTINGS)
    {
        gSUPER_SETTINGS.temat_rebin_x = 2;
        gSUPER_SETTINGS.temat_rebin_y = 2;
        gSUPER_SETTINGS.use_gaussian  = true;
        gSUPER_SETTINGS.valid_only    = true;
        // gSUPER_SETTINGS.interpolator_type      = "akima";
        gSUPER_SETTINGS.interpolator_type      = "";
        gSUPER_SETTINGS.interpolator_smoothing = false;
        gSUPER_SETTINGS.smoother_type          = TEC::SmootherType::NONE;
        gSUPER_SETTINGS.smoother_par           = -1;
        gSUPER_SETTINGS.cost                   = std::numeric_limits<double>::quiet_NaN();

        run_ccm_super_settings(TEMAT_original, gSUPER_SETTINGS,
                               get_conffilename(gDIR, gRUN, gCRYSTAL));
        if (gREFERENCE_RUN != -1) { run_chained_runs(gSUPER_SETTINGS); }
        return 0;
    }

    std::cout << "\nTesting following settings: " << std::endl;
    ccm_settings::print_header(std::cout);

    auto start = high_resolution_clock::now();

    auto result =
        ccm_optimizer_global(std::shared_ptr<TH2>(TEMAT_original), [&](TH1 *histo) {
            return get_fwfm(histo, gFIT_PEAK.at(0), gFIT_PEAK.at(1), gFIT_PEAK.at(2));
        });

    std::sort(result.begin(), result.end(),
              [](const auto &a, const auto &b) { return a.cost < b.cost; });

    auto stop     = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration: " << duration.count() << " ms" << std::endl;

    ccm_settings::print_header(std::cout);
    for (const auto &r : result) { r.print_values(std::cout); }

    std::string minimization_file = gDIR + "/diagnostics/" + "CCMconf_r" +
                                    std::to_string(gRUN) + "_" + gCRYSTAL + ".txt";
    std::fstream out_file(minimization_file.c_str(), std::ios::out);
    if (!out_file)
    {
        std::cerr << "Error opening CCM conf file for writing." << std::endl;
    }
    ccm_settings::print_header(out_file);
    for (const auto &r : result) { r.print_values(out_file); }

    gSUPER_SETTINGS = result.front();
    run_ccm_super_settings(TEMAT_original, gSUPER_SETTINGS,
                           get_conffilename(gDIR, gRUN, gCRYSTAL));
    run_chained_runs(gSUPER_SETTINGS);

    return 0;
}
