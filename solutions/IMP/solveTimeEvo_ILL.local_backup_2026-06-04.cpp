#include <algorithm>
// #include <array>
#include <chrono> //measure time
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
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
#include "TheuerkaufPeak.h"
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
const std::vector<uint> gRebinX{1};
const std::vector<uint> gRebinY{1};
// const std::vector<double> gSmooth_param_lowess{.2, .4, .6, .8, 1.0};
// const std::vector<double> gSmooth_param_others{1, 2, 5, 10, 20, 50, 100, 200};
const std::vector<double> gSmooth_param_others{1, 2, 5};

// Global variables for input parameters
std::vector<std::vector<float>> gROIarrs;
std::vector<float> gREFERENCE_TIME;
std::vector<float> gFIT_PEAK;
bool               gUSE_SUPER_SETTINGS{false};
bool               gUSE_GAIN_WITH_OFFSET{false};
std::string        gROOTFILE      = "";
int                gSTART_RUN     = -1;
int                gEND_RUN       = -1;
int                gDETECTOR      = -1;
std::string        gROOT_DIR_NAME = "";
std::string        gMATRIX_NAME   = "";
std::vector<std::vector<float>> gREFERENCE_VECTORS;
std::vector<std::pair<int, int>> gCHAIN_RUN_RANGES;

std::string get_ill_directory_name() { return "RunVsEnergy_" + std::to_string(gSTART_RUN) + "_" + std::to_string(gEND_RUN); }

std::string get_ill_matrix_name() { return "RunVsEnergy_det" + std::to_string(gDETECTOR); }

std::filesystem::path get_ill_output_base()
{
    const std::filesystem::path root_path(gROOTFILE);
    auto                        parent = root_path.parent_path();
    if (parent.empty()) { parent = "."; }

    return parent / (get_ill_directory_name() + "_det" + std::to_string(gDETECTOR));
}

std::string get_output_conffilename() { return (get_ill_output_base().string() + "_TimeEvoCC.conf"); }

std::string get_output_diagnostic_filename() { return (get_ill_output_base().string() + "_correctedTimeEvo.root"); }

std::string get_output_minimization_filename() { return (get_ill_output_base().string() + "_CCMconf.txt"); }

std::string get_correction_formula()
{
    return gUSE_GAIN_WITH_OFFSET ? "[0] + [1]*x" : "[0]*x";
}

std::shared_ptr<TH2> get_ill_matrix_from_file(TFile *matfile)
{
    if (!matfile || matfile->IsZombie()) { throw std::runtime_error("Error! could not open/find the " + gROOTFILE + " file"); }

    gROOT_DIR_NAME              = get_ill_directory_name();
    gMATRIX_NAME                = get_ill_matrix_name();
    const std::string matrix_path = gROOT_DIR_NAME + "/" + gMATRIX_NAME;
    std::shared_ptr<TH2> TEMAT((TH2 *)matfile->Get(matrix_path.c_str()));
    if (!TEMAT) { throw std::runtime_error("Error! could not open/find the " + matrix_path + " matrix"); }

    return TEMAT;
}

std::vector<RegionOfInterest> make_rois(const std::shared_ptr<TH2> &matrix)
{
    std::vector<RegionOfInterest> ROIs;
    ROIs.reserve(gROIarrs.size());
    for (const auto &roi : gROIarrs)
    {
        ROIs.emplace_back(RegionOfInterest(matrix, roi.at(1), roi.at(2), roi.at(3), roi.at(4), roi.at(0)));
    }
    return ROIs;
}

std::pair<double, double> get_diagnostic_energy_range()
{
    if (gROIarrs.empty()) { throw std::runtime_error("No ROI configured"); }

    double low  = gROIarrs.front().at(1) + gROIarrs.front().at(3);
    double high = gROIarrs.front().at(2) + gROIarrs.front().at(4);
    for (const auto &roi : gROIarrs)
    {
        low  = std::min(low, static_cast<double>(roi.at(1) + roi.at(3)));
        high = std::max(high, static_cast<double>(roi.at(2) + roi.at(4)));
    }
    return {low, high};
}

void adjust_peak_energy(std::shared_ptr<TH2> TEMAT, std::vector<float> &peak_array)

{
    TEMAT->GetXaxis()->SetRange(1, TEMAT->GetXaxis()->GetNbins());
    TEMAT->GetYaxis()->SetRange(1, TEMAT->GetYaxis()->GetNbins());
    if (peak_array.size() < 3)
    {
        std::cerr << "Error: Peak array must contain at least 3 elements (position, min, max)." << std::endl;
        return;
    }
    std::unique_ptr<TH1> proj(TEMAT->ProjectionY(Form("adjust_roi_energy_%f", static_cast<double>(peak_array.at(0))),
                                                 TEMAT->GetXaxis()->FindBin(gREFERENCE_TIME.at(0)),
                                                 TEMAT->GetXaxis()->FindBin(gREFERENCE_TIME.at(1))));
    proj->SetDirectory(0);
    proj->GetXaxis()->SetRangeUser(peak_array.at(1), peak_array.at(2));
    auto mean = proj->GetMean();

    // std::cout << peak_array.at(1) << " " << peak_array.at(2) << " " << mean <<
    // std::endl;

    TheuerkaufFitter fitter(peak_array.at(1), peak_array.at(2));
    fitter.AddPeak(mean, true, false, false);

    fitter.Fit(proj.get(), "OUTPUT_NONE");
    peak_array.at(0) = static_cast<float>(fitter.GetPeak(0)->GetPos());
}

void write_timeevo_conf_format(std::shared_ptr<CCM> corrections, std::string fname = "TimeEvoCC.conf")
{
    assert(!fname.empty());
    std::ofstream file(fname);
    if (!file.is_open())
    {
        std::cerr << "Error opening output file: " << fname << std::endl;
        return;
    }

    std::cout << "Corrections are being written to: " << fname << std::endl;

    auto   matrix = corrections->GetInputMatrix();
    double run;

    // Write the header
    if (gUSE_GAIN_WITH_OFFSET)
    {
        file << "#" << std::setw(15) << "run" << std::setw(22) << "offset" << std::setw(22) << "gain" << "\n";
    }
    else
    {
        file << "#" << std::setw(15) << "run" << std::setw(22) << "gain" << "\n";
    }

    for (int bin = 1; bin < matrix->GetXaxis()->GetNbins(); bin++)
    {
        run = matrix->GetXaxis()->GetBinCenter(bin);

        const auto fit = corrections->GetCorrectionFit(run);

        file << std::fixed << std::setprecision(0) << std::setw(16) << static_cast<Long64_t>(std::llround(run));
        if (gUSE_GAIN_WITH_OFFSET)
        {
            const auto offset = fit.coef.size() == 2 ? fit.coef.at(0) : 0.0;
            const auto gain   = fit.coef.size() == 2 ? fit.coef.at(1) : (fit.coef.size() == 1 ? fit.coef.front() : 0.0);
            file << std::fixed << std::setprecision(10) << std::setw(22) << offset << std::setw(22) << gain << "\n";
        }
        else
        {
            const auto gain = fit.coef.size() == 1 ? fit.coef.front() : 0.0;
            file << std::fixed << std::setprecision(10) << std::setw(22) << gain << "\n";
        }
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
        os << std::setw(10) << "Cost" << std::setw(8) << "RebX" << std::setw(8) << "RebY" << std::setw(12) << "Gaussian" << std::setw(10) << "Valid"
           << std::setw(20) << "InterpType" << std::setw(12) << "Smoothing" << std::setw(15) << "SmootherType" << std::setw(12) << "SmoothPar"
           << std::setw(12) << std::endl;
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

        os << std::setw(10) << std::setprecision(6) << this->cost << std::setw(8) << this->temat_rebin_x << std::setw(8) << this->temat_rebin_y
           << std::setw(12) << std::boolalpha << this->use_gaussian << std::setw(10) << this->valid_only << std::setw(20) << this->interpolator_type
           << std::setw(12) << this->interpolator_smoothing << std::setw(15) << smoother_type_str << std::setw(12) << std::setprecision(5)
           << this->smoother_par << std::endl;
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

void run_ccm_super_settings(std::shared_ptr<TH2> TEMAT, const ccm_settings &settings, std::string output_conffilename)
{
    std::cout << "Running final corrections for runs " << gSTART_RUN << "-" << gEND_RUN << ", detector " << gDETECTOR << " with super settings..."
              << std::endl;
    settings.print_header(std::cout);
    settings.print_values(std::cout);

    auto mr = TEMAT->Rebin2D(static_cast<int>(settings.temat_rebin_x), static_cast<int>(settings.temat_rebin_y),
                             Form("%s_rebin_%ux_%uy", TEMAT->GetName(), settings.temat_rebin_x, settings.temat_rebin_y));

    if (!mr) { throw std::runtime_error("Error: Rebinning TEMAT failed!"); }
    std::shared_ptr<TH2> rTEMAT = std::shared_ptr<TH2>(mr);

    auto ROIs = make_rois(rTEMAT);

    std::shared_ptr<CCM> ccm_fix = nullptr;

    if (!gREFERENCE_VECTORS.empty())
    {
        // we need to "invent" the reference time or it may throw error if the time is
        // outside of this matrix range
        float stupid_ref_start = static_cast<float>(rTEMAT->GetXaxis()->GetBinLowEdge(1));
        float stupid_ref_end   = static_cast<float>(rTEMAT->GetXaxis()->GetBinUpEdge(rTEMAT->GetXaxis()->GetNbins()));
        if (ROIs.size() != gREFERENCE_VECTORS.size())
        {
            throw std::runtime_error("Number of shared reference vectors does not match number of ROIs in solveTimeEvo_ILL");
        }
        ccm_fix                = std::make_shared<CCM>(rTEMAT, ROIs, stupid_ref_start, stupid_ref_end);
        for (size_t roi_index = 0; roi_index < gREFERENCE_VECTORS.size(); ++roi_index)
        {
            ccm_fix->SetReferenceVector(static_cast<unsigned int>(roi_index), gREFERENCE_VECTORS.at(roi_index));
        }
    }
    else
    {
        ccm_fix = std::make_shared<CCM>(rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));
    }

    // CCM ccm_fix(rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));

    std::string addressStr = "gain_fcn_" + get_pointer_string(ccm_fix.get());
    TF1         fcn(addressStr.c_str(), get_correction_formula().c_str(), 0, 32000);

    ccm_fix->SetCorrectionFunction(fcn, "");
    ccm_fix->CalculateEnergyShifts(8);

    if (settings.use_gaussian) { ccm_fix->UseGaussianResult(); }
    else
    {
        ccm_fix->UsePolynomialResult();
    }

    if (settings.interpolator_type.empty() && !settings.interpolator_smoothing) { ccm_fix->DisableInterpolation(); }
    if (!settings.interpolator_type.empty() && !settings.interpolator_smoothing)
    {
        ccm_fix->EnableInterpolation();
        ccm_fix->ConfigureShiftInterpolator(settings.interpolator_type, settings.valid_only);
    }

    if (settings.interpolator_smoothing) { ccm_fix->SmoothShifts(settings.smoother_type, settings.smoother_par); }

    write_timeevo_conf_format(ccm_fix, output_conffilename);

    {
        std::string diagnostic_file_name = get_output_diagnostic_filename();
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
            const auto [roi_range_low, roi_range_high] = get_diagnostic_energy_range();

            TEMAT_fixed->GetYaxis()->SetRangeUser(roi_range_low, roi_range_high);
            TEMAT->GetYaxis()->SetRangeUser(roi_range_low, roi_range_high);

            TEMAT->Write();
            proj->Write();
            TEMAT_fixed->Write();
            for (size_t roi_index = 0; roi_index < ccm_fix->GetNumberOfROIs(); ++roi_index)
            {
                auto shifts  = ccm_fix->GetROIShifts(roi_index);
                auto profile = ccm_fix->GetInterpolationGraph(roi_index, static_cast<int>(settings.temat_rebin_x), true);
                shifts->Write();
                profile->Write();
            }
        }
    }
}

void set_reference_vector()
{
    TFile *matfile = TFile::Open(gROOTFILE.c_str(), "READ");
    auto   TEMAT   = get_ill_matrix_from_file(matfile);
    auto   ROIs    = make_rois(TEMAT);
    CCM    ccm_fix(TEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));
    gREFERENCE_VECTORS.clear();
    gREFERENCE_VECTORS.reserve(ROIs.size());
    for (size_t roi_index = 0; roi_index < ROIs.size(); ++roi_index)
    {
        gREFERENCE_VECTORS.emplace_back(ccm_fix.GetReferenceVector(roi_index));
    }
}

void run_chained_ranges(const ccm_settings &optimal_settings)
{
    if (gCHAIN_RUN_RANGES.empty()) { return; }

    const int reference_start = gSTART_RUN;
    const int reference_end   = gEND_RUN;
    set_reference_vector();

    for (const auto &[chain_start, chain_end] : gCHAIN_RUN_RANGES)
    {
        gSTART_RUN = chain_start;
        gEND_RUN   = chain_end;

        TFile *matfile = TFile::Open(gROOTFILE.c_str(), "READ");
        auto   TEMAT   = get_ill_matrix_from_file(matfile);

        const auto conf_filename = get_output_conffilename();
        if (!can_create_file(conf_filename)) throw std::runtime_error("Problem with creating output configuration file\n");

        std::cout << std::endl
                  << "Running chained ILL range: " << gSTART_RUN << "-" << gEND_RUN
                  << " using reference range " << reference_start << "-" << reference_end << std::endl;
        run_ccm_super_settings(TEMAT, optimal_settings, conf_filename);
    }
}

std::vector<ccm_settings> ccm_local_optimizer(const std::shared_ptr<TH2>          original_TEMAT,
                                              std::shared_ptr<CCM>                ccm_fix,
                                              ccm_settings                        settings,
                                              const std::function<double(TH1 *)> &costFcn)
{
    std::vector<ccm_settings> results;
    // get fcn unique name by setting it to the address of the CCM
    // ojbect
    std::string addressStr = "gain_fcn_" + get_pointer_string(ccm_fix.get());
    TF1         fcn("gain_fcn", get_correction_formula().c_str(), 0, 32000);

    ccm_fix->SetCorrectionFunction(fcn, "");
    ccm_fix->CalculateEnergyShifts(1);

    // now, loop over all the settings for given CCM
    for (const bool use_gaussian : {true, false})
    {
        settings.use_gaussian = use_gaussian;
        if (use_gaussian) { ccm_fix->UseGaussianResult(); }
        else
        {
            ccm_fix->UsePolynomialResult();
        }

        // {
        //     settings.interpolator_smoothing = false;
        //     settings.interpolator_type      = "";
        //     settings.smoother_type          = TEC::SmootherType::NONE;
        //     settings.smoother_par           = -1.;
        //     ccm_fix->DisableInterpolation();

        //     auto        TEMAT_fixed = ccm_fix->FixMatrix(original_TEMAT.get());
        //     std::string proj_name   = "projY_" +
        //     get_pointer_string(TEMAT_fixed.get()); TH1        *proj        =
        //     TEMAT_fixed->ProjectionY(proj_name.c_str()); proj->SetDirectory(0);
        //     double cost   = costFcn(proj);
        //     settings.cost = cost;
        //     settings.print_values(std::cout);
        //     results.emplace_back(settings);
        //     delete proj;
        // }

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

        // for (const auto smoother : {TEC::SmootherType::KERNEL,
        // TEC::SmootherType::LOWESS,
        //                             TEC::SmootherType::SUPER})
        for (const auto smoother : {TEC::SmootherType::SUPER})
        {
            std::vector<double> smooth_param;
            if (smoother == TEC::SmootherType::LOWESS)
            {
                // smooth_param = gSmooth_param_lowess;
            }
            else
            {
                smooth_param = gSmooth_param_others;
            }

            for (const auto par : smooth_param)
            {
                settings.interpolator_smoothing = true;
                settings.smoother_type          = smoother;
                settings.smoother_par           = par;
                ccm_fix->SmoothShifts(smoother, par);

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
        }
    }
    return results;
}

std::vector<ccm_settings> ccm_optimizer_global(const std::shared_ptr<TH2> TEMAT, const std::function<double(TH1 *)> &costFcn)
{
    std::vector<ccm_settings> global_results;

    std::list<std::future<std::vector<ccm_settings>>> futures;

    ccm_settings s;
    s.valid_only = true;
    // rebinX loop
    for (const auto rebinX : gRebinX)
    {
        s.temat_rebin_x = rebinX;
        // rebinY loop
        for (const auto rebinY : gRebinY)
        {
            s.temat_rebin_y = rebinY;
            std::shared_ptr<TH2> rTEMAT;
            if (rebinX == 1 && rebinY == 1) { rTEMAT = TEMAT; }
            else
            {
                rTEMAT = std::shared_ptr<TH2>(TEMAT->Rebin2D(static_cast<int>(s.temat_rebin_x), static_cast<int>(s.temat_rebin_y),
                                                             Form("%s_rebin_%ux_%uy", TEMAT->GetName(), s.temat_rebin_x, s.temat_rebin_y)));
            }

            if (rTEMAT.get() == nullptr) { throw std::runtime_error("Error: Rebinning TEMAT failed!"); }

            auto ROIs = make_rois(rTEMAT);

            std::shared_ptr<CCM> ccm_fix = std::make_shared<CCM>(rTEMAT, ROIs, gREFERENCE_TIME.at(0), gREFERENCE_TIME.at(1));

            auto res = ccm_local_optimizer(rTEMAT, ccm_fix, s, costFcn);

            global_results.insert(global_results.end(), res.begin(), res.end());
        }
    }

    return global_results;
}

// Function to print help message
void print_help()
{
    std::cout << "Usage: solveTimeEvo_ILL [OPTIONS]\n\n"
              << "Options:\n";
    std::cout << "  --help                     Display this help message "
                 "and exit.\n";
    std::cout << "  --rootfile <1>             Direct path to the ROOT file including its name\n";
    std::cout << "  --start_run <1>            START run number; used in RunVsEnergy_START_END\n";
    std::cout << "  --end_run <1>              END run number; used in RunVsEnergy_START_END\n";
    std::cout << "  --detector <1>             Detector number; reads RunVsEnergy_detNUMBER\n";
    std::cout << "  --ROI <1> <2> <3> <4> <5>  Specify a Region of "
                 "Interest (ROI). Can be used multiple times:\n"
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

    std::cout << "  --ROIsource <1>            Add ROI for calibration sources. "
                 "Currently "
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
    std::cout << "  --super_settings           Run corrections with hardcoded "
                 "parameters \n";
    std::cout << "  --gain_with_offset         Use [0] + [1]*x correction instead of [0]*x. "
                 "Requires at least two ROIs.\n";
    std::cout << "  --chain_ranges <s> <e> [...] Run additional ILL run ranges using the "
                 "reference vector from --start_run/--end_run. Values must be "
                 "start/end pairs.\n";
    std::cout << std::endl << std::endl;
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
        else if (arg == "--fit_peak") { gFIT_PEAK = parse_space_separated_floats(i, argc, argv, 3); }
        else if (arg == "--rootfile" || arg == "--rfile")
        {
            if (i + 1 < argc) { gROOTFILE = argv[++i]; }
            else
            {
                throw std::runtime_error("Missing value for --rootfile");
            }
        }
        else if (arg == "--start_run" || arg == "--start")
        {
            if (i + 1 < argc)
            {
                try
                {
                    gSTART_RUN = std::stoi(argv[++i]);
                }
                catch (const std::invalid_argument &)
                {
                    throw std::runtime_error("Invalid integer value for --start_run");
                }
            }
            else
            {
                throw std::runtime_error("Missing value for --start_run");
            }
        }
        else if (arg == "--end_run" || arg == "--end")
        {
            if (i + 1 < argc)
            {
                try
                {
                    gEND_RUN = std::stoi(argv[++i]);
                }
                catch (const std::invalid_argument &)
                {
                    throw std::runtime_error("Invalid integer value for --end_run");
                }
            }
            else
            {
                throw std::runtime_error("Missing value for --end_run");
            }
        }
        else if (arg == "--detector" || arg == "--det")
        {
            if (i + 1 < argc)
            {
                try
                {
                    gDETECTOR = std::stoi(argv[++i]);
                }
                catch (const std::invalid_argument &)
                {
                    throw std::runtime_error("Invalid integer value for --detector");
                }
            }
            else
            {
                throw std::runtime_error("Missing value for --detector");
            }
        }
        else if (arg == "--ROIsource")
        {
            std::vector<float> roi;
            std::vector<float> peak;
            if (i + 1 < argc)
            {
                parse_ROI_source(argv[++i], roi, peak);
                gROIarrs.emplace_back(std::move(roi));
                if (gFIT_PEAK.empty()) { gFIT_PEAK = peak; }
                else
                {
                    std::cerr << "Warning: --fit_peak is overwriting peak defined by "
                                 "--ROIsource\n";
                }
            }
            else
            {
                throw std::runtime_error("Missing value for --ROIsource");
            }
        }
        else if (arg == "--ROI")
        {
            gROIarrs.emplace_back(parse_space_separated_floats(i, argc, argv, 5)); // Expecting 5 floats
        }
        else if (arg == "--ref_time")
        {
            gREFERENCE_TIME = parse_space_separated_floats(i, argc, argv, 2); // Expecting 2 floats
        }
        else if (arg == "--chain_ranges")
        {
            std::vector<int> chain_values;
            while (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
            {
                try
                {
                    chain_values.emplace_back(std::stoi(argv[++i]));
                }
                catch (const std::invalid_argument &)
                {
                    throw std::runtime_error("Invalid integer value for --chain_ranges");
                }
            }
            if (chain_values.empty() || chain_values.size() % 2 != 0)
            {
                throw std::runtime_error("--chain_ranges must be followed by one or more start/end integer pairs");
            }
            for (size_t range_index = 0; range_index < chain_values.size(); range_index += 2)
            {
                gCHAIN_RUN_RANGES.emplace_back(chain_values.at(range_index), chain_values.at(range_index + 1));
            }
        }
        else if (arg == "--gain_with_offset") { gUSE_GAIN_WITH_OFFSET = true; }
        else if (arg == "--super_settings") { gUSE_SUPER_SETTINGS = true; }
        else
        {
            print_help();
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (gROOTFILE.empty())
    {
        print_help();
        throw std::runtime_error("--rootfile is required");
    }

    if (gSTART_RUN == -1)
    {
        print_help();
        throw std::runtime_error("--start_run is required");
    }

    if (gEND_RUN == -1)
    {
        print_help();
        throw std::runtime_error("--end_run is required");
    }

    if (gDETECTOR == -1)
    {
        print_help();
        throw std::runtime_error("--detector is required");
    }
    if (gROIarrs.empty())
    {
        print_help();
        throw std::runtime_error("At least one --ROI or --ROIsource must be defined");
    }
    for (const auto &roi : gROIarrs)
    {
        if (roi.size() != 5)
        {
            print_help();
            throw std::runtime_error("--ROI must have exactly 5 float values, but one ROI has " + std::to_string(roi.size()));
        }
    }
    if (gREFERENCE_TIME.size() != 2)
    {
        print_help();
        throw std::runtime_error("--ref_time must have exactly 2 float values");
    }
    if (gUSE_GAIN_WITH_OFFSET && gROIarrs.size() < 2)
    {
        print_help();
        throw std::runtime_error("--gain_with_offset requires at least two ROIs");
    }

    if (gFIT_PEAK.size() != 3 && !gUSE_SUPER_SETTINGS)
    {
        print_help();
        throw std::runtime_error("if running without --super_settings flag, grid search for optimal "
                                 "parameters is used and thus --fit_peak must be defined");
    }

    gROOT_DIR_NAME = get_ill_directory_name();
    gMATRIX_NAME   = get_ill_matrix_name();

    // Print all parsed input parameters
    std::cout << "Parsed Input Parameters:" << std::endl;
    std::cout << "  Root File: " << gROOTFILE << std::endl;
    std::cout << "  ROOT Directory: " << gROOT_DIR_NAME << std::endl;
    std::cout << "  Start Run: " << gSTART_RUN << std::endl;
    std::cout << "  End Run: " << gEND_RUN << std::endl;
    std::cout << "  Detector: " << gDETECTOR << std::endl;
    std::cout << "  Matrix Name: " << gMATRIX_NAME << std::endl;
    std::cout << "  ROIs: " << gROIarrs.size() << std::endl;
    for (size_t roi_index = 0; roi_index < gROIarrs.size(); ++roi_index)
    {
        std::cout << "    ROI " << roi_index << ": ";
        for (const auto &val : gROIarrs.at(roi_index)) { std::cout << val << " "; }
        std::cout << std::endl;
    }
    std::cout << "  Reference Time: ";
    for (const auto &val : gREFERENCE_TIME) { std::cout << val << " "; }
    std::cout << std::endl;
    std::cout << "  Fit Peak: ";
    for (const auto &val : gFIT_PEAK) { std::cout << val << " "; }
    std::cout << std::endl;
    std::cout << "  Chained Ranges: ";
    for (const auto &[chain_start, chain_end] : gCHAIN_RUN_RANGES) { std::cout << chain_start << "-" << chain_end << " "; }
    std::cout << std::endl;
    std::cout << "  Gain With Offset: " << std::boolalpha << gUSE_GAIN_WITH_OFFSET << std::endl;
    std::cout << "  Use Super Settings: " << std::boolalpha << gUSE_SUPER_SETTINGS << std::endl;
    auto conffile = get_output_conffilename();

    std::cout << "Output configuration file will be saved to: " << conffile << std::endl;
    if (!can_create_file(conffile)) throw std::runtime_error("Problem with creating output configuration file\n");
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    // this makes it slower!!!
    // ROOT::EnableImplicitMT();
    // ROOT::EnableThreadSafety();

    TFile *matfile = TFile::Open(gROOTFILE.c_str(), "READ");
    auto   TEMAT_original = get_ill_matrix_from_file(matfile);

    for (size_t roi_index = 0; roi_index < gROIarrs.size(); ++roi_index)
    {
        adjust_peak_energy(TEMAT_original, gROIarrs.at(roi_index));
        std::cout << "Adjusted ROI " << roi_index << " energy to: " << gROIarrs.at(roi_index).at(0) << std::endl;
    }
    if (gFIT_PEAK.size() == 3)
    {
        adjust_peak_energy(TEMAT_original, gFIT_PEAK);
        std::cout << "Adjusted cost peak energy to: " << gFIT_PEAK.at(0) << std::endl;
    }

    if (gUSE_SUPER_SETTINGS)
    {
        gSUPER_SETTINGS.temat_rebin_x = 1;
        gSUPER_SETTINGS.temat_rebin_y = 1;
        gSUPER_SETTINGS.use_gaussian  = false;
        gSUPER_SETTINGS.valid_only    = true;
        // gSUPER_SETTINGS.interpolator_type      = "akima";
        gSUPER_SETTINGS.interpolator_type      = "akima";
        gSUPER_SETTINGS.interpolator_smoothing = false;
        gSUPER_SETTINGS.smoother_type          = TEC::SmootherType::NONE;
        gSUPER_SETTINGS.smoother_par           = -1;
        gSUPER_SETTINGS.cost                   = std::numeric_limits<double>::quiet_NaN();

        run_ccm_super_settings(TEMAT_original, gSUPER_SETTINGS, get_output_conffilename());
        run_chained_ranges(gSUPER_SETTINGS);
        return 0;
    }

    std::cout << "\nTesting following settings: " << std::endl;
    ccm_settings::print_header(std::cout);

    auto start = high_resolution_clock::now();

    auto result = ccm_optimizer_global(std::shared_ptr<TH2>(TEMAT_original),
                                       [&](TH1 *histo) { return get_fwfm(histo, gFIT_PEAK.at(0), gFIT_PEAK.at(1), gFIT_PEAK.at(2)); });

    std::sort(result.begin(), result.end(), [](const auto &a, const auto &b) { return a.cost < b.cost; });

    auto stop     = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration: " << duration.count() << " ms" << std::endl;

    ccm_settings::print_header(std::cout);
    for (const auto &r : result) { r.print_values(std::cout); }

    std::string  minimization_file = get_output_minimization_filename();
    std::fstream out_file(minimization_file.c_str(), std::ios::out);
    if (!out_file) { std::cerr << "Error opening CCM conf file: " << minimization_file << "  for writing" << std::endl; }
    ccm_settings::print_header(out_file);
    for (const auto &r : result) { r.print_values(out_file); }

    gSUPER_SETTINGS = result.front();
    run_ccm_super_settings(TEMAT_original, gSUPER_SETTINGS, get_output_conffilename());
    run_chained_ranges(gSUPER_SETTINGS);

    return 0;
}
