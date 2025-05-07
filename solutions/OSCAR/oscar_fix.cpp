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
#include <TRandom3.h>
#include <TTree.h>

#include "CCM.h"

using namespace std::chrono;
using namespace TEC;

#ifndef DATA_PATH
#define DATA_PATH ""
#endif

// int    bins       = 1000;
// double time_start = 22E12;
// double time_end   = 50E12;

const std::string INPUT_ROOT_FILE =
    // "/home/mbalogh/data/56ni/temat_merged_compressed.root";
    "/home/mbalogh/data/56ni/temat_merged.root";

// const std::string MATRIX_NAME = "strips/summed_compressed_temat_strip_0";
const std::string MATRIX_NAME = "strips/summed_temat_strip_0";

struct ccm_cost
{
    double ks_distance{std::numeric_limits<double>::max()};
    double chi2_cumulative{std::numeric_limits<double>::max()};
    double chi2_normalized{std::numeric_limits<double>::max()};
};

struct ccm_settings
{
    ccm_cost cost{};
    uint     temat_rebin_x{1};
    uint     temat_rebin_y{1};
    // std::vector<RegionOfInterest> ROIs{};
    bool        use_gaussian{true}; // gaussian or polynomial dot product fit
    bool        valid_only{true};
    std::string interpolator_type{"AKIMA"};
    bool        interpolator_smoothing{false};
    std::string smoother_type{""};
    double      smoother_par0{0};
    double      smoother_par1{0};
    double      smoother_par2{0};

    // ROIs
    // correction function
    // reference times

    void print_header(std::ostream &os)
    {
        os << "cost_KS, cost_chi2, cost_chi2_cumul, temat_rebin_x, temat_rebin_y, ROIs, "
              "use_gaussian, valid_only, "
              "interpolator_type, "
              "interpolator_smoothing, smoother_type, smoother_par0, smoother_par1, "
              "smoother_par2 \n ";
    }

    ccm_settings() = default;
    ccm_settings(const uint         rebinx,
                 const uint         rebiny,
                 const bool         gaussian,
                 const bool         valid,
                 const std::string &interpolator)
        : temat_rebin_x(rebinx), temat_rebin_y(rebiny), use_gaussian(gaussian),
          valid_only(valid){};

    ccm_settings(const ccm_settings &)            = default;
    ccm_settings &operator=(const ccm_settings &) = default;
    ~ccm_settings()                               = default;

    void print_values(std::ostream &os)
    {
        os << this->cost.ks_distance << " " << this->cost.chi2_normalized << " "
           << this->cost.chi2_cumulative << " " << this->temat_rebin_x << " "
           << this->temat_rebin_y << " " << std::boolalpha << this->use_gaussian << " "
           << this->valid_only << " " << this->interpolator_type << " "
           << this->interpolator_smoothing << " " << this->smoother_type << " "
           << this->smoother_par0 << " " << this->smoother_par1 << " "
           << this->smoother_par2 << "\n";
    }
};

std::vector<ccm_settings> CCM_OPTIONS;

double calculate_chi2(double *a2, double *a1, int size)
{
    double norm1 = 0;
    for (auto i = 0; i < size; ++i) { norm1 += a1[i]; }
    for (auto i = 0; i < size; ++i) { a1[i] = a1[i] / norm1; }
    double norm2 = 0;
    for (auto i = 0; i < size; ++i) { norm2 += a2[i]; }
    for (auto i = 0; i < size; ++i) { a2[i] = a2[i] / norm2; }

    double chi2 = 0;
    for (auto i = 0; i < size; i++)
    {
        if (a2[i] == 0) continue;
        chi2 += pow(a1[i] - a2[i], 2) / a2[i];
    }
    return chi2;
}

double calculate_cumulative_chi2(double *a2, double *a1, int size)
{
    double chi2 = 0;
    double norm = 0;
    for (auto i = 0; i < size; ++i) { norm += a1[i]; }
    for (auto i = 0; i < size; ++i) { a1[i] = a1[i] / norm; }
    norm = 0;
    for (auto i = 0; i < size; ++i) { norm += a2[i]; }
    for (auto i = 0; i < size; ++i) { a2[i] = a2[i] / norm; }

    double cumul_a1{0}, cumul_a2{0};
    for (auto i = 0; i < size; i++)
    {
        cumul_a1 += a1[i];
        cumul_a2 += a2[i];
        if (cumul_a2 == 0) continue;
        chi2 += pow(cumul_a1 - cumul_a2, 2) / cumul_a2;
    }
    return chi2;
}

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

#include <numeric> // for std::accumulate
#include <vector>

double calculateCenterOfMass(const double *array, int size)
{
    double weightedSum  = 0.0;
    double sumOfWeights = 0.0;

    for (size_t i = 0; i < size; ++i)
    {
        weightedSum += i * array[i];
        sumOfWeights += array[i];
    }

    if (sumOfWeights == 0)
    {
        throw std::runtime_error(
            "Sum of weights is zero, cannot calculate center of mass.");
    }

    return weightedSum / sumOfWeights;
}

void test_dp_methods()
{
    std::unique_ptr<TH2D> TEMAT(get_matrix());
    TH2D                 *TEMAT_copy = (TH2D *)TEMAT->Rebin2D(20, 8, "TEMAT_rebinned");

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    // calculate center of mass of a roi
    const double reference_time_bgn = 4.85e12;
    const double reference_time_end = 5.35e12;

    std::vector<RegionOfInterest> ROIs;
    // ROIs.emplace_back(RegionOfInterest(*TEMAT_copy, 5500, 8000, -100, 1500, 6166));
    {
        int bin_start = TEMAT_copy->GetXaxis()->FindBin(reference_time_bgn);
        int bin_end   = TEMAT_copy->GetXaxis()->FindBin(reference_time_end);

        auto projection = TEMAT_copy->ProjectionY("proj", bin_start, bin_end);

        int  bin_en_start = projection->FindBin(5500);
        int  bin_en_end   = projection->FindBin(8000);
        auto center       = calculateCenterOfMass(&projection->GetArray()[bin_en_start],
                                                  bin_en_end - bin_en_start);
        center            = projection->GetBinCenter(bin_en_start + center);
        std::cout << "center of mass: " << center << std::endl;

        auto center_bin_value =
            projection->GetBinCenter((bin_en_end - bin_en_start) / 2 + bin_en_start);
        std::cout << "center bin value " << center_bin_value << std::endl;
        ROIs.emplace_back(
            RegionOfInterest(*TEMAT_copy, 5500, 8000, -100, 1500, center_bin_value));
    }

    CCM fix(TEMAT_copy, ROIs, reference_time_bgn, reference_time_end);

    fix.SetCorrectionFunction(fcn, "");
    fix.CalculateEnergyShifts(6);
    fix.UseGaussianResult();

    std::unique_ptr<TH1D> proj_old(
        TEMAT->ProjectionY("proj_old", 1030, TEMAT->GetXaxis()->GetNbins()));

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::unique_ptr<TH2D> TEMAT_fixed(fix.FixMatrix(TEMAT.get()));
    std::unique_ptr<TH1D> proj_gauss(TEMAT_fixed->ProjectionY(
        "proj_gauss", 1030, TEMAT_fixed->GetXaxis()->GetNbins()));
    {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "time " << duration / 1e6 << " s                           "
                  << std::endl;
        t1 = t2;
    }

    fix.UsePolynomialResult();
    TEMAT_fixed.reset((fix.FixMatrix(TEMAT.get())));
    std::unique_ptr<TH1D> proj_poly2(TEMAT_fixed->ProjectionY(
        "proj_poly2", 1030, TEMAT_fixed->GetXaxis()->GetNbins()));
    {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "time " << duration / 1e6 << " s                           "
                  << std::endl;
        t1 = t2;
    }

    fix.UseMaxDPResult();
    TEMAT_fixed.reset((fix.FixMatrix(TEMAT.get())));
    std::unique_ptr<TH1D> proj_maxdp(TEMAT_fixed->ProjectionY(
        "proj_maxdp", 1030, TEMAT_fixed->GetXaxis()->GetNbins()));
    {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "time " << duration / 1e6 << " s                           "
                  << std::endl;
        t1 = t2;
    }

    std::unique_ptr<TH1D> model(
        get_model_histogram(TEMAT.get(), reference_time_bgn, reference_time_end));
    model->Scale(proj_old->Integral() / model->Integral());

    TApplication app("app", 0, 0);
    proj_gauss->SetLineColor(kRed);
    proj_poly2->SetLineColor(kBlue);
    proj_maxdp->SetLineColor(kOrange);
    proj_old->SetLineColor(kBlack);
    model->SetLineColor(kCyan);
    proj_gauss->Draw("hist");
    proj_poly2->Draw("same hist");
    proj_maxdp->Draw("same hist");
    proj_old->Draw("same hist");
    model->Draw("same hist");

    double chi2_test = model->Chi2Test(proj_gauss.get(), "UU");
    std::cout << "chi2 test: " << chi2_test << std::endl;

    chi2_test = model->Chi2Test(proj_poly2.get(), "UU");
    std::cout << "chi2 test: " << chi2_test << std::endl;

    chi2_test = model->Chi2Test(proj_maxdp.get(), "UU");
    std::cout << "chi2 test: " << chi2_test << std::endl;

    for (int i = 1; i <= proj_gauss->FindBin(2200); i++)
    {
        proj_gauss->SetBinContent(i, 0);
        proj_poly2->SetBinContent(i, 0);
        proj_maxdp->SetBinContent(i, 0);
        proj_old->SetBinContent(i, 0);
        model->SetBinContent(i, 0);
    }

    for (int i = proj_gauss->FindBin(9000); i <= proj_gauss->GetXaxis()->GetNbins(); i++)
    {
        proj_gauss->SetBinContent(i, 0);
        proj_poly2->SetBinContent(i, 0);
        proj_maxdp->SetBinContent(i, 0);
        proj_old->SetBinContent(i, 0);
        model->SetBinContent(i, 0);
    }

    {
        double ks_distance = model->KolmogorovTest(proj_gauss.get(), "M");
        std::cout << std::setprecision(8) << "KS distance gauss:  " << ks_distance
                  << std::endl;

        int  bin_start = model->FindBin(2200);
        int  bin_end   = model->FindBin(9000);
        auto val       = calculate_cumulative_chi2(&model->GetArray()[bin_start],
                                                   &proj_gauss->GetArray()[bin_start],
                                                   bin_end - bin_start);
        std::cout << std::setprecision(8) << "Cumulative chi2: " << val << std::endl;

        auto val2 =
            calculate_chi2(&model->GetArray()[bin_start],
                           &proj_gauss->GetArray()[bin_start], bin_end - bin_start);
        // std::cout << "chi2 " << val2 << std::endl;

        std::cout << std::setprecision(8) << "chi2: " << val2 << std::endl;
    }
    {
        double ks_distance = model->KolmogorovTest(proj_poly2.get(), "M");
        std::cout << "KS distance poly2:  " << ks_distance << std::endl;

        int  bin_start = model->FindBin(2200);
        int  bin_end   = model->FindBin(9000);
        auto val       = calculate_cumulative_chi2(&model->GetArray()[bin_start],
                                                   &proj_poly2->GetArray()[bin_start],
                                                   bin_end - bin_start);
        std::cout << std::setprecision(8) << "Cumulative chi2: " << val << std::endl;

        auto val2 =
            calculate_chi2(&model->GetArray()[bin_start],
                           &proj_poly2->GetArray()[bin_start], bin_end - bin_start);
        // std::cout << "chi2 " << val2 << std::endl;

        std::cout << std::setprecision(8) << "chi2: " << val2 << std::endl;
    }
    {
        double ks_distance = model->KolmogorovTest(proj_maxdp.get(), "M");
        std::cout << "KS distance maxdp: " << ks_distance << std::endl;

        int  bin_start = model->FindBin(2200);
        int  bin_end   = model->FindBin(9000);
        auto val       = calculate_cumulative_chi2(&model->GetArray()[bin_start],
                                                   &proj_maxdp->GetArray()[bin_start],
                                                   bin_end - bin_start);
        std::cout << std::setprecision(8) << "Cumulative chi2: " << val << std::endl;

        auto val2 =
            calculate_chi2(&model->GetArray()[bin_start],
                           &proj_maxdp->GetArray()[bin_start], bin_end - bin_start);
        // std::cout << "chi2 " << val2 << std::endl;

        std::cout << std::setprecision(8) << "chi2: " << val2 << std::endl;
    }
    // auto interpol = fix.GetInterpolationGraph(0, 20, true);
    // interpol->SetLineColor(kRed);
    // new TCanvas();
    // fix.GetInterpolationGraph(0, 10, true)->Draw("AP");

    app.Run();
}

void test_interpolation()
{
    std::unique_ptr<TH2D> TEMAT(get_matrix());
    TH2D                 *TEMAT_copy = (TH2D *)TEMAT->Rebin2D(20, 2, "TEMAT_rebinned");

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    std::vector<RegionOfInterest> ROIs;
    ROIs.emplace_back(RegionOfInterest(*TEMAT_copy, 5500, 8000, -100, 1500, 6166));

    const double reference_time_bgn = 4.9e12;
    const double reference_time_end = 5.1e12;

    CCM fix(TEMAT_copy, ROIs, reference_time_bgn, reference_time_end);

    fix.SetCorrectionFunction(fcn, "");
    fix.CalculateEnergyShifts(6);
    fix.DisableInterpolation();
    // fix.CalculateCorrectionFits();

    auto no_interpol = fix.GetInterpolationGraph(0, 20, true);
    no_interpol->SetLineColor(kGreen);

    std::unique_ptr<TH2D> TEMAT_fixed(fix.FixMatrix(TEMAT.get()));

    std::unique_ptr<TH1D> proj_old(
        TEMAT->ProjectionY("proj_old", 1030, TEMAT->GetXaxis()->GetNbins()));

    std::unique_ptr<TH1D> proj_fixed_NI(TEMAT_fixed->ProjectionY(
        "proj_fixed_NI", 1030, TEMAT_fixed->GetXaxis()->GetNbins()));

    fix.EnableInterpolation();
    fix.ConfigureShiftInterpolator(ROOT::Math::Interpolation::Type::kCSPLINE);

    std::unique_ptr<TH2D> TEMAT_fixed_2(fix.FixMatrix(TEMAT.get()));

    std::unique_ptr<TH1D> proj_fixed_I(TEMAT_fixed_2->ProjectionY(
        "proj_fixed_I", 1030, TEMAT_fixed_2->GetXaxis()->GetNbins()));

    TApplication app("app", 0, 0);
    proj_fixed_NI->SetLineColor(kRed);
    proj_fixed_I->SetLineColor(kBlue);
    proj_old->SetLineColor(kBlack);
    proj_fixed_NI->Draw("hist");
    proj_fixed_I->Draw("same hist");
    proj_old->Draw("same hist");

    auto interpol = fix.GetInterpolationGraph(0, 20, true);
    interpol->SetLineColor(kRed);
    // new TCanvas();
    // fix.GetInterpolationGraph(0, 10, true)->Draw("AP");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(interpol);
    mg->Add(no_interpol);
    mg->SetTitle("shifts");

    new TCanvas();
    mg->Draw("AL");

    new TCanvas();
    interpol->Draw("AP");

    app.Run();
}

int main(int argc, char **argv)
{
    test_dp_methods();
    test_interpolation();

    CCM_OPTIONS.emplace_back(1, 1, false, true, "");
    CCM_OPTIONS.emplace_back(1, 1, false, true, "cspline");
    CCM_OPTIONS.emplace_back(1, 1, false, true, "polynomial");
    CCM_OPTIONS.emplace_back(1, 1, false, true, "akima");

    CCM_OPTIONS.emplace_back(2, 1, false, true, "");
    CCM_OPTIONS.emplace_back(2, 1, false, true, "cspline");
    CCM_OPTIONS.emplace_back(2, 1, false, true, "polynomial");
    CCM_OPTIONS.emplace_back(2, 1, false, true, "akima");

    CCM_OPTIONS.emplace_back(3, 1, false, true, "");
    CCM_OPTIONS.emplace_back(3, 1, false, true, "cspline");
    CCM_OPTIONS.emplace_back(3, 1, false, true, "polynomial");
    CCM_OPTIONS.emplace_back(3, 1, false, true, "akima");

    CCM_OPTIONS.emplace_back(5, 1, false, true, "");
    CCM_OPTIONS.emplace_back(5, 1, false, true, "cspline");
    CCM_OPTIONS.emplace_back(5, 1, false, true, "polynomial");
    CCM_OPTIONS.emplace_back(5, 1, false, true, "akima");

    std::unique_ptr<TH2D> TEMAT(get_matrix());

    TF1 fcn("gain_fcn", "[0]*x", 0, 1);

    for (auto option : CCM_OPTIONS)
    {

        std::unique_ptr<TH2D> TEMAT_copy((TH2D *)TEMAT->Clone("TEMAT_copy"));
        TEMAT_copy->SetDirectory(0); // claim ownership of the object
        TEMAT_copy->Rebin2D(option.temat_rebin_x, option.temat_rebin_y);

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        std::vector<RegionOfInterest> ROIs;
        ROIs.emplace_back(RegionOfInterest(*TEMAT_copy, 5500, 8000, -100, 1500, 6166));
        // ROIs.emplace_back(RegionOfInterest(*TEMAT, 3000, 4300, -500, 1500, 3900));

        const double reference_time_bgn = 4.9e12;
        const double reference_time_end = 5.1e12;

        // create CCM object

        CCM fix(TEMAT_copy, ROIs, reference_time_bgn, reference_time_end);

        fix.SetCorrectionFunction(fcn, "");
        fix.CalculateEnergyShifts(6);
        if (option.interpolator_type.empty()) fix.DisableInterpolation(0);
        else
            fix.EnableInterpolation(0);
        if (option.use_gaussian) fix.UseGaussianResult();
        fix.CalculateCorrectionFits();

        std::unique_ptr<TH2D> TEMAT_fixed(fix.FixMatrix(TEMAT.get()));
        TEMAT_fixed->SetDirectory(0);

        std::unique_ptr<TH1D> proj_old(
            TEMAT->ProjectionY("proj_old", 1030, TEMAT->GetXaxis()->GetNbins()));

        std::unique_ptr<TH1D> proj_fixed(TEMAT_fixed->ProjectionY(
            "proj_fixed", 1030, TEMAT_fixed->GetXaxis()->GetNbins()));

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "TOTAL DURATION OF " << duration / (double)1E6 << " seconds"
                  << std::endl;

        std::unique_ptr<TH1D> model(
            get_model_histogram(TEMAT.get(), reference_time_bgn, reference_time_end));
        model->Scale(proj_fixed->Integral() / model->Integral());

        // evaluate cost
        double ks_distance = model->KolmogorovTest(proj_fixed.get(), "M");
        // std::cout << std::setprecision(9) << "KS-distance: " << ks_distance <<
        // std::endl;
        option.cost.ks_distance = ks_distance;
        {
            int  bin_start = model->FindBin(2200);
            int  bin_end   = model->FindBin(9000);
            auto val       = calculate_cumulative_chi2(&model->GetArray()[bin_start],
                                                       &proj_fixed->GetArray()[bin_start],
                                                       bin_end - bin_start);
            // std::cout << "Cumulative chi2 test: " << val << std::endl;

            auto val2 =
                calculate_chi2(&model->GetArray()[bin_start],
                               &proj_fixed->GetArray()[bin_start], bin_end - bin_start);
            // std::cout << "chi2 " << val2 << std::endl;

            option.cost.chi2_cumulative = val;
            option.cost.chi2_normalized = val2;
            option.print_values(std::cout);
            std::cout << "KS-distance: " << ks_distance << std::endl;
            std::cout << "chi2: " << val2 << std::endl;
            std::cout << "Cumulative chi2: " << val << std::endl;
        }

        // TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        // legend->AddEntry(proj_old, "Original", "l");
        // legend->AddEntry(proj_fixed, "Corrected", "l");
        // legend->AddEntry(model, "Model", "l");
        // legend->Draw();

        // for (int bin = model->GetXaxis()->FindBin(3000); bin > 0; --bin)
        // {
        //     model->SetBinContent(bin, 0);
        //     proj_fixed->SetBinContent(bin, 0);
        // }

        // for (int bin = model->GetXaxis()->FindBin(8000); bin <=
        // model->GetXaxis()->GetNbins();
        //      bin++)
        // {
        //     model->SetBinContent(bin, 0);
        //     proj_fixed->SetBinContent(bin, 0);
        // }

        // new TCanvas();
        // TEMAT->Rebin2D(6, 4);
        // TEMAT->Draw("COLZ");

        // new TCanvas();
        // TEMAT_fixed->Draw("COLZ");

        // new TCanvas();

        // fix.GetROIShifts(0, true)->Draw("AP");

        // new TCanvas();
        // model->Draw();
        // proj_fixed->Scale(1. / proj_fixed->Integral());
        // model->Scale(1. / model->Integral());
        // proj_fixed->Draw("SAME HIST");

        // model->ResetStats();
        // proj_fixed->ResetStats();

        // double ks_distance = model->Chi2Test(proj_fixed, "UU P ");
        // std::cout << std::setprecision(9) << "p-value: " << ks_distance << std::endl;

        // ks_distance = model->KolmogorovTest(proj_old, "M");
        // std::cout << std::setprecision(9) << "KS-distance: " << ks_distance <<
        // std::endl;

        // fix.SaveToRootFile();
    }

    TApplication app("app", &argc, argv);

    TGraph g_ks;
    TGraph g_chi2;
    TGraph g_chi2_cumul;

    int iterator = 0;

    std::array<double, 3> norm{0, 0, 0};
    for (const auto option : CCM_OPTIONS)
    {
        norm[0] = std::max(norm[0], option.cost.ks_distance);
        norm[1] = std::max(norm[1], option.cost.chi2_normalized);
        norm[2] = std::max(norm[2], option.cost.chi2_cumulative);
    }

    for (const auto option : CCM_OPTIONS)
    {
        g_ks.AddPoint(iterator, option.cost.ks_distance);
        g_chi2.AddPoint(iterator, option.cost.chi2_normalized);
        g_chi2_cumul.AddPoint(iterator, option.cost.chi2_cumulative);
        iterator++;
    }

    g_ks.SetLineColor(kBlue);
    g_chi2.SetLineColor(kRed);
    g_chi2_cumul.SetLineColor(kGreen);
    g_ks.SetLineWidth(2);
    g_chi2.SetLineWidth(2);
    g_chi2_cumul.SetLineWidth(2);
    g_ks.SetMarkerStyle(20);
    g_chi2.SetMarkerStyle(20);
    g_chi2_cumul.SetMarkerStyle(20);

    // Create a TMultiGraph and add the TGraph objects to it
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(&g_ks);
    mg->Add(&g_chi2);
    mg->Add(&g_chi2_cumul);
    mg->SetTitle("Results");

    // Create a canvas and draw the TMultiGraph
    TCanvas *c1 = new TCanvas("c1", "TMultiGraph Example", 800, 600);
    mg->Draw("ALP");

    new TCanvas();
    g_ks.Draw("ALP");
    g_chi2.Draw("ALP");
    g_chi2_cumul.Draw("ALP");

    // Run the application
    app.Run();

    // app.Run();
    return 0;
}