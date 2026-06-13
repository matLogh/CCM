#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TString.h>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <vector>

// Global pointers so ScaleY(...) can modify the currently loaded detector
TFile*   gFile_timeevo = nullptr;
TCanvas* gCanvas_timeevo = nullptr;
TH2D*    h_raw = nullptr;
TH2D*    h_fixed = nullptr;
TGraph*  gBadRunGraph = nullptr;
int      gCurrentDet = -1;

void SetNicePad()
{
    gPad->SetCrosshair();
    gPad->SetLogz();

    // Reduced margins to remove extra white space
    gPad->SetLeftMargin(0.025);
    gPad->SetRightMargin(0.035);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.08);
}

void DrawBadRunOverlay()
{
    if (gBadRunGraph) {
        gBadRunGraph->Draw("P SAME");
    }
}

void LoadDet(int det,
             bool show_gain_only = false,
             bool compact_width = false,
             int custom_canvas_width = 0)
{
    gCurrentDet = det;

    TString file_name =
        Form("output/combined_offset_and_gain_timeevo/"
             "RunVsEnergy_all_ranges_det%d_correctedTimeEvo.root",
             det);

    TString top_name = show_gain_only
        ? Form("RunVsEnergy_det%d_fixed_all_ranges", det)
        : Form("RunVsEnergy_det%d_all_ranges", det);

    TString fixed_name = show_gain_only
        ? Form("RunVsEnergy_det%d_gain_only_all_ranges", det)
        : Form("RunVsEnergy_det%d_fixed_all_ranges", det);

    TString top_label = show_gain_only ? "offset+gain" : "raw";
    TString fixed_label = show_gain_only ? "gain-only" : "fixed";

    // Close previous file if one is open
    if (gFile_timeevo) {
        gFile_timeevo->Close();
        delete gFile_timeevo;
        gFile_timeevo = nullptr;
    }

    gFile_timeevo = TFile::Open(file_name, "READ");

    if (!gFile_timeevo || gFile_timeevo->IsZombie()) {
        std::cerr << "Cannot open file: " << file_name << std::endl;
        return;
    }

    h_raw = (TH2D*)gFile_timeevo->Get(top_name);
    h_fixed = (TH2D*)gFile_timeevo->Get(fixed_name);

    if (!h_raw) {
        std::cerr << "Cannot find histogram: " << top_name << std::endl;
        return;
    }

    if (!h_fixed) {
        std::cerr << "Cannot find histogram: " << fixed_name << std::endl;
        return;
    }

    h_raw->SetStats(0);
    h_fixed->SetStats(0);

    // Delete old canvas if it exists
    if (gCanvas_timeevo) {
        delete gCanvas_timeevo;
        gCanvas_timeevo = nullptr;
    }

    if (gBadRunGraph) {
        delete gBadRunGraph;
        gBadRunGraph = nullptr;
    }

    int canvas_width = custom_canvas_width > 0
        ? custom_canvas_width
        : (compact_width ? 1512 : 3840);
    int canvas_height = compact_width ? 820 : 864;

    // Full width spans two 4K monitors; compact width fits better on a laptop display.
    gCanvas_timeevo =
        new TCanvas(Form("c_det%d", det),
                    Form("%s and %s det%d",
                         top_label.Data(),
                         fixed_label.Data(),
                         det),
                    0, 0, canvas_width, canvas_height);

    gCanvas_timeevo->SetWindowPosition(0, 0);
    gCanvas_timeevo->SetWindowSize(canvas_width, canvas_height);
    gCanvas_timeevo->SetCrosshair();

    // Small gaps between pads
    gCanvas_timeevo->Divide(1, 2, 0.001, 0.001);

    gCanvas_timeevo->cd(1);
    SetNicePad();
    h_raw->Draw("colz");

    gCanvas_timeevo->cd(2);
    SetNicePad();
    h_fixed->Draw("colz");
    DrawBadRunOverlay();

    gCanvas_timeevo->Update();

    std::cout << "Loaded detector " << det << std::endl;
    std::cout << top_label << ": " << top_name << std::endl;
    std::cout << fixed_label << ": " << fixed_name << std::endl;
}

void LoadDetWidth(int det, bool show_gain_only, int canvas_width)
{
    LoadDet(det, show_gain_only, false, canvas_width);
}

void OverlayBadRuns(const char* bad_runs_file,
                    const std::vector<double>& marker_y_values)
{
    if (!h_fixed || !gCanvas_timeevo) {
        std::cerr << "No detector loaded. First run LoadDet(det)." << std::endl;
        return;
    }

    if (marker_y_values.empty()) {
        std::cerr << "No y values were provided for the bad-run overlay." << std::endl;
        return;
    }

    std::ifstream in(bad_runs_file);
    if (!in) {
        std::cerr << "Cannot open bad-runs file: " << bad_runs_file << std::endl;
        return;
    }

    std::vector<int> bad_runs;
    std::string line;

    while (std::getline(in, line)) {
        std::stringstream ss(line);
        int run_number = 0;

        if (ss >> run_number) {
            bad_runs.push_back(run_number);
        }
    }

    if (bad_runs.empty()) {
        std::cerr << "No bad run numbers found in: " << bad_runs_file << std::endl;
        return;
    }

    std::vector<double> run_numbers;
    std::vector<double> y_values;
    run_numbers.reserve(bad_runs.size() * marker_y_values.size());
    y_values.reserve(bad_runs.size() * marker_y_values.size());

    for (double y_value : marker_y_values) {
        for (int run_number : bad_runs) {
            run_numbers.push_back(run_number);
            y_values.push_back(y_value);
        }
    }

    if (gBadRunGraph) {
        delete gBadRunGraph;
        gBadRunGraph = nullptr;
    }

    gBadRunGraph = new TGraph(run_numbers.size(),
                              run_numbers.data(),
                              y_values.data());
    gBadRunGraph->SetName("gBadRunGraph");
    gBadRunGraph->SetTitle("bad run number");
    gBadRunGraph->SetMarkerColor(kRed);
    gBadRunGraph->SetMarkerStyle(5);
    gBadRunGraph->SetMarkerSize(1.0);

    gCanvas_timeevo->cd(2);
    DrawBadRunOverlay();
    gPad->Modified();
    gCanvas_timeevo->Modified();
    gCanvas_timeevo->Update();

    std::cout << "Overlayed " << run_numbers.size()
              << " bad run markers from " << bad_runs_file
              << " using " << marker_y_values.size()
              << " y value(s)" << std::endl;
}

void OverlayBadRuns(const char* bad_runs_file,
                    std::initializer_list<double> marker_y_values)
{
    OverlayBadRuns(bad_runs_file, std::vector<double>(marker_y_values));
}

template <typename... YValues>
void OverlayBadRuns(const char* bad_runs_file,
                    double first_y_value,
                    YValues... more_y_values)
{
    std::vector<double> marker_y_values =
        {first_y_value, static_cast<double>(more_y_values)...};
    OverlayBadRuns(bad_runs_file, marker_y_values);
}

void OverlayBadRuns(double y_value)
{
    OverlayBadRuns("../bad_runs/det7_badruns.txt", y_value);
}

void OverlayBadRuns(std::initializer_list<double> marker_y_values)
{
    OverlayBadRuns("../bad_runs/det7_badruns.txt", marker_y_values);
}

template <typename... YValues>
void OverlayBadRuns(double first_y_value, YValues... more_y_values)
{
    std::vector<double> marker_y_values =
        {first_y_value, static_cast<double>(more_y_values)...};
    OverlayBadRuns("../bad_runs/det7_badruns.txt", marker_y_values);
}

void ScaleY(double y_low, double y_high)
{
    if (!h_raw || !h_fixed || !gCanvas_timeevo) {
        std::cerr << "No detector loaded. First run LoadDet(det)." << std::endl;
        return;
    }

    h_raw->GetYaxis()->SetRangeUser(y_low, y_high);
    h_fixed->GetYaxis()->SetRangeUser(y_low, y_high);

    gCanvas_timeevo->cd(1);
    h_raw->Draw("colz");

    gCanvas_timeevo->cd(2);
    h_fixed->Draw("colz");
    DrawBadRunOverlay();

    gCanvas_timeevo->Modified();
    gCanvas_timeevo->Update();

    std::cout << "Set Y range for det " << gCurrentDet
              << " to " << y_low << " - " << y_high << std::endl;
}

void ScaleX(double x_low, double x_high)
{
    if (!h_raw || !h_fixed || !gCanvas_timeevo) {
        std::cerr << "No detector loaded. First run LoadDet(det)." << std::endl;
        return;
    }

    h_raw->GetXaxis()->SetRangeUser(x_low, x_high);
    h_fixed->GetXaxis()->SetRangeUser(x_low, x_high);

    gCanvas_timeevo->cd(1);
    h_raw->Draw("colz");

    gCanvas_timeevo->cd(2);
    h_fixed->Draw("colz");
    DrawBadRunOverlay();

    gCanvas_timeevo->Modified();
    gCanvas_timeevo->Update();

    std::cout << "Set X range for det " << gCurrentDet
              << " to " << x_low << " - " << x_high << std::endl;
}

void ResetY()
{
    if (!h_raw || !h_fixed || !gCanvas_timeevo) {
        std::cerr << "No detector loaded. First run LoadDet(det)." << std::endl;
        return;
    }

    h_raw->GetYaxis()->UnZoom();
    h_fixed->GetYaxis()->UnZoom();

    gCanvas_timeevo->cd(1);
    h_raw->Draw("colz");

    gCanvas_timeevo->cd(2);
    h_fixed->Draw("colz");
    DrawBadRunOverlay();

    gCanvas_timeevo->Modified();
    gCanvas_timeevo->Update();

    std::cout << "Reset Y range for det " << gCurrentDet << std::endl;
}
