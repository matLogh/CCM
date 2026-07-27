#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.cpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TRint.h>

std::vector<int>         gRUNLIST;
std::vector<std::string> gCRYSTALLIST;
// std::string gCONFDIR{""};
std::string gDIR{""};
// minimum number of holes is by default 2 to account for starting and ending padding of the temat
const int gMINIMUM_NUMBER_OF_HOLES_TO_IGNORE = 2;
bool      gMODIFY_TIMEEVO_CONF               = false;
bool      gDRAW_CUTS                         = true;
double    gTHRESHOLD                         = 0.7;
bool      gPRINTDEAD                         = false;
bool      gSAVE_RESULT                       = false;
std::string gSAVE_FILENAME                   = "";
bool      gUSE_MATRIX_ZOOM                   = false;
double    gMATRIX_ZOOM_LOW                   = 0.0;
double    gMATRIX_ZOOM_HIGH                  = 0.0;

std::vector<std::shared_ptr<TObject>> gROOTOBJECTS;

const Long64_t MINUTE_TO_TIMESTAMPS       = 60l * 100000000l;
const double   GAIN_FACTOR_TO_REMOVE_DATA = -1.;

struct ValidationCheckResult
{
    std::vector<std::pair<double, double>> cut_times;
    double                                average_integral   = 0.0;
    double                                threshold_integral = 0.0;
};

Long64_t get_first_timestamp(const std::shared_ptr<TH2> temat)
{
    const int nbinsx = temat->GetXaxis()->GetNbins();
    const int nbinsy = temat->GetYaxis()->GetNbins();

    for (int binx = 1; binx < nbinsx; binx++)
    {
        if (temat->Integral(binx, binx, 1, nbinsy) > 0) { return static_cast<Long64_t>(temat->GetXaxis()->GetBinLowEdge(binx)); }
    }
    throw std::runtime_error("empty matrix?");
}

std::shared_ptr<TH1> make_events_per_time_plot(const std::shared_ptr<TH2> temat)
{
    if (!temat) { throw std::runtime_error("Cannot make events-per-time plot from a null matrix"); }

    std::shared_ptr<TH1> events_per_time(
        temat->ProjectionX(Form("EventsPerTime_%s", temat->GetName()), 1, temat->GetYaxis()->GetNbins(), "e"));
    events_per_time->SetDirectory(nullptr);
    events_per_time->SetTitle(Form("Events per time for %s", temat->GetName()));

    const char *time_axis_title = temat->GetXaxis()->GetTitle();
    events_per_time->GetXaxis()->SetTitle(time_axis_title && std::string(time_axis_title).size() > 0 ? time_axis_title : "Time");
    events_per_time->GetYaxis()->SetTitle("Events");
    events_per_time->SetLineColor(kBlue + 1);
    events_per_time->SetLineWidth(2);

    return events_per_time;
}

bool modify_conffile(const std::string &filename, const std::vector<std::pair<double, double>> &cut_ranges)
{
    // std::cout << "modifying conf file " << filename << " with " << cut_ranges.size() << " cut ranges" << std::endl;
    // for (const auto &[start, end] : cut_ranges)
    // {
    //     std::cout << "   from " << start * MINUTE_TO_TIMESTAMPS << " to " << end * MINUTE_TO_TIMESTAMPS << std::endl;
    // }
    if (cut_ranges.empty()) { return false; }

    // convert cut ranges to timestamp format and sort them by start time
    std::vector<std::pair<Long64_t, Long64_t>> cut_ranges_ts;
    cut_ranges_ts.reserve(cut_ranges.size());
    for (const auto &cut : cut_ranges)
    {
        cut_ranges_ts.push_back({static_cast<Long64_t>(cut.first * static_cast<double>(MINUTE_TO_TIMESTAMPS)),
                                 static_cast<Long64_t>(cut.second * static_cast<double>(MINUTE_TO_TIMESTAMPS))});
        // std::cout << "TS from " << cut_ranges_ts.back().first << " to " << cut_ranges_ts.back().second << std::endl;
        // std::cout << "   minutes from " << cut_ranges_ts.back().first / MINUTE_TO_TIMESTAMPS << " to "
        //           << cut_ranges_ts.back().second / MINUTE_TO_TIMESTAMPS << std::endl;
    }

    std::sort(cut_ranges_ts.begin(), cut_ranges_ts.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    cut_ranges_ts.erase(std::remove_if(cut_ranges_ts.begin(), cut_ranges_ts.end(), [](const auto &cut) { return cut.second <= cut.first; }),
                        cut_ranges_ts.end());
    if (cut_ranges_ts.empty()) { return false; }

    std::vector<std::pair<Long64_t, Long64_t>> merged_cut_ranges_ts;
    merged_cut_ranges_ts.reserve(cut_ranges_ts.size());
    for (const auto &cut : cut_ranges_ts)
    {
        if (merged_cut_ranges_ts.empty() || cut.first > merged_cut_ranges_ts.back().second) { merged_cut_ranges_ts.push_back(cut); }
        else
        {
            merged_cut_ranges_ts.back().second = std::max(merged_cut_ranges_ts.back().second, cut.second);
        }
    }
    cut_ranges_ts = std::move(merged_cut_ranges_ts);

    // helper function to format a line for the conf file
    auto format_line = [](Long64_t ts_start, Long64_t ts_end, double gain) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(0) << std::setw(22) << ts_start << std::setw(22) << ts_end << std::fixed << std::setprecision(10)
            << std::setw(22) << gain;
        return out.str();
    };

    const std::filesystem::path conf_path(filename);
    if (conf_path.has_parent_path())
    {
        std::error_code ec;
        std::filesystem::create_directories(conf_path.parent_path(), ec);
        if (ec)
        {
            std::cerr << "error creating conf directory " << conf_path.parent_path() << ": " << ec.message() << std::endl;
            return false;
        }
    }

    // if file exist, create a backup
    const bool file_exists = std::filesystem::exists(filename);
    if (file_exists)
    {
        const std::string backup_filename = filename + ".backup_before_validation";
        std::error_code   ec;
        std::filesystem::copy_file(filename, backup_filename, std::filesystem::copy_options::overwrite_existing, ec);
        if (ec)
        {
            std::cerr << "error creating backup file " << backup_filename << " from " << filename << ": " << ec.message() << std::endl;
            return false;
        }
    }

    std::vector<std::string> lines;
    bool                     modified = false;

    if (!file_exists)
    {
        Long64_t cursor = 0;
        lines.push_back("#             TS_start                TS_end                  gain");
        for (const auto &cut : cut_ranges_ts)
        {
            if (cut.first > cursor) { lines.push_back(format_line(cursor, cut.first, 1.0)); }
            if (cut.second > cut.first) { lines.push_back(format_line(cut.first, cut.second, GAIN_FACTOR_TO_REMOVE_DATA)); }
            cursor = std::max(cursor, cut.second);
        }
        if (cursor < std::numeric_limits<Long64_t>::max()) { lines.push_back(format_line(cursor, std::numeric_limits<Long64_t>::max(), 1.0)); }
        modified = true;
    }
    else
    {
        std::ifstream conf_file(filename);
        if (!conf_file.is_open())
        {
            std::cerr << "error opening conf file " << filename << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(conf_file, line))
        {
            std::istringstream ls(line);
            Long64_t           ts_start = 0;
            Long64_t           ts_end   = 0;
            double             gain     = 0.0;

            if (!(ls >> ts_start >> ts_end >> gain))
            {
                lines.push_back(line);
                continue;
            }

            if (ts_end <= ts_start)
            {
                lines.push_back(line);
                continue;
            }

            bool     line_modified = false;
            Long64_t cursor        = ts_start;
            for (const auto &cut : cut_ranges_ts)
            {
                if (cut.second <= cursor) { continue; }
                if (cut.first >= ts_end) { break; }

                if (cursor < cut.first)
                {
                    lines.push_back(format_line(cursor, std::min(cut.first, ts_end), gain));
                    cursor = std::min(cut.first, ts_end);
                }

                const Long64_t dead_start = std::max(cursor, cut.first);
                const Long64_t dead_end   = std::min(ts_end, cut.second);
                if (dead_start < dead_end)
                {
                    lines.push_back(format_line(dead_start, dead_end, GAIN_FACTOR_TO_REMOVE_DATA));
                    cursor        = dead_end;
                    line_modified = true;
                }
            }

            if (line_modified)
            {
                if (cursor < ts_end) { lines.push_back(format_line(cursor, ts_end, gain)); }
                modified = true;
            }
            else
            {
                lines.push_back(line);
            }
        }
        conf_file.close();
    }

    if (!modified) { return false; }

    const std::string tmp_filename = filename + ".tmp_validation_write";
    std::ofstream     out_file(tmp_filename, std::ios::trunc);
    if (!out_file.is_open())
    {
        std::cerr << "error writing temp conf file " << tmp_filename << std::endl;
        return false;
    }
    for (const auto &ln : lines) { out_file << ln << "\n"; }
    out_file.close();

    std::error_code ec;
    std::filesystem::rename(tmp_filename, filename, ec);
    if (ec)
    {
        std::filesystem::remove(filename, ec);
        ec.clear();
        std::filesystem::rename(tmp_filename, filename, ec);
        if (ec)
        {
            std::cerr << "error replacing conf file " << filename << " with " << tmp_filename << ": " << ec.message() << std::endl;
            return false;
        }
    }

    return true;
}

std::string get_validation_output_filename(const int run, const std::string &crystal)
{
    if (!gSAVE_FILENAME.empty()) { return gSAVE_FILENAME; }
    return "vCheck_run" + fourCharInt(run) + "_crys" + crystal + ".root";
}

void draw_cuts(std::shared_ptr<TH2> temat,
               const ValidationCheckResult &validation,
               const std::string          &save_filename = "")
{
    gROOTOBJECTS.push_back(temat);
    auto events_per_time = make_events_per_time_plot(temat);
    gROOTOBJECTS.push_back(events_per_time);

    std::shared_ptr<TCanvas> c(new TCanvas(Form("Cuts_%s", temat->GetName()), Form("Cuts_%s", temat->GetName()), 1400, 1000));
    gROOTOBJECTS.push_back(c);
    c->Divide(1, 2);
    c->cd(1);
    c->SetGrid();
    c->SetCrosshair(1);
    if (gUSE_MATRIX_ZOOM) { temat->GetYaxis()->SetRangeUser(gMATRIX_ZOOM_LOW, gMATRIX_ZOOM_HIGH); }
    temat->Draw("colz");

    for (const auto &cut : validation.cut_times)
    {
        std::shared_ptr<TGraph> gr(new TGraph());
        gROOTOBJECTS.push_back(gr);
        gr->SetName(Form("CutGraph_%f_%f", cut.first, cut.second));
        auto ymin = gUSE_MATRIX_ZOOM ? gMATRIX_ZOOM_LOW : temat->GetYaxis()->GetXmin();
        auto ymax = gUSE_MATRIX_ZOOM ? gMATRIX_ZOOM_HIGH : temat->GetYaxis()->GetXmax();
        gr->AddPoint(cut.first, ymin);
        gr->AddPoint(cut.first, ymax);
        gr->AddPoint(cut.second, ymax);
        gr->AddPoint(cut.second, ymin);
        gr->SetFillColor(kRed);
        gr->SetFillStyle(3013);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(2);
        gr->SetBit(TObject::kCannotPick);

        gr->Draw("same F");
    }

    c->cd(2);
    gPad->SetGrid();
    gPad->SetCrosshair(1);
    const double reference_max =
        std::max({events_per_time->GetMaximum(), validation.average_integral, validation.threshold_integral, 1.0});
    events_per_time->SetMaximum(reference_max * 1.05);
    events_per_time->Draw("hist");
    const double ymax = events_per_time->GetMaximum();

    for (const auto &cut : validation.cut_times)
    {
        std::shared_ptr<TGraph> gr(new TGraph());
        gROOTOBJECTS.push_back(gr);
        gr->SetName(Form("EventsPerTimeCutGraph_%f_%f", cut.first, cut.second));
        gr->AddPoint(cut.first, 0.0);
        gr->AddPoint(cut.first, ymax);
        gr->AddPoint(cut.second, ymax);
        gr->AddPoint(cut.second, 0.0);
        gr->SetFillColor(kRed);
        gr->SetFillStyle(3013);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(2);
        gr->SetBit(TObject::kCannotPick);

        gr->Draw("same F");
    }

    const double xmin = events_per_time->GetXaxis()->GetXmin();
    const double xmax = events_per_time->GetXaxis()->GetXmax();
    std::shared_ptr<TLine> average_line(new TLine(xmin, validation.average_integral, xmax, validation.average_integral));
    gROOTOBJECTS.push_back(average_line);
    average_line->SetLineColor(kGreen + 2);
    average_line->SetLineWidth(2);
    average_line->Draw("same");

    std::shared_ptr<TLine> threshold_line(new TLine(xmin, validation.threshold_integral, xmax, validation.threshold_integral));
    gROOTOBJECTS.push_back(threshold_line);
    threshold_line->SetLineColor(kRed + 1);
    threshold_line->SetLineStyle(2);
    threshold_line->SetLineWidth(2);
    threshold_line->Draw("same");

    std::shared_ptr<TLegend> legend(new TLegend(0.65, 0.72, 0.92, 0.90));
    gROOTOBJECTS.push_back(legend);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(events_per_time.get(), "Events per time bin", "l");
    legend->AddEntry(average_line.get(), Form("Average = %.1f", validation.average_integral), "l");
    legend->AddEntry(threshold_line.get(), Form("Threshold = %.1f", validation.threshold_integral), "l");
    legend->Draw("same");

    c->Update();

    if (!save_filename.empty())
    {
        if (!can_create_file(save_filename))
        {
            std::cerr << "Cannot create output file: " << save_filename << std::endl;
            return;
        }
        TFile out_file(save_filename.c_str(), "recreate");
        if (!out_file.IsOpen())
        {
            std::cerr << "Error opening output file: " << save_filename << std::endl;
            return;
        }
        out_file.cd();
        c->Write();
        temat->Write();
        events_per_time->Write();
        std::cout << "Saved validation check canvas to: " << save_filename << std::endl;
    }
}

double calculate_total_duration(const std::vector<std::pair<double, double>> &cut_times)
{
    double total_duration = 0.0;
    for (const auto &cut : cut_times) { total_duration += (cut.second - cut.first); }
    return total_duration;
}

double calculate_run_duration(const std::shared_ptr<TH2> temat)
{
    const int nbinsx = temat->GetXaxis()->GetNbins();
    if (nbinsx <= 1) return 0.0;

    double first_time = temat->GetXaxis()->GetBinLowEdge(1);
    double last_time  = temat->GetXaxis()->GetBinUpEdge(nbinsx);

    return last_time - first_time;
}

double calculate_cut_padding(const std::vector<std::pair<double, double>> &cut_times)
{
    if (cut_times.size() < 2) return 0.0;
    double padding_front = 0.0;
    double padding_back  = 0.0;
    padding_front        = (cut_times.front().second - cut_times.front().first);
    padding_back         = (cut_times.back().second - cut_times.back().first);

    // There are at least 2 holes due to the padding at the beginning and at the end, however this calculation could be wrong if validation hole
    // overlaps with the start or end of the run - let's assume that at least one is ok
    return std::min(padding_front, padding_back);
}

/// @brief check input matrix for missing validation windows. Missing validation window is identified if number of events contained in a y-projection
/// slice is less then median multiplied by average_threshold parameter
/// @param temat
/// @param average_threshold
/// @return
ValidationCheckResult getMissingValidation(const std::shared_ptr<TH2> temat, const double average_threshold = 0.7)
{
    assert(average_threshold > 0.0 && average_threshold <= 2.0);

    const int nbinsx = temat->GetXaxis()->GetNbins();
    const int nbinsy = temat->GetYaxis()->GetNbins();

    std::vector<double> integrals;
    std::vector<double> times;
    integrals.reserve(static_cast<size_t>(nbinsx));
    times.reserve(static_cast<size_t>(nbinsx));

    for (int binx = 1; binx < nbinsx; binx++)
    {
        times.emplace_back(temat->GetXaxis()->GetBinCenter(binx));
        integrals.emplace_back(temat->Integral(binx, binx, 1, nbinsy));
    }

    ValidationCheckResult result;
    {
        auto tmp = integrals;
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::remove(tmp.begin(), tmp.end(), 0), tmp.end());

        // average
        result.average_integral = std::accumulate(tmp.begin(), tmp.end(), 0.0) / static_cast<double>(tmp.size());
        // double median_integral = tmp[tmp.size() / 2];
        result.threshold_integral = result.average_integral * average_threshold;
    }

    TGraph gr(static_cast<int>(integrals.size()), times.data(), integrals.data());
    gr.SetName("IntegralGraph");
    gr.SetTitle("Integral of time slices");
    gr.GetXaxis()->SetTitle("Time [s]");
    gr.GetYaxis()->SetTitle("Integral [counts]");

    for (size_t i = 0; i < integrals.size(); i++)
    {
        if (integrals[i] < result.threshold_integral)
        {
            int time_bin = temat->GetXaxis()->FindBin(times[i]);
            if (result.cut_times.size() == 0)
            {
                result.cut_times.emplace_back(temat->GetXaxis()->GetBinLowEdge(time_bin), temat->GetXaxis()->GetBinUpEdge(time_bin));
                continue;
            }
            if (result.cut_times.back().second == temat->GetXaxis()->GetBinLowEdge(time_bin))
            {
                // extend the last interval
                result.cut_times.back().second = temat->GetXaxis()->GetBinUpEdge(time_bin);
            }
            else
            {
                // add a new interval
                result.cut_times.emplace_back(temat->GetXaxis()->GetBinLowEdge(time_bin), temat->GetXaxis()->GetBinUpEdge(time_bin));
            }
        }
    }

    return result;
}

void print_help()
{
    std::cout << "Usage: getMissingValidation [options]\n"
              << "Options:\n"
              << "  --run <run_number>       Specify run number(s) (can be repeated)\n"
              << "  --crys <crystal>         Specify crystal(s) (can be repeated)\n"
              << "  --allcrys                Use all crystals\n"
              << "  --dir <directory>        Specify data directory\n"
              << "  --threshold <value>      Valid bin content threshold calculated using threshold factor multiplied average integral for given "
                 "detector (default: 0.7)\n"
              << "  --nodraw                 Do not draw missing-validation cut windows\n"
              << "  --save [filename]        Save each validation canvas to a ROOT file. If filename is omitted, use "
                 "vCheck_runXXXX_crysYYY.root\n"
              << "  --zoom <emin> <emax>     Zoom the matrix Y axis to the specified energy range\n"
              << "  --printdead              If validation holes are found, print their duration in seconds starting from first non-zero bin\n"
              << "  --modify-conf            If holes are found, put zero gain in overlapping "
                 "TimeEvoCC.conf bins\n"
              << "  --help, -h               Show this help message\n";

    exit(0);
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
                while (i + 1 < argc && std::isdigit(argv[i + 1][0])) { gRUNLIST.emplace_back(std::stoi(argv[++i])); }
                if (gRUNLIST.empty()) { throw std::invalid_argument("--run must be followed by at least one integer value"); }
            }
            else
            {
                throw std::invalid_argument("Missing value for --run");
            }
        }
        else if (arg == "--crys" || arg == "--crystal" || arg == "--crystals")
        {
            auto _crys = parse_space_separated_crystals(i, argc, argv);
            // avoid duplicates
            for (const auto &cry : _crys)
            {
                if (std::find(gCRYSTALLIST.begin(), gCRYSTALLIST.end(), cry) == gCRYSTALLIST.end()) { gCRYSTALLIST.emplace_back(cry); }
            }
        }
        else if (arg == "--allcrys")
        {
            std::vector<std::string> _c = {"00A", "00B", "00C", "01A", "01C", "02A", "02B", "02C", "04A", "04B", "04C",
                                           "05B", "05C", "06A", "06B", "06C", "07A", "07B", "08A", "08B", "09A", "09B",
                                           "09C", "10A", "10B", "10C", "11A", "11B", "11C", "14A", "14B", "14C"};
            for (const auto &cry : _c)
            {
                if (std::find(gCRYSTALLIST.begin(), gCRYSTALLIST.end(), cry) == gCRYSTALLIST.end()) { gCRYSTALLIST.emplace_back(cry); }
            }
        }
        else if (arg == "--dir")
        {
            if (i + 1 < argc) { gDIR = argv[++i]; }
            else
            {
                throw std::invalid_argument("Missing value for --dir");
            }
        }
        else if (arg == "--threshold")
        {
            if (i + 1 < argc)
            {
                gTHRESHOLD = std::stod(argv[++i]);
                if (!(gTHRESHOLD > 0.0 && gTHRESHOLD <= 2.0)) { throw std::invalid_argument("--threshold must be in the open interval (0, 2)"); }
            }
            else
            {
                throw std::invalid_argument("Missing value for --threshold");
            }
        }
        else if (arg == "--nodraw") { gDRAW_CUTS = false; }
        else if (arg == "--save" || arg == "--save-result")
        {
            gSAVE_RESULT = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') { gSAVE_FILENAME = argv[++i]; }
        }
        else if (arg == "--zoom")
        {
            if (i + 2 < argc)
            {
                gMATRIX_ZOOM_LOW  = std::stod(argv[++i]);
                gMATRIX_ZOOM_HIGH = std::stod(argv[++i]);
                if (!(gMATRIX_ZOOM_LOW < gMATRIX_ZOOM_HIGH)) { throw std::invalid_argument("--zoom requires emin < emax"); }
                gUSE_MATRIX_ZOOM = true;
            }
            else
            {
                throw std::invalid_argument("--zoom must be followed by emin and emax");
            }
        }
        else if (arg == "--modify-conf") { gMODIFY_TIMEEVO_CONF = true; }
        else if (arg == "--printdead") { gPRINTDEAD = true; }

        else
        {
            print_help();
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }

    if (gRUNLIST.empty()) { throw std::invalid_argument("--run must be provided with at least one run number"); }
    if (gCRYSTALLIST.empty()) { throw std::invalid_argument("--crystal must be provided with at least one crystal"); }
    if (gSAVE_RESULT && !gSAVE_FILENAME.empty() && gRUNLIST.size() * gCRYSTALLIST.size() > 1)
    {
        throw std::invalid_argument("--save filename can only be used with one run/crystal pair. Omit filename to use default per-crystal names.");
    }

    std::cout << "Parameters used are:" << std::endl;
    std::cout << "Run number(s):    ";
    for (const auto &run : gRUNLIST) { std::cout << run << " "; }
    std::cout << std::endl;
    std::cout << "Crystals:         ";
    for (const auto &cry : gCRYSTALLIST) { std::cout << cry << " "; }
    std::cout << std::endl;
    std::cout << "Data directory:        " << gDIR << std::endl;
    std::cout << "Threshold:             " << gTHRESHOLD << std::endl;
    std::cout << "Draw cuts:             " << std::boolalpha << gDRAW_CUTS << std::endl;
    std::cout << "Save result:           " << std::boolalpha << gSAVE_RESULT << std::endl;
    if (gSAVE_RESULT && !gSAVE_FILENAME.empty()) { std::cout << "Save filename:         " << gSAVE_FILENAME << std::endl; }
    if (gUSE_MATRIX_ZOOM) { std::cout << "Matrix Y zoom:         " << gMATRIX_ZOOM_LOW << " " << gMATRIX_ZOOM_HIGH << std::endl; }
    std::cout << "Modify TimeEvo conf:   " << std::boolalpha << gMODIFY_TIMEEVO_CONF << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
}

/// @brief Construct merged intervals of missed validation windows across crystals in the same run (specified in a single input vector (detector)
/// vector of intervals) and apply them to TimeEvoCC.conf if requested. Intervals are merged if they overlap or touch each other (i.e. end of one
/// interval is the same as the start of the other).
/// @param input
/// @return
std::vector<std::pair<double, double>> mergeIntervals(std::vector<std::vector<std::pair<double, double>>> input)
{
    std::vector<std::pair<double, double>> return_intervals;
    if (input.size() == 0) return return_intervals;
    // sort by initial time

    return_intervals.reserve(input.size() * input[0].size()); // reserve roughly enough space to avoid multiple allocations
    // merge all intervals into 1 vector
    for (size_t i = 0; i < input.size(); i++) { return_intervals.insert(return_intervals.end(), input[i].begin(), input[i].end()); }
    // sort them by start time
    std::sort(return_intervals.begin(), return_intervals.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    for (uint i = 1; i < return_intervals.size(); i++)
    {
        auto &previous = return_intervals[i - 1];
        auto &current  = return_intervals[i];
        // check overlap with previous interval
        if (current.first <= previous.second)
        {
            // merge intervals
            previous.second = std::max(previous.second, current.second);
            // remove current interval
            return_intervals.erase(return_intervals.begin() + i);
            i--; // stay at the same index to check with the next interval
        }
    }
    return return_intervals;
}

int main(int argc, char **argv)
{

    parse_args(argc, argv);
    // TRint app("app", &argc, argv);
    TApplication app("app", 0, 0);

    Long64_t ts_offset = std::numeric_limits<Long64_t>::max();

    for (const auto &run : gRUNLIST)
    {
        std::vector<std::vector<std::pair<double, double>>> cut_times_in_run;
        for (const auto &crystal : gCRYSTALLIST)
        {
            std::string rootfilename = get_rootfilename(gDIR, run, crystal);
            TFile      *matfile      = TFile::Open(rootfilename.c_str(), "READ");
            if (!matfile || matfile->IsZombie())
            {
                std::cerr << "Error! could not open/find the " << rootfilename << " file" << std::endl;
                continue;
            }
            std::string          matrixname = get_matrixname(crystal);
            std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(matrixname.c_str()));
            if (!TEMAT_original)
            {
                std::cerr << "Error! could not open/find the " << matrixname << " matrix" << std::endl;
                continue;
            }
            auto validation = getMissingValidation(TEMAT_original, gTHRESHOLD);

            if (validation.cut_times.size() > gMINIMUM_NUMBER_OF_HOLES_TO_IGNORE)
            {
                cut_times_in_run.emplace_back(validation.cut_times);

                double missing_duration   = calculate_total_duration(validation.cut_times);
                double total_run_duration = calculate_run_duration(TEMAT_original);
                double padding_duration   = calculate_cut_padding(validation.cut_times);
                missing_duration -= padding_duration;
                total_run_duration -= padding_duration;

                // to seconds
                missing_duration *= 60.0;
                total_run_duration *= 60.0;

                double percentage = (total_run_duration > 0) ? (missing_duration / total_run_duration) * 100.0 : 0.0;

                std::cout << "Lost due to validation: run " << run << " cry " << crystal << " " << std::fixed << std::setprecision(2)
                          << missing_duration << "s out of total run duration " << total_run_duration << "s " << percentage << "%" << std::endl;
                std::cout << "Found validation " << validation.cut_times.size() << " holes for run " << run << " crystal " << crystal << ": "
                          << std::endl;
            }
            if (gDRAW_CUTS || gSAVE_RESULT)
            {
                draw_cuts(TEMAT_original, validation, gSAVE_RESULT ? get_validation_output_filename(run, crystal) : "");
            }
            ts_offset = std::min(ts_offset, get_first_timestamp(TEMAT_original));
        }
        auto merged_cut_times = mergeIntervals(cut_times_in_run);

        std::cout << "After merging holes, run " << run << " has " << merged_cut_times.size() << " cut intervals " << std::endl;

        if (gMODIFY_TIMEEVO_CONF && !merged_cut_times.empty())
        {
            std::cout << "Applying " << merged_cut_times.size() << " merged cut interval(s) to TimeEvoCC.conf for run " << run << std::endl;
            for (const auto &crystal : gCRYSTALLIST)
            {
                const std::string conf_filename = get_conffilename(gDIR, run, crystal);
                if (modify_conffile(conf_filename, merged_cut_times)) { std::cout << "Updated config: " << conf_filename << std::endl; }
                else
                {
                    std::cout << "No overlapping bins in config: " << conf_filename << std::endl;
                }
            }
        }
        if (gPRINTDEAD && !merged_cut_times.empty())
        {
            std::cout << "\n\n********** Printing list of dead intervals in seconds since start of the run **********\n";
            std::cout << "# Dead intervals run " << run << " by merging crystals ";
            for (const auto &cry : gCRYSTALLIST) { std::cout << cry << " "; }
            std::cout << "\n";
            for (const auto [ts_start, ts_stop] : merged_cut_times)
            {
                double sec_start = (ts_start - static_cast<double>(ts_offset)) * 60.;
                double sec_stop  = (ts_stop - static_cast<double>(ts_offset)) * 60.;
                if (sec_start < 0) sec_start = 0;
                if (sec_stop < 0) continue;
                std::cout << "DEAD_INTERVAL                   " << std::setprecision(1) << std::setw(10) << sec_start << " " << std::setw(10)
                          << sec_stop << "\n";
            }
        }
        std::cout << "-------------------------DONE----------------------------------" << std::endl;
    }
    if (gDRAW_CUTS && !gROOTOBJECTS.empty()) { app.Run(); }
    return 0;
}
