#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.cpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>

std::vector<int>         gRUNLIST;
std::vector<std::string> gCRYSTALLIST;
// std::string gCONFDIR{""};
std::string gDIR{""};
const int   gMINIMUM_NUMBER_OF_HOLES_TO_IGNORE = 2;
bool        gSUMMARIZE_DURATION                = false;

std::vector<std::shared_ptr<TObject>> gROOTOBJECTS;

const double MINUTE_TO_TIMESTAMPS = 6.e9;

void bckp_conf_file(const std::string &infile, std::string &outfile)
{
    if (infile.empty() || outfile.empty())
    {
        std::cerr << "Error: infile or outfile is empty!" << std::endl;
        return;
    }

    std::ifstream src(infile);
    std::ofstream dst(outfile);
    if (!src.is_open())
    {
        std::cerr << "Error opening source file: " << infile << std::endl;
        return;
    }
    if (!dst.is_open())
    {
        std::cerr << "Error opening destination file: " << outfile << std::endl;
        return;
    }

    dst << src.rdbuf();
    src.close();
    dst.close();
}

void modify_conffile(const std::string                            &filename,
                     const std::vector<std::pair<double, double>> &to_be_cut)
{
    (void)to_be_cut; // currently not used, but we might want to add the cut times to the
                     // conf file in the future
    std::ifstream conf_file(filename);
    if (!conf_file.is_open())
    {
        std::cerr << "error opening conf file " << filename << std::endl;
        exit(1);
    }

    std::string              line;
    std::vector<std::string> lines;
    while (std::getline(conf_file, line)) { lines.push_back(line); }
}

void draw_cuts(std::shared_ptr<TH2>                          temat,
               const std::vector<std::pair<double, double>> &to_be_cut)
{
    if (to_be_cut.empty()) { return; }

    gROOTOBJECTS.push_back(temat);

    std::shared_ptr<TCanvas> c(new TCanvas("Cuts", "Cuts", 1400, 1000));
    gROOTOBJECTS.push_back(c);
    c->SetGrid();
    c->SetCrosshair(1);
    temat->Draw("colz");

    for (const auto &cut : to_be_cut)
    {
        std::shared_ptr<TGraph> gr(new TGraph());
        gROOTOBJECTS.push_back(gr);
        gr->SetName(Form("CutGraph_%f_%f", cut.first, cut.second));
        auto ymin = temat->GetYaxis()->GetXmin();
        auto ymax = temat->GetYaxis()->GetXmax();
        gr->AddPoint(cut.first, ymin);
        gr->AddPoint(cut.first, ymax);
        gr->AddPoint(cut.second, ymax);
        gr->AddPoint(cut.second, ymin);
        gr->SetFillColor(kRed);
        gr->SetFillStyle(3013);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(2);

        gr->Draw("same F");
    }

    c->Update();
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
    double padding = 0.0;
    padding += (cut_times.front().second - cut_times.front().first);
    padding += (cut_times.back().second - cut_times.back().first);

    return padding;
}

std::vector<std::pair<double, double>> checkValidation(
    const std::shared_ptr<TH2> temat, const double integral_threshold = 0.5)
{
    assert(integral_threshold > 0.0 && integral_threshold < 1.);

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

    double avg_integral = std::accumulate(integrals.begin(), integrals.end(), 0.0) /
                          static_cast<double>(integrals.size());

    TGraph gr(static_cast<int>(integrals.size()), times.data(), integrals.data());
    gr.SetName("IntegralGraph");
    gr.SetTitle("Integral of time slices");
    gr.GetXaxis()->SetTitle("Time [s]");
    gr.GetYaxis()->SetTitle("Integral [counts]");

    std::vector<std::pair<double, double>> to_be_cut;
    for (size_t i = 0; i < integrals.size(); i++)
    {
        if (integrals[i] < avg_integral * integral_threshold)
        {
            int time_bin = temat->GetXaxis()->FindBin(times[i]);
            if (to_be_cut.size() == 0)
            {
                to_be_cut.emplace_back(temat->GetXaxis()->GetBinLowEdge(time_bin),
                                       temat->GetXaxis()->GetBinUpEdge(time_bin));
                continue;
            }
            if (to_be_cut.back().second == temat->GetXaxis()->GetBinLowEdge(time_bin))
            {
                // extend the last interval
                to_be_cut.back().second = temat->GetXaxis()->GetBinUpEdge(time_bin);
            }
            else
            {
                // add a new interval
                to_be_cut.emplace_back(temat->GetXaxis()->GetBinLowEdge(time_bin),
                                       temat->GetXaxis()->GetBinUpEdge(time_bin));
            }
        }
    }

    return to_be_cut;
}

void print_help()
{
    std::cout << "Usage: checkValidation [options]\n"
              << "Options:\n"
              << "  --run <run_number>       Specify run number(s) (can be repeated)\n"
              << "  --crys <crystal>         Specify crystal(s) (can be repeated)\n"
              << "  --allcrys                Use all crystals\n"
              << "  --dir <directory>        Specify data directory\n"
              << "  --summarize              Show only total duration of missing "
                 "validation windows\n"
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
            else
            {
                throw std::invalid_argument("Missing value for --dir");
            }
        }
        else if (arg == "--summarize") { gSUMMARIZE_DURATION = true; }

        else
        {
            print_help();
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }

    if (gRUNLIST.empty())
    {
        throw std::invalid_argument(
            "--run must be provided with at least one run number");
    }
    if (gCRYSTALLIST.empty())
    {
        throw std::invalid_argument(
            "--crystal must be provided with at least one crystal");
    }

    std::cout << "Parameters used are:" << std::endl;
    std::cout << "Run number(s):    ";
    for (const auto &run : gRUNLIST) { std::cout << run << " "; }
    std::cout << std::endl;
    std::cout << "Crystals:         ";
    for (const auto &cry : gCRYSTALLIST) { std::cout << cry << " "; }
    std::cout << std::endl;
    std::cout << "Data directory:        " << gDIR << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
}

int main(int argc, char **argv)
{

    parse_args(argc, argv);
    TApplication app("app", 0, 0);

    for (const auto &run : gRUNLIST)
    {
        for (const auto &crystal : gCRYSTALLIST)
        {
            std::string rootfilename = get_rootfilename(gDIR, run, crystal);
            TFile      *matfile      = TFile::Open(rootfilename.c_str(), "READ");
            if (!matfile || matfile->IsZombie())
            {
                std::cerr << "Error! could not open/find the " << rootfilename << " file"
                          << std::endl;
                continue;
            }
            std::string          matrixname = get_matrixname(crystal);
            std::shared_ptr<TH2> TEMAT_original((TH2 *)matfile->Get(matrixname.c_str()));
            if (!TEMAT_original)
            {
                std::cerr << "Error! could not open/find the " << matrixname << " matrix"
                          << std::endl;
                continue;
            }
            auto cut_times = checkValidation(TEMAT_original, 0.5);
            if (cut_times.size() > gMINIMUM_NUMBER_OF_HOLES_TO_IGNORE)
            {
                if (gSUMMARIZE_DURATION)
                {
                    double missing_duration   = calculate_total_duration(cut_times);
                    double total_run_duration = calculate_run_duration(TEMAT_original);
                    double padding_duration   = calculate_cut_padding(cut_times);
                    missing_duration -= padding_duration;
                    total_run_duration -= padding_duration;

                    // to seconds
                    missing_duration *= 60.0;
                    total_run_duration *= 60.0;

                    double percentage =
                        (total_run_duration > 0)
                            ? (missing_duration / total_run_duration) * 100.0
                            : 0.0;

                    std::cout << "Run " << run << " cry " << crystal << " " << std::fixed
                              << std::setprecision(2) << missing_duration << " s "
                              << total_run_duration << " s " << percentage << "%"
                              << std::endl;
                }
                else
                {
                    std::cout << "Found validation " << cut_times.size()
                              << " holes for run " << run << " crystal " << crystal
                              << ": " << std::endl;
                    draw_cuts(TEMAT_original, cut_times);
                }
            }
        }
    }
    if (!gROOTOBJECTS.empty() && !gSUMMARIZE_DURATION) { app.Run(); }
    return 0;
}