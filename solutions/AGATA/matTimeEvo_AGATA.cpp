#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common.cpp"

using namespace std;

struct user_config
{
    int                      run;
    Long64_t                 maxentries;
    std::vector<std::string> crystals;
    int                      number_of_seconds_per_bin;
    std::string              outDir;
    std::string              replayDir;
    std::vector<double>      energy_binning{};
};

std::vector<double> gENERGY_BINING{};

const int gTEMAT_PADDING_BINS = 20;

#include <filesystem>
#include <iostream>

void createDirectoryIfNotExists(const std::string &path)
{
    std::filesystem::path dirPath(path);

    // Check if the directory exists
    if (!std::filesystem::exists(dirPath))
    {
        // Create the directory
        if (!std::filesystem::create_directories(dirPath))
        {
            std::cerr << "Failed to create directory: " << path << std::endl;
        }
    }
    // else { std::cout << "Directory already exists: " << path << std::endl; }
}

int SegmentsTimeEvo(int              runNr,
                    std::string      crystal,
                    std::vector<int> segments,
                    int              seconds_per_bin,
                    Long64_t         maxEntries = 0,
                    std::string      outDir     = "",
                    std::string      replayDir  = "")
{

    string inFilePattern;
    if (replayDir.empty())
        inFilePattern = "run_" + fourCharInt(runNr) + "/Out/Analysis" + "/Tree_";
    else
        inFilePattern = replayDir + "/Tree_";
    // string outDirName    = "run_" + fourCharInt(runNr) + "/TimeEvo";
    string outDirName = "timeEvo";
    if (!outDir.empty()) outDirName = outDir;
    if (outDirName.back() == '/') outDirName.pop_back();

    createDirectoryIfNotExists(outDirName);

    int crystalId = get_crystal_id(crystal);
    std::cout << "Crystal: " << crystal << "\tID: " << crystalId << std::endl;

    std::cout << "\n";
    std::cout << "input file pattern: " << inFilePattern << std::endl;
    std::cout << "output directory: " << outDirName << std::endl;

    TChain *tree = new TChain(("TreeMaster"));
    tree->Add((inFilePattern + "*.root").c_str());

    Long64_t TotalNumberOfEntries = tree->GetEntries();
    if (maxEntries != 0) TotalNumberOfEntries = maxEntries;

    ULong64_t minTS = 0;
    ULong64_t maxTS = 0;

    // Get the first and last TS
    std::array<ULong64_t, 100> hitTS;
    std::array<Float_t, 100>   hitE;
    std::array<int, 100>       hitId;
    std::array<int, 100>       hitSg;
    int                        nbhits;

    tree->SetBranchAddress("TSHit", hitTS.data());
    tree->SetBranchAddress("hitE", hitE.data());
    tree->SetBranchAddress("hitId", hitId.data());
    tree->SetBranchAddress("hitSg", hitSg.data());
    tree->SetBranchAddress("nbHits", &nbhits);

    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("TSHit", true);
    tree->SetBranchStatus("hitE", true);
    tree->SetBranchStatus("hitId", true);
    tree->SetBranchStatus("hitSg", true);
    tree->SetBranchStatus("nbHits", true);

    // find min/max TS
    Long64_t index = 0;
    while (minTS == 0 && index < TotalNumberOfEntries)
    {
        tree->GetEntry(index);
        index++;
        if (nbhits == 0) continue;
        minTS = hitTS[0];
    }
    if (minTS == 0)
    {
        std::cout << "\nRun is empty\n" << std::endl;
        return 1;
    }

    index = TotalNumberOfEntries - 1;
    while (maxTS == 0 && index >= 0)
    {
        tree->GetEntry(index);
        index--;
        if (nbhits == 0) continue;
        maxTS = hitTS[0];
    }
    if (maxTS == 0)
    {
        std::cout << "\nRun is empty\n" << std::endl;
        return -1;
    }

    std::cout << "minTS: " << minTS << "\n";
    std::cout << "maxTS: " << maxTS << std::endl;

    // allocate matrices
    // TS is in 10 ns units; matrix time axis is in minutes.
    const double time_bin_width_minutes = static_cast<double>(seconds_per_bin) / 60.0;
    const double first_time             = static_cast<double>(minTS) * 1.e-8 / 60.0;
    const double last_time              = static_cast<double>(maxTS) * 1.e-8 / 60.0;
    const int    data_time_bins         = static_cast<int>((last_time - first_time) / time_bin_width_minutes) + 1;
    const double minTime                = first_time - static_cast<double>(gTEMAT_PADDING_BINS) * time_bin_width_minutes;
    const int    nTimeBins              = data_time_bins + 2 * gTEMAT_PADDING_BINS;
    const double maxTime                = minTime + static_cast<double>(nTimeBins) * time_bin_width_minutes;

    std::cout << "time binning: " << nTimeBins << "\trange: " << minTime << " " << maxTime
              << std::endl;

    std::vector<std::shared_ptr<TFile>> root_files;

    std::vector<std::shared_ptr<TH2F>> timeEvoMatrices;
    for (const auto &seg : segments)
    {
        string outFileName = outDirName + "/temat_" + fourCharInt(runNr) + "_" + crystal +
                             "_seg_" + twoCharInt(seg) + ".root";
        root_files.emplace_back(std::make_shared<TFile>(outFileName.c_str(), "recreate"));

        std::string mat_name = "hE0_TS_" + crystal + "_seg_" + twoCharInt(seg);
        timeEvoMatrices.emplace_back(std::make_shared<TH2F>(
            mat_name.c_str(), mat_name.c_str(), nTimeBins, minTime, maxTime,
            gENERGY_BINING.at(0), gENERGY_BINING.at(1), gENERGY_BINING.at(2)));

        timeEvoMatrices.back()->SetXTitle("Time [min]");
        timeEvoMatrices.back()->SetYTitle("Energy [keV]");
    }

    // sort data
    for (Long64_t entry = 0; entry < TotalNumberOfEntries; entry++)
    {
        // Start timing the loop
        static auto start_time = std::chrono::steady_clock::now();

        // Process the entry
        tree->GetEntry(entry);
        if (nbhits == 0) continue;

        for (uint n = 0; n < static_cast<uint>(nbhits); n++)
        {
            if (hitId[n] != crystalId) continue;
            auto it = std::find(segments.begin(), segments.end(), hitSg[n]);
            if (it != segments.end())
            {
                uint _index = static_cast<uint>(std::distance(segments.begin(), it));
                timeEvoMatrices[_index]->Fill(
                    static_cast<double>(hitTS[n]) * 1.e-8 / 60.0, hitE[n]);
            }
        }

        // Print progress every X entries
        if (entry % 1000000 == 0)
        {
            auto current_time    = std::chrono::steady_clock::now();
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
                                       current_time - start_time)
                                       .count();
            double progress =
                static_cast<double>(entry) / static_cast<double>(TotalNumberOfEntries);
            double estimated_total_time = static_cast<double>(elapsed_seconds) / progress;
            auto   remaining_time =
                static_cast<int64_t>(estimated_total_time) - elapsed_seconds;

            std::cout << "\rProcessed " << entry << " / " << TotalNumberOfEntries << " ("
                      << (progress * 100.0) << "%)"
                      << " | Elapsed: " << elapsed_seconds << "s"
                      << " | Remaining: " << remaining_time
                      << "s                                         " << std::flush;
        }
    }
    for (uint i = 0; i < root_files.size(); i++)
    {
        std::cout << "\nWriting matrix for crystal " << crystal << " segment "
                  << segments[i] << std::endl;

        auto file   = root_files[i];
        auto hE0_TS = timeEvoMatrices[i];
        file->cd();
        hE0_TS->Write();
    }
    return 0;
}
int CoresTimeEvo(int                 runNr,
                 std::vector<string> crystals,
                 int                 seconds_per_bin,
                 Long64_t            maxEntries = 0,
                 std::string         outDir     = "",
                 std::string         replayDir  = "")
{

    string inFilePattern;
    if (replayDir.empty())
        inFilePattern = "run_" + fourCharInt(runNr) + "/Out/Analysis" + "/Tree_";
    else
        inFilePattern = replayDir + "/Tree_";
    // string outDirName    = "run_" + fourCharInt(runNr) + "/TimeEvo";
    string outDirName = "timeEvo";
    if (!outDir.empty()) outDirName = outDir;
    if (outDirName.back() == '/') outDirName.pop_back();

    createDirectoryIfNotExists(outDirName);

    std::vector<Int_t> crystalIds;
    for (const auto &cry : crystals)
    {
        crystalIds.emplace_back(get_crystal_id(cry));
        std::cout << "Crystal: " << cry << "\tID: " << crystalIds.back() << std::endl;
    }

    std::cout << "\n";
    std::cout << "input file pattern: " << inFilePattern << std::endl;
    std::cout << "output directory: " << outDirName << std::endl;

    TChain *tree = new TChain(("TreeMaster"));
    tree->Add((inFilePattern + "*.root").c_str());

    Long64_t TotalNumberOfEntries = tree->GetEntries();
    if (maxEntries != 0) TotalNumberOfEntries = maxEntries;

    ULong64_t minTS = 0;
    ULong64_t maxTS = 0;

    // Get the first and last TS
    std::array<ULong64_t, 100> coreTS;
    std::array<Float_t, 100>   coreE0;
    std::array<int, 100>       coreId;
    int                        nbcores;

    tree->SetBranchAddress("coreTS", coreTS.data());
    tree->SetBranchAddress("coreE0", coreE0.data());
    tree->SetBranchAddress("coreId", coreId.data());
    tree->SetBranchAddress("nbCores", &nbcores);

    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("coreTS", true);
    tree->SetBranchStatus("coreE0", true);
    tree->SetBranchStatus("coreId", true);
    tree->SetBranchStatus("nbCores", true);

    Long64_t index = 0;
    while (minTS == 0 && index < TotalNumberOfEntries)
    {
        tree->GetEntry(index);
        index++;
        if (nbcores == 0) continue;
        minTS = coreTS[0];
    }
    if (minTS == 0)
    {
        std::cout << "\nRun is empty\n" << std::endl;
        return 1;
    }

    index = TotalNumberOfEntries - 1;
    while (maxTS == 0 && index >= 0)
    {
        tree->GetEntry(index);
        index--;
        if (nbcores == 0) continue;
        maxTS = coreTS[0];
    }
    if (maxTS == 0)
    {
        std::cout << "\nRun is empty\n" << std::endl;
        return -1;
    }

    std::cout << "minTS: " << minTS << "\n";
    std::cout << "maxTS: " << maxTS << std::endl;

    // TS is in 10 ns units; matrix time axis is in minutes.
    const double time_bin_width_minutes = static_cast<double>(seconds_per_bin) / 60.0;
    const double first_time             = static_cast<double>(minTS) * 1.e-8 / 60.0;
    const double last_time              = static_cast<double>(maxTS) * 1.e-8 / 60.0;
    const int    data_time_bins         = static_cast<int>((last_time - first_time) / time_bin_width_minutes) + 1;
    const double minTime                = first_time - static_cast<double>(gTEMAT_PADDING_BINS) * time_bin_width_minutes;
    const int    nTimeBins              = data_time_bins + 2 * gTEMAT_PADDING_BINS;
    const double maxTime                = minTime + static_cast<double>(nTimeBins) * time_bin_width_minutes;

    std::cout << "time binning: " << nTimeBins << "\trange: " << minTime << " " << maxTime
              << std::endl;

    std::vector<std::shared_ptr<TFile>> root_files;

    std::vector<std::shared_ptr<TH2F>> timeEvoMatrices;
    for (const auto &cry : crystals)
    {
        string outFileName =
            outDirName + "/temat_" + fourCharInt(runNr) + "_" + cry + ".root";
        root_files.emplace_back(std::make_shared<TFile>(outFileName.c_str(), "recreate"));

        timeEvoMatrices.emplace_back(std::make_shared<TH2F>(
            ("hE0_TS_" + cry).c_str(), ("hE0_TS_" + cry).c_str(), nTimeBins, minTime,
            maxTime, gENERGY_BINING.at(0), gENERGY_BINING.at(1), gENERGY_BINING.at(2)));

        timeEvoMatrices.back()->SetXTitle("Time [min]");
        timeEvoMatrices.back()->SetYTitle("Energy [keV]");
    }

    for (Long64_t entry = 0; entry < TotalNumberOfEntries; entry++)
    {
        // Start timing the loop
        static auto start_time = std::chrono::steady_clock::now();

        // Process the entry
        tree->GetEntry(entry);
        if (nbcores == 0) continue;

        for (uint n = 0; n < static_cast<uint>(nbcores); n++)
        {
            auto it = std::find(crystalIds.begin(), crystalIds.end(), coreId[n]);
            if (it != crystalIds.end())
            {
                uint _index = static_cast<uint>(std::distance(crystalIds.begin(), it));
                timeEvoMatrices[_index]->Fill(
                    static_cast<double>(coreTS[n]) * 1.e-8 / 60.0, coreE0[n]);
            }
        }

        // Print progress every X entries
        if (entry % 1000000 == 0)
        {
            auto current_time    = std::chrono::steady_clock::now();
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
                                       current_time - start_time)
                                       .count();
            double progress =
                static_cast<double>(entry) / static_cast<double>(TotalNumberOfEntries);
            double estimated_total_time = static_cast<double>(elapsed_seconds) / progress;
            auto   remaining_time =
                static_cast<int64_t>(estimated_total_time) - elapsed_seconds;

            std::cout << "\rProcessed " << entry << " / " << TotalNumberOfEntries << " ("
                      << (progress * 100.0) << "%)"
                      << " | Elapsed: " << elapsed_seconds << "s"
                      << " | Remaining: " << remaining_time
                      << "s                                         " << std::flush;
        }
    }
    for (uint i = 0; i < root_files.size(); i++)
    {
        std::cout << "\nWriting matrix for crystal " << crystals[i] << std::endl;

        auto file   = root_files[i];
        auto hE0_TS = timeEvoMatrices[i];
        file->cd();
        hE0_TS->Write();
    }
    return 0;
}

void printHelp()
{
    std::cout
        << "To use the code, you should be in the directory where you ran replays\n\n";
    std::cout << "Usage: program [OPTIONS]\n";
    std::cout << "Options:\n";
    std::cout << "  --help                    Display this help message\n";
    std::cout << "  --run <integer>           Specify the run number (required)\n";
    std::cout << "  --crys <3-letter strings> Specify crystals (can be multiple "
                 "3-character strings)\n";
    std::cout
        << "  --seg <integers>          Specify segments (can be multiple integers)\n";
    std::cout << "                            Only possible to use for 1 crystal!\n";
    std::cout
        << "  --maxentries <integer>    Set the maximum number of entries (optional)\n";
    std::cout << "  --allcrys                 Run for all crystals of EXP_035\n";
    std::cout
        << "  --Tbinning <integer>      Set number of seconds per bin (default 30)\n";
    std::cout << "  --Ebinning <1> <2> <3>    Set energy binning as: \n"
              << "                                <1> number of bins (default 32 000)\n"
              << "                                <2> min energy (default 0) \n"
              << "                                <3> max energy (default 8 000)\n";
    std::cout << "  --outdir <string>       Specify output directory (default: "
                 "TimeEvo/)\n";
    std::cout << "  --replaydir <string>    Specify replay directory that contain ROOT "
                 "trees, default is run_XXXX/Out/Analysis \n";

    std::cout << std::endl << std::endl;
}

void parseArguments(int                       argc,
                    char                    **argv,
                    int                      &binning,
                    int                      &run,
                    Long64_t                 &maxEntries,
                    std::vector<std::string> &crystals,
                    std::vector<int>         &segments,
                    std::string              &outDir,
                    std::string              &replayDir)
{
    replayDir = "";
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--help")
        {
            printHelp();
            exit(0);
        }
        else if (arg == "--Tbinning")
        {
            if (i + 1 < argc) { binning = std::stoi(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --binning");
            }
        }
        else if (arg == "--Ebinning")
        {
            if (i + 3 < argc)
            {
                gENERGY_BINING.emplace_back(std::stoi(argv[++i]));
                gENERGY_BINING.emplace_back(std::stof(argv[++i]));
                gENERGY_BINING.emplace_back(std::stof(argv[++i]));
            }
            else
            {
                throw std::invalid_argument("Missing values for --Ebinning");
            }
        }
        else if (arg == "--run")
        {
            if (i + 1 < argc) { run = std::stoi(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --run");
            }
        }
        else if (arg == "--maxentries")
        {
            if (i + 1 < argc) { maxEntries = std::stol(argv[++i]); }
            else
            {
                throw std::invalid_argument("Missing value for --maxentries");
            }
        }
        else if (arg == "--crys" || arg == "--crystal" || arg == "--crystals")
        {
            auto _crys = parse_space_separated_crystals(i, argc, argv);
            // avoid duplicates
            for (const auto &cry : _crys)
            {
                if (std::find(crystals.begin(), crystals.end(), cry) == crystals.end())
                {
                    crystals.emplace_back(cry);
                }
            }
        }
        else if (arg == "--seg" || arg == "--segment" || arg == "--segments")
        {
            auto _segs = parse_space_separated_ints(i, argc, argv);
            for (const auto &seg : _segs) { segments.emplace_back(seg); }
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
                if (std::find(crystals.begin(), crystals.end(), cry) == crystals.end())
                {
                    crystals.emplace_back(cry);
                }
            }
        }
        else if (arg == "--outdir")
        {
            if (i + 1 < argc) { outDir = argv[++i]; }
            else
            {
                throw std::invalid_argument("Missing value for --outdir");
            }
        }
        else if (arg == "--replaydir")
        {
            if (i + 1 < argc) { replayDir = argv[++i]; }
            else
            {
                throw std::invalid_argument("Missing value for --replaydir");
            }
            if (replayDir.back() == '/') { replayDir.pop_back(); }
        }
        else
        {
            printHelp();
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }
    if (gENERGY_BINING.size() == 0) { gENERGY_BINING = {32000, 0, 8000}; }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "No arguments provided." << std::endl << std::endl;
        printHelp();
        return 1;
    }

    int                      run;
    Long64_t                 maxentries = 0;
    std::vector<std::string> crystals{};
    std::vector<int>         segments{};
    int                      number_of_seconds_per_bin = 30;
    std::string              outDir{};
    std::string              replayDir{};
    int                      segmentId = -1; // -1 = core

    (void)segmentId; // currently not used, but can be used in the future to decide
                     // whether to run core or segment time evolution

    parseArguments(argc, argv, number_of_seconds_per_bin, run, maxentries, crystals,
                   segments, outDir, replayDir);

    std::cout << "Parameters used are:" << std::endl;
    std::cout << "Run number:       " << run << std::endl;
    std::cout << "Seconds per bin:  " << number_of_seconds_per_bin << std::endl;
    std::cout << "Max entries:      " << maxentries << std::endl;
    std::cout << "Output directory: " << outDir << std::endl;
    std::cout << "Replay directory: " << replayDir << std::endl;

    for (const auto &cry : crystals) { std::cout << "   crystal: " << cry << std::endl; }

    if (crystals.empty())
    {
        std::cerr << "No crystals specified. Use --crys option." << std::endl;
        return 1;
    }
    if (run == 0)
    {
        std::cerr << "No run number specified. Use --run option." << std::endl;
        return 1;
    }
    if (segments.size() == 0)
    {
        return CoresTimeEvo(run, crystals, number_of_seconds_per_bin, maxentries, outDir,
                            replayDir);
    }
    else
    {
        if (crystals.size() > 1)
        {
            std::cerr << "Multiple crystals specified. Segments can only be used for 1 "
                         "crystal. Please specify only 1 crystal when using --seg option."
                      << std::endl;
            return 1;
        }
        return SegmentsTimeEvo(run, crystals[0], segments, number_of_seconds_per_bin,
                               maxentries, outDir, replayDir);
    }
}
