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

std::vector<double> gENERGY_BINING{};

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

int CoresTimeEvo(int                 runNr,
                 std::vector<string> crystals,
                 int                 seconds_per_bin,
                 ULong64_t           maxEntries = 0,
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

    ULong64_t TotalNumberOfEntries = tree->GetEntries();
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

    // 1 bin = 5 mins
    Long64_t minTime   = Long64_t(minTS * 1e-8) / 60. - 10;
    Long64_t maxTime   = Long64_t(maxTS * 1e-8) / 60. + 10;
    int      nTimeBins = (maxTime - minTime) * 60 / seconds_per_bin; // 30 seconds per bin

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

    for (ULong64_t entry = 0; entry < TotalNumberOfEntries; entry++)
    {
        // Start timing the loop
        static auto start_time = std::chrono::steady_clock::now();

        // Process the entry
        tree->GetEntry(entry);
        if (nbcores == 0) continue;

        for (int n = 0; n < nbcores; n++)
        {
            auto it = std::find(crystalIds.begin(), crystalIds.end(), coreId[n]);
            if (it != crystalIds.end())
            {
                int index = std::distance(crystalIds.begin(), it);
                timeEvoMatrices[index]->Fill(coreTS[n] * 1.e-8 / 60.0, coreE0[n]);
            }
        }

        // Print progress every X entries
        if (entry % 1000000 == 0)
        {
            auto   current_time    = std::chrono::steady_clock::now();
            double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
                                         current_time - start_time)
                                         .count();
            double progress = static_cast<double>(entry) / TotalNumberOfEntries;
            double estimated_total_time = elapsed_seconds / progress;
            int    remaining_time       = estimated_total_time - elapsed_seconds;

            std::cout << "\rProcessed " << entry << " / " << TotalNumberOfEntries << " ("
                      << (progress * 100.0) << "%)"
                      << " | Elapsed: " << elapsed_seconds << "s"
                      << " | Remaining: " << remaining_time
                      << "s                                         " << std::flush;
        }
    }
    for (int i = 0; i < root_files.size(); i++)
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
            else { throw std::invalid_argument("Missing value for --binning"); }
        }
        else if (arg == "--Ebinning")
        {
            if (i + 3 < argc)
            {
                gENERGY_BINING.emplace_back(std::stoi(argv[++i]));
                gENERGY_BINING.emplace_back(std::stof(argv[++i]));
                gENERGY_BINING.emplace_back(std::stof(argv[++i]));
            }
            else { throw std::invalid_argument("Missing values for --Ebinning"); }
        }
        else if (arg == "--run")
        {
            if (i + 1 < argc) { run = std::stoi(argv[++i]); }
            else { throw std::invalid_argument("Missing value for --run"); }
        }
        else if (arg == "--maxentries")
        {
            if (i + 1 < argc) { maxEntries = std::stol(argv[++i]); }
            else { throw std::invalid_argument("Missing value for --maxentries"); }
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
            else { throw std::invalid_argument("Missing value for --outdir"); }
        }
        else if (arg == "--replaydir")
        {
            if (i + 1 < argc) { replayDir = argv[++i]; }
            else { throw std::invalid_argument("Missing value for --replaydir"); }
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
    std::vector<std::string> crystals;
    int                      number_of_seconds_per_bin = 30;
    std::string              outDir{};
    std::string              replayDir{};

    parseArguments(argc, argv, number_of_seconds_per_bin, run, maxentries, crystals,
                   outDir, replayDir);

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
    return CoresTimeEvo(run, crystals, number_of_seconds_per_bin, maxentries, outDir,
                        replayDir);
}