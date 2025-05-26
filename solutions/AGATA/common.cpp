#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <atomic>
#include <csignal>

std::atomic<bool> gTerminate_program{false};

void signal_handler(int signal)
{
    if (signal == SIGINT)
    {
        std::cout << "\nCtrl+C detected. Terminating gracefully..." << std::endl;
        gTerminate_program = true;
    }
}

void setup_signal_handling() { std::signal(SIGINT, signal_handler); }

std::string fourCharInt(int I)
{
    std::stringstream ID;
    ID << std::setfill('0') << std::setw(4) << I;
    return ID.str();
}

bool can_create_file(const std::string &path)
{
    try
    {
        // Extract the directory from the path
        std::filesystem::path file_path(path);
        std::filesystem::path dir_path = file_path.parent_path();

        // Check if the directory exists
        if (!dir_path.empty() && !std::filesystem::exists(dir_path))
        {
            if (!std::filesystem::create_directories(dir_path))
            {
                std::cerr << "Error: Unable to create directory: " << dir_path
                          << std::endl;
                return false;
            }
        }

        // Try to create and write to the file
        std::ofstream file(path);
        if (!file.is_open())
        {
            std::cerr << "Error: Unable to create file at " << path << std::endl;
            return false;
        }

        // Close the file and delete it (cleanup)
        file.close();
        std::filesystem::remove(path);

        return true;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return false;
    }
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

std::string get_pointer_string(void *address)
{
    std::ostringstream oss;
    oss << address;
    return oss.str();
}

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

std::vector<std::string> parse_space_separated_crystals(int &i, int argc, char **argv)
{
    std::vector<std::string> crystals;
    while (i + 1 < argc && std::string(argv[i + 1]).size() == 3)
    {
        std::string cry = argv[i + 1];
        if (!(cry.size() == 3 && std::isdigit(cry[0]) && std::isdigit(cry[1]) &&
              (cry[2] == 'A' || cry[2] == 'B' || cry[2] == 'C')))
        {
            break;
        }
        if (std::find(crystals.begin(), crystals.end(), argv[i + 1]) == crystals.end())
        {
            crystals.emplace_back(argv[++i]);
        }
        else
        {
            std::cerr << "Crystal " << argv[i + 1] << " already specified. Skipping."
                      << std::endl;
            continue;
        }
    }
    return crystals;
}

void parse_ROI_source(char *argv, std::vector<float> &ROI, std::vector<float> &fit_peak)
{
    std::string source = argv;
    std::transform(source.begin(), source.end(), source.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (source.compare("60co") == 0 || source.compare("co60") == 0)
    {
        ROI      = {1332.492, 1300., 1370., -20, 20}; // Example values for 60Co
        fit_peak = {1173.228, 1165., 1185.};
    }
    else if (source.compare("133ba") == 0 || source.compare("ba133") == 0)
    {
        ROI      = {356.012, 340., 370., -10, 10};
        fit_peak = {302.85, 290., 310.};
    }
    else if (source.compare("152eu") == 0 || source.compare("eu152") == 0)
    {
        throw std::runtime_error("152Eu source is not implemented yet");
    }
    else if (source.compare("226ra") == 0 || source.compare("ra226") == 0)
    {
        // decay of 214Bi
        ROI      = {1764.491, 1720, 1780, -50, 50};
        fit_peak = {2204.1, 2150., 2250.};
    }
    else if (source.compare("66ga") == 0 || source.compare("ga66") == 0)
    {
        ROI      = {2751.835, 2720., 2780., -50, 50};
        fit_peak = {4295.187, 4220., 4360.};
    }
    else if (source.compare("56co") == 0 || source.compare("co56") == 0)
    {
        throw std::runtime_error("56Co source is not implemented yet");
    }
    else if (source.compare("22na") == 0 || source.compare("na22") == 0)
    {
        throw std::runtime_error("Na-22 source is not implemented yet");
    }
    else if (source.compare("cs137") == 0 || source.compare("137cs") == 0)
    {
        throw std::runtime_error("Cs-137 source is not implemented yet");
    }
    else { throw std::runtime_error("Unknown source: " + source); }
}

std::string get_conffilename(const std::string &dir,
                             const int          run,
                             const std::string &crystal)
{
    std::string conffile =
        dir + "/" + "run_" + fourCharInt(run) + "/Conf/" + crystal + "/TimeEvoCC.conf";
    if (!can_create_file(conffile))
        throw std::runtime_error("Problem with creating output configuration file\n");
    return conffile;
}

std::string get_rootfilename(const std::string &dir,
                             const int          run,
                             const std::string &crystal)
{
    return dir + "/temat_" + fourCharInt(run) + "_" + crystal + ".root";
}