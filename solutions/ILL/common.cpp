#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

bool can_create_file(const std::string &path)
{
    try
    {
        std::filesystem::path file_path(path);
        std::filesystem::path dir_path = file_path.parent_path();

        if (std::filesystem::exists(file_path)) { return true; }

        if (!dir_path.empty() && !std::filesystem::exists(dir_path))
        {
            if (!std::filesystem::create_directories(dir_path))
            {
                std::cerr << "Error: Unable to create directory: " << dir_path << std::endl;
                return false;
            }
        }

        std::ofstream file(path);
        if (!file.is_open())
        {
            std::cerr << "Error: Unable to create file at " << path << std::endl;
            return false;
        }

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
        else
        {
            throw std::runtime_error("Missing float value for parameter");
        }
    }
    return result;
}

void parse_ROI_source(char *argv, std::vector<float> &ROI, std::vector<float> &fit_peak)
{
    std::string source = argv;
    std::transform(source.begin(), source.end(), source.begin(), [](unsigned char c) { return std::tolower(c); });
    if (source.compare("60co") == 0 || source.compare("co60") == 0)
    {
        ROI      = {1332.492, 1300., 1370., -20, 20};
        fit_peak = {1173.228, 1165., 1185.};
    }
    else if (source.compare("133ba") == 0 || source.compare("ba133") == 0)
    {
        ROI      = {383.95, 365, 400., -15, 15};
        fit_peak = {2614, 2550., 2650};
    }
    else if (source.compare("152eu") == 0 || source.compare("eu152") == 0)
    {
        ROI      = {1408.013, 1378., 1438., -40, 40};
        fit_peak = {1528.10, 1498, 1558};
    }
    else if (source.compare("226ra") == 0 || source.compare("ra226") == 0)
    {
        ROI      = {2204.1, 2150., 2250., -50, 50};
        fit_peak = {2447.69, 2400, 2500};
    }
    else if (source.compare("66ga") == 0 || source.compare("ga66") == 0)
    {
        ROI      = {2751.835, 2720., 2780., -90, 50};
        fit_peak = {4295.187, 4220., 4360.};
    }
    else if (source.compare("56co") == 0 || source.compare("co56") == 0) { throw std::runtime_error("56Co source is not implemented yet"); }
    else if (source.compare("22na") == 0 || source.compare("na22") == 0) { throw std::runtime_error("Na-22 source is not implemented yet"); }
    else if (source.compare("cs137") == 0 || source.compare("137cs") == 0) { throw std::runtime_error("Cs-137 source is not implemented yet"); }
    else
    {
        throw std::runtime_error("Unknown source: " + source);
    }
}
