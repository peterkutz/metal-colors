// This program computes F82-Tint parameters for real metals.
// For more details and usage instructions, see the README.txt file in the same directory as this file.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

struct Color
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Color() = default;
    Color(double arg_x, double arg_y, double arg_z)
    : x(arg_x), y(arg_y), z(arg_z) {}
};

Color operator*(Color a, Color b)
{
    return Color(a.x * b.x, a.y * b.y, a.z * b.z);
}

Color operator/(Color a, Color b)
{
    return Color(a.x / b.x, a.y / b.y, a.z / b.z);
}

Color operator+(Color a, Color b)
{
    return Color(a.x + b.x, a.y + b.y, a.z + b.z);
}

Color operator-(Color a, Color b)
{
    return Color(a.x - b.x, a.y - b.y, a.z - b.z);
}

Color operator*(Color a, double b)
{
    return Color(a.x * b, a.y * b, a.z * b);
}

Color operator/(Color a, double b)
{
    return Color(a.x / b, a.y / b, a.z / b);
}

std::ostream& operator<<(std::ostream& os, Color color)
{
  return os << color.x << ", " << color.y << ", " << color.z;
}

void set_printing_precision(int digits)
{
    std::cout << std::setprecision(digits);
}

template<typename T>
void print(const T& x)
{
    std::cout << x;
}

template<typename THead, typename... TTail>
void print(const THead& head, const TTail&... tail)
{
    print(head);
    print(tail...);
}

void print_line()
{
    std::cout << std::endl;
}

template<typename THead, typename... TTail>
void print_line(const THead& head, const TTail&... tail)
{
    print(head);
    print_line(tail...);
}

static constexpr bool verbose_printing = false;

template<typename... TArgs>
void print_verbose(const TArgs&... args)
{
    if (verbose_printing)
    {
        print(args...);
    }
}

template<typename... TArgs>
void print_line_verbose(const TArgs&... args)
{
    if (verbose_printing)
    {
        print_line(args...);
    }
}

void check_assumption(bool condition)
{
    if (!condition)
    {
        throw std::runtime_error("ASSERTION FAILED");
    }
}

template <typename T>
constexpr T square(const T& x)
{
    return x * x;
}

template <typename T>
constexpr T cube(const T& x)
{
    return x * x * x;
}

template <typename T>
constexpr T fourth_power(const T& x)
{
    return square(square(x));
}

template <typename T>
constexpr T fifth_power(const T& x)
{
    return fourth_power(x) * x;
}

template <typename T>
constexpr T sixth_power(const T& x)
{
    return cube(square(x));
}

template <typename T>
constexpr T lerp(const T& x, const T& y, double t)
{
    return x * (1 - t) + y * t;
}

// Complex unpolarized Fresnel reflectivity (for metals).
double metal_fresnel(double n, double k, double cos_theta)
{
    double theta = std::acos(cos_theta);
    double two_a_squared = std::sqrt(square(square(n) - square(k) - square(std::sin(theta))) + 4.0 * square(n) * square(k)) + (square(n) - square(k) - square(std::sin(theta)));
    double two_b_squared = std::sqrt(square(square(n) - square(k) - square(std::sin(theta))) + 4.0 * square(n) * square(k)) - (square(n) - square(k) - square(std::sin(theta)));
    double a = std::sqrt(two_a_squared / 2.0);
    double b = std::sqrt(two_b_squared / 2.0);
    double Fs_theta = (square(a) + square(b) - 2.0 * a * cos_theta + square(cos_theta))
                    / (square(a) + square(b) + 2.0 * a * cos_theta + square(cos_theta));
    double Fp_theta = Fs_theta * (square(a) + square(b) - 2.0 * a * std::sin(theta) * std::tan(theta) + square(std::sin(theta)) * square(std::tan(theta)))
                               / (square(a) + square(b) + 2.0 * a * std::sin(theta) * std::tan(theta) + square(std::sin(theta)) * square(std::tan(theta)));
    double F_theta = (Fs_theta + Fp_theta) / 2.0;
    return F_theta;
}

Color metal_schlick_with_F82_tint(Color F0, Color F82_tint, double cos_theta)
{
    check_assumption(cos_theta >= 0.0);
    constexpr double cos_theta_max = 1.0 / 7.0;
    Color r = F0; // switch to terminology from slides
    Color t = F82_tint; // custom terminology: "t" (for "tint") = h / (r + (1 - r) * (1 - cos_theta_max)^5)
    constexpr double one_minus_cos_theta_max_to_the_fifth = fifth_power(1.0 - cos_theta_max);
    constexpr double one_minus_cos_theta_max_to_the_sixth = sixth_power(1.0 - cos_theta_max);
    Color white_minus_r = Color(1.0, 1.0, 1.0) - r;
    Color b_numerator = (r + white_minus_r * one_minus_cos_theta_max_to_the_fifth) * (Color(1.0, 1.0, 1.0) - t); // note "* (1 - t)" instead of "- h"
    constexpr double b_denominator = cos_theta_max * one_minus_cos_theta_max_to_the_sixth;
    constexpr double b_denominator_reciprocal = 1.0 / b_denominator;
    Color b = b_numerator * b_denominator_reciprocal; // analogous to "a" in slides
    double one_minus_cos_theta = 1.0 - cos_theta;
    Color offset_from_r = (white_minus_r - b * cos_theta * one_minus_cos_theta) * fifth_power(one_minus_cos_theta);
    Color F_theta = r + offset_from_r;
    return F_theta;
}

Color convert_F82_to_F82_tint(Color F0, Color F82)
{
    constexpr double cos_theta_max = 1.0 / 7.0;
    Color r = F0;
    Color h = F82;
    Color white_minus_r = Color(1.0, 1.0, 1.0) - r;
    Color t = h / (r + white_minus_r * fifth_power(1.0 - cos_theta_max));
    return t;
}

double sum_components(Color c)
{
    return c.x + c.y + c.z;
}

Color convert_XYZ_to_linear_Rec709_RGB(Color xyz)
{
    return Color(3.24096994 * xyz.x + -1.53738318 * xyz.y + -0.49861076 * xyz.z,
                 -0.96924364 * xyz.x + 1.8759675 * xyz.y + 0.04155506 * xyz.z,
                 0.05563008 * xyz.x + -0.20397696 * xyz.y + 1.05697151 * xyz.z);
}

Color convert_XYZ_to_linear_p3d65_RGB(Color xyz)
{
    return Color(2.493497 * xyz.x + -0.931384 * xyz.y + -0.402711 * xyz.z,
                 -0.829489 * xyz.x + 1.762664 * xyz.y + 0.023625 * xyz.z,
                 0.035846 * xyz.x + -0.076172 * xyz.y + 0.956885 * xyz.z);
}

Color convert_XYZ_to_linear_adobe_RGB(Color xyz)
{
    return Color(2.041588 * xyz.x + -0.565007 * xyz.y + -0.344731 * xyz.z,
                 -0.969244 * xyz.x + 1.875968 * xyz.y + 0.041555 * xyz.z,
                 0.013444 * xyz.x + -0.118362 * xyz.y + 1.015175 * xyz.z);
}

Color convert_XYZ_to_linear_rec2020_RGB(Color xyz)
{
    return Color(1.71666343 * xyz.x + -0.35567332 * xyz.y + -0.25336809 * xyz.z,
                 -0.66667384 * xyz.x + 1.61645574 * xyz.y + 0.0157683 * xyz.z,
                 0.01764248 * xyz.x + -0.04277698 * xyz.y + 0.94224328 * xyz.z);
}

double convert_linear_Rec709_RGB_to_sRGB_for_one_channel(double u)
{
    if (u <= 0.0031308)
    {
        return 12.92 * u;
    }
    else
    {
        return 1.055 * std::pow(u, 1.0 / 2.4) - 0.055;
    }
}

Color convert_linear_Rec709_RGB_to_sRGB(Color xyz)
{
    return Color(convert_linear_Rec709_RGB_to_sRGB_for_one_channel(xyz.x),
                 convert_linear_Rec709_RGB_to_sRGB_for_one_channel(xyz.y),
                 convert_linear_Rec709_RGB_to_sRGB_for_one_channel(xyz.z));
}

Color clip_color(Color c)
{
    return Color(std::clamp(c.x, 0.0, 1.0),
                 std::clamp(c.y, 0.0, 1.0), 
                 std::clamp(c.z, 0.0, 1.0));
}

Color convert_to_8_bit(Color c)
{
    return Color(std::round(c.x * 255.0),
                 std::round(c.y * 255.0),
                 std::round(c.z * 255.0));
}

struct LoadedData
{
    std::vector<double> CMF_wavelength_nm;
    std::vector<Color> CMF_XYZ;

    std::vector<double> illuminant_wavelength_nm;
    std::vector<double> illuminant_value;

    std::vector<std::string> metal_symbol;
    std::vector<std::vector<double>> metal_n_wavelength_nm;
    std::vector<std::vector<double>> metal_n;
    std::vector<std::vector<double>> metal_k_wavelength_nm;
    std::vector<std::vector<double>> metal_k;
};

LoadedData gather_data()
{
    LoadedData loaded_data;

    set_printing_precision(9);

    auto is_blank = [](const std::string &s) {
        return s.empty() || std::all_of(s.begin(), s.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        });
    };
    
    namespace fs = std::filesystem;
    
    const fs::path CMF_directory{ "ColorMatchingFunctions" };
    bool first_and_only_CMF_file = true;
    for (const auto& directory_entry : fs::directory_iterator(CMF_directory))
    {
        check_assumption(directory_entry.is_regular_file());
        if (directory_entry.path().extension().string() != ".csv")
        {
            continue;
        }

        check_assumption(first_and_only_CMF_file);
        first_and_only_CMF_file = false;

        std::string data_file_path = directory_entry.path().string();
        std::string data_file_stem = directory_entry.path().filename().stem().string();
        print_line_verbose(data_file_stem);
        
        std::ifstream data_file{ data_file_path };
        check_assumption(static_cast<bool>(data_file));
        std::string line;
        while (std::getline(data_file, line))
        {
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }

            if (line.empty() || is_blank(line))
            {
                print_line_verbose("[EMPTY LINE]");
            }
            else
            {
                int start = 0;
                int end = line.find(',');
                double wavelength_nm = std::stod(line.substr(start, end));
                start = end + 1;
                end = line.find(',', start);
                double x = std::stod(line.substr(start, end));
                start = end + 1;
                end = line.find(',', start);
                double y = std::stod(line.substr(start, end));
                start = end + 1;
                double z = std::stod(line.substr(start));

                print_verbose("wavelength_nm: ", wavelength_nm);
                print_verbose(" x: ", x);
                print_verbose(" y: ", y);
                print_line_verbose(" z: ", z);

                loaded_data.CMF_wavelength_nm.emplace_back(wavelength_nm);
                loaded_data.CMF_XYZ.emplace_back(x, y, z);
            }
        }
    }

    const fs::path illuminant_directory{ "Illuminants" };
    bool first_and_only_illuminant_file = true;
    for (const auto& directory_entry : fs::directory_iterator(illuminant_directory))
    {
        check_assumption(directory_entry.is_regular_file());
        if (directory_entry.path().extension().string() != ".csv")
        {
            continue;
        }

        check_assumption(first_and_only_illuminant_file);
        first_and_only_illuminant_file = false;

        std::string data_file_path = directory_entry.path().string();
        std::string data_file_stem = directory_entry.path().filename().stem().string();
        print_line_verbose(data_file_stem);
        
        std::ifstream data_file{ data_file_path };
        check_assumption(static_cast<bool>(data_file));
        bool odd = false;
        std::string line;
        while (std::getline(data_file, line))
        {
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }

            if (line.empty() || is_blank(line))
            {
                print_line_verbose("[EMPTY LINE]");
                continue;
            }
            else
            {
                // Split on the first comma: "wavelength,value"
                auto comma = line.find(',');
                double wavelength_nm = std::stod(line.substr(0, comma));
                double value         = std::stod(line.substr(comma + 1));

                print_line_verbose("wavelength_nm: ", wavelength_nm);
                print_line_verbose("value: ", value);

                loaded_data.illuminant_wavelength_nm.emplace_back(wavelength_nm);
                loaded_data.illuminant_value.emplace_back(value);
            }
        }

        check_assumption(loaded_data.illuminant_wavelength_nm.size() == loaded_data.illuminant_value.size());
    }

    std::vector<fs::path> ior_file_paths;

    const fs::path ior_directory{ "ComplexRefractiveIndexes" };
    for (const auto& directory_entry : fs::directory_iterator(ior_directory))
    {
        check_assumption(directory_entry.is_regular_file());
        if (directory_entry.path().extension().string() != ".csv")
        {
            continue;
        }

        ior_file_paths.emplace_back(directory_entry.path());
    }

    std::sort(ior_file_paths.begin(), ior_file_paths.end(),
        [](const auto& lhs, const auto& rhs) {
            return lhs.filename().string() < rhs.filename().string();
        });

    for (const auto& ior_file_path : ior_file_paths)
    {
        std::string data_file_path = ior_file_path.string();
        std::string data_file_stem = ior_file_path.filename().stem().string();
        print_line_verbose(data_file_stem);

        loaded_data.metal_symbol.emplace_back(data_file_stem);
        loaded_data.metal_n_wavelength_nm.emplace_back();
        loaded_data.metal_n.emplace_back();
        loaded_data.metal_k_wavelength_nm.emplace_back();
        loaded_data.metal_k.emplace_back();
        
        std::ifstream data_file{ data_file_path };
        check_assumption(static_cast<bool>(data_file));
        bool collecting_n = false;
        bool collecting_k = false;
        std::string line;
        while (std::getline(data_file, line))
        {
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }

            if (line.empty() || is_blank(line))
            {
                print_line_verbose("[EMPTY LINE]");

                collecting_n = false;
                collecting_k = false;
            }
            else if (line == "wl,n")
            {
                print_line_verbose("[", line, "]" );

                collecting_n = true;
                collecting_k = false;
            }
            else if (line == "wl,k")
            {
                print_line_verbose("[", line, "]" );

                collecting_n = false;
                collecting_k = true;
            }
            else if (collecting_n)
            {
                auto comma = line.find(',');
                double wavelength_nm = std::stod(line.substr(0, comma)) * 1000.0;  // convert from micrometers to nanometers
                double n             = std::stod(line.substr(comma + 1));

                print_verbose("wavelength_nm: ", wavelength_nm);
                print_line_verbose(" n: ", n);

                loaded_data.metal_n_wavelength_nm.back().emplace_back(wavelength_nm);
                loaded_data.metal_n.back().emplace_back(n);
            }
            else if (collecting_k)
            {
                auto comma = line.find(',');
                double wavelength_nm = std::stod(line.substr(0, comma)) * 1000.0;  // convert from micrometers to nanometers
                double k             = std::stod(line.substr(comma + 1));

                print_verbose("wavelength_nm: ", wavelength_nm);
                print_line_verbose(" k: ", k);

                loaded_data.metal_k_wavelength_nm.back().emplace_back(wavelength_nm);
                loaded_data.metal_k.back().emplace_back(k);
            }
            else
            {
                print_line_verbose("[UNEXPECTED]");
                check_assumption(false);
            }
        }

        check_assumption(loaded_data.metal_n_wavelength_nm.back().size() == loaded_data.metal_k_wavelength_nm.back().size());
    }

    return loaded_data;
}

struct Results
{
    std::vector<Color> base_color;
    std::vector<Color> edge_color;
};

Results process_data(const LoadedData& loaded_data)
{
    Results results;

    for (int index_metal = 0; index_metal < loaded_data.metal_symbol.size(); ++index_metal)
    {
        constexpr double cos_theta0 = 1.0;
        constexpr double cos_theta82 = 1.0 / 7.0;
        constexpr double cos_theta90 = 0.0;

        Color F0_XYZ(0.0, 0.0, 0.0);
        Color F82_XYZ(0.0, 0.0, 0.0);
        Color F90_XYZ(0.0, 0.0, 0.0);
        Color D65_XYZ(0.0, 0.0, 0.0);

        for (int index_CMF = 0; index_CMF < loaded_data.CMF_wavelength_nm.size(); ++index_CMF)
        {
            double target_wavelength_nm = loaded_data.CMF_wavelength_nm[index_CMF];

            double D65 = 0.0;
            for (int index_illuminant = 0; index_illuminant < loaded_data.illuminant_wavelength_nm.size(); ++index_illuminant)
            {
                if (loaded_data.illuminant_wavelength_nm[index_illuminant] == target_wavelength_nm)
                {
                    D65 = loaded_data.illuminant_value[index_illuminant];
                    break;
                }
            }
            check_assumption(D65 > 0.0);

            double n = 0.0;
            double k = 0.0;
            for (int index_current = 0; index_current < loaded_data.metal_n_wavelength_nm[index_metal].size(); ++index_current)
            {
                check_assumption(loaded_data.metal_n_wavelength_nm[index_metal][index_current] == loaded_data.metal_k_wavelength_nm[index_metal][index_current]);
                double current_metal_wavelength_nm = loaded_data.metal_n_wavelength_nm[index_metal][index_current];

                if (current_metal_wavelength_nm == target_wavelength_nm)
                {
                    n = loaded_data.metal_n[index_metal][index_current];
                    k = loaded_data.metal_k[index_metal][index_current];

                    break;
                }
                else if (index_current == 0
                         && current_metal_wavelength_nm > target_wavelength_nm)
                {
                    print_line_verbose(loaded_data.metal_symbol[index_metal], ": EXTRAPOLATING DOWN FROM INSUFFICIENT DATA!");

                    n = loaded_data.metal_n[index_metal][index_current];
                    k = loaded_data.metal_k[index_metal][index_current];

                    break;
                }
                else if (index_current == loaded_data.metal_n_wavelength_nm[index_metal].size() - 1
                         && current_metal_wavelength_nm < target_wavelength_nm)
                {
                    print_line_verbose(loaded_data.metal_symbol[index_metal], ": EXTRAPOLATING UP FROM INSUFFICIENT DATA!");

                    n = loaded_data.metal_n[index_metal][index_current];
                    k = loaded_data.metal_k[index_metal][index_current];

                    break;
                }
                else if (current_metal_wavelength_nm > target_wavelength_nm)
                {
                    int index_previous = index_current - 1;
                    check_assumption(index_previous >= 0);
                    double previous_metal_wavelength_nm = loaded_data.metal_n_wavelength_nm[index_metal][index_previous];

                    double progress = (target_wavelength_nm - previous_metal_wavelength_nm) / (current_metal_wavelength_nm - previous_metal_wavelength_nm);

                    double previous_n = loaded_data.metal_n[index_metal][index_previous];
                    double current_n = loaded_data.metal_n[index_metal][index_current];

                    double previous_k = loaded_data.metal_k[index_metal][index_previous];
                    double current_k = loaded_data.metal_k[index_metal][index_current];

                    n = lerp(previous_n, current_n, progress);
                    k = lerp(previous_k, current_k, progress);

                    break;
                }
            }
            check_assumption(n > 0.0 && k > 0.0);

            double F0 = metal_fresnel(n, k, cos_theta0);
            double F82 = metal_fresnel(n, k, cos_theta82);
            double F90 = metal_fresnel(n, k, cos_theta90);

            Color CMF_XYZ_times_D65 = loaded_data.CMF_XYZ[index_CMF] * D65;
            F0_XYZ = F0_XYZ + (CMF_XYZ_times_D65 * F0);
            F82_XYZ = F82_XYZ + (CMF_XYZ_times_D65 * F82);
            F90_XYZ = F90_XYZ + (CMF_XYZ_times_D65 * F90);
            D65_XYZ = D65_XYZ + (CMF_XYZ_times_D65);
        }

        print_line(loaded_data.metal_symbol[index_metal]);
        print_line();

        print_line("F0 (non-Y-normalized XYZ): ", F0_XYZ);
        print_line("F82 (non-Y-normalized XYZ): ", F82_XYZ);
        print_line("F90 (non-Y-normalized XYZ): ", F90_XYZ);
        print_line("D65 (non-Y-normalized XYZ): ", D65_XYZ);

        Color F0_XYZ_normalized = F0_XYZ / D65_XYZ.y;
        Color F82_XYZ_normalized = F82_XYZ / D65_XYZ.y;
        Color F90_XYZ_normalized = F90_XYZ / D65_XYZ.y;
        Color D65_XYZ_normalized = D65_XYZ / D65_XYZ.y;

        print_line("F0 (XYZ): ", F0_XYZ_normalized);
        print_line("F82 (XYZ): ", F82_XYZ_normalized);
        print_line("F90 (XYZ): ", F90_XYZ_normalized);
        print_line("D65 (XYZ): ", D65_XYZ_normalized);

        print_line("F0 (xyz chromaticity): ", F0_XYZ_normalized / sum_components(F0_XYZ_normalized));
        print_line("F82 (xyz chromaticity): ", F82_XYZ_normalized / sum_components(F82_XYZ_normalized));
        print_line("F90 (xyz chromaticity): ", F90_XYZ_normalized / sum_components(F90_XYZ_normalized));
        print_line("D65 (xyz chromaticity): ", D65_XYZ_normalized / sum_components(D65_XYZ_normalized));

        Color F0_linear_RGB = convert_XYZ_to_linear_Rec709_RGB(F0_XYZ_normalized);
        Color F82_linear_RGB = convert_XYZ_to_linear_Rec709_RGB(F82_XYZ_normalized);
        Color F90_linear_RGB = convert_XYZ_to_linear_Rec709_RGB(F90_XYZ_normalized);
        Color D65_linear_RGB = convert_XYZ_to_linear_Rec709_RGB(D65_XYZ_normalized);

        print_line("F0 (non-D65-relative linear Rec. 709 RGB): ", F0_linear_RGB);
        print_line("F82 (non-D65-relative linear Rec. 709 RGB): ", F82_linear_RGB);
        print_line("F90 (non-D65-relative linear Rec. 709 RGB): ", F90_linear_RGB);
        print_line("D65 (non-D65-relative linear Rec. 709 RGB): ", D65_linear_RGB);

        Color F0_linear_RGB_relative = F0_linear_RGB / D65_linear_RGB;
        Color F82_linear_RGB_relative = F82_linear_RGB / D65_linear_RGB;
        Color F90_linear_RGB_relative = F90_linear_RGB / D65_linear_RGB;
        Color D65_linear_RGB_relative = D65_linear_RGB / D65_linear_RGB;

        print_line("F0 (linear Rec. 709 RGB): ", F0_linear_RGB_relative);
        print_line("F82 (linear Rec. 709 RGB): ", F82_linear_RGB_relative);
        print_line("F90 (linear Rec. 709 RGB): ", F90_linear_RGB_relative);
        print_line("D65 (linear Rec. 709 RGB): ", D65_linear_RGB_relative);

        print_line("F0 (linear P3-D65 RGB): ", convert_XYZ_to_linear_p3d65_RGB(F0_XYZ_normalized) / convert_XYZ_to_linear_p3d65_RGB(D65_XYZ_normalized));
        print_line("F82 (linear P3-D65 RGB): ", convert_XYZ_to_linear_p3d65_RGB(F82_XYZ_normalized) / convert_XYZ_to_linear_p3d65_RGB(D65_XYZ_normalized));
        print_line("F90 (linear P3-D65 RGB): ", convert_XYZ_to_linear_p3d65_RGB(F90_XYZ_normalized) / convert_XYZ_to_linear_p3d65_RGB(D65_XYZ_normalized));
        print_line("D65 (linear P3-D65 RGB): ", convert_XYZ_to_linear_p3d65_RGB(D65_XYZ_normalized) / convert_XYZ_to_linear_p3d65_RGB(D65_XYZ_normalized));

        print_line("F0 (linear Adobe RGB): ", convert_XYZ_to_linear_adobe_RGB(F0_XYZ_normalized) / convert_XYZ_to_linear_adobe_RGB(D65_XYZ_normalized));
        print_line("F82 (linear Adobe RGB): ", convert_XYZ_to_linear_adobe_RGB(F82_XYZ_normalized) / convert_XYZ_to_linear_adobe_RGB(D65_XYZ_normalized));
        print_line("F90 (linear Adobe RGB): ", convert_XYZ_to_linear_adobe_RGB(F90_XYZ_normalized) / convert_XYZ_to_linear_adobe_RGB(D65_XYZ_normalized));
        print_line("D65 (linear Adobe RGB): ", convert_XYZ_to_linear_adobe_RGB(D65_XYZ_normalized) / convert_XYZ_to_linear_adobe_RGB(D65_XYZ_normalized));

        print_line("F0 (linear Rec. 2020 RGB): ", convert_XYZ_to_linear_rec2020_RGB(F0_XYZ_normalized) / convert_XYZ_to_linear_rec2020_RGB(D65_XYZ_normalized));
        print_line("F82 (linear Rec. 2020 RGB): ", convert_XYZ_to_linear_rec2020_RGB(F82_XYZ_normalized) / convert_XYZ_to_linear_rec2020_RGB(D65_XYZ_normalized));
        print_line("F90 (linear Rec. 2020 RGB): ", convert_XYZ_to_linear_rec2020_RGB(F90_XYZ_normalized) / convert_XYZ_to_linear_rec2020_RGB(D65_XYZ_normalized));
        print_line("D65 (linear Rec. 2020 RGB): ", convert_XYZ_to_linear_rec2020_RGB(D65_XYZ_normalized) / convert_XYZ_to_linear_rec2020_RGB(D65_XYZ_normalized));

        print_line("F0 (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(F0_linear_RGB_relative));
        print_line("F82 (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(F82_linear_RGB_relative));
        print_line("F90 (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(F90_linear_RGB_relative));
        print_line("D65 (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(D65_linear_RGB_relative));

        Color base_color_linear_RGB = F0_linear_RGB_relative;
        Color edge_color_linear_RGB = convert_F82_to_F82_tint(F0_linear_RGB_relative, F82_linear_RGB_relative);

        print_line("Base Color (linear Rec. 709 RGB): ", base_color_linear_RGB);
        print_line("Specular Edge Color (linear Rec. 709 RGB): ", edge_color_linear_RGB);
        
        print_line("Base Color (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(base_color_linear_RGB));
        print_line("Specular Edge Color (sRGB): ", convert_linear_Rec709_RGB_to_sRGB(edge_color_linear_RGB));

        print_line("Base Color (clipped 8-bit sRGB): ", convert_to_8_bit(clip_color(convert_linear_Rec709_RGB_to_sRGB(base_color_linear_RGB))));
        print_line("Specular Edge Color (clipped 8-bit sRGB): ", convert_to_8_bit(clip_color(convert_linear_Rec709_RGB_to_sRGB(edge_color_linear_RGB))));
    
        double scale;
        double max_component = std::max(std::max(F0_linear_RGB_relative.x, F0_linear_RGB_relative.y), F0_linear_RGB_relative.z);
        if (max_component > 1.0)
        {
            scale = 1.0 / max_component;
        }
        else
        {
            scale = 1.0;
        }

        Color base_color_scaled_linear_RGB = F0_linear_RGB_relative * scale;
        Color edge_color_scaled_linear_RGB = convert_F82_to_F82_tint(F0_linear_RGB_relative * scale, F82_linear_RGB_relative * scale);

        print_line("Base Color (scaled linear Rec. 709 RGB): ", base_color_scaled_linear_RGB);
        print_line("Specular Edge Color (scaled linear Rec. 709 RGB): ", edge_color_scaled_linear_RGB);
        
        print_line("Base Color (scaled sRGB): ", convert_linear_Rec709_RGB_to_sRGB(base_color_scaled_linear_RGB));
        print_line("Specular Edge Color (scaled sRGB): ", convert_linear_Rec709_RGB_to_sRGB(edge_color_scaled_linear_RGB));

        Color final_base_color = convert_to_8_bit(clip_color(convert_linear_Rec709_RGB_to_sRGB(base_color_scaled_linear_RGB)));
        Color final_edge_color = convert_to_8_bit(clip_color(convert_linear_Rec709_RGB_to_sRGB(edge_color_scaled_linear_RGB)));

        print_line("Base Color (clipped scaled 8-bit sRGB): ", final_base_color);
        print_line("Specular Edge Color (clipped scaled 8-bit sRGB): ", final_edge_color);

        results.base_color.emplace_back(final_base_color);
        results.edge_color.emplace_back(final_edge_color);

        print_line("F0 relative error: ", (metal_schlick_with_F82_tint(base_color_scaled_linear_RGB, edge_color_scaled_linear_RGB, cos_theta0) - F0_linear_RGB_relative) / F0_linear_RGB_relative);
        print_line("F82 relative error: ", (metal_schlick_with_F82_tint(base_color_scaled_linear_RGB, edge_color_scaled_linear_RGB, cos_theta82) - F82_linear_RGB_relative) / F82_linear_RGB_relative);
        print_line("F90 relative error: ", (metal_schlick_with_F82_tint(base_color_scaled_linear_RGB, edge_color_scaled_linear_RGB, cos_theta90) - F90_linear_RGB_relative) / F90_linear_RGB_relative);
    
        print_line();
    }

    return results;
}

void print_results(const LoadedData& loaded_data, const Results& results)
{    
    for (int index_metal = 0; index_metal < loaded_data.metal_symbol.size(); ++index_metal)
    {
        print_line(loaded_data.metal_symbol[index_metal]);
        print_line();

        print_line("Base Color (8-bit sRGB): ", results.base_color[index_metal]);
        print_line("Specular Edge Color (8-bit sRGB): ", results.edge_color[index_metal]);

        print_line();
    }
}

class MetalColors
{
public:
    MetalColors()
    {
        LoadedData loaded_data = gather_data();

        print_line("==============================");
        print_line("DETAILED INFORMATION");
        print_line("==============================");
        print_line();

        Results results = process_data(loaded_data);

        print_line("==============================");
        print_line("SUMMARY");
        print_line("==============================");
        print_line();

        print_results(loaded_data, results);

        print_line("==============================");
    }
};

int main()
{
    MetalColors metal_colors;
    
    return 0;
}
