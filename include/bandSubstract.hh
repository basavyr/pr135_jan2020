#ifndef BANDSUBSTRACT_HH
#define BANDSUBSTRACT_HH

#include <vector>
#include <chrono>
#include <ctime>
#include <string>
#include <utility>
#include <memory>
#include <fstream>

class Data_ENSDF
{
public:
    double substractValue;
    //class constructor and destructor
public:
    //class constructor and destructor
    // Data_ENSDF();

    //class constructor and destructor
    ~Data_ENSDF();

    Data_ENSDF(double substraction);
    //containers to hold the experimental data for pr 135
    //containers hold the normalized energies
    //first band is excluded from the fitting procedure!
public:
    std::vector<double> spin1 = {7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    std::vector<double> spin2 = {8.5, 10.5, 12.5, 14.5, 16.5};
    std::vector<double> yrastExp = {358.06, 730.83, 1390.9, 2244.9, 3244.8, 4319.8, 5163.8, 5997.8, 6879.8, 7802.8, 8716};
    std::vector<double> wobbExp = {1433.3, 2158.5, 3001.0, 3957.0, 5069.0};
};

class Data_MATTA
{
    //containers to hold the experimental data for pr 135
    //containers hold the normalized energies
    //first band is excluded from the fitting procedure
public:
    std::vector<double> spin1 = {};
    std::vector<double> spin2 = {};
    std::vector<double> yrastExp = {};
    std::vector<double> wobbExp = {};
};

class BandSubstract
{

    //the useful constants
public:
    static constexpr double j_singleParticle = 5.5;
    static constexpr double PI = 3.14159265358979;

    //the useful formulae for calculus
public:
    static void testApp_Substraction();
    static void smartPointerTest();
    template <typename T>
    static void arrayPrinter(std::vector<T> &array)
    {
        for (int i = 0; i < array.size(); ++i)
        {
            if (i == array.size() - 1)
                std::cout << array.at(i) << "\n";
            else
            {
                std::cout << array.at(i) << " , ";
            }
        }
    }
    template <typename T>
    static void generatePlotData(const std::string &filename, std::vector<T> &spin, std::vector<T> &energy, std::vector<T> &energyAdjusted)
    {
        std::ofstream fout(filename);
        for (int i = 0; i < spin.size(); ++i)
        {
            fout << spin.at(i) << " " << energy.at(i) << " " << energyAdjusted.at(i) << "\n";
        }
    }
};

#endif // BANDSUBSTRACT_HH
