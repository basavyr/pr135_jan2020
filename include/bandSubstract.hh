#ifndef BANDSUBSTRACT_HH
#define BANDSUBSTRACT_HH

#include <vector>
#include <chrono>
#include <ctime>
#include <string>
#include <utility>
#include <memory>
#include <cmath>
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
    std::vector<double> yrastExp = {0.37277, 1.03284, 1.88684, 2.88674, 3.96174, 4.80574, 5.63974,
                                    6.52174, 7.44474, 8.35794};
    std::vector<double> wobbExp = {1.07524, 1.80044, 2.64294, 3.59894, 4.71094};

    //containers to hold the substracted energies
    std::vector<double> yrastExp_Sub;
    std::vector<double> wobExp_Sub;

    //containers for storing the theoretical energies (obtained after the best root mean square)
    std::vector<double> yrastTh;
    std::vector<double> wobblingTh;
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

    //containers to hold the substracted energies
    std::vector<double> yrastExp_Sub;
    std::vector<double> wobExp_Sub;

    //containers for storing the theoretical energies (obtained after the best root mean square)
    std::vector<double> yrastTh;
    std::vector<double> wobblingTh;
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
    static void smartPointerTest(double);

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
    
    template <typename T>
    static void bandSubstracter(std::vector<T> &object, std::vector<T> &newObject, double subtractor)
    {
        auto halfSize = static_cast<size_t>(object.size() / 2);
        for (int i = 0; i < object.size(); ++i)
        {
            if (i < halfSize)
            {
                newObject.emplace_back(object.at(i) - subtractor);
            }
            else
            {
                newObject.emplace_back(object.at(i));
            }
        }
        // std::cout << "the half band is of length= " << halfSize;
        // std::cout << "\n";
    }
};

class EnergyFormulae
{
public:
    //calculates the obbling frequency of the nucleus
    static double omega(double spin, double i1, double i2, double i3, double theta);

    //calculates the wobbling frequency of the nucleus without the single particle term in the total hamiltonian
    static double omegaPrime(double spin, double i1, double i2, double i3, double theta);

    //returns the energy calculation of an arbitrary spin state from a band
    static double energyExpression(int N, double spin, double i2, double i2, double i3, double theta);

    //tranforms the moments of inertia into inertia factors A_k
    static double inertiaFactor(double moi);

    //computes the components of the single particle angular momentum
    static double j_Component(int k, double theta);
}

#endif // BANDSUBSTRACT_HH
