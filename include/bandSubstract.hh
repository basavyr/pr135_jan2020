#ifndef BANDSUBSTRACT_HH
#define BANDSUBSTRACT_HH

#include <vector>
#include <chrono>
#include <ctime>
#include <string>
#include <utility>
#include <memory>

class Data_ENSDF
{
    //class constructor and destructor
public:
    //class constructor and destructor
    Data_ENSDF();

    //class constructor and destructor
    ~Data_ENSDF();

    //containers to hold the experimental data for pr 135
    //containers hold the normalized energies
    //first band is excluded from the fitting procedure!
public:
    std::vector<double> spin1 = {};
    std::vector<double> spin2 = {};
    std::vector<double> yrastExp = {};
    std::vector<double> wobbExp = {};
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
};

#endif // BANDSUBSTRACT_HH
