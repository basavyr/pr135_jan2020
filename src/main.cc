#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "../include/pr135.hh"
#include "../include/data.hh"
#include "../include/energyExpressions.hh"

template <typename T>
void defaultPrinter(std::vector<T> &vec)
{
    for (auto &&n : vec)
    {
        std::cout << n << " ";
    }
    std::cout << "\n";
}

//GENERATING CONTAINER WITH DATA FOR BAND 1
template <typename T>
std::vector<T> band1(Pr135Experimental &obj, EnergyFormula &energy, double param)
{
    std::vector<T> retval;
    energy.normalize(obj.band1);
    energy.kevTOmevBand<T>(obj.band1);
    for (int i = 0; i < obj.band1.size(); ++i)
    {
        // retval.emplace_back(static_cast<T>(EnergyFormula::yrastBand(obj.band1.at(i).spin, param)));
    }
    return retval;
}

//GENERATING CONTAINER WITH DATA FOR BAND 2
template <typename T>
std::vector<T> band2(Pr135Experimental &obj, EnergyFormula &energy, double param)
{
    std::vector<T> retval;
    energy.normalize(obj.band2);
    energy.kevTOmevBand<T>(obj.band2);
    for (int i = 0; i < obj.band2.size(); ++i)
    {
        // retval.emplace_back(static_cast<T>(EnergyFormula::wobblingBand(obj.band2.at(i).spin, param)));
    }
    return retval;
}

//ensdf data initialization
void _init(Pr135Experimental &nucleus, EnergyFormula &formulas)
{
    //generate the two experimental bands from the input data
    nucleus.init_ENSDF(nucleus);

    //band1
    formulas.normalize(nucleus.band1);
    formulas.kevTOmevBand<double>(nucleus.band1);

    //band2
    formulas.normalize(nucleus.band2);
    formulas.kevTOmevBand<double>(nucleus.band2);

    //print the energies
    // nucleus.printer(nucleus.band1);
    // nucleus.printer(nucleus.band2);

    // nucleus.mathPrinter(nucleus.band1);
    // nucleus.mathPrinter(nucleus.band2);
}

int main()
{
    //ENSDF DATA
    Pr135Experimental *nucleus = new Pr135Experimental;
    ChiSquared *chisquared = new ChiSquared;
    EnergyFormula *formulas = new EnergyFormula;
    MinimumValueParameter *paramSet = new MinimumValueParameter;

    _init(*nucleus, *formulas);
    // std::cout << paramSet->calculateMinimumValue<double>(*nucleus, *formulas, *chisquared);

    for (auto &&n : nucleus->band1)
    {
        auto I1 = 100;
        auto I2 = 10;
        auto I3 = 40;
        auto A1 = EnergyFormula::inertiaFactor(static_cast<double>(I1));
        auto A2 = EnergyFormula::inertiaFactor(static_cast<double>(I2));
        auto A3 = EnergyFormula::inertiaFactor(static_cast<double>(I3));
        auto theta = 26.0;
        std::cout << A1 << " " << A2 << " " << A3 << " " << formulas->omega(n.spin, I1, I2, I3, theta) << " " << formulas->yrastBand(n.spin, I1, I2, I3, theta);
        std::cout << "\n";
    }
    std::cout << "\n";

    for (auto &&n : nucleus->band2)
    {
        auto I1 = 100;
        auto I2 = 10;
        auto I3 = 40;
        auto A1 = EnergyFormula::inertiaFactor(static_cast<double>(I1));
        auto A2 = EnergyFormula::inertiaFactor(static_cast<double>(I2));
        auto A3 = EnergyFormula::inertiaFactor(static_cast<double>(I3));
        auto theta = 26.0;
        std::cout << A1 << " " << A2 << " " << A3 << " " << formulas->omega(n.spin, I1, I2, I3, theta) << " " << formulas->wobblingBand(n.spin, I1, I2, I3, theta);
        std::cout << "\n";
    }

    delete nucleus;
    delete chisquared;
    delete formulas;
    delete paramSet;
}