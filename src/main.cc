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
        retval.emplace_back(static_cast<T>(EnergyFormula::yrastBand(obj.band1.at(i).spin, param)));
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
        retval.emplace_back(static_cast<T>(EnergyFormula::wobblingBand(obj.band2.at(i).spin, param)));
    }
    return retval;
}

//ensdf data initialization
void _init(Pr135Experimental &nucleus, EnergyFormula &formulas)
{
    nucleus.init_ENSDF(nucleus);

    //band1
    formulas.normalize(nucleus.band1);
    formulas.kevTOmevBand<double>(nucleus.band1);

    //band2
    formulas.normalize(nucleus.band2);
    formulas.kevTOmevBand<double>(nucleus.band2);

    nucleus.printer(nucleus.band1);
    nucleus.printer(nucleus.band2);
}

int main()
{
    //ENSDF DATA
    Pr135Experimental *nucleus = new Pr135Experimental;
    ChiSquared *chisquared = new ChiSquared;
    EnergyFormula *formulas = new EnergyFormula;
    MinimumValueParameter *paramSet = new MinimumValueParameter;

    _init(*nucleus, *formulas);
    std::cout << paramSet->calculateMinimumValue<double>(*nucleus, *formulas, *chisquared);

    delete nucleus;
    delete chisquared;
    delete formulas;
    delete paramSet;
}