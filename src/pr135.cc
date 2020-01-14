#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "../include/pr135.hh"
#include "../include/data.hh"

void Pr135Experimental::init_ENSDF(Pr135Experimental &obj)
{
    //ENSDF DATA
    //band1 experimental values
    for (int i = 0; i < Constants::dim1; ++i)
    {
        obj.exp1.emplace_back(band());
        exp1.at(i).spin = static_cast<double>(ExperimentalData_ENSDF::spin1[i]);
        exp1.at(i).energy = static_cast<double>(ExperimentalData_ENSDF::energy1[i]);
    }

    //ENSDF DATA
    //band2 experimental values
    for (int i = 0; i < Constants::dim2; ++i)
    {
        obj.exp2.emplace_back(band());
        exp2.at(i).spin = static_cast<double>(ExperimentalData_ENSDF::spin2[i]);
        exp2.at(i).energy = static_cast<double>(ExperimentalData_ENSDF::energy2[i]);
    }
}

void Pr135Experimental::init_MATTA(Pr135Experimental &obj)
{
    //MATTA DATA
    //band1 experimental values
    for (int i = 0; i < Constants::dim1; ++i)
    {
        obj.exp1.emplace_back(band());
        exp1.at(i).spin = static_cast<double>(ExperimentalData_MATTA::spin1[i]);
        exp1.at(i).energy = static_cast<double>(ExperimentalData_MATTA::energy1[i]);
    }

    //MATTA DATA
    //band2 experimental values
    for (int i = 0; i < Constants::dim2; ++i)
    {
        obj.exp2.emplace_back(band());
        exp2.at(i).spin = static_cast<double>(ExperimentalData_MATTA::spin2[i]);
        exp2.at(i).energy = static_cast<double>(ExperimentalData_MATTA::energy2[i]);
    }
}

void Pr135Experimental::newLine()
{
    std::cout << "\n";
}

void Pr135Experimental::printer(std::vector<Pr135Experimental::band> &array)
{
    for (auto &&n : array)
    {
        std::cout << n.spin << " " << n.energy;
        newLine();
    }
}

std::vector<Pr135Experimental::band> Pr135Experimental::data1Exp(Pr135Experimental &obj)
{
    std::vector<Pr135Experimental::band> retval;
    for (auto &&n : obj.exp1)
    {
        retval.emplace_back(n);
    }
    return retval;
}

std::vector<Pr135Experimental::band> Pr135Experimental::data2Exp(Pr135Experimental &obj)
{
    std::vector<Pr135Experimental::band> retval;
    for (auto &&n : obj.exp2)
    {
        retval.emplace_back(n);
    }
    return retval;
}