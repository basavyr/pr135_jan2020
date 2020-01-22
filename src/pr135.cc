#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "../include/pr135.hh"
#include "../include/data.hh"

//fetches the experimental data from C arrays into container object BAND
void Pr135Experimental::init_ENSDF(Pr135Experimental &obj)
{
    //ENSDF DATA
    //band1 experimental values
    for (int i = 0; i < Constants::dim1; ++i)
    {
        obj.band1.emplace_back(band());
        //determine R component of the total spin I
        // auto spinCorrection = static_cast<double>(ExperimentalData_ENSDF::spin1[i] - Constants::oddSpin);
        //energy formulas work with the spin correction I=R+j
        // band1.at(i).spin = static_cast<double>(spinCorrection);
        //energy formulas with the standard I spin from the experimental data
        band1.at(i).spin = static_cast<double>(ExperimentalData_ENSDF::spin1[i]);
        band1.at(i).energy = static_cast<double>(ExperimentalData_ENSDF::energy1[i]);
    }

    //ENSDF DATA
    //band2 experimental values
    for (int i = 0; i < Constants::dim2; ++i)
    {
        obj.band2.emplace_back(band());
        //determine R component of the total spin I
        // auto spinCorrection = static_cast<double>(ExperimentalData_ENSDF::spin2[i] - Constants::oddSpin);
        //energy formulas work with the spin correction I=R+j
        // band2.at(i).spin = static_cast<double>(spinCorrection);
        //energy formulas with the standard I spin from the experimental data
        band2.at(i).spin = static_cast<double>(ExperimentalData_ENSDF::spin2[i]);
        band2.at(i).energy = static_cast<double>(ExperimentalData_ENSDF::energy2[i]);
    }
}

//fetches the experimental data from C arrays into container object BAND
void Pr135Experimental::init_MATTA(Pr135Experimental &obj)
{
    //MATTA DATA
    //band1 experimental values
    for (int i = 0; i < Constants::dim1; ++i)
    {
        obj.band1.emplace_back(band());
        band1.at(i).spin = static_cast<double>(ExperimentalData_MATTA::spin1[i]);
        band1.at(i).energy = static_cast<double>(ExperimentalData_MATTA::energy1[i]);
    }

    //MATTA DATA
    //band2 experimental values
    for (int i = 0; i < Constants::dim2; ++i)
    {
        obj.band2.emplace_back(band());
        band2.at(i).spin = static_cast<double>(ExperimentalData_MATTA::spin2[i]);
        band2.at(i).energy = static_cast<double>(ExperimentalData_MATTA::energy2[i]);
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

void Pr135Experimental::mathPrinter(std::vector<Pr135Experimental::band> &array)
{
    std::cout << "{ ";
    for (int i = 0; i < array.size(); ++i)
    {
        if (i == array.size() - 1)
        {
            std::cout << array.at(i).energy;
            break;
        }
        std::cout << array.at(i).energy << " , ";
    }
    std::cout << "};\n";
}


std::vector<Pr135Experimental::band> Pr135Experimental::data1Exp(Pr135Experimental &obj)
{
    std::vector<Pr135Experimental::band> retval;
    for (auto &&n : obj.band1)
    {
        retval.emplace_back(n);
    }
    return retval;
}

std::vector<Pr135Experimental::band> Pr135Experimental::data2Exp(Pr135Experimental &obj)
{
    std::vector<Pr135Experimental::band> retval;
    for (auto &&n : obj.band2)
    {
        retval.emplace_back(n);
    }
    return retval;
}

void Pr135Experimental::cleanYrastBand(Pr135Experimental &obj)
{
    auto firstElement = obj.band1.at(0);
    auto size = obj.band1.size();
    // std::cout << "the first element of the yrast band is :" << firstElement.energy << " " << firstElement.spin << " " << size << "\n";
    obj.band1.erase(obj.band1.begin(), obj.band1.begin() + 1);
    firstElement = obj.band1.at(0);
    size = obj.band1.size();
    // std::cout << "the first element of the yrast band is :" << firstElement.energy << " " << firstElement.spin << " " << size << "\n";
}

void Pr135Experimental::cleanWobblingBand(Pr135Experimental &obj)
{
    auto firstElement = obj.band2.at(0);
    auto size = obj.band2.size();
    obj.band2.erase(obj.band2.begin(), obj.band2.begin() + 1);
    firstElement = obj.band2.at(0);
    size = obj.band2.size();
}