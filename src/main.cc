#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "../include/pr135.hh"

void Pr135Experimental::init_ENSDF(Pr135Experimental &obj)
{
    //ENSDF DATA
    double spin1[Pr135Experimental::dim1] = {5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    double energy1[Pr135Experimental::dim1] = {358.06, 730.83, 1390.9, 2244.9, 3244.8, 4319.8, 5163.8, 5997.8, 6879.8, 7802.8, 8716};
    //band1 experimental values
    for (int i = 0; i < Pr135Experimental::dim1; ++i)
    {
        obj.exp1.emplace_back(band());
        exp1.at(i).spin = static_cast<double>(spin1[i]);
        exp1.at(i).energy = static_cast<double>(energy1[i]);
    }

    //ENSDF DATA
    double spin2[Pr135Experimental::dim2] = {8.5, 10.5, 12.5, 14.5, 16.5};
    double energy2[Pr135Experimental::dim2] = {1433.3, 2158.5, 3001.0, 3957.0, 5069.0};
    //band2 experimental values
    for (int i = 0; i < Pr135Experimental::dim2; ++i)
    {
        obj.exp2.emplace_back(band());
        exp2.at(i).spin = static_cast<double>(spin2[i]);
        exp2.at(i).energy = static_cast<double>(energy2[i]);
    }
}

void Pr135Experimental::init_MATTA(Pr135Experimental &obj)
{
    //ENSDF DATA
    double spin1[Pr135Experimental::dim1] = {5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    double energy1[Pr135Experimental::dim1] = {358.06, 730.83, 1390.9, 2244.9, 3244.8, 4319.8, 5163.8, 5997.8, 6879.8, 7802.8, 8716};
    //band1 experimental values
    for (int i = 0; i < Pr135Experimental::dim1; ++i)
    {
        obj.exp1.emplace_back(band());
        exp1.at(i).spin = static_cast<double>(spin1[i] - 1);
        exp1.at(i).energy = static_cast<double>(energy1[i] - 1);
    }

    //ENSDF DATA
    double spin2[Pr135Experimental::dim2] = {8.5, 10.5, 12.5, 14.5, 16.5};
    double energy2[Pr135Experimental::dim2] = {1433.3, 2158.5, 3001.0, 3957.0, 5069.0};
    //band2 experimental values
    for (int i = 0; i < Pr135Experimental::dim2; ++i)
    {
        obj.exp2.emplace_back(band());
        exp2.at(i).spin = static_cast<double>(spin2[i] - 1);
        exp2.at(i).energy = static_cast<double>(energy2[i] - 1);
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
        Pr135Experimental::newLine();
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

int main()
{
    Pr135Experimental *nucleu = new Pr135Experimental;

    nucleu->init_ENSDF(*nucleu);
    std::vector<Pr135Experimental::band> v1 = nucleu->data1Exp(*nucleu);
    Pr135Experimental::printer(v1);

    std::cout << "***************";
    Pr135Experimental::newLine();

    std::vector<Pr135Experimental::band> v2 = nucleu->data2Exp(*nucleu);
    Pr135Experimental::printer(v2);
    delete nucleu;
    v1.clear();
    v2.clear();

    nucleu = new Pr135Experimental;
    nucleu->init_MATTA(*nucleu);
    v1 = nucleu->data1Exp(*nucleu);
    Pr135Experimental::printer(v1);

    std::cout << "***************";
    Pr135Experimental::newLine();

    v2 = nucleu->data2Exp(*nucleu);
    Pr135Experimental::printer(v2);
    delete nucleu;
    return 0;
}