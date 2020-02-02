#include <iostream>

#include "../include/bandSubstract.hh"

void BandSubstract::testApp_Substraction()
{
    auto now = std::chrono::system_clock::now();
    std::time_t datenow = std::chrono::system_clock::to_time_t(now);
    std::cout << ctime(&datenow);
    std::cout << "App works";
}

//testing the smart pointers

void BandSubstract::smartPointerTest()
{
    auto smartPtr = std::make_unique<Data_ENSDF>(0.1);

    auto yrast = smartPtr->yrastExp;
    auto wobbling = smartPtr->wobbExp;
    auto e0 = static_cast<double>(yrast.at(0) / 1000.0);

    //STEP1: transform keV to MeV
    for (int i = 0; i < yrast.size(); ++i)
    {
        // auto x = yrast.at(i);
        yrast.at(i) = static_cast<double>(yrast.at(i) / 1000.0);
    }
    for (int i = 0; i < wobbling.size(); ++i)
    {
        // auto x = wobbling.at(i);
        wobbling.at(i) = static_cast<double>(wobbling.at(i) / 1000.0);
    }
    //normalize the containers to the first energy band
    for (int i = 0; i < yrast.size(); ++i)
    {
        // auto x = yrast.at(i);
        yrast.at(i) = yrast.at(i) - e0;
        // x = x - e0;
    }
    for (int i = 0; i < wobbling.size(); ++i)
    {
        // auto x = wobbling.at(i);
        wobbling.at(i) = wobbling.at(i) - e0;
        // x = x - e0;
    }
    arrayPrinter(yrast);
    arrayPrinter(wobbling);

    std::string file1 = "../output/plot.dat";
    generatePlotData(file1, smartPtr->spin1, yrast, yrast);
}

// Data_ENSDF::Data_ENSDF()
// {
//     std::string constructorMessage = "The class for storing the ENSDF experimental data is succesfully created";
//     std::cout << constructorMessage;
//     std::cout << std::endl;
// }

Data_ENSDF::Data_ENSDF(double substraction)
{
    std::string part1Message = "Class container s constructed with the const value ";
    this->substractValue = static_cast<double>(substraction);
    std::string part2Message = " as a band substracter";
    std::cout << part1Message << substraction << part2Message;
    std::cout << std::endl;
}

Data_ENSDF::~Data_ENSDF()
{
    std::string destructorMessage = "The class container has been succsesfully destroyed";
    std::cout << destructorMessage;
    std::cout << std::endl;
}