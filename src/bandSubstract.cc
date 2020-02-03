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

void BandSubstract::smartPointerTest(double subtractor)
{
    auto smartPtr = std::make_unique<Data_ENSDF>(subtractor);

    auto yrast = smartPtr->yrastExp;
    auto wobbling = smartPtr->wobbExp;
    auto e0 = static_cast<double>(yrast.at(0));

    //########################################################################
    //########################################################################
    // //STEP1: transform keV to MeV
    // for (int i = 0; i < yrast.size(); ++i)
    // {
    //     // auto x = yrast.at(i);
    //     yrast.at(i) = static_cast<double>(yrast.at(i) / 1000.0);
    // }
    // for (int i = 0; i < wobbling.size(); ++i)
    // {
    //     // auto x = wobbling.at(i);
    //     wobbling.at(i) = static_cast<double>(wobbling.at(i) / 1000.0);
    // }
    // //normalize the containers to the first energy band
    // for (int i = 0; i < yrast.size(); ++i)
    // {
    //     // auto x = yrast.at(i);
    //     yrast.at(i) = yrast.at(i) - e0;
    //     // x = x - e0;
    // }
    // for (int i = 0; i < wobbling.size(); ++i)
    // {
    //     // auto x = wobbling.at(i);
    //     wobbling.at(i) = wobbling.at(i) - e0;
    //     // x = x - e0;
    // }
    //########################################################################
    //########################################################################
    // arrayPrinter(yrast);
    // arrayPrinter(wobbling);

    std::string yrastFile = "../output/yrastBand.dat";
    std::string wobblingFile = "../output/wobblingBand.dat";
    bandSubstracter<double>(yrast, smartPtr->yrastExp_Sub, subtractor);
    bandSubstracter<double>(wobbling, smartPtr->wobExp_Sub, subtractor);
    generatePlotData(yrastFile, smartPtr->spin1, yrast, smartPtr->yrastExp_Sub);
    generatePlotData(wobblingFile, smartPtr->spin2, wobbling, smartPtr->wobExp_Sub);
    arrayPrinter(smartPtr->yrastExp_Sub);
    arrayPrinter(smartPtr->wobExp_Sub);
}

// Data_ENSDF::Data_ENSDF()
// {
//     std::string constructorMessage = "The class for storing the ENSDF experimental data is succesfully created";
//     std::cout << constructorMessage;
//     std::cout << std::endl;
// }

Data_ENSDF::Data_ENSDF(double substraction)
{
    std::string part1Message = "Class containers constructed with the const value ";
    this->substractValue = static_cast<double>(substraction);
    std::string part2Message = " as a band substracter";
    std::cout << part1Message << substraction << part2Message;
    std::cout << std::endl;
}

Data_ENSDF::~Data_ENSDF()
{
    std::string destructorMessage = "Class containers has been succsesfully destroyed";
    std::cout << destructorMessage;
    std::cout << std::endl;
}

double EnergyFormulae::inertiaFactor(double moi)
{
    return static_cast<double>(1.0 / (2.0 * moi));
}

double EnergyFormulae::j_Component(int k, double theta)
{
    if (k == 1)
        return BandSubstract::j_singleParticle * cos(theta * BandSubstract::PI / 180.00);
    return BandSubstract::j_singleParticle * sin(theta * BandSubstract::PI / 180);
}

double EnergyFormulae::omega(double spin, double i1, double i2, double i3, double theta)
{
    //generate the inertia factors
    auto a1 = inertiaFactor(i1);
    auto a2 = inertiaFactor(i2);
    auto a3 = inertiaFactor(i3);
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);
    auto term1 = (2.0 * spin + 1.0) * (a2 - a1 - (a2 * j2 / spin)) - 2.0 * a1 * j1;
    auto term2 = (2.0 * spin + 1.0) * (a3 - a1) - 2.0 * a1 * j1;
    auto term3 = (a3 - a1) * (a2 - a1 - (a2 * j2 / spin));
    auto result = sqrt(term1 * term2 - term3);
    return static_cast<double>(result);
}

double EnergyFormulae::omegaPrime(double spin, double i1, double i2, double i3, double theta)
{
    //generate the inertia factors
    auto a1 = inertiaFactor(i1);
    auto a2 = inertiaFactor(i2);
    auto a3 = inertiaFactor(i3);
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);
    auto term1 = (2.0 * spin + 1.0) * (a2 - a1 - (a2 * j2 / spin)) + 2.0 * a1 * j1;
    auto term2 = (2.0 * spin + 1.0) * (a3 - a1) + 2.0 * a1 * j1;
    auto term3 = (a3 - a1) * (a2 - a1 - (a2 * j2 / spin));
    auto result = sqrt(term1 * term2 - term3);
    return static_cast<double>(result);
}