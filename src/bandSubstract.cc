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
    //create the subtracted energy set
    auto smartPtr = std::make_unique<Data_ENSDF>(subtractor);
    //create the set for storing best fit params
    auto bestParams = std::make_unique<RMS_Calculus::minParamSet>();

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
    // generatePlotData(yrastFile, smartPtr->spin1, yrast, smartPtr->yrastExp_Sub);
    // generatePlotData(wobblingFile, smartPtr->spin2, wobbling, smartPtr->wobExp_Sub);
    // arrayPrinter(smartPtr->yrastExp_Sub);
    // arrayPrinter(smartPtr->wobExp_Sub);
    // bandGeneration(*smartPtr, 1, 1, 1, 1);
    RMS_Calculus::searchMinimum<Data_ENSDF>(*smartPtr, *bestParams);
    RMS_Calculus::paramPrinter<RMS_Calculus::minParamSet>(*bestParams);
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
    std::cout << "\n";
    std::cout << "*********************************************";
    std::cout << std::endl;
}

BandSubstract::rms_Tuple BandSubstract::bandGeneration(Data_ENSDF &object, double i1, double i2, double i3, double theta)
{
    auto yrast = [&](auto spin) { return EnergyFormulae::energyExpression(0, spin, i1, i2, i3, theta); };
    auto wobbling = [&](auto spin) { return EnergyFormulae::energyExpression(0, spin, i1, i2, i3, theta); };
    std::vector<double> dataExp;
    std::vector<double> dataExp_pure;
    std::vector<double> dataTh;
    for (int i = 0; i < object.yrastExp_Sub.size(); ++i)
    {
        dataExp.emplace_back(object.yrastExp_Sub.at(i));
        dataExp_pure.emplace_back(object.yrastExp.at(i));
        auto currentValue = yrast(object.spin1.at(i));
        if (currentValue == 6969)
        {
            /* THIS IS NOT OK*/
            //should never fail!
        }
        else
        {
            dataTh.emplace_back(static_cast<double>(currentValue));
        }
    }
    for (int i = 0; i < object.wobExp_Sub.size(); ++i)
    {
        dataExp.emplace_back(object.wobExp_Sub.at(i));
        dataExp_pure.emplace_back(object.wobbExp.at(i));
        auto currentValue = wobbling(object.spin2.at(i));
        if (currentValue == 6969)
        {
            /* THIS IS NOT OK*/
            //should never fail!
        }
        else
        {
            dataTh.emplace_back(static_cast<double>(currentValue));
        }
    }
    // dataExp.clear();
    // dataTh.clear();
    auto rms_subtracted = RMS_Calculus::rootMeanSquare(dataExp, dataTh);
    auto rms_pure = RMS_Calculus::rootMeanSquare(dataExp_pure, dataTh);
    rms_Tuple result;
    result.rms_substracted = rms_subtracted;
    result.rms_pure = rms_pure;
    return result;
    // arrayPrinter(dataExp);
    // arrayPrinter(dataTh);
    // std::cout << std::endl;
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
    //generate the single particle a.m. components
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);
    auto term1 = (2.0 * spin + 1.0) * (a2 - a1 - (a2 * j2 / spin)) - 2.0 * a1 * j1;
    auto term2 = (2.0 * spin + 1.0) * (a3 - a1) - 2.0 * a1 * j1;
    auto term3 = (a3 - a1) * (a2 - a1 - (a2 * j2 / spin));
    auto result = sqrt(term1 * term2 - term3);
    if (isnan(result))
        return 6969;
    return static_cast<double>(result);
}

double EnergyFormulae::omegaPrime(double spin, double i1, double i2, double i3, double theta)
{
    //generate the inertia factors
    auto a1 = inertiaFactor(i1);
    auto a2 = inertiaFactor(i2);
    auto a3 = inertiaFactor(i3);
    //generate the single particle a.m. components
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);
    auto term1 = (2.0 * spin + 1.0) * (a2 - a1 - (a2 * j2 / spin)) + 2.0 * a1 * j1;
    auto term2 = (2.0 * spin + 1.0) * (a3 - a1) + 2.0 * a1 * j1;
    auto term3 = (a3 - a1) * (a2 - a1 - (a2 * j2 / spin));
    auto result = sqrt(term1 * term2 - term3);
    if (isnan(result))
        return 6969;
    return static_cast<double>(result);
}

double EnergyFormulae::minHamiltonian(double spin, double i1, double i2, double i3, double theta)
{
    //generate the inertia factors
    auto a1 = inertiaFactor(i1);
    auto a2 = inertiaFactor(i2);
    auto a3 = inertiaFactor(i3);
    //generate the single particle a.m. components
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);

    auto singleParticleSum = a1 * (j1 * j1) + a2 * (j2 * j2);
    auto minTerm = a1 * spin * spin + (2.0 * spin + 1.0) * a1 * j1 - spin * a2 * j2;
    auto result = minTerm + singleParticleSum;
    if (isnan(result))
        return 6969;
    return static_cast<double>(result);
}

double EnergyFormulae::energyExpression(int N, double spin, double i1, double i2, double i3, double theta)
{
    auto omega = EnergyFormulae::omega(spin, i1, i2, i3, theta);
    //stop immediately if the wobbling frequency is not real
    if (omega == 6969)
        return 6969;
    //generate the inertia factors
    auto a1 = inertiaFactor(i1);
    auto a2 = inertiaFactor(i2);
    auto a3 = inertiaFactor(i3);
    auto result = static_cast<double>(minHamiltonian(spin, i1, i2, i3, theta) + (N + 0.5) * EnergyFormulae::omega(spin, i1, i2, i3, theta));
    if (isnan(result))
        return 6969;
    return static_cast<double>(result);
}