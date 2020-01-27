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
    formulas.normalizeBands(nucleus);
    // formulas.normalize(nucleus.band1);
    formulas.kevTOmevBand<double>(nucleus.band1);

    //band2
    // formulas.normalize(nucleus.band2);
    formulas.kevTOmevBand<double>(nucleus.band2);
    nucleus.cleanYrastBand(nucleus);
    // nucleus.cleanWobblingBand(nucleus);

    //print the energies
    // nucleus.printer(nucleus.band1);
    // nucleus.printer(nucleus.band2);

    // nucleus.mathPrinter(nucleus.band1);
    // nucleus.mathPrinter(nucleus.band2);
}

//matta data initialization
void _init_matta(Pr135Experimental &nucleus, EnergyFormula &formulas)
{
    //generate the two experimental bands from the input data
    nucleus.init_MATTA(nucleus);

    //band1
    formulas.normalizeBands(nucleus);
    // formulas.normalize(nucleus.band1);
    formulas.kevTOmevBand<double>(nucleus.band1);

    //band2
    // formulas.normalize(nucleus.band2);
    formulas.kevTOmevBand<double>(nucleus.band2);
    nucleus.cleanYrastBand(nucleus);
    // nucleus.cleanWobblingBand(nucleus);

    //print the energies
    nucleus.printer(nucleus.band1);
    nucleus.printer(nucleus.band2);

    // nucleus.mathPrinter(nucleus.band1);
    // nucleus.mathPrinter(nucleus.band2);
}

void testForComplexValue_yrast()
{
    std::cout << "*********"
              << "\n";
    std::cout << "YRAST..."
              << "\n";
    EnergyFormula::ParameterSet params;
    int noIterations = 0;
    int noFails = 0;
    for (double i1 = params.A_left; i1 <= params.A_right; i1 += params.A_step)
    {
        auto a1 = EnergyFormula::inertiaFactor(i1);
        for (double i2 = params.A_left; i2 <= params.A_right; i2 += params.A_step)
        {
            auto a2 = EnergyFormula::inertiaFactor(i2);
            for (double i3 = params.A_left; i3 <= params.A_right; i3 += params.A_step)
            {
                auto a3 = EnergyFormula::inertiaFactor(i3);
                for (double theta = params.theta_left; theta <= params.theta_right; theta += params.theta_step)
                {

                    //this condition is avoided with the new updated starting param for A
                    /*  if (!a1 || !a2 || !a3)
                    {
                        std::cout << "CAN'T COMPUTE INERTIA FACTOR for I_k= " << i1 << " " << i2 << " " << i3 << " ..MOI IS ZERO"
                                  << "\n ";
                        break;
                    } */
                    {
                        /* code */
                        for (int k = 0; k < 11; ++k)
                        {
                            auto energy = EnergyFormula::yrastBand(5.5 + 2 * k, i1, i2, i3, theta);
                            auto omega = EnergyFormula::omega(5.5 + 2 * k, i1, i2, i3, theta);
                            //check the complex wobbling frequency
                            if (omega == 6969)
                            {
                                // std::cout << "It failed for...E= " << energy << " with spin I= " << 5.5 + 2 * k << " and params { " << i1 << " " << i2 << " " << i3 << " " << theta << "\n";
                                // i1 = params.A_right + params.A_step;
                                // i2 = params.A_right + params.A_step;
                                // i3 = params.A_right + params.A_step;
                                // theta = params.theta_right + params.theta_step;
                                noFails++;
                                // break;
                            }
                            noIterations++;
                        }
                        // std::cout << omega << " " << energy << "\n";
                        /* 
                        auto omega = EnergyFormula::omega(5.5, a1, a2, a3, theta);
                        if (omega == 6969)
                        {
                            std::cout << "FOUND COMPLEX FREQUENCY..." << i1 << " " << i2 << " " << i3 << " " << omega
                                      << "\n";
                            if (noIterations == 215793)
                                std::cout << typeid(omega).name() << "\n";
                            break;
                        } */
                        // std::cout << i1 << " " << i2 << " " << i3 << " " << a1 << " " << a2 << " " << a3 << omega << "\n";
                    }
                }
            }
        }
    }
    std::cout << "TEST PASSED..."
              << "\n";
    std::cout << "TOTAL N.O. ITERATIONS... " << noIterations << "\n";
    std::cout << "TOTAL N.O. FAILS... " << noFails << "\n";
    std::cout << "Fail rate = " << static_cast<double>((static_cast<double>(noFails) / static_cast<double>(noIterations)) * 100) << "%\n";
}
void testForComplexValue_wobbling()
{
    std::cout << "*********"
              << "\n";
    std::cout << "WOBBLING..."
              << "\n";
    EnergyFormula::ParameterSet params;
    int noIterations = 0;
    int noFails = 0;
    for (double i1 = params.A_left; i1 <= params.A_right; i1 += params.A_step)
    {
        auto a1 = EnergyFormula::inertiaFactor(i1);
        for (double i2 = params.A_left; i2 <= params.A_right; i2 += params.A_step)
        {
            auto a2 = EnergyFormula::inertiaFactor(i2);
            for (double i3 = params.A_left; i3 <= params.A_right; i3 += params.A_step)
            {
                auto a3 = EnergyFormula::inertiaFactor(i3);
                for (double theta = params.theta_left; theta <= params.theta_right; theta += params.theta_step)
                {

                    //this condition is avoided with the new updated starting param for A
                    /*  if (!a1 || !a2 || !a3)
                    {
                        std::cout << "CAN'T COMPUTE INERTIA FACTOR for I_k= " << i1 << " " << i2 << " " << i3 << " ..MOI IS ZERO"
                                  << "\n ";
                        break;
                    } */
                    {
                        /* code */
                        for (int k = 0; k < 5; ++k)
                        {
                            auto omega = EnergyFormula::omega(8.5 + 2 * k, i1, i2, i3, theta);
                            auto energy = EnergyFormula::wobblingBand(8.5 + 2 * k, i1, i2, i3, theta);
                            //check the complex wobbling frequency
                            if (omega == 6969)
                            {
                                // std::cout << "It failed for..." << energy<< " with spin I= " << 8.5 + 2 * k << " "<< " and params { " << i1 << " " << i2 << " " << i3 << " " << theta << "\n";
                                // i1 = params.A_right + params.A_step;
                                // i2 = params.A_right + params.A_step;
                                // i3 = params.A_right + params.A_step;
                                // theta = params.theta_right + params.theta_step;
                                noFails++;
                                // break;
                            }
                            noIterations++;
                        }
                        // std::cout << omega << " " << energy << "\n";
                        /* 
                        auto omega = EnergyFormula::omega(5.5, a1, a2, a3, theta);
                        if (omega == 6969)
                        {
                            std::cout << "FOUND COMPLEX FREQUENCY..." << i1 << " " << i2 << " " << i3 << " " << omega
                                      << "\n";
                            if (noIterations == 215793)
                                std::cout << typeid(omega).name() << "\n";
                            break;
                        } */
                        // std::cout << i1 << " " << i2 << " " << i3 << " " << a1 << " " << a2 << " " << a3 << omega << "\n";
                    }
                }
            }
        }
    }
    std::cout << "TEST PASSED..."
              << "\n";
    std::cout << "TOTAL N.O. ITERATIONS... " << noIterations << "\n";
    std::cout << "TOTAL N.O. FAILS... " << noFails << "\n";
    std::cout << "Fail rate = " << static_cast<double>((static_cast<double>(noFails) / static_cast<double>(noIterations)) * 100) << "%\n";
}

int main()
{
    //ENSDF DATA
    Pr135Experimental *nucleus = new Pr135Experimental;
    ChiSquared *chisquared = new ChiSquared;
    EnergyFormula *formulas = new EnergyFormula;
    MinimumValueParameter *paramSet = new MinimumValueParameter;

    double paramsImported[5] = {129.144 , 9.91615 , 0.364585 , -1.98};

    //initialize the containers with the experimental data from ENSDF
    _init_matta(*nucleus, *formulas);
    // std::cout << paramSet->calculateMinimumValue<double>(*nucleus, *formulas, *chisquared);

    //initialize the containers with the experimental data from MATTA
    // _init_matta(*nucleus, *formulas);
    // std::cout << paramSet->calculateMinimumValue<double>(*nucleus, *formulas, *chisquared);
    // std::cout << chisquared->applyEnergies<double>(*nucleus, EnergyFormula::inertiaFactor(paramsImported[0]), EnergyFormula::inertiaFactor(paramsImported[1]), EnergyFormula::inertiaFactor(paramsImported[2]), EnergyFormula::inertiaFactor(paramsImported[3]));
    std::cout << chisquared->applyEnergies<double>(*nucleus, paramsImported[0], paramsImported[1], paramsImported[2], paramsImported[3]);

    for (auto &&n : nucleus->band1)
    {
        auto I1 = 100;
        auto I2 = 10;
        auto I3 = 40;
        auto A1 = EnergyFormula::inertiaFactor(static_cast<double>(I1));
        auto A2 = EnergyFormula::inertiaFactor(static_cast<double>(I2));
        auto A3 = EnergyFormula::inertiaFactor(static_cast<double>(I3));
        auto theta = 26.0;
        // std::cout << A1 << " " << A2 << " " << A3 << " " << formulas->omega(n.spin, I1, I2, I3, theta) << " " << formulas->yrastBand(n.spin, I1, I2, I3, theta) << "\n";
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
        // std::cout << A1 << " " << A2 << " " << A3 << " " << formulas->omega(n.spin, I1, I2, I3, theta) << " " << formulas->wobblingBand(n.spin, I1, I2, I3, theta)<< "\n";
    }

    // testForComplexValue_yrast();

    // testForComplexValue_wobbling();

    // std::cout << EnergyFormula::omega(5.5, 27, 14, 9, 0) << "\n";

    // std::cout << chisquared->applyEnergies<double>(*nucleus, 1, 1, 1, 60) << "\n";
    // // std::cout << formulas->yrastBand(5.5, 116, 1, 1, 85) << "\n";
    // // std::cout << formulas->yrastBand(7.5, 116, 1, 1, 86) << "\n";
    // std::cout << "**********\n";
    // // std::cout << formulas->wobblingBand(8.5, 101, 1, 66, 0) << "\n";
    // // std::cout << formulas->wobblingBand(10.5, 101, 1, 66, 1) << "\n";
    // std::cout << chisquared->applyEnergies<double>(*nucleus, 6, 91, 1, 70) << "\n";

    /* for (int i = 0; i < nucleus->band1.size(); ++i)
    {
        auto spin = nucleus->band1.at(i).spin;
        std::cout << spin << " " << EnergyFormula::yrastBand(spin, 116, 1, 1, 85) << "\n";
    }
    for (int i = 0; i < nucleus->band2.size(); ++i)
    {
        auto spin = nucleus->band2.at(i).spin;
        std::cout << spin << " " << EnergyFormula::wobblingBand(spin, 116, 1, 76, 0) << "\n";
    }
 */
    delete nucleus;
    delete chisquared;
    delete formulas;
    delete paramSet;
}