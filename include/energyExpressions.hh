#ifndef ENERGYEXPRESSIONS_HH
#define ENERGYEXPRESSIONS_HH

#include <vector>
#include <cmath>
#include "../include/pr135.hh"
#include "../include/data.hh"
#include <iostream>
#include <chrono>
#include <utility>

class EnergyFormula
{
public:
    struct triplet
    {
        double firstSumComponent;
        double secondSumComponent;
        double totalsum;
        double totalSum()
        {
            return firstSumComponent + secondSumComponent;
        }
    };

public:
    //N=0 wobbling band
    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double yrastBand(double spin, double a1, double a2, double a3, double theta);

    //N=1 wobbling band
    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double wobblingBand(double spin, double a1, double a2, double a3, double theta);

    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double omega(double spin, double a1, double a2, double a3, double theta);

    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double energyExpression(int N, double spin, double a1, double a2, double a3, double theta);

    static double j_Component(int n, double theta);

    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double inertiaFactor(double a);

    //print the energies (exp and th) as a function of the best fit parameters
    static void parametricPrinter(Pr135Experimental &nucleus, double i1, double i2, double i3, double theta);

    //print as mathematica format
    static void mathmematicaPrinter(Pr135Experimental &nucleus, std::vector<double> &vec);

    //print the spins as mathematica format
    static void mathematicaSpinPrinter(Pr135Experimental &nucleus);

    //the minimum hamiltonian. ingore the wobbling escitations and the particle alignment
    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double energy_minHamiltonian(double spin, double a1, double a2, double a3, double theta);

    //energy expression for theta equal to zero => no single particle alignment
    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double energy_ZeroTheta(int n, double spin, double a1, double a2, double a3, double theta);

    //energy expression for theta equal to zero => no single particle alignment
    //a1,a2,a3 are the MOMENTS OF INERTIA (function depends on the inertia factor though!)
    static double energy_justOmega(int n, double spin, double a1, double a2, double a3, double theta);

    //checks the single particle sum from the total Hamiltonian
    static triplet singleParticleSum(double a1, double a2, double a3, double theta);

    //energy expression for

    //normalize to first state in yrast band
    template <typename T>
    std::vector<T> normalize(std::vector<Pr135Experimental::band> &vec)
    {
        std::vector<T> retval;
        // auto zero = NULL;
        // retval.emplace_back(static_cast<T>(zero));
        for (int i = 0; i < vec.size(); ++i)
        {
            auto diff = vec.at(i).energy - vec.at(0).energy;
            retval.emplace_back(static_cast<T>(diff));
        }
        return retval;
    }

    //**************************
    //WRONG NORMALIZE: must substract only E0_yrast
    //**************************
    void normalize(std::vector<Pr135Experimental::band> &vec)
    {
        auto e0 = vec.at(0).energy;
        for (int i = 0; i < vec.size(); ++i)
        {
            vec.at(i).energy = vec.at(i).energy - e0;
        }
    }

    void normalizeBands(Pr135Experimental &obj)
    {
        auto size = obj.band1.size();
        auto e0 = obj.band1.at(0).energy;
        for (int i = 0; i < size; ++i)
        {
            obj.band1.at(i).energy = obj.band1.at(i).energy - e0;
        }
        size = obj.band2.size();
        for (int i = 0; i < size; ++i)
        {
            obj.band2.at(i).energy = obj.band2.at(i).energy - e0;
        }
    }

    template <typename T>
    void kevTOmev(std::vector<T> &vec)
    {
        auto mev = 1000;
        for (int i = 0; i < vec.size(); ++i)
        {
            vec.at(i) = static_cast<T>(vec.at(i) / mev);
        }
    }
    template <typename T>
    void kevTOmevBand(std::vector<Pr135Experimental::band> &vec)
    {
        auto mev = static_cast<T>(1000);
        for (int i = 0; i < vec.size(); ++i)
        {
            vec.at(i).energy = static_cast<T>(vec.at(i).energy / mev);
        }
    }

    template <typename T>
    static T testFunction(T spin, T param1)
    {
        T retval;
        retval = 0.5 * spin + param1 * spin;
        return retval;
    }
    struct ParameterSet
    {
        // double A1, A2, A3;
        const double A_left = 1.0;
        const double A_right = 120.0;
        const double A_step = 1;
        const double theta_left = 1.0;
        const double theta_right = 180.0;
        const double theta_step = 1.0;
        // double theta;
    };
};

class ChiSquared
{
public:
    template <typename T>
    T meanSquaredError(std::vector<T> &exp, std::vector<T> &th)
    {
        if (exp.size() != th.size())
            return 987654321;
        auto n = exp.size();
        T retval = static_cast<T>(0);
        for (int i = 0; i < exp.size(); ++i)
        {
            auto elem = pow(exp.at(i) - th.at(i), 2);
            if (isValid(elem))
            {
                retval += elem;
                n--;
            }
        }
        if (!n)
            return sqrt(retval / exp.size());
        return 987654321;
    }

    //generate the two containers for computing the RMS.
    //experimental and theoretical containers with the excitation energies
    //returns the RMS value provided by the rms function given the two energy containers
    template <typename T>
    T applyEnergies(Pr135Experimental &nucleus, double a1, double a2, double a3, double theta)
    {
        std::vector<T> dataExp;
        std::vector<T> dataTh;
        auto yrast = [&](auto spin) {
            return EnergyFormula::yrastBand(static_cast<double>(spin), a1, a2, a3, theta);
        };
        auto wobbling = [&](auto spin) {
            return EnergyFormula::wobblingBand(static_cast<double>(spin), a1, a2, a3, theta);
        };

        //add the exp data and th data for first band
        for (int i = 0; i < nucleus.band1.size(); ++i)
        {
            dataExp.emplace_back(static_cast<T>(nucleus.band1.at(i).energy));
            auto theoreticalEnergy = yrast(nucleus.band1.at(i).spin);
            if (theoreticalEnergy == 6969)
            {
                // std::cout << "energy (yrast) value not valid!..." << nucleus.band1.at(i).spin << " " << a1 << " " << a2 << " " << a3 << " " << theta << " } \n";
                // break;
            }
            else // (theoreticalEnergy != 6969)
            {
                dataTh.emplace_back(static_cast<T>(theoreticalEnergy));
            }
        }

        //add the exp data and th data for the second
        for (int i = 0; i < nucleus.band2.size(); ++i)
        {
            dataExp.emplace_back(static_cast<T>(nucleus.band2.at(i).energy));
            auto theoreticalEnergy = wobbling(nucleus.band2.at(i).spin);
            if (theoreticalEnergy == 6969)
            {
                // std::cout << "energy (wobbling) value not valid!..." << nucleus.band2.at(i).spin << " " << a1 << " " << a2 << " " << a3 << " " << theta << " } \n";
                // break;
            }
            else //(theoreticalEnergy != 6969)
            {
                dataTh.emplace_back(static_cast<T>(theoreticalEnergy));
            }
        }

        auto chiVal = meanSquaredError(dataExp, dataTh);
        return static_cast<T>(chiVal);
    }

    //generates a pair of containers with experimental and theoretical data for pr135
    //uses only the pure energies (un-normalized)
    //combines both bands into one dataset
    //returns the RMS value provided by the rms function given the two energy containers
    template <typename T>
    T applyRawEnergies(double i1, double i2, double i3, double theta)
    {
        std::vector<T> dataExp;
        std::vector<T> dataTh;
        //populates the containers with anydata
        auto filler = [](std::vector<T> &v, auto y) {
            v.emplace_back(std::move(y));
        };

        for (int i = 0; i < Constants::dim1; ++i)
        {
            filler(dataExp, static_cast<T>(ExperimentalData_ENSDF::energy1[i] / 1000.0));
            auto theoreticalEnergy = EnergyFormula::energyExpression(0, ExperimentalData_ENSDF::spin1[i], i1, i2, i3, theta);
            if (theoreticalEnergy != 6969)
                filler(dataTh, theoreticalEnergy);
        }
        for (int i = 0; i < Constants::dim2; ++i)
        {
            filler(dataExp, static_cast<T>(ExperimentalData_ENSDF::energy2[i] / 1000.0));
            auto theoreticalEnergy = EnergyFormula::energyExpression(0, ExperimentalData_ENSDF::spin2[i], i1, i2, i3, theta);
            if (theoreticalEnergy != 6969)
                filler(dataTh, theoreticalEnergy);
        }
        auto rmsValue = meanSquaredError(dataExp, dataTh);
        return static_cast<T>(rmsValue);
    }

    //generates a pair of containers with experimental and theoretical data for pr135
    //uses only the pure energies (un-normalized)
    //combines both bands into one dataset
    //returns the RMS value provided by the rms function given the two energy containers
    template <typename T>
    T applyRawEnergies_matta(double i1, double i2, double i3, double theta)
    {
        std::vector<T> dataExp;
        std::vector<T> dataTh;
        //populates the containers with anydata
        auto filler = [](std::vector<T> &v, auto y) {
            v.emplace_back(std::move(y));
        };

        for (int i = 0; i < Constants::dim1; ++i)
        {
            filler(dataExp, static_cast<T>(ExperimentalData_MATTA::energy1[i] / 1000.0));
            auto theoreticalEnergy = EnergyFormula::energyExpression(0, ExperimentalData_MATTA::spin1[i], i1, i2, i3, theta);
            if (theoreticalEnergy != 6969)
                filler(dataTh, theoreticalEnergy);
        }
        for (int i = 0; i < Constants::dim2; ++i)
        {
            filler(dataExp, static_cast<T>(ExperimentalData_MATTA::energy2[i] / 1000.0));
            auto theoreticalEnergy = EnergyFormula::energyExpression(0, ExperimentalData_MATTA::spin2[i], i1, i2, i3, theta);
            if (theoreticalEnergy != 6969)
                filler(dataTh, theoreticalEnergy);
        }
        auto rmsValue = meanSquaredError(dataExp, dataTh);
        return static_cast<T>(rmsValue);
    }

public:
    template <typename T>
    bool isValid(T arg)
    {
        if (!isnan(arg))
            return 1;
        return 10;
    }
};

class MinimumValueParameter
{
public:
    //the minimal set of params for the energy fit
    struct minimumSetOfParams
    {
        double A1_min;
        double A2_min;
        double A3_min;
        double theta_min;
        double I1_min;
        double I2_min;
        double I3_min;
    };

    template <typename T>
    T calculateMinimumValue(Pr135Experimental &nucleus, EnergyFormula &formulas, ChiSquared &chi)
    {
        auto start = std::chrono::system_clock::now();
        std::vector<T> chi_Stack;
        bool ok = 1;
        T minvalue;
        minimumSetOfParams min_set;
        EnergyFormula::ParameterSet params;
        std::vector<minimumSetOfParams> minSetOfParamsArray;
        auto minArrayIndex = 0;

        for (double I1 = params.A_left; I1 <= params.A_right && ok; I1 += params.A_step)
        {
            for (double I2 = params.A_left; I2 <= params.A_right && ok; I2 += params.A_step)
            {
                for (double I3 = params.A_left; I3 <= params.A_right && ok; I3 += params.A_step)
                {
                    {
                        for (double theta = params.theta_left; theta <= params.theta_right && ok; theta += params.theta_step)
                        {
                            //perfornm the minimum chi2 calculus using the normalized energies
                            // auto currentChi = chi.applyEnergies<T>(nucleus, I1, I2, I3, theta);

                            //perform the minimum chi2 calculus using the raw (un-normalized) energies provided by ENSF
                            // auto currentChi = chi.applyRawEnergies<T>(I1, I2, I3, theta);

                            //perform the minimum chi2 calculus using the raw (un-normalized) energies provided by Matta et al
                            auto currentChi = chi.applyRawEnergies<T>(I1, I2, I3, theta);
                            if (isnan(currentChi))
                            {
                                break;
                                ok = 0;
                            }
                            else
                            {
                                chi_Stack.emplace_back(currentChi);
                                minSetOfParamsArray.emplace_back(minimumSetOfParams());
                                minSetOfParamsArray.at(minArrayIndex).I1_min = I1;
                                minSetOfParamsArray.at(minArrayIndex).I2_min = I2;
                                minSetOfParamsArray.at(minArrayIndex).I3_min = I3;
                                minSetOfParamsArray.at(minArrayIndex).theta_min = theta;
                                minArrayIndex++;
                            }
                        }
                    }
                }
            }
        }

        // }
        // for (auto &&i : chi_Stack)
        // {
        //     std::cout << i << " ";
        // }
        // std::cout << chi_Stack.size();
        // std::cout << std::endl;
        auto index = std::distance(chi_Stack.begin(), std::min_element(chi_Stack.begin(), chi_Stack.end()));
        minvalue = chi_Stack.at(static_cast<int>(index));
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "minvalue evaluation took " << static_cast<double>(duration / 1000.0) << " seconds...";
        std::cout << "\n";
        std::cout << "index of the minimum in the stack is " << index << "\n";
        std::cout << "Minimumm set of parameters are: I1= " << minSetOfParamsArray.at(index).I1_min << " I2= " << minSetOfParamsArray.at(index).I2_min << " I3= " << minSetOfParamsArray.at(index).I3_min << " theta= " << minSetOfParamsArray.at(index).theta_min << "\n";
        std::cout << "******** DATA ********"
                  << "\n";
        EnergyFormula::parametricPrinter(nucleus, minSetOfParamsArray.at(index).I1_min, minSetOfParamsArray.at(index).I2_min, minSetOfParamsArray.at(index).I3_min, minSetOfParamsArray.at(index).theta_min);

        return minvalue;
    }
};
#endif // ENERGYEXPRESSIONS_HH
