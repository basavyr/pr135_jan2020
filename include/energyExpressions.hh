#ifndef ENERGYEXPRESSIONS_HH
#define ENERGYEXPRESSIONS_HH

#include <vector>
#include <cmath>
#include "../include/pr135.hh"
#include <iostream>
#include <chrono>

class EnergyFormula
{
public:
    //N=0 wobbling band
    static double yrastBand(double, double);

    //N=1 wobbling band
    static double wobblingBand(double, double);

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

    void normalize(std::vector<Pr135Experimental::band> &vec)
    {
        auto e0 = vec.at(0).energy;
        for (int i = 0; i < vec.size(); ++i)
        {
            vec.at(i).energy = vec.at(i).energy - e0;
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
};

class ChiSquared
{
public:
    template <typename T>
    T meanSquaredError(std::vector<T> &exp, std::vector<T> &th)
    {
        if (exp.size() != th.size())
            return static_cast<T>(0);
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

    template <typename T>
    T meanSquaredError(std::vector<Pr135Experimental::band> &exp, std::vector<Pr135Theoretical::band> &th)
    {
        return 0;
    }

    //sum implementation

    template <typename T>
    T apply(Pr135Experimental &nucleus, double param)
    {
        std::vector<T> dataExp;
        std::vector<T> dataTh;
        auto yrast = [&param](auto spin) {
            return EnergyFormula::yrastBand(static_cast<double>(spin), static_cast<double>(param));
        };
        auto wobbling = [&param](auto spin) {
            return EnergyFormula::wobblingBand(static_cast<double>(spin), static_cast<double>(param));
        };
        auto squaredDiff = [](auto x, auto y) {
            return pow((x - y), 2);
        };

        //add the exp data and th data for first band
        for (int i = 0; i < nucleus.band1.size(); ++i)
        {
            dataExp.emplace_back(static_cast<T>(nucleus.band1.at(i).energy));
            dataTh.emplace_back(static_cast<T>(yrast(nucleus.band1.at(i).spin)));
        }

        //add the exp data and th data for the second
        for (int i = 0; i < nucleus.band2.size(); ++i)
        {
            dataExp.emplace_back(static_cast<T>(nucleus.band2.at(i).energy));
            dataTh.emplace_back(static_cast<T>(wobbling(nucleus.band2.at(i).spin)));
        }

        auto chiVal = meanSquaredError(dataExp, dataTh);
        return static_cast<T>(chiVal);
    }

    // template <typename T>
    // T apply_Wobbling(Pr135Experimental &nucleus)
    // {
    //     auto wobbling = [](auto spin, auto param) {
    //         return EnergyFormula::wobblingBand(static_cast<double>(spin), static_cast<double>(param));
    //     };
    //     return wobbling(nucleus.band2.at(0).spin, 1);
    // }

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
    template <typename T>
    T calculateMinimumValue(Pr135Experimental &nucleus, EnergyFormula &formulas, ChiSquared &chi)
    {
        auto start = std::chrono::system_clock::now();
        std::vector<T>
            chi_Stack;
        bool ok = 1;
        T minvalue;

        //set of params for the energy fit
        struct minimumSetOfParams
        {
            double param1;
            double param1_left = -5.0;
            double param1_right = 5.0;
            double param1_step = 0.13;
        };
        minimumSetOfParams min_set;

        for (double param1 = min_set.param1_left; param1 <= min_set.param1_right && ok; param1 += min_set.param1_step)
        {
            auto currentChi = chi.apply<T>(nucleus, param1);
            if (isnan(currentChi))
            {
                break;
                ok = 0;
            }
            else
            {
                chi_Stack.emplace_back(currentChi);
            }
        } /* 
        for (auto &&i : chi_Stack)
        {
            std::cout << i << " ";
        } */
        // std::cout << std::endl;
        auto index = std::distance(chi_Stack.begin(), std::min_element(chi_Stack.begin(), chi_Stack.end()));
        minvalue = chi_Stack.at(static_cast<int>(index));
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "minvalue evaluation took " << static_cast<double>(duration / 1000) << " seconds...";
        std::cout << "\n";
        std::cout << "index of the minimum in the stack is " << index << "\n";
        std::cout << "the value of param1 is " << static_cast<double>(min_set.param1_left + (index * min_set.param1_step)) << "\n ";
        return minvalue;
    }
};

#endif // ENERGYEXPRESSIONS_HH
