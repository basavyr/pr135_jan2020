#ifndef ENERGYEXPRESSIONS_HH
#define ENERGYEXPRESSIONS_HH

#include <vector>
#include <cmath>
#include "../include/pr135.hh"

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
        T retval = static_cast<T>(0);
        for (int i = 0; i < exp.size(); ++i)
        {
            auto elem = pow(exp.at(i) - th.at(i), 2);
            if (isValid(elem))
            {
                retval += elem;
            }
        }
        return sqrt(retval / exp.size());
    }

    template <typename T>
    T meanSquaredError(std::vector<Pr135Experimental::band> &exp, std::vector<Pr135Theoretical::band> &th)
    {
        return 0;
    }

    template <typename T>
    T apply_Yrast(Pr135Experimental &nucleus)
    {
        auto yrast = [](auto spin, auto param) {
            return EnergyFormula::yrastBand(static_cast<double>(spin), static_cast<double>(param));
        };
        return yrast(nucleus.band1.at(0).spin, 1);
    }

    template <typename T>
    T apply_Wobbling(Pr135Experimental &nucleus)
    {
        auto wobbling = [](auto spin, auto param) {
            return EnergyFormula::wobblingBand(static_cast<double>(spin), static_cast<double>(param));
        };
        return wobbling(nucleus.band2.at(0).spin, 1);
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
    template <typename T>
    T calculateMinimumValue(Pr135Experimental &nucleus, EnergyFormula &formulas)
    {
        auto retval = 0;
        return retval;
    }
};

#endif // ENERGYEXPRESSIONS_HH
