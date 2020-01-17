#include "../include/energyExpressions.hh"

#include <cmath>

double EnergyFormula::yrastBand(double spin, double params)
{
    return params * spin * spin + spin;
}

double EnergyFormula::wobblingBand(double spin, double params)
{
    return params * spin * spin + spin;
}

// template <typename T>
// void EnergyFormula::kevTOmev(std::vector<T> &vec)
// {
//     auto mev = 1000;
//     for (int i = 0; i < vec.size(); ++i)
//     {
//         vec.at(i) = static_cast<T>(vec.at(i) / mev);
//     }
// }

// std::vector<Pr135Experimental::band> EnergyFormula::normalize(std::vector<Pr135Experimental::band> &vec)
// {
//     std::vector<Pr135Experimental::band> retval;
//     for (int i = 1; i < vec.size(); ++i)
//     {
//         retval.emplace_back(Pr135Experimental::band());
//     }
//     return retval;
// }

// template <typename T>
// bool ChiSquared::isValid(T arg)
// {
//     if (!isnan(arg))
//         return 1;
//     return 0;
//     // return !isnan(arg) ? 1 : 0;
// }

// template <typename T>
// T ChiSquared::meanSquaredError(std::vector<T> &exp, std::vector<T> &th)
// {
//     if (exp.size() != th.size())
//         return 0;
//     T retval;
//     for (int i = 0; i < exp.size(); ++i)
//     {
//         auto elem = pow(exp.at(i) - th.at(i), 2);
//         if (isValid(elem))
//         {
//             retval += elem;
//         }
//     }
//     return retval;
// }