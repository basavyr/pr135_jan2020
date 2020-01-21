#include "../include/energyExpressions.hh"
#include "../include/data.hh"
#include <cmath>

double EnergyFormula::yrastBand(double spin, double a1, double a2, double a3, double theta)
{
    auto retval = energyExpression(0, spin, a1, a2, a3, theta);
    if (retval)
        return retval;
    return 0;
}

double EnergyFormula::wobblingBand(double spin, double a1, double a2, double a3, double theta)
{
    auto retval = energyExpression(1, spin, a1, a2, a3, theta);
    if (retval)
        return retval;
    return 0;
}

double EnergyFormula::j_Component(int n, double theta)
{
    if (n == 1)
        return Constants::oddSpin * cos(theta * Constants::PI / 180.0);
    return Constants::oddSpin * sin(theta * Constants::PI / 180.0);
}

double EnergyFormula::energyExpression(int N, double spin, double a1, double a2, double a3, double theta)
{

    //stop calculus if the moments of inerta are NULL
    if (!a1 || !a2 || !a3)
        return 6969;
    if (omega(spin, a1, a2, a3, theta) == 6969)
        return 6969;

    auto j = Constants::oddSpin;
    // auto j1 = j * cos(theta * Constants::PI / 180.0);
    // auto j2 = j * sin(theta * Constants::PI / 180.0);

    //j component functions
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);

    //get the INERTIA FACTORS FROM THE MOMENTS OF INERTIA
    auto A1 = EnergyFormula::inertiaFactor(a1);
    auto A2 = EnergyFormula::inertiaFactor(a2);
    auto A3 = EnergyFormula::inertiaFactor(a3);

    auto term1 = A1 * spin * spin + (2.0 * spin + 1.0) * A1 * j1 - spin * A2 * j2;
    auto term2 = omega(spin, a1, a2, a3, theta) * (N + 0.5);

    //checking for complex numers
    if (!term2)
        std::cout << "found complex number at..." << spin << " and params { " << a1 << " " << a2 << " " << a3 << " " << theta << " } "
                  << "\n";

    auto sum = A1 * j1 * j1 + A2 * j2 * j2;
    auto retval = static_cast<double>(term1 + term2 + sum);
    if (!isnan(retval))
        return retval;
    return 0;
}

double EnergyFormula::inertiaFactor(double a)
{
    if (a)
        return 1.0 / (2.0 * a);
    return 0;
}

double EnergyFormula::omega(double spin, double a1, double a2, double a3, double theta)
{
    auto j = Constants::oddSpin;

    // auto j1 = j * sin(theta * Constants::PI / 180.0);
    // auto j2 = j * cos(theta * Constants::PI / 180.0);

    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);

    //get the INERTIA FACTORS FROM THE MOMENTS OF INERTIA
    auto A1 = EnergyFormula::inertiaFactor(a1);
    auto A2 = EnergyFormula::inertiaFactor(a2);
    auto A3 = EnergyFormula::inertiaFactor(a3);

    auto term1 = (2.0 * spin + 1.0) * (A2 - A1 - static_cast<double>((A2 * j2) / spin)) - 2.0 * A1 * j1;
    auto term2 = (2.0 * spin + 1.0) * (A3 - A1) - 2.0 * A1 * j1;
    auto term3 = (A3 - A1) * (A2 - A1 - static_cast<double>((A2 * j2) / spin));

    // std::cout << "\n";
    // std::cout << "in omega... " << term1 << " " << term2 << " " << term3 << "\n";

    auto retval = static_cast<double>(sqrt(term1 * term2 - term3));
    //check if the return value is actually real or not wobbling frequency
    if (!isnan(retval))
        return retval;
    return 6969;
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