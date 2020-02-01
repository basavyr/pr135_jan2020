#include "../include/energyExpressions.hh"
#include "../include/data.hh"
#include <cmath>

double EnergyFormula::yrastBand(double spin, double a1, double a2, double a3, double theta)
{
    auto retval = energyExpression(0, spin, a1, a2, a3, theta);
    auto e0 = energyExpression(0, 5.5, a1, a2, a3, theta);
    if (retval && e0 != 6969)
        return retval - e0;
    //avoid showing energy which containes complex number (e0 can be complex!)
    return 6969;
}

double EnergyFormula::wobblingBand(double spin, double a1, double a2, double a3, double theta)
{
    //wobbling band has also n=0
    auto retval = energyExpression(0, spin, a1, a2, a3, theta);
    auto e0 = energyExpression(0, 5.5, a1, a2, a3, theta);
    if (retval && e0 != 6969)
        return retval - e0;
    //avoid showing energy which containes complex number (e0 can be complex!)
    return 6969;
}

double EnergyFormula::j_Component(int n, double theta)
{
    if (n == 1)
        return Constants::oddSpin * cos(theta * Constants::PI / 180.0);
    return Constants::oddSpin * sin(theta * Constants::PI / 180.0);
}

double EnergyFormula::energyExpression(int N, double spin, double a1, double a2, double a3, double theta)
{

    //stop calculus if the moments of inertia are NULL
    if (!a1 || !a2 || !a3)
        return 6969;

    //checking for complex numers
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

    auto case1 = 6969 * 0.5;
    auto case2 = 6969 * 1.5;
    if (term2 == case1 || term2 == case2)
    {
        std::cout << "found complex number at..." << spin << " and params { " << a1 << " " << a2 << " " << a3 << " " << theta << " } "
                  << "\n";
        return 6969;
    }

    auto sum = A1 * j1 * j1 + A2 * j2 * j2;
    auto retval = static_cast<double>(term1 + term2 + sum);
    if (!isnan(retval))
        return retval;
    return 6969;
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

void EnergyFormula::parametricPrinter(Pr135Experimental &nucleus, double i1, double i2, double i3, double theta)
{
    std::vector<double> theoreticalEnergies;
    auto size = nucleus.band1.size();
    for (int i = 0; i < size; ++i)
    {
        auto spin = nucleus.band1.at(i).spin;
        theoreticalEnergies.emplace_back(EnergyFormula::yrastBand(spin, i1, i2, i3, theta));
    }
    size = nucleus.band2.size();
    for (int i = 0; i < size; ++i)
    {
        auto spin = nucleus.band2.at(i).spin;
        theoreticalEnergies.emplace_back(EnergyFormula::wobblingBand(spin, i1, i2, i3, theta));
    }
    size = nucleus.band1.size();
    auto index = 0;

    //classic printing
    /* for (int i = 0; i < size; ++i)
    {
        std::cout << nucleus.band1.at(i).spin << " " << nucleus.band1.at(i).energy << " " << theoreticalEnergies.at(index);
        std::cout << "\n";
        index++;
    }
    size = nucleus.band2.size();
    for (int i = 0; i < size; ++i)
    {
        std::cout << nucleus.band2.at(i).spin << " " << nucleus.band2.at(i).energy << " " << theoreticalEnergies.at(index);
        std::cout << "\n";
        index++;
    } */
    mathmematicaPrinter(nucleus, theoreticalEnergies);
    std::cout << std::endl;
}

double EnergyFormula::energy_ZeroTheta(int n, double spin, double a1, double a2, double a3, double theta)
{
    theta = static_cast<double>(0);

    // ######################################
    //don't have to transform the moments of inertia since the energy function takes care of that by default

    //taking care of the moments of inertia
    //convert the moments of inertia in inertia factors by the conformal transformation
    // a1 = inertiaFactor(a1);
    // a2 = inertiaFactor(a2);
    // a3 = inertiaFactor(a3);

    //#########################################
    // don't have to generate the compoents of the ingle particle angular momentun since the energy function takes care of that by default

    //generate the components of the single particle angular momentun
    // auto j1 = j_Component(1, theta);
    // auto j2 = j_Component(2, theta);

    auto result = energyExpression(n, spin, a1, a2, a3, theta);
    return static_cast<double>(result);
}

double EnergyFormula::energy_minHamiltonian(double spin, double a1, double a2, double a3, double theta)
{

    //taking care of the moments of inertia
    //convert the moments of inertia in inertia factors by the conformal transformation
    a1 = inertiaFactor(a1);
    a2 = inertiaFactor(a2);
    a3 = inertiaFactor(a3);

    //generate the components of the single particle angular momentun
    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);

    auto minimalTerm = a1 * spin * spin + (2.0 * spin + 1.0) * a1 * j1 - a2 * j2 * spin;
    auto result = minimalTerm;
    return static_cast<double>(result);
}

double EnergyFormula::energy_justOmega(int n, double spin, double a1, double a2, double a3, double theta)
{
    auto result = 1;
    return static_cast<double>(result);
}

EnergyFormula::triplet EnergyFormula::singleParticleSum(double theta, double a1, double a2, double a3)
{
    // auto x = new triplet;
    triplet result;
    a1 = inertiaFactor(a1);
    a2 = inertiaFactor(a2);
    a3 = inertiaFactor(a3);

    auto j1 = j_Component(1, theta);
    auto j2 = j_Component(2, theta);

    result.firstSumComponent = (a1 * j1 * j1);
    result.secondSumComponent = (a2 * j2 * j2);
    result.totalsum = result.totalSum();
    return result;
}

void EnergyFormula::mathmematicaPrinter(Pr135Experimental &nucleus, std::vector<double> &vec)
{
    auto index = 0;
    //spins
    // mathematicaSpinPrinter(nucleus);
    //energies
    /*  std::cout << "exp1= { ";
    for (int i = 0; i < nucleus.band1.size(); ++i)
    {
        if (i == nucleus.band1.size() - 1)
        {
            std::cout << nucleus.band1.at(i).energy;
        }
        else
        {
            std::cout << nucleus.band1.at(i).energy << " , ";
        }
    }
    std::cout << " };"
              << "\n"; */
    std::cout << "th1= { ";
    for (int i = 0; i < nucleus.band1.size(); ++i)
    {
        if (i == nucleus.band1.size() - 1)
        {
            std::cout << vec.at(index);
            index++;
        }
        else
        {
            std::cout << vec.at(index) << " , ";
            index++;
        }
    }
    std::cout << " };"
              << "\n";
    /*  std::cout << "exp2= { ";
    for (int i = 0; i < nucleus.band2.size(); ++i)
    {
        if (i == nucleus.band2.size() - 1)
        {
            std::cout << nucleus.band2.at(i).energy;
        }
        else
        {
            std::cout << nucleus.band2.at(i).energy << " , ";
        }
    }
    std::cout << " };"
              << "\n"; */
    std::cout << "th2= { ";
    for (int i = 0; i < nucleus.band2.size(); ++i)
    {
        if (i == nucleus.band2.size() - 1)
        {
            std::cout << vec.at(index);
            index++;
        }
        else
        {
            std::cout << vec.at(index) << " , ";
            index++;
        }
    }
    std::cout << " };"
              << "\n";
}

void EnergyFormula::mathematicaSpinPrinter(Pr135Experimental &nucleus)
{
    std::cout << "spin1= { ";
    for (auto i = 0; i < nucleus.band1.size(); ++i)
    {
        if (i == nucleus.band1.size() - 1)
        {
            std::cout << nucleus.band1.at(i).spin;
        }
        else
        {
            std::cout << nucleus.band1.at(i).spin << " , ";
        }
    }
    std::cout << " };"
              << "\n";
    std::cout << "spin2= { ";
    for (auto i = 0; i < nucleus.band2.size(); ++i)
    {
        if (i == nucleus.band2.size() - 1)
        {
            std::cout << nucleus.band2.at(i).spin;
        }
        else
        {
            std::cout << nucleus.band2.at(i).spin << " , ";
        }
    }
    std::cout << " };"
              << "\n";
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