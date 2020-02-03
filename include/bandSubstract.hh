#ifndef BANDSUBSTRACT_HH
#define BANDSUBSTRACT_HH

#include <vector>
#include <chrono>
#include <ctime>
#include <string>
#include <utility>
#include <memory>
#include <cmath>
#include <fstream>

class Data_ENSDF
{
public:
    double substractValue;
    //class constructor and destructor
public:
    //class constructor and destructor
    // Data_ENSDF();

    //class constructor and destructor
    ~Data_ENSDF();

    Data_ENSDF(double substraction);
    //containers to hold the experimental data for pr 135
    //containers hold the normalized energies
    //first band is excluded from the fitting procedure!
public:
    std::vector<double> spin1 = {7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    std::vector<double> spin2 = {8.5, 10.5, 12.5, 14.5, 16.5};
    std::vector<double> yrastExp = {0.37277, 1.03284, 1.88684, 2.88674, 3.96174, 4.80574, 5.63974,
                                    6.52174, 7.44474, 8.35794};
    std::vector<double> wobbExp = {1.07524, 1.80044, 2.64294, 3.59894, 4.71094};

    //containers to hold the substracted energies
    std::vector<double> yrastExp_Sub;
    std::vector<double> wobExp_Sub;

    //containers for storing the theoretical energies (obtained after the best root mean square)
    std::vector<double> yrastTh;
    std::vector<double> wobblingTh;
};

class Data_MATTA
{
    //containers to hold the experimental data for pr 135
    //containers hold the normalized energies
    //first band is excluded from the fitting procedure
public:
    std::vector<double> spin1 = {};
    std::vector<double> spin2 = {};
    std::vector<double> yrastExp = {};
    std::vector<double> wobbExp = {};

    //containers to hold the substracted energies
    std::vector<double> yrastExp_Sub;
    std::vector<double> wobExp_Sub;

    //containers for storing the theoretical energies (obtained after the best root mean square)
    std::vector<double> yrastTh;
    std::vector<double> wobblingTh;
};

class BandSubstract
{

    //the useful constants
public:
    static constexpr double j_singleParticle = 5.5;
    static constexpr double PI = 3.14159265358979;

    //the useful formulae for calculus
public:
    static void testApp_Substraction();
    static void smartPointerTest(double, std::ofstream &);

    template <typename T>
    static void arrayPrinter(std::vector<T> &array)
    {
        for (int i = 0; i < array.size(); ++i)
        {
            if (i == array.size() - 1)
                std::cout << array.at(i) << "\n";
            else
            {
                std::cout << array.at(i) << " , ";
            }
        }
    }

    template <typename T>
    static void generatePlotData(const std::string &filename, std::vector<T> &spin, std::vector<T> &energy, std::vector<T> &energyAdjusted)
    {
        std::ofstream fout(filename);
        for (int i = 0; i < spin.size(); ++i)
        {
            fout << spin.at(i) << " " << energy.at(i) << " " << energyAdjusted.at(i) << "\n";
        }
    }

    template <typename T>
    static void bandSubstracter(std::vector<T> &object, std::vector<T> &newObject, double subtractor)
    {
        auto halfSize = static_cast<size_t>(object.size() / 2);
        for (int i = 0; i < object.size(); ++i)
        {
            if (i < halfSize)
            {
                newObject.emplace_back(object.at(i) - subtractor);
            }
            else
            {
                newObject.emplace_back(object.at(i));
            }
        }

        // std::cout << "the half band is of length= " << halfSize;
        // std::cout << "\n";
    }

    //generate the full data set for RMS calculation
    //prepare the data for FIT
    struct rms_Tuple
    {
        double rms_substracted;
        double rms_pure;
    };
    static rms_Tuple bandGeneration(Data_ENSDF &, double, double, double, double);
};

class EnergyFormulae
{
public:
    //calculates the obbling frequency of the nucleus
    static double omega(double spin, double i1, double i2, double i3, double theta);

    //calculates the wobbling frequency of the nucleus without the single particle term in the total hamiltonian
    static double omegaPrime(double spin, double i1, double i2, double i3, double theta);

    //returns the energy calculation of an arbitrary spin state from a band
    static double energyExpression(int N, double spin, double i1, double i2, double i3, double theta);

    //tranforms the moments of inertia into inertia factors A_k
    static double inertiaFactor(double moi);

    //computes the components of the single particle angular momentum
    static double j_Component(int k, double theta);

    //returns the minimal term from the hamiltonian: term independent on the wobbling frequency
    static double minHamiltonian(double spin, double i1, double i2, double i3, double theta);
};

class RMS_Calculus
{
public:
    //the set of parameters
    struct minParamSet
    {
        double i1_min, i2_min, i3_min;
        double theta_min;
        //the rms calculated with the substracted energies
        double rms;
        //rms calculated with the pure experimental energies
        double pure_rms;
    };

    //the set of limits for fitting params (moments of inertia and the single particle angle)
    struct paramsLimits
    {
        const double i_left = 1.0;
        const double i_right = 120.0;
        const double i_step = 1.0;
        const double theta_left = 0.0;
        const double theta_right = 180.0;
        const double theta_step = 1.0;
    };
    std::vector<minParamSet> minimalParameters;

public:
    //rms calculation for any arbitrary energy containers
    template <typename T>
    static double pure_RootMeanSquare(std::vector<T> &exp, std::vector<T> &th)
    {
        return rootMeanSquare<T>(&exp, &th);
    }

    //root mean square error calculation for the SUBSTRACTED energy containers
    template <typename T>
    static double rootMeanSquare(std::vector<T> &exp, std::vector<T> &th)
    {
        if (exp.size() != th.size())
            return 6969;
        double sum = 0;
        int dimCheck = 0;
        for (int i = 0; i < exp.size(); ++i)
        {
            auto currentValue = pow(exp.at(i) - th.at(i), 2);
            if (!isnan(currentValue))
            {
                dimCheck++;
                sum += currentValue;
            }
        }
        auto rms = static_cast<double>(sqrt(1.0 / exp.size() * sum));
        if (!isnan(rms) && dimCheck == exp.size())
            return rms;
        return 6969;
    }

    //search for the minimal set of parameters after the substraction of the bands with a fixed quantity.
    template <typename T>
    static void searchMinimum(T &object, minParamSet &bestParams)
    {
        auto limits = std::make_unique<paramsLimits>();
        std::vector<minParamSet> minSetOfParams;
        std::vector<double> RMS_stack;
        int index = 0;
        for (auto I1 = limits->i_left; I1 < limits->i_right; I1 += limits->i_step)
        {
            for (auto I2 = limits->i_left; I2 < limits->i_right; I2 += limits->i_step)
            {
                for (auto I3 = limits->i_left; I3 < limits->i_right; I3 += limits->i_step)
                {
                    for (auto theta = limits->theta_left; theta < limits->theta_right; theta += limits->theta_step)
                    {
                        auto currentRMS = BandSubstract::bandGeneration(object, I1, I2, I3, theta);
                        // if (!isnan(currentRMS) && currentRMS != 6969)
                        {
                            minSetOfParams.emplace_back(minParamSet());
                            minSetOfParams.at(index).i1_min = I1;
                            minSetOfParams.at(index).i2_min = I2;
                            minSetOfParams.at(index).i3_min = I3;
                            minSetOfParams.at(index).theta_min = theta;
                            minSetOfParams.at(index).pure_rms = currentRMS.rms_pure;
                            RMS_stack.emplace_back(currentRMS.rms_substracted);
                            index++;
                        }
                    }
                }
            }
        }
        auto minIndex = std::distance(RMS_stack.begin(), std::min_element(RMS_stack.begin(), RMS_stack.end()));
        bestParams.i1_min = minSetOfParams.at(minIndex).i1_min;
        bestParams.i2_min = minSetOfParams.at(minIndex).i2_min;
        bestParams.i3_min = minSetOfParams.at(minIndex).i3_min;
        bestParams.theta_min = minSetOfParams.at(minIndex).theta_min;
        bestParams.pure_rms = minSetOfParams.at(minIndex).pure_rms;
        bestParams.rms = RMS_stack.at(minIndex);
    }

    template <typename T>
    static void paramPrinter(std::ofstream &out, T &params)
    {
        std::cout << "The best fit parameters are"
                  << "\n";
        out << "The best fit parameters are"
            << "\n";
        std::cout << "I1= " << params.i1_min << " I2= " << params.i2_min << " I3= " << params.i3_min << " th= " << params.theta_min << "\n";
        out << "I1= " << params.i1_min << " I2= " << params.i2_min << " I3= " << params.i3_min << " th= " << params.theta_min << "\n";
        std::cout << "E_RMS= " << params.rms << "\n";
        out << "E_RMS= " << params.rms << "\n";
        std::cout << "E_RMS_pure= " << params.pure_rms << "\n";
        out << "E_RMS_pure= " << params.pure_rms << "\n";
    }
};

#endif // BANDSUBSTRACT_HH
