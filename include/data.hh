#ifndef DATA_HH
#define DATA_HH

//constants for Pr135 nucleus (experimental data)
class Constants
{
public:
    static constexpr double beta = 0.17;
    //number of states of the yrast band (nw=0)
    static constexpr int dim1 = 11;
    //number of states of the first wobbling band (nw=1)
    static constexpr int dim2 = 5;
    static constexpr double gamma = 26.0;
    //number of states that have exp values for transition probabilities
    static constexpr int noTransitions = 3;
    static constexpr double oddSpin = 5.5;
    static constexpr double PI = 3.14159265358979;
};

//experimental data for the Pr135 energies (taken from ENSDF data)
class ExperimentalData_ENSDF
{
public:
    static constexpr double spin1[Constants::dim1] = {5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    static constexpr double spin2[Constants::dim2] = {8.5, 10.5, 12.5, 14.5, 16.5};
    static constexpr double energy1[Constants::dim1] = {358.06, 730.83, 1390.9, 2244.9, 3244.8, 4319.8, 5163.8, 5997.8, 6879.8, 7802.8, 8716};
    static constexpr double energy2[Constants::dim2] = {1433.3, 2158.5, 3001.0, 3957.0, 5069.0};
};

//experimental data for the Pr135 energies (taken from MATTA et al 2015 paper)
class ExperimentalData_MATTA
{
public:
    static constexpr double spin1[Constants::dim1] = {5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5};
    static constexpr double spin2[Constants::dim2] = {8.5, 10.5, 12.5, 14.5, 16.5};
    static constexpr double energy1[Constants::dim1] = {358., 730.8, 1391, 2245.2, 3245.3, 4320.8, 5164.8, 5998.7, 6880.4, 7803.3, 8760.3};
    static constexpr double energy2[Constants::dim2] = {1477.8, 2203.8, 2999.8, 3955.4, 5067.4};
    static void show();
};

class EM_Transitions
{
public:
    static constexpr double spin[Constants::noTransitions] = {10.5, 12.5, 14.5};
    static constexpr double BE2out_BE2in[Constants::noTransitions] = {0.843, 0.500, 0.261};
    static constexpr double BM1out_BE2in[Constants::noTransitions] = {0.164, 0.035, 0.016};
    static constexpr double delta_BE2out_BE2in[Constants::noTransitions] = {0.032, 0.025, 0.014};
    static constexpr double delta_BM1out_BE2in[Constants::noTransitions] = {0.014, 0.009, 0.004};
};

#endif // DATA_HH
