#ifndef PR135_HH
#define PR135_HH

#include <vector>

class Pr135Experimental
{
public:
    //create (init) the experimental containers with ENSDF values
    void init_ENSDF(Pr135Experimental &);

    //create (init) the experimental containers with Matta article's values
    void init_MATTA(Pr135Experimental &);

public:
    //data structure to be created
    struct band
    {
        double spin, energy;
    };

    struct reducedProbabilities
    {
        double BE2out_BE2in;
        double BM1OUT_BE2in;
        double delta;
    };

    //containers to store the experimental data (spin and energy)
    std::vector<band> band1;
    //containers to store the experimental data (spin and energy)
    std::vector<band> band2;
    //containers to store the experimental data for electromagnetic transitions;
    std::vector<reducedProbabilities> transitions;

    std::vector<band> data1Exp(Pr135Experimental &);
    std::vector<band> data2Exp(Pr135Experimental &);

    //additional methods:
public:
    static void printer(std::vector<band> &);
    static void newLine();
};

class Pr135Theoretical
{
public:
    struct band
    {
        double spin, energy;
    };
    std::vector<band> band1;
    std::vector<band> band2;
};

#endif // PR135_HH
