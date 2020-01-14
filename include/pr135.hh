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

    //containers to store the experimental data (spin and energy)
    std::vector<band> exp1;
    //containers to store the experimental data (spin and energy)
    std::vector<band> exp2;

    std::vector<band> data1Exp(Pr135Experimental &);
    std::vector<band> data2Exp(Pr135Experimental &);

    //additional methods:
public:
    static void printer(std::vector<band> &);
    static void newLine();
};

#endif // PR135_HH
