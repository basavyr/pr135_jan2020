#ifndef PR135_HH
#define PR135_HH

#include <vector>

class Pr135Experimental
{
public:
    void init_ENSDF(Pr135Experimental &);
    void init_MATTA(Pr135Experimental &);
    static const int dim1 = 11;
    static const int dim2 = 5;

public:
    struct band
    {
        double spin, energy;
    };
    std::vector<band> exp1;
    std::vector<band> exp2;
    static void printer(std::vector<band> &);
    static void newLine();
    std::vector<band> data1Exp(Pr135Experimental &);
    std::vector<band> data2Exp(Pr135Experimental &);
};

#endif // PR135_HH
