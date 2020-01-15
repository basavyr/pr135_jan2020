#ifndef ENERGYEXPRESSIONS_HH
#define ENERGYEXPRESSIONS_HH

#include <vector>


class EnergyFormula
{
public:
    //N=0 wobbling band
    double yrastBand(double params);

    //N=1 wobbling band
    double wobblingBand(double params);
};

#endif // ENERGYEXPRESSIONS_HH
