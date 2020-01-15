#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "../include/pr135.hh"
#include "../include/data.hh"
#include "../include/energyExpressions.hh"
int main()
{
    //ENSDF DATA
    Pr135Experimental *nucleu = new Pr135Experimental;

    nucleu->init_ENSDF(*nucleu);
    std::vector<Pr135Experimental::band> v1 = nucleu->data1Exp(*nucleu);
    Pr135Experimental::printer(v1);

    std::cout << "***************";
    Pr135Experimental::newLine();

    std::vector<Pr135Experimental::band> v2 = nucleu->data2Exp(*nucleu);
    Pr135Experimental::printer(v2);
    delete nucleu;
    v1.clear();
    v2.clear();
    return 0;
}