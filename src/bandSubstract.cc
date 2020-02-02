#include <iostream>

#include "../include/bandSubstract.hh"

void BandSubstract::testApp_Substraction()
{
    auto now = std::chrono::system_clock::now();
    std::time_t datenow = std::chrono::system_clock::to_time_t(now);
    std::cout << ctime(&datenow);
    std::cout << "App works";
}

//testing the smart pointers

void BandSubstract::smartPointerTest()
{
    auto smartPtr = std::make_unique<Data_ENSDF>();
}

Data_ENSDF::Data_ENSDF()
{
    std::string constructorMessage = "The class for storing the ENSDF experimental data is succesfully created";
    std::cout << constructorMessage;
    std::cout << std::endl;
}

Data_ENSDF::~Data_ENSDF()
{
    std::string destructorMessage = "The class container has been succsesfully destroyed";
    std::cout << destructorMessage;
    std::cout << std::endl;
}