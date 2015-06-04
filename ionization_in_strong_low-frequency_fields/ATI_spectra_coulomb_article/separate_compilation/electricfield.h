#ifndef ELECTRICFIELD_H
#define ELECTRICFIELD_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

//Define characteristics and properties of an eletric field
class ElectricField
{

public:

//Field parameters
double waveLenght;
double fieldAmpl;
double pulsation;
double opticalCycle;
double pulseDuration;
double phase;
double ellipticity;

//Constructor
ElectricField();

//Methods which return the force induced by the electric field on the electron on each cartesian components 
double operator()(char component, const double& t);

//Method which returns the vector potential
double vectPot(char component, const double& t); 

};

#endif
