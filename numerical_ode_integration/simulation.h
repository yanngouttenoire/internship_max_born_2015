#ifndef SIMULATION_H
#define SIMULATION_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string>
#include"solve.h" 
#include"ode.h" 
#include"coulomb.h" 
#include"ic.h"
#include"plot.h"
#include"simulation.h" 

class Simulation
{
public:
Simulation(Ode* m_ode);
~Simulation();

//We create a method which set the ode
void setOde(Ode* m_ode);

//We create a method which run the numerical integration
void run(std::string algorithm, double step, unsigned int ptsNbr);

double* q;
double* t;
double step;
unsigned int ptsNbr;

//We create pointers for all the others classes we use in Simulation -> That is the spirit of object-oriented programming 
//A big advantage: it is easier for the second user to understand the architecture of the code because the only read of .h can be sufficient
Ode* m_ode;
Ic* m_ic;
Solve* m_solve;
Plot* m_plot;

};

#endif
