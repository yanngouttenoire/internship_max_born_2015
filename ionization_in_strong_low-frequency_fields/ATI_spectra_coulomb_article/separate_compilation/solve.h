#ifndef SOLVE_H
#define SOLVE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"

//Forward declaration
template<typename state_type>
class System;

template<typename state_type>
class Solve
{

public:

//We implement a method which solves an ode system using Runge-Kutta of order 5 with adaptive step-size algorithm 
//The adaptive step-size algorithm uses Runge-Kutta of order 4 embedded in RK5
//For more details: http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c16-2.pdf

void controlledRK5(System<state_type> &system, state_type& x, double &t, double &dt, double &error, double desiredErrorMin, double desiredErrorMax);


};



#endif
