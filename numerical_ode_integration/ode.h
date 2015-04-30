#ifndef ODE_H
#define ODE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

//Abstract class which can contains Ordinary differential system of equations

class Ode
{

public:
Ode();
~Ode();

//We define a virtual method which will contains the ode
virtual void diff(double* q, double t)=0;

//We define the array which will contains the first derivative of the ode
double *qp;

};

#endif
