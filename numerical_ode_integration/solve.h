#ifndef SOLVE_H
#define SOLVE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"ode.h"

//Class which can solves ordinary differential system of equations

class Solve
{
public:
Solve(double* q, double* t, Ode* m_ode);

//For instance, we create a method which solve the ode using Runge-Kutta-4
void rk4(double dt);

//We create a method which set the ode
void setOde(Ode* m_ode);

double* q;
double* t;

//We declare the ode we wants to solve
Ode *m_ode;

};

#endif
