#ifndef COULOMB_H
#define COULOMB_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"ode.h"

//This class implements Coulomb ode and inherits from Ode

class Coulomb : public Ode
{

public:
Coulomb(int au);

void diff(double* q, double t);


double K;
int au;

};

#endif
