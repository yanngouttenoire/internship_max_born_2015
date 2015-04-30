#ifndef IC_H
#define IC_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>


class Ic
{

public:
Ic(double* q);
void setIc(double x0, double y0, double z0, double vx0, double vy0, double vz0, double* t);

double* q;

};

#endif
