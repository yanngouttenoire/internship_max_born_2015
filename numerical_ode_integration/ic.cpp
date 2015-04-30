#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"ic.h" 

Ic::Ic(double* a_q)
{
q=a_q;
}

void Ic::setIc(double x0, double y0, double z0, double vx0, double vy0, double vz0, double* t)
{
//We initialise the position of the orbit in the phase space
q[0]=x0;
q[1]=y0;
q[2]=z0;
q[3]=vx0;
q[4]=vy0;
q[5]=vz0;

//We initialise the time variable
*t=0.;
}
