#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"coulomb.h" 

Coulomb::Coulomb(int a_au)
: Ode()
{
au=a_au;
double Z=1.; //Positive charge of the nucleus
double e=1.6E-19;
double mass=9.11E-31;
double eps0=8.85E-12;

if(au==0)
K=e*e*Z/(4*M_PI*eps0*mass);
else
K=Z;

}

void Coulomb::diff(double* q, double t)
{
double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],3./2.);
qp[0]=q[3];
qp[1]=q[4];
qp[2]=q[5];
qp[3]=-K*q[0]/a;
qp[4]=-K*q[1]/a;
qp[5]=-K*q[2]/a;
}
