#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"ode.h" 

Ode::Ode()
{
qp=new double[6];
}

Ode::~Ode()
{
delete[] qp;
}
