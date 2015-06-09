#ifndef ELECTROSTATICPOTENTIAL_H
#define ELECTROSTATICPOTENTIAL_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>
/*
template<typename state_type> class A
{
public:
double a;
double softParameter;
double radialDistanceCube;
A()
{
a=0.5;
softParameter=0.;
}
virtual void preparePotential(const state_type &x);
};
template<typename state_type>
void A<state_type>::preparePotential(const state_type &x)
{
radialDistanceCube=pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.); 
}

*/

//Define characteristics and properties of an eletric field
template<typename state_type> class ElectrostaticPotential
{

public:

//Electrostatic potential parameters
double IP;
double softParameter;


//Constructor
ElectrostaticPotential()
{
IP=0.5;
softParameter=0.;
}

virtual void preparePotential(const state_type &x);

//Methods which return the coulomb force on each components 
virtual double operator()(char component, const state_type &x);

//Method which returns the electrostatic potential energy
virtual double potentialEnergy(const state_type &x);

};


template<typename state_type>
void ElectrostaticPotential<state_type>::preparePotential(const state_type &x)
{

}



template<typename state_type>
double ElectrostaticPotential<state_type>::operator()(char component, const state_type &x)
{
double radialDistanceCube=pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.); 
switch(component)
{
case 'X' :
return -x[0]/radialDistanceCube;

case 'Y' :
return -x[1]/radialDistanceCube;

case 'Z' :
return -x[2]/radialDistanceCube;
}
}

template<typename state_type>
double ElectrostaticPotential<state_type>::potentialEnergy(const state_type& x)
{
return -1/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,1./2.); 
}





#endif
