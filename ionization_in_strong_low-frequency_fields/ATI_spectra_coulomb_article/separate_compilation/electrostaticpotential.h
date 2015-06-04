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

//Define characteristics and properties of an eletric field
template<typename state_type>
class ElectrostaticPotential
{

public:

//Electrostatic potential parameters
double IP;
double charge; //Positive charge of the nucleus
double softParameter;

//radial distance
double radialDistanceCube;

//Constructor
ElectrostaticPotential();

//Method which sets the radial distance  (it is used to compute radialDistanceCube in advance and then decrease the number of operations)
void preparePotential(const state_type &x);

//Methods which return the coulomb force on each components 
double operator()(char component, const state_type &x);

//Method which returns the electrostatic potential energy
double potentialEnergy(const state_type &x);

};

template<typename state_type>
ElectrostaticPotential<state_type>::ElectrostaticPotential()
{

//CONSTANTS
double uaEnergy=27.211608;

//Electrostatic potential parameters
IP=13.605804/uaEnergy;
charge=1.; //Positive charge of the nucleus
softParameter=0.;

}


template<typename state_type>
void ElectrostaticPotential<state_type>::preparePotential(const state_type &x)
{
radialDistanceCube=pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.); 
}



template<typename state_type>
double ElectrostaticPotential<state_type>::operator()(char component, const state_type &x)
{
switch(component)
{
case 'X' :
return -charge*x[0]/radialDistanceCube;

case 'Y' :
return -charge*x[1]/radialDistanceCube;

case 'Z' :
return -charge*x[2]/radialDistanceCube;
}
}

template<typename state_type>
double ElectrostaticPotential<state_type>::potentialEnergy(const state_type& x)
{
return -charge/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.); 
}




#endif
