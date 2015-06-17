#ifndef HYDROGEN_H
#define HYDROGEN_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

#include"electrostaticpotential.h"

template<typename state_type> class Hydrogen : public ElectrostaticPotential<state_type>
{

public:

//Electrostatic potential parameters
double charge;
double inverseRadialDistanceCube;

//Constructor
Hydrogen();

//We set the ionization potential
void setIP(double a_IP);

//Method which compute some quantities in advance
void preparePotential(const state_type &x);

//Methods which return the coulomb force on each components 
double operator()(char component, const state_type &x);

//Method which returns the electrostatic potential energy
double potentialEnergy(const state_type &x);

};

template<typename state_type>
Hydrogen<state_type>::Hydrogen()
{

//Nuclei charges
charge=1.;
//Ionization potential
this->IP=0.5;

}

template<typename state_type>
void Hydrogen<state_type>::setIP(double a_IP)
{
this->IP=a_IP;
}



template<typename state_type>
void Hydrogen<state_type>::preparePotential(const state_type &x)
{

inverseRadialDistanceCube = charge/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,3./2.) ; 

}


template<typename state_type>
double Hydrogen<state_type>::operator()(char component, const state_type &x)
{

switch(component)
{
case 'X' :
return -x[0]*inverseRadialDistanceCube;

case 'Y' :
return -x[1]*inverseRadialDistanceCube;

case 'Z' :
return -x[2]*inverseRadialDistanceCube;
}


}

template<typename state_type>
double Hydrogen<state_type>::potentialEnergy(const state_type& x)
{

return -charge/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,1./2.); 

}



#endif
