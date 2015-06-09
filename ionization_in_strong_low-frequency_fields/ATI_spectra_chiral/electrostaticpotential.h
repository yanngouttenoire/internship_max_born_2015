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
template<typename state_type> class ElectrostaticPotential
{

public:

//Electrostatic potential parameters
double IP;
double softParameter;
double charge;
double radialDistanceCube;


//Constructor
ElectrostaticPotential(): IP(0.5), charge(1.), softParameter(0) {}

virtual void preparePotential(const state_type &x);

//Methods which return the coulomb force on each components 
virtual double operator()(char component, const state_type &x);

//Method which returns the electrostatic potential energy
virtual double potentialEnergy(const state_type &x);

};


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


//gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg


template<typename state_type> class Molecule : public ElectrostaticPotential<state_type>
{

public:

double IP;
double softParameter;

//Electrostatic potential parameters
double charge[4];
double bondLength[3]; 
double inverseRadialDistanceCube[4];
double sumInverseRadialDistanceCube;

//Constructor
Molecule()  //: ElectrostaticPotential()
{

//Nuclei charges
charge[0]=0.4;
charge[1]=0.2;
charge[2]=0.2;
charge[3]=0.2;
//Covalent bond length
bondLength[0]=2.0;
bondLength[1]=1.0;
bondLength[2]=3.0;
}


//Method which compute some quantities in advance
void preparePotential(const state_type &x);

//Methods which return the coulomb force on each components 
double operator()(char component, const state_type &x);

//Method which returns the electrostatic potential energy
double potentialEnergy(const state_type &x);

};
/*
template<typename state_type>
Molecule<state_type>::Molecule(): ElectrostaticPotential(0.5)
{


//CONSTANTS
double uaEnergy=27.211608;

//Ionization potential
IP=13.6/uaEnergy;
//Nuclei charges
charge[0]=0.4;
charge[1]=0.2;
charge[2]=0.2;
charge[3]=0.2;
//Covalent bond length
bondLength[0]=2.0;
bondLength[1]=1.0;
bondLength[2]=3.0;


}
*/



template<typename state_type>
void Molecule<state_type>::preparePotential(const state_type &x)
{

inverseRadialDistanceCube[0] = charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.) ;
inverseRadialDistanceCube[1] = charge[1]/pow((x[0]-bondLength[0])*(x[0]-bondLength[0])+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.) ;
inverseRadialDistanceCube[2] = charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[1])*(x[1]-bondLength[1])+x[2]*x[2]+softParameter*softParameter,3./2.) ;
inverseRadialDistanceCube[3] = charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[2])*(x[2]-bondLength[2])+softParameter*softParameter,3./2.) ;
sumInverseRadialDistanceCube = inverseRadialDistanceCube[0]+inverseRadialDistanceCube[1]+inverseRadialDistanceCube[2]+inverseRadialDistanceCube[3];

}



template<typename state_type>
double Molecule<state_type>::operator()(char component, const state_type &x)
{

switch(component)
{
case 'X' :
return -x[0]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[1]*bondLength[0];

case 'Y' :
return -x[1]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[2]*bondLength[1];

case 'Z' :
return -x[2]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[3]*bondLength[2];
}


}


template<typename state_type>
double Molecule<state_type>::potentialEnergy(const state_type& x)
{

return -charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,1./2.)-charge[1]/pow((x[0]-bondLength[0])*(x[0]-bondLength[0])+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,1./2.)-charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[1])*(x[1]-bondLength[1])+x[2]*x[2]+softParameter*softParameter,1./2.)-charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[2])*(x[2]-bondLength[2])+softParameter*softParameter,1./2.);

}





#endif
