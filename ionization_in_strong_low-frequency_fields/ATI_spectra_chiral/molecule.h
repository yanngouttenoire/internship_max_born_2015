#ifndef MOLECULE_H
#define MOLECULE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

#include"electrostaticpotential.h"


//Define characteristics and properties of an electrostatic potential for a given molecule
template<typename state_type> class Molecule : public ElectrostaticPotential<state_type>
{

 public:


  //Electrostatic potential parameters
  double charge[4];
  double bondLength[3]; 
  double inverseRadialDistanceCube[4];
  double sumInverseRadialDistanceCube;

  //Constructor
  Molecule();
  
  //Necessary method for copy procedure when using polymophism (abstract parent class)
  Molecule* Clone()
  {
  return new Molecule(*this);
  }

  //Method which compute some quantities in advance
  void preparePotential(const state_type &x);

  //Methods which return the coulomb force on each components 
  double operator()(char component, const state_type &x);

  //Method which returns the electrostatic potential energy
  double potentialEnergy(const state_type &x);

};

template<typename state_type>
Molecule<state_type>::Molecule()
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
  //Ionization potential
  this->IP=0.32772257;
}




template<typename state_type>
void Molecule<state_type>::preparePotential(const state_type &x)
{

  inverseRadialDistanceCube[0] = charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,3./2.);
  inverseRadialDistanceCube[1] = charge[1]/pow((x[0]-bondLength[0])*(x[0]-bondLength[0])+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,3./2.) ;
  inverseRadialDistanceCube[2] = charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[1])*(x[1]-bondLength[1])+x[2]*x[2]+this->softParameter*this->softParameter,3./2.) ;
  inverseRadialDistanceCube[3] = charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[2])*(x[2]-bondLength[2])+this->softParameter*this->softParameter,3./2.) ;
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

  return -charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[1]/pow((x[0]-bondLength[0])*(x[0]-bondLength[0])+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[1])*(x[1]-bondLength[1])+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[2])*(x[2]-bondLength[2])+this->softParameter*this->softParameter,1./2.);

}



#endif
