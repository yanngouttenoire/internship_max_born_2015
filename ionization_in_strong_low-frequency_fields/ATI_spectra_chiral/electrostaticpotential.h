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


//Define characteristics and properties of an electrostatic potential for an hydrogenoid atom
template<typename state_type> class ElectrostaticPotential
{

 public:

  //Electrostatic potential parameters
  double IP;
  double softParameter;


  //Constructor
  ElectrostaticPotential()
    {
      //Parameter for soft coulomb potential
      softParameter=0.;
    }
    
   //Necessary method for copy procedure when using polymophism (abstract parent class)
  virtual ElectrostaticPotential* Clone()=0;
  
  virtual void preparePotential(const state_type &x)=0;

  //Methods which return the coulomb force on each components 
  virtual double operator()(char component, const state_type &x)=0;

  //Method which returns the electrostatic potential energy
  virtual double potentialEnergy(const state_type &x)=0;

};





#endif
