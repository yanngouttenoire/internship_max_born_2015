#ifndef SYSTEM_H
#define SYSTEM_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

#include"electricfield.h"
#include"electrostaticpotential.h"
#include"hydrogen.h"

//Ordinary differential equations of the dynamic of an electron in coulomb potential 
template<typename state_type>
class System
{
  //We declare an object of type ElectrostaticPotential for electrostatic potential properties
  //ElectrostaticPotential<state_type> *myPotential;	
  Hydrogen<state_type> myPotential;

  //We declare field variables
  ElectricField myField;		

 public:

  //Contructor 
  //System(ElectrostaticPotential<state_type> myPotential, ElectricField &myField);
  System(Hydrogen<state_type> myPotential, ElectricField &myField);
  
  //Copy constructor (in order to use System as firstprivate in OPEN MP)
 // System(const System & source);
  
  //Copy assignement operator
 // System & operator= (const System & source);
  
  //Destructor
  ~System();
  
  //We implement the method which contains dynamic equations of motion
  //We overload the operator() for convenience interest
  void operator() (const state_type &x , state_type &dxdt , const double&  t);

};

#endif
