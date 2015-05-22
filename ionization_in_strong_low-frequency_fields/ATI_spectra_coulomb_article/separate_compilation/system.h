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

//Ordinary differential equations of the dynamic of an electron in coulomb potential 
template<typename state_type>
class System
{
  //We declare the coulomb constant
  double K;
  //We declare a softening parameter in order to implement an eventual soft coulomb potential
  double softParameter=0.;

  //We declare field variables
  double field[3];
  std::string whichField;

  public:

  //Contructor 
  System(std::string whichField);
  
  //We implement the method which contains dynamic equations of motion
  //We overload the operator() for convenience interest
  void operator() (const state_type &x , state_type &dxdt , const double  t);

};

#endif
