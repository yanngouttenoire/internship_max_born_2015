#ifndef IC_H
#define IC_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"electrostaticpotential.h"
#include"electricfield.h"
 
template<typename state_type>
class IC
{

public:


//We declare an object of type Field in order to access the field properties
ElectricField myField;

//We declare an object of type ElectrostaticPotential in order to access electrostatic potential properties
ElectrostaticPotential<state_type> myPotential;

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare a variable for the initial perpendicular velocity
double vPerpBirth;

//We declare a variable for the probability of ionization with a given fieldBirth and a given vPerpBirth
double weightIonization;

//We declare electron radial position after tunneling
double rhoBirth;

//We declare polar coordinnates of the electron after tunneling
double ZBirth;
double XBirth;



//We declare constructor
IC(ElectrostaticPotential<state_type> myPotential, ElectricField myField, double tBirth=0.871, double vPerpBirth=6E-5);
 
//we set the initial ionization time
void setTBirth(int iFieldBirth, int nFieldBirth);

//We set the initial field value
void setFieldBirth();

//We set the initial perpendicular velocity
void setVPerpBirth(int iVPerpBirth, int nVPerpBirth);

//We compute the ionization rate with a given field and a given transverse velocity
void setWeightIonization();

//We set the electron position after tunneling
void setRhoBirth();

//We set the electron polar position after tunneling
void setPolarCoordBirth();

//We build the function which sets the initial conditions
void setIC(state_type &x, double& t);


};


//We declare constructor
template<typename state_type>
IC<state_type>::IC(ElectrostaticPotential<state_type> myPotential, ElectricField myField, double tBirth, double vPerpBirth) : tBirth(tBirth), vPerpBirth(vPerpBirth) 
{
setFieldBirth(); 
setRhoBirth();
}


//we set the initial ionization time tBirth
template<typename state_type>
void IC<state_type>::setTBirth(int iFieldBirth, int nFieldBirth)
{ 
//We want to scan field phase values from -pi/4 to pi/4
tBirth=-M_PI/4./myField.pulsation+double(iFieldBirth)/double(nFieldBirth)*M_PI/2./myField.pulsation;
}


//We set the initial field value
template<typename state_type>
void IC<state_type>::setFieldBirth()
{
//We use tBirth
 // double  fieldBirth=myField.fieldAmpl*sqrt(pow(cos(myField.pulsation*tBirth+myField.phase),2)+pow(myField.ellipticity*sin(myField.pulsation*tBirth+myField.phase),2));
double  fieldBirth=myField.fieldAmpl;
  fieldBirth=fabs(fieldBirth);
}


//We set the initial perpendicular velocity
template<typename state_type>
void IC<state_type>::setVPerpBirth(int iVPerpBirth, int nVPerpBirth)
{
  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldBirth/sqrt(2.*myPotential.IP));

  //Perpendicular velocity after tunneling
  vPerpBirth=2.*double(iVPerpBirth)/double(nVPerpBirth)*sigma_V;
  vPerpBirth=2.*vPerpBirth/4.;

weightIonization=4./fieldBirth*exp(-2.*pow(2.*myPotential.IP,3./2.)/3./fieldBirth)*fabs(vPerpBirth)/fieldBirth/M_PI*exp(-pow(2.*myPotential.IP,1./2.)*vPerpBirth*vPerpBirth/fieldBirth);
}


//Ionization rate with a given field and a given transverse velocity according ADK distribution
//REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
template<typename state_type>
void IC<state_type>::setWeightIonization()
{
  weightIonization=4./fieldBirth*exp(-2.*pow(2.*myPotential.IP,3./2.)/3./fieldBirth)*fabs(vPerpBirth)/fieldBirth/M_PI*exp(-pow(2.*myPotential.IP,1./2.)*vPerpBirth*vPerpBirth/fieldBirth);
}


//We set the electron radial position after tunneling
template<typename state_type>
void IC<state_type>::setRhoBirth()
{
  //We approximate the radial position of the electron after tunneling like:
  rhoBirth=myPotential.IP/fieldBirth;	
}


//We set the electron polar position after tunneling
template<typename state_type>
void IC<state_type>::setPolarCoordBirth()
{
 /*ZBirth=-rhoBirth*myField('Z',tBirth)/fieldBirth;
 XBirth=-rhoBirth*myField('X',tBirth)/fieldBirth;*/
}


//We build the function which sets the initial conditions
template<typename state_type>
void IC<state_type>::setIC(state_type &x, double& t) 
{

  //We set the initial time
   t=tBirth; 

  //We set the initial position in phase space
  x[0]=rhoBirth;
  x[1]=0.;
  x[2]=0.;
  x[3]=vPerpBirth;
  x[4]=0.;
  x[5]=0.;
 
}



#endif
