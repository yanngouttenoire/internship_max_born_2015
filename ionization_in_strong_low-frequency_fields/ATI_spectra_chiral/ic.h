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
ElectrostaticPotential<state_type> *myPotential;

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare the norm of the initial perpendicular velocity
double vPerpBirth;

//We declare the initial perpendicular velocity along Y
double vYPerpBirth;

//We declare the initial perpendicular velocity along Z', a vector of the rotating frame
double vZPrimPerpBirth;

//We declare the initial perpendicular velocity along X and Z which depend on the one along Z'
double vXPerpBirth;
double vZPerpBirth;

//We declare the probability of ionization with a given fieldBirth and a given vPerpBirth
double weightIonization;

//We declare electron radial position after tunneling
double rhoBirth;

//We declare polar coordinnates of the electron after tunneling
double ZBirth;
double XBirth;

//We declare constructor
IC(ElectrostaticPotential<state_type> *myPotential, ElectricField myField);
 
//we set the initial ionization time
void setTBirth(int iFieldBirth, int nFieldBirth);

//We set the initial field value
void setFieldBirth();

//The following method gives a values to the initial perpendicular velocity
double getVPerpBirth(int iVPerpBirth, int nVPerpBirth);

//We set the initial perpendicular velocity along Y
void setVYPerpBirth(int iVYPerpBirth, int nVYPerpBirth);

//We set the initial perpendicular velocity along X and Z
void setVXZPerpBirth(int iVZPrimPerpBirth, int nVZPrimPerpBirth);

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
IC<state_type>::IC(ElectrostaticPotential<state_type> *myPotential, ElectricField myField) : myPotential(myPotential), myField(myField){}


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
  fieldBirth=myField.fieldAmpl*sqrt(pow(cos(myField.pulsation*tBirth+myField.phase),2)+pow(myField.ellipticity*sin(myField.pulsation*tBirth+myField.phase),2));
  fieldBirth=fabs(fieldBirth);
}

//The following method gives a values to the initial perpendicular velocity
template<typename state_type>
double IC<state_type>::getVPerpBirth(int iVPerpBirth, int nVPerpBirth)
{
//if nVPerpBirth equals to 1, we set VPerpBirth equals to 0
//Thus, we can choose to put VYPerpBirth or VZPrimPerpBirth equals to 0 for all the simulation
if(nVPerpBirth==1)
return 0.0;
else
{
  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldBirth/sqrt(2.*myPotential->IP));

  //Perpendicular velocity after tunneling
  double vPerpBirth=4.*double(iVPerpBirth)/double(nVPerpBirth)*sigma_V-2*sigma_V;
  return vPerpBirth;
}
 
}


//We set the initial perpendicular velocity along Y
template<typename state_type>
void IC<state_type>::setVYPerpBirth(int iVYPerpBirth, int nVYPerpBirth)
{
vYPerpBirth=getVPerpBirth(iVYPerpBirth, nVYPerpBirth);
}

//We set the initial perpendicular velocity along X and Z
template<typename state_type>
void IC<state_type>::setVXZPerpBirth(int iVZPrimPerpBirth, int nVZPrimPerpBirth)
{
vZPrimPerpBirth=getVPerpBirth(iVZPrimPerpBirth, nVZPrimPerpBirth);

vZPerpBirth=+vZPrimPerpBirth*myField('X',tBirth)/fieldBirth;
vXPerpBirth=-vZPrimPerpBirth*myField('Z',tBirth)/fieldBirth;

}

//Ionization rate with a given field and a given transverse velocity according ADK distribution
//REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
template<typename state_type>
void IC<state_type>::setWeightIonization()
{
  vPerpBirth=pow(vXPerpBirth*vXPerpBirth+vYPerpBirth*vYPerpBirth+vZPerpBirth*vZPerpBirth, 1./2.);

  weightIonization=4./fieldBirth*exp(-2.*pow(2.*myPotential->IP,3./2.)/3./fieldBirth)*fabs(vPerpBirth)/fieldBirth/M_PI*exp(-pow(2.*myPotential->IP,1./2.)*vPerpBirth*vPerpBirth/fieldBirth);
}


//We set the electron radial position after tunneling
template<typename state_type>
void IC<state_type>::setRhoBirth()
{
  //We approximate the radial position of the electron after tunneling like:
  rhoBirth=myPotential->IP/fieldBirth;	
}


//We set the electron polar position after tunneling
template<typename state_type>
void IC<state_type>::setPolarCoordBirth()
{
   ZBirth=-rhoBirth*myField('Z',tBirth)/fieldBirth;
   XBirth=-rhoBirth*myField('X',tBirth)/fieldBirth;
}


//We build the function which sets the initial conditions
template<typename state_type>
void IC<state_type>::setIC(state_type &x, double& t) 
{

  //We set the initial time
   t=tBirth; 

  //We set the initial position in phase space
  x[0]=XBirth;
  x[1]=0.;
  x[2]=ZBirth;
  x[3]=vXPerpBirth;
  x[4]=vYPerpBirth;
  x[5]=vZPerpBirth;
 
}



#endif
