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

//Return the ionization rate with a given field and a given transverse velocity
double getWeightIonization(double t, double vPerp, int i);

//Compute the maximum ionization rate
double getMaxWeightIonization(int k=2);

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
//We want to scan field phase values from -200° to 20°

double anglei=-200.*M_PI/180.;
double anglef=20.*M_PI/180.;
tBirth=anglei/myField.pulsation+double(iFieldBirth)/double(nFieldBirth)*(anglef-anglei)/myField.pulsation;
}


//We set the initial field value
template<typename state_type>
void IC<state_type>::setFieldBirth()
{
//We use tBirth
  fieldBirth=myField.getInstFieldAmpl(tBirth);
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
  double sigma_V=0.3;

  //Perpendicular velocity after tunneling
  return double(iVPerpBirth)/double(nVPerpBirth)*sigma_V;
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
double IC<state_type>::getWeightIonization(double t, double vPerp, int i)
{
double field=myField.getInstFieldAmpl(t);
if(i==1)
return 4./field*exp(-2.*pow(2.*myPotential->IP,3./2.)/3./field)*fabs(vPerp)/field*exp(-pow(2.*myPotential->IP,1./2.)*vPerp*vPerp/field);
if(i==2)
return 1./field/field*exp(-2.*pow((2.*myPotential->IP+vPerp*vPerp),3./2.)/3/field)*vPerp/sqrt(1+vPerp*vPerp/2./myPotential->IP);
}

//Compute the maximum ionization rate
template<typename state_type>
double IC<state_type>::getMaxWeightIonization(int k)
{
int ni=100;
int nj=100;
double vPerpMax=0.3;
double ti;
double vPerp;
double weight;
double max=0.;

for(int i=0; i<ni; i++)
{
for(int j=0; j<nj; j++)
{
 ti=-M_PI/myField.pulsation/90.*25+i/double(ni)*2.*M_PI/myField.pulsation/90.*25;
 vPerp=0.01+j*vPerpMax/nj;

weight=getWeightIonization(ti, vPerp,k);
if(weight>max)
max=weight;

}
}
return max;
}


//We set the ionization rate
template<typename state_type>
void IC<state_type>::setWeightIonization()
{
  vPerpBirth=pow(vXPerpBirth*vXPerpBirth+vYPerpBirth*vYPerpBirth+vZPerpBirth*vZPerpBirth, 1./2.);
  weightIonization=getWeightIonization(tBirth, vPerpBirth,2);
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
