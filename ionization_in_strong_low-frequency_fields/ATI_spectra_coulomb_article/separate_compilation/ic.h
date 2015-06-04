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

//We declare a variable for the electron position after tunneling
double rhoBirth;



//We declare constructor
IC(ElectrostaticPotential<state_type> myPotential, ElectricField myField, double tBirth=0.871, double vPerpBirth=6E-5);
 
//we set the initial ionization time
void setTBirth(int iFieldBirth, int nFieldBirth);

//We set the initial field value
void setFieldBirth();

//We set the initial perpendicular velocity
void setVPerpBirth(int iVPerpBirth, int nVPerpBirth);

//We set the electron position after tunneling according the article
void setRhoBirth();

//We build the function which sets the initial conditions
void setIC(state_type &x, double t);


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
fieldBirth=myField.fieldAmpl*cos(myField.pulsation*tBirth+myField.phase);
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


  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
  weightIonization=4./fieldBirth*exp(-2.*pow(2.*myPotential.IP,3./2.)/3./fieldBirth)*fabs(vPerpBirth)/fieldBirth/M_PI*exp(-pow(2.*myPotential.IP,1./2.)*vPerpBirth*vPerpBirth/fieldBirth);

}


//We set the electron position after tunneling according the article
template<typename state_type>
void IC<state_type>::setRhoBirth()
{

  //We set the initial field value
  double  fieldBirth=myField.fieldAmpl*cos(myField.pulsation*tBirth+myField.phase);
  fieldBirth=fabs(fieldBirth);

  // we find the root of the function which gives us the position of the electron after tunneling
  double x=100.;

  for(int k=0; k<50; k++)                       
    {
       x=x-(-1./4./x-1./8./x/x-1./8.*fabs(fieldBirth)*x+1./8.)/(+1./4./x/x+1./4./x/x/x-fabs(fieldBirth)/8.);
    }
 
  rhoBirth=-x/2.*fabs(fieldBirth)/fieldBirth;

}


//We build the function which sets the initial conditions
template<typename state_type>
void IC<state_type>::setIC(state_type &x, double t) 
{

  //We set the initial time
   t=tBirth; 

  //We set the initial position in phase space
  x[0]=0.;
  x[1]=0.;
  x[2]=rhoBirth;
  x[3]=vPerpBirth;
  x[4]=0.;
  x[5]=0.;
 
}



#endif
