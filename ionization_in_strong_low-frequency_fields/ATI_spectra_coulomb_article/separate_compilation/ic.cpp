#include"ic.h"

using namespace std;

//We declare constructor
IC::IC(double tBirth, double vPerpBirth) : tBirth(tBirth), vPerpBirth(vPerpBirth) 
{
setFieldBirth(); 
setRhoBirth();
}


//we set the initial ionization time tBirth
void IC::setTBirth(double iFieldBirth, double nFieldBirth)
{ 
//We want to scan field phase values from -pi/4 to pi/4
tBirth=-M_PI/4./pulsation+double(iFieldBirth)/double(nFieldBirth)*M_PI/2./pulsation;
}


//We set the initial field value
void IC::setFieldBirth()
{
//We use tBirth
fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
fieldBirth=fabs(fieldBirth);
}


//We set the initial perpendicular velocity
void IC::setVPerpBirth(double iVPerpBirth, double nVPerpBirth)
{

  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldBirth/sqrt(2.*IP));

  //Perpendicular velocity after tunneling
  vPerpBirth=2.*double(iVPerpBirth)/double(nVPerpBirth)*sigma_V;
  vPerpBirth=2.*vPerpBirth/4.;


  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
  weightIonization=4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerpBirth)/fieldBirth/M_PI*exp(-pow(2.*IP,1./2.)*vPerpBirth*vPerpBirth/fieldBirth);

}


//We set the electron position after tunneling according the article
void IC::setRhoBirth()
{

  //We set the initial field value
  double  fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
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
void IC::setIC(state_type &x, double t) 
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



template void IC::setIC<vector<double> >;
template void IC::setIC<double* >;
//For more details about templates, see https://www.cs.umd.edu/class/fall2002/cmsc214/Projects/P2/proj2.temp.html
