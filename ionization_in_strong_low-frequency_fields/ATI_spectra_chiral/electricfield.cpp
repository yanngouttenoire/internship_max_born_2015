#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"electricfield.h"

using namespace std;

ElectricField::ElectricField(double ellipticity)
  : ellipticity(ellipticity)
{

  //CONSTANTS
  double uaIntensity=3.5094451E16;
  double uaTime=24.1888421562712E-18;
  double lightSpeed=2.99792458E8;

  //Field parameters
  pulsation=0.05;
  fieldAmpl=0.04;
  cyclesNbr=4;
  phase=0.;
  waveLenght=2.*M_PI*lightSpeed/pulsation*uaTime;
  opticalCycle=2.*M_PI/pulsation;
}

//The electric field is a pulse whose full duration at half maximum contains 2 optical cycles
//The envelope can be
//either exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))
//or pow(cos(M_PI*t/2./opticalCycle/cyclesNbr),2)
//with if(fabs(t)<opticalCycle*cyclesNbr) return [..] else return 0
double ElectricField::operator()(char component, const double& t)
{
  switch(component)
    {
    case 'X' :
      return exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))*ellipticity*fieldAmpl*sin(pulsation*t+phase);

    case 'Y' :
      return 0.;

    case 'Z' :
      return exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))*fieldAmpl*cos(pulsation*t+phase);
    }

}

//We build a method which returns the instantaneous field amplitude for a given time
double ElectricField::getInstFieldAmpl(const double& t)
{
  return sqrt(pow(operator()('X',t),2)+pow(operator()('Y',t),2)+pow(operator()('Z',t),2));
}

//We build the asssociated vector potential
double ElectricField::vectPot(char component, const double& t)
{

  switch(component)
    {
      //specific to the article
      return 0;

    case 'X' :
      return exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))*ellipticity*fieldAmpl/pulsation*cos(pulsation*t+phase);

    case 'Y' :
      return 0.;

    case 'Z' :
      return -exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))*fieldAmpl/pulsation*sin(pulsation*t+phase);
    }

}
