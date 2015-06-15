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
waveLenght=2E-6;
fieldAmpl=0.0534;
cyclesNbr=20;
phase=0.;
pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
opticalCycle=2.*M_PI/pulsation;
}

//The electric field is a pulse whose full duration at half maximum contains 2 optical cycles
double ElectricField::operator()(char component, const double& t)
{
if(fabs(t)<opticalCycle*cyclesNbr)
{
switch(component)
{
case 'X' :
return pow(cos(M_PI*t/2./opticalCycle/cyclesNbr),2)*ellipticity*fieldAmpl*sin(pulsation*t+phase);

case 'Y' :
return 0.;

case 'Z' :
return pow(cos(M_PI*t/2./opticalCycle/cyclesNbr),2)*fieldAmpl*cos(pulsation*t+phase);
}
}
else
return 0.;
}

//We build the asssociated vector potential
double ElectricField::vectPot(char component, const double& t)
{
if(fabs(t)<opticalCycle*cyclesNbr)
{
switch(component)
{
case 'X' :
return pow(cos(M_PI*t/2./opticalCycle/cyclesNbr),2)*ellipticity*fieldAmpl/pulsation*cos(pulsation*t+phase);

case 'Y' :
return 0.;

case 'Z' :
return -pow(cos(M_PI*t/2./opticalCycle/cyclesNbr),2)*fieldAmpl/pulsation*sin(pulsation*t+phase);
}
}
else
return 0.;
}
