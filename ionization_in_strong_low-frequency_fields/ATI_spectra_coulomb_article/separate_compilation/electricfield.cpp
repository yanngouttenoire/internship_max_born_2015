#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"electricfield.h"

using namespace std;

ElectricField::ElectricField()
{

//CONSTANTS
double uaIntensity=3.5094451E16;
double uaTime=24.1888421562712E-18;
double lightSpeed=2.99792458E8;

//Field parameters
waveLenght=1064E-9;
fieldAmpl=sqrt(13*(1E13)/uaIntensity);
pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
opticalCycle=2.*M_PI/pulsation;
pulseDuration=3.*opticalCycle;
phase=0.;
ellipticity=0.1;

}

double ElectricField::operator()(char component, const double& t)
{
switch(component)
{
case 'X' :
return 0.;

case 'Y' :
return 0.;

case 'Z' :
return -fieldAmpl*cos(pulsation*t+phase);
}
}

double ElectricField::vectPot(char component, const double& t)
{
switch(component)
{
case 'X' :
return 0.;

case 'Y' :
return 0.;

case 'Z' :
return -fieldAmpl/pulsation*sin(pulsation*t+phase);
}
}
