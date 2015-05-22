#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

//CONSTANTS
double uaIntensity=3.5094451E16;
double uaTime=24.1888421562712E-18;
double uaEnergy=27.211608;
double lightSpeed=2.99792458E8;

//Field parameters
double waveLenght=1064E-9;
double fieldAmpl=sqrt(13*(1E13)/uaIntensity);
double pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
double opticalCycle=2.*M_PI/pulsation;
double pulseDuration=3.*opticalCycle;
double phase=0;
double ellipticity=0.1;

//Coulomb potential parameters
double IP=13.605804/uaEnergy;
double Z=1.; //Positive charge of the nucleus

//Numbers of computed points
int nFieldBirth=100, nVPerpBirth=100;

//Desired accuracy
double dtMin=1E-20;
double desiredErrorMax=1E-13;
double desiredErrorMin=desiredErrorMax/10.;

//Some observers
int dtMinNbr=0;
int spectraPointsNbr=0;
string elapsedTime="0";
