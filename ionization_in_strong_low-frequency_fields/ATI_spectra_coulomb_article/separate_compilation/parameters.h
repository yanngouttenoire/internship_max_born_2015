#ifndef PARAMETERS_H
#define PARAMETERS_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

extern double waveLenght;
extern double fieldAmpl;
extern double pulsation;
extern double opticalCycle;
extern double pulseDuration;
extern double phase;
extern double ellipticity;

extern double IP;
extern double Z;

extern int nFieldBirth=100, nVPerpBirth=100;

extern double dtMin=1E-20;
extern double desiredErrorMax=1E-13;
extern double desiredErrorMin=desiredErrorMax/10.;

extern int dtMinNbr=0;
extern int spectraPointsNbr=0;
extern std::string elapsedTime="0";


#endif
