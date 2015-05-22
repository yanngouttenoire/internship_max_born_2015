#ifndef IC_H
#define IC_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"parameters.h"
 
class IC
{

public:


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
IC(double tBirth=0.871, double vPerpBirth=6E-5);
 
//we set the initial ionization time
void setTBirth(double iFieldBirth, double nFieldBirth);

//We set the initial field value
void setFieldBirth();

//We set the initial perpendicular velocity
void setVPerpBirth(double iVPerpBirth, double nVPerpBirth);

//We set the electron position after tunneling according the article
void setRhoBirth();

//We build the function which sets the initial conditions
template<typename state_type>
void setIC(state_type &x, double t);


};

#endif
