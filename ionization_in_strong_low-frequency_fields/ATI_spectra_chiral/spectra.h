#ifndef SPECTRA_H
#define SPECTRA_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"electricfield.h"
#include"electrostaticpotential.h"

//We implement a class for computing the photo-electron spectrum as an histogram of asymptotic Energy according probability of ionization (weightIonization)
template<typename state_type>
class Spectra
{

public:

double binsWidth;

//We declare a counter for the number of events stored
int spectraPointsNbr;

//We declare a counter for the number of events not stored because the electron always trapped by the atom after the pulse
int trappedElectronNbr;

//We declare an object of type ElectrostaticPotential
ElectrostaticPotential<state_type> *myPotential;

//We declare an object of type ElectricField
ElectricField myField;


//Constructor
Spectra(ElectrostaticPotential<state_type> *myPotential, ElectricField myField, double binsWidth=0.1);

//We compute the asymptotic energy
double asymptoticEnergy(const state_type& x, const double& t);

//We store asymptotic velocities in containers of map type
void writeData(const state_type& x, const double& t, double tBirth, double fieldBirth, double vPerpBirth, double closestApproach, double weightIonization, std::fstream& dataFile, const bool& inexpectedStop);

};


//We set the histogram intervals width
template<typename state_type>
Spectra<state_type>::Spectra(ElectrostaticPotential<state_type> *myPotential, ElectricField myField, double binsWidth) : myPotential(myPotential), myField(myField), binsWidth(binsWidth) 
{
spectraPointsNbr=0;
trappedElectronNbr=0;
}

//We store asymptotic energies in containers of map type with a view to make a data binning
template<typename state_type>
void Spectra<state_type>::writeData(const state_type& x, const double& t, double tBirth, double fieldBirth, double vPerpBirth, double closestApproach, double weightIonization, std::fstream& dataFile, const bool& unexpectedStop)
{

  //If the computation of the trajectory has been stopped unexpectedly, we do not consider the event
  if(unexpectedStop) return;
  
   double energy=asymptoticEnergy(x,t);

  //If the energy of the electron is negative, the electron is not free and we do not consider the event
  //If the range equals zero we do not consider this neither (esthetic choice)
  if(energy<=0) 
  {
  trappedElectronNbr++;
  return;
  }

//Angular resolution
//We consider electrons differently depending if they are detected along the polarization of the field or not
if(fabs(atan(sqrt(x[3]*x[3]+x[4]*x[4])/x[5]))*180./M_PI>5)
return;
//if(x[2]*myField('Z',myIC->tBirth)>=0)


  //If all is ok, we increment the trajectory ID
  spectraPointsNbr++;
  
//We write data in a file
 dataFile<<spectraPointsNbr<<" "<<fieldBirth<<" "<<vPerpBirth<<" "<<closestApproach<<" "<<weightIonization<<" "<<energy*37.3<<" "<<tBirth<<std::endl;

}


//We compute asymptotic energy
template<typename state_type>
double Spectra<state_type>::asymptoticEnergy(const state_type &x, const double& t)
{
  double Vsq=(x[3]-myField.vectPot('X',t))*(x[3]-myField.vectPot('X',t))+(x[4]-myField.vectPot('Y',t))*(x[4]-myField.vectPot('Y',t))+(x[5]-myField.vectPot('Z',t))*(x[5]-myField.vectPot('Z',t));
   double E=Vsq/2-myPotential->potentialEnergy(x);

  return E;

}



#endif
