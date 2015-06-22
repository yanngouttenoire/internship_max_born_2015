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

//We declare two containers of map type, they can contain pairs which are couple of objects: first the range, second the ionization probability value

//We declare one for electrons which velocity vector is oriented along positive values of y
std::map<int,double> asymptEnergyUp;

//And one for electrons which velocity vector is oriented along negative values of y
std::map<int,double> asymptEnergyDown;

//We declare containers of map type for the closest approach, the asymptotic energy and the initial conditions
std::vector<double> asymptEnergy;
std::vector<double> closestApproach;
std::vector<double> fieldBirth;
std::vector<double> velocityBirth;
std::vector<double> weightIonization;

//We declare a variable for the bin interval width
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
void storeDataBinning(const state_type& x, const double& t, double a_fieldBirth, double a_velocityBirth, double a_closestApproach, double a_weightIonization, const bool& inexpectedStop);

//The following method is called by storeDataBinning and insert element <range,weightIonization> in map asymptEnergy
void insertInMap(std::map<int,double>& asymptEnergy, const int& range, const double& weightIonization);

//Finally we write all the data binning in a file
void writeDataBinning(std::fstream& dataFile);

//The following method is called by writeDataBinning and write all the elements range and weightIonization of map asymptEnergy in a file
void getFromMap(std::fstream& dataFile, std::map<int,double>& asymptEnergy);

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
void Spectra<state_type>::storeDataBinning(const state_type& x, const double& t, double a_fieldBirth, double a_velocityBirth, double a_closestApproach, double a_weightIonization, const bool& unexpectedStop)
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

  //If all is ok, we increment the trajectory ID
  spectraPointsNbr++;
  
//we insert values associated with the new event in the corresponding container according if the electron propagated in y>0 or y<0
 fieldBirth.push_back(a_fieldBirth);
 velocityBirth.push_back(a_velocityBirth); 
 closestApproach.push_back(a_closestApproach);
 asymptEnergy.push_back(energy);
 weightIonization.push_back(a_weightIonization);
}


  //Finally we write all the data binning in a file
template<typename state_type>
void Spectra<state_type>::writeDataBinning(std::fstream& dataFile)
{

  std::vector<double>::iterator it= asymptEnergy.begin();
  //We write the data in a file	
  for(int k=0; it!=asymptEnergy.end(); it++, k++)
    {
      dataFile<<k<<" "<<fieldBirth[k]<<" "<<velocityBirth[k]<<" "<<closestApproach[k]<<" "<<weightIonization[k]<<" "<<(*it)*37.3<<std::endl;
    }
 
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
