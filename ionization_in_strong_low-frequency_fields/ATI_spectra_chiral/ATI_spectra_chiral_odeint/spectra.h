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

//We declare a variable for the bin interval width
double binsWidth;

//We declare a counter for the number of events stored
int spectraPointsNbr;

//We declare an object of type ElectrostaticPotential
ElectrostaticPotential<state_type> myPotential;

//We declare an object of type ElectricField
ElectricField myField;


//Constructor
Spectra(ElectrostaticPotential<state_type> myPotential, ElectricField myField, double binsWidth=0.1);

//We compute the asymptotic energy
double asymptoticEnergy(const state_type& x, const double& t);

//We store asymptotic velocities in containers of map type
void storeDataBinning(const state_type& x, const double& t, const double& weightIonization, const bool& inexpectedStop);

//The following method is called by storeDataBinning and insert element <range,weightIonization> in map asymptEnergy
void insertInMap(std::map<int,double>& asymptEnergy, const int& range, const double& weightIonization);

//Finally we write all the data binning in a file
void writeDataBinning(std::fstream& dataFile);

//The following method is called by writeDataBinning and write all the elements range and weightIonization of map asymptEnergy in a file
void getFromMap(std::fstream& dataFile, std::map<int,double>& asymptEnergy);

};


//We set the histogram intervals width
template<typename state_type>
Spectra<state_type>::Spectra(ElectrostaticPotential<state_type> myPotential, ElectricField myField, double binsWidth) : myPotential(myPotential), myField(myField), binsWidth(binsWidth) {spectraPointsNbr=0;}


//We store asymptotic energies in containers of map type with a view to make a data binning
template<typename state_type>
void Spectra<state_type>::storeDataBinning(const state_type& x, const double& t, const double& weightIonization, const bool& unexpectedStop)
{

  //If the computation of the trajectory has been stopped unexpectedly, we do not consider the event
  if(unexpectedStop) return;
  
  //We compute the x value of the new point in the histogram
 int  range=int(asymptoticEnergy(x,t)/binsWidth);
  
//we insert values associated with the new event in the corresponding container
if(x[1]>=0)
 insertInMap(asymptEnergyUp, range, weightIonization);
else
 insertInMap(asymptEnergyDown, range, weightIonization);
}

//The following method is called by storeDataBinning and insert element <range,weightIonization> in map asympEnergy
template<typename state_type>
void Spectra<state_type>::insertInMap(std::map<int,double>&  asymptEnergy, const int& range, const double& weightIonization)
{

//we declare iterators which can iterate through differents elements of a container of type map
std::map<int,double>::iterator findRange=asymptEnergy.find(range);

  spectraPointsNbr++;
    
  if(findRange==asymptEnergy.end())
    {
      asymptEnergy[range]=weightIonization;
    }
  else
    {
      asymptEnergy[range]=asymptEnergy[range]+weightIonization;
    }
}
  

  //Finally we write all the data binning in a file
template<typename state_type>
void Spectra<state_type>::writeDataBinning(std::fstream& dataFile)
{
  
//we write values contained in asymptEnergyUp in file "dataFile"
 getFromMap(dataFile, asymptEnergyUp);
//we leave two lines break
dataFile<<" "<<std::endl;
dataFile<<" "<<std::endl;
//we write values contained in asymptEnergyDown in file "dataFile"
 getFromMap(dataFile, asymptEnergyDown);
 
}

//The following method is called by writeDataBinning and write all the elements range and weightIonization of map asymptEnergy in a file
template<typename state_type>
void Spectra<state_type>::getFromMap(std::fstream& dataFile, std::map<int,double>& asymptEnergy)
{
 std::map<int,double>::iterator it= asymptEnergy.begin();

  //We normalize the distribution law
  double sum=0.;
  for(1; it!=asymptEnergy.end(); it++)
    {
          sum=sum+it->second;
    }

  //Remove the following line if you wan a normalized distribution law
  sum=1.;

  //We write the histogram in a file	
  for(it=asymptEnergy.begin(); it!=asymptEnergy.end(); it++)
    {
      dataFile<<(it->first)*binsWidth<<" "<<log10( (it->second)/sum )<<std::endl;
    }

}


//We compute asymptotic energy
template<typename state_type>
double Spectra<state_type>::asymptoticEnergy(const state_type &x, const double& t)
{
  double Vsq=(x[3]-myField.vectPot('X',t))*(x[3]-myField.vectPot('X',t))+(x[4]-myField.vectPot('Y',t))*(x[4]-myField.vectPot('Y',t))+(x[5]-myField.vectPot('Z',t))*(x[5]-myField.vectPot('Z',t));
   double E=Vsq/2-myPotential.potentialEnergy(x);

  return E;

}



#endif
