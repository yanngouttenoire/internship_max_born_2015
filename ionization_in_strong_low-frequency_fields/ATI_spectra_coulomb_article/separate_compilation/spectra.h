#ifndef SPECTRA_H
#define SPECTRA_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

//We implement a class for computing the photo-electron spectrum as an histogram of asymptotic Energy according probability of ionization (weightIonization)
template<typename state_type>
class Spectra
{

//We declare a container which will contain the electron spectrum
std::map<int,double> asymptEnergy;

//We declare a variable for the bin interval width
double binsWidth;

public:

//Constructor
Spectra(double binsWidth);

//We compute the asymptotic energy
void asymptoticEnergy(const state_type& x);

//We set the histogram intervals width
void setBinsWidth(int binsNumber=100);

//We store asymptotic velocities in containers of map type
void storeDataBinning(const state_type& x, const double& weightIonization, const bool& inexpectedStop);

//Finally we write all the data binning in a file
void writeDataBinning(std::fstream& dataFile);


};

#endif
