#ifndef SPECTRA_H
#define SPECTRA_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<algorithm>

#include"electricfield.h"
#include"electrostaticpotential.h"
#include"ic.h"

//We introduce to different way of detection for the electron: UP and DOWN
enum detectionType {UP, DOWN};

//We implement a class for computing the photo-electron spectrum as an histogram of asymptotic Energy according probability of ionization (weightIonization)
template<typename state_type>
class Spectra
{

 public:

  //We declare two containers of map type, they can contain pairs which are couple of objects: first the rangeEnergy, second the ionization probability value

  //We declare one for electrons which velocity vector is oriented along positive values of y
  std::map<int,double> asympEnergyUP;

  //And one for electrons which velocity vector is oriented along negative values of y
  std::map<int,double> asympEnergyDOWN;
  
   //And one for initial phase of the field 
  std::map<int,double> initialPhase;
  
  //And one for initial velocity component along XZ
  std::map<int,double> initialVXZ;
  
  //And one for initial velocity component along Y
  std::map<int,double> initialVY;

  //We declare the horizontal variable in the histogram
  int rangeEnergy;
  
  //We declare the asymptotic energy
  double asympEnergy;

  //We declare a variable for the bin interval width
  double binsWidthEnergy;
  double binsWidthPhase;
  double binsWidthVelocity;
    
    //We declare the conversion ratio from au to the desired unit
   int energyUnit;
  
  //We declare the angle between the vector velocity and the laser polarization within which we detect electrons
  double angleDetection;

  //We declare a counter for the number of events stored
  int electronsDetectedNbr;
  
  //We declare a counter for the number of electrons detected UP
  int electronsDetectedUPNbr;
  
  //We declare a counter for the number of electrons detected DOWN
  int electronsDetectedDOWNNbr;

  //We declare a counter for the number of events not stored because the electron always trapped by the atom after the pulse
  int trappedElectronNbr;

  //We declare a counter which counts how many events have not been accepted because the angle between the velocity vector and the laser polarization was too large
  int angleTooLargeNbr;
  

  //We declare an object of type ElectrostaticPotential
  ElectrostaticPotential<state_type> *myPotential;

  //We declare an object of type ElectricField
  ElectricField myField;

  //We declare an object of type IC
  IC<state_type> *myIC;


  //Constructor
  Spectra(ElectrostaticPotential<state_type> *myPotential, ElectricField myField, IC<state_type> *myIC, double angleDetection=180., double binsWidthEnergy=0.1);

  //We compute the asymptotic energy
  double asymptoticEnergy(const state_type& x, const double& t);
  
  //We return a bool which tell us if a trajectory has a good profile
  bool hasTrajectoryGoodProfile(const state_type& x, const double& t, const bool& unexpectedStop);
  
  //We return a int (0 or 1) which tell us which is the profile of the trajectory
  detectionType whichProfile(const state_type& x, const double& t);

  //We store asymptotic velocities and weight ionization in containers of map type
  void storeDataBinning(const state_type& x, const double& t, const double& weightIonization, const double& phase, const double& vXZPerp, const double& vYPerp, const bool& unexpectedStop);
    
  //We store PES
  void storePES(const state_type& x, const double& t, const double& weightIonization);
  
    //We store initial conditions
  void storeIC(const state_type& x, const double& t, const double& weightIonization, const double& phase, const double& vXZPerp, const double& vYPerp);
  
  //The following method is called by storeDataBinning and insert element <range,weightIonization> in map asympEnergy
  void insertInMap(std::map<int,double>& asympEnergy, const int& range, const double& weightIonization);

  //The following method gather together all the data binning in one
  void mergeSpectra(std::vector<Spectra<state_type> >& mySpectra);

  //Finally we write all the data binning in a file
  void writeDataBinning(std::fstream& PESFile,std::fstream& ICFile);

  //The following method is called by writeDataBinning and write all the elements rangeEnergy and weightIonization of map asympEnergy in a file
  void getFromMap(std::fstream& DataFile, std::map<int,double>& asympEnergy,const double& binsWidth);
};


//We set the histogram intervals width
template<typename state_type>
Spectra<state_type>::Spectra(ElectrostaticPotential<state_type> *myPotential, ElectricField myField, IC<state_type> *a_myIC, double angleDetection, double binsWidthEnergy) : myPotential(myPotential), myField(myField), myIC(a_myIC), angleDetection(angleDetection), binsWidthEnergy(binsWidthEnergy) 
{
  electronsDetectedNbr=0;
  electronsDetectedUPNbr=0;
  electronsDetectedDOWNNbr=0;    
  trappedElectronNbr=0;
  angleTooLargeNbr=0;
  energyUnit=27.2;
  
  binsWidthPhase=0.1;
  binsWidthVelocity=0.01;
}

  //We return a bool which tell us if a trajectory has a good profile
template<typename state_type>
bool Spectra<state_type>::hasTrajectoryGoodProfile(const state_type& x, const double& t, const bool& unexpectedStop)
{
  //If the computation of the trajectory has been stopped unexpectedly, we do not consider the event
  if(unexpectedStop) return false;

  //If the energy of the electron is negative, the electron is not free and we do not consider the event
  if(asympEnergy<0) 
    {
      trappedElectronNbr++;
      return false;
    }
  
  //we insert values associated with the new event in the corresponding container according if the electron propagated in y>0 or y<0
  //We consider electron differently depending if they are detected along the polarization of the field or not
  if(fabs(atan(sqrt(x[3]*x[3]+x[4]*x[4])/x[5]))*180./M_PI<=angleDetection)
    {
    return true;
    }
   else 
   {
    return false;
   angleTooLargeNbr++;
   }
}

  //We return a int (0 or 1) which tell us which is the profile of the trajectory
template<typename state_type>
detectionType Spectra<state_type>::whichProfile(const state_type& x, const double& t)
{
 //if(x[2]*myField('Z',myIC->tBirth)>=0)
 if(x[1]>=0)
 return UP;
 else
 return DOWN;
}

//We store asymptotic energies in containers of map type with a view to make a data binning
template<typename state_type>
void Spectra<state_type>::storeDataBinning(const state_type& x, const double& t, const double& weightIonization, const double& phase, const double& vXZPerp, const double& vYPerp, const bool& unexpectedStop)
{
	asymptoticEnergy(x,t);

	if(hasTrajectoryGoodProfile(x, t, unexpectedStop)==true)
	 {
         electronsDetectedNbr++; 
        storePES(x,t,weightIonization);
 	storeIC(x,t,weightIonization,phase, vXZPerp, vYPerp);      	        
	 }
   
}

  //We store PES
template<typename state_type>
void Spectra<state_type>::storePES(const state_type& x, const double& t, const double& weightIonization)
{
	rangeEnergy=int(asympEnergy*energyUnit/binsWidthEnergy);
	
	 if(whichProfile(x,t)==UP)
     	        {
		insertInMap(asympEnergyUP, rangeEnergy, weightIonization);
      	        electronsDetectedUPNbr++;
      	        }
      	 else
     	        {
		insertInMap(asympEnergyDOWN, rangeEnergy, weightIonization);
      	        electronsDetectedDOWNNbr++;		      	        
      	        }

}

//We store initial conditions
template<typename state_type>
void Spectra<state_type>::storeIC(const state_type& x, const double& t, const double& weightIonization, const double& phase, const double& vXZPerp, const double& vYPerp)
{
    double energy=asympEnergy*27.2;
  //Angle between the velocity vector of the electron and the electric vector force
  int rangePhase=int(phase/binsWidthPhase);
  int rangeVXZ=int(vXZPerp/binsWidthVelocity);
  int rangeVY=int(vYPerp/binsWidthVelocity);
    
if(energy>0.15 && energy<0.80)
{
  insertInMap(initialPhase, rangePhase, weightIonization);
  insertInMap(initialVXZ, rangeVXZ, weightIonization);
  insertInMap(initialVY, rangeVY, weightIonization);
}
  
}
  
//The following method is called by storeDataBinning and insert element <range,weightIonization> in map asympEnergy
template<typename state_type>
void Spectra<state_type>::insertInMap(std::map<int,double>&  asympEnergy, const int& range, const double& weightIonization)
{

  //we declare iterators which can iterate through differents elements of a container of type map
  std::map<int,double>::iterator findRange=asympEnergy.find(range);

  if(findRange==asympEnergy.end())
    {
      asympEnergy[range]=weightIonization;
    }
  else
    {
      asympEnergy[range]=asympEnergy[range]+weightIonization;
    }

}

 //The following method gather together all the data binning in one
template<typename state_type>
void Spectra<state_type>::mergeSpectra(std::vector<Spectra<state_type> > &mySpectra)
{

/******************************We merge the data binning******************************/
  int range;
  double weightIonization;
  std::map<int,double>::iterator itmap;
  typename std::vector<Spectra<state_type> >::iterator itSpectra;

//We start from 1 and not from 0 because the current object is the one owned by the master thread
   for(itSpectra=mySpectra.begin(); itSpectra!=mySpectra.end(); itSpectra++)  
    {
     for(itmap=(itSpectra->asympEnergyUP).begin(); itmap!=(itSpectra->asympEnergyUP).end(); itmap++)
     {
       range=itmap->first;
       weightIonization=itmap->second;    
       insertInMap(asympEnergyUP, range, weightIonization);
     }  
   
      for(itmap=(itSpectra->asympEnergyDOWN).begin(); itmap!=(itSpectra->asympEnergyDOWN).end(); itmap++)
     {
       range=itmap->first;
       weightIonization=itmap->second;    
       insertInMap(asympEnergyDOWN, range, weightIonization);
     }  

     for(itmap=(itSpectra->initialPhase).begin(); itmap!=(itSpectra->initialPhase).end(); itmap++)
     {
       range=itmap->first;
       weightIonization=itmap->second;    
       insertInMap(initialPhase, range, weightIonization);
     }  
      for(itmap=(itSpectra->initialVXZ).begin(); itmap!=(itSpectra->initialVXZ).end(); itmap++)
     {
       range=itmap->first;
       weightIonization=itmap->second;    
       insertInMap(initialVXZ, range, weightIonization);
     }  
      for(itmap=(itSpectra->initialVY).begin(); itmap!=(itSpectra->initialVY).end(); itmap++)
     {
       range=itmap->first;
       weightIonization=itmap->second;    
       insertInMap(initialVY, range, weightIonization);
     }  
      
   }
   
/******************************We merge the counters******************************/

   for(itSpectra=mySpectra.begin(); itSpectra!=mySpectra.end(); itSpectra++)  
    {
     electronsDetectedNbr+=itSpectra->electronsDetectedNbr;
     electronsDetectedUPNbr+=itSpectra->electronsDetectedUPNbr;
     electronsDetectedDOWNNbr+=itSpectra->electronsDetectedDOWNNbr; 
     trappedElectronNbr+=itSpectra->trappedElectronNbr;
     angleTooLargeNbr+=itSpectra->angleTooLargeNbr;
    }

}


//Finally we write all the data binning in a file
template<typename state_type>
void Spectra<state_type>::writeDataBinning(std::fstream& PESFile,std::fstream& ICFile)
{
  
  //we write values contained in asympEnergyUP in file "PESFile"
  getFromMap(PESFile, asympEnergyUP,binsWidthEnergy);
  //we leave two lines break
  PESFile<<" "<<std::endl;
  PESFile<<" "<<std::endl;
  //we write values contained in asympEnergyDOWN in file "PESFile"
  getFromMap(PESFile, asympEnergyDOWN,binsWidthEnergy);
  
  //we write values contained in asympEnergyUP in file "PESFile"
  getFromMap(ICFile, initialPhase,binsWidthPhase);
  ICFile<<" "<<std::endl;
  ICFile<<" "<<std::endl;
  //we write values contained in asympEnergyDOWN in file "PESFile"
  getFromMap(ICFile, initialVXZ,binsWidthVelocity);
  ICFile<<" "<<std::endl;
  ICFile<<" "<<std::endl;
  //we write values contained in asympEnergyDOWN in file "PESFile"
  getFromMap(ICFile, initialVY,binsWidthVelocity);
  ICFile<<" "<<std::endl;
  ICFile<<" "<<std::endl;
 
}


//The following method is called by writeDataBinning and write all the elements rangeEnergy and weightIonization of map asympEnergy in a file
template<typename state_type>
void Spectra<state_type>::getFromMap(std::fstream& dataFile, std::map<int,double>& asympEnergy, const double& binsWidth)
{
  std::map<int,double>::iterator it= asympEnergy.begin();

  //We normalize the distribution law
  double sum=0.;
  for(1; it!=asympEnergy.end(); it++)
    {
      sum=sum+it->second;
    }

  //Remove the following line if you wan a normalized distribution law
  sum=1.;

  //We write the histogram in a file	
  for(it=asympEnergy.begin(); it!=asympEnergy.end(); it++)
    {
      dataFile<<(it->first)*binsWidth<<" "<<(it->second)/sum<<std::endl;
    }

}


//We compute asymptotic energy
template<typename state_type>
double Spectra<state_type>::asymptoticEnergy(const state_type &x, const double& t)
{
  double Vsq=(x[3]-myField.vectPot('X',t))*(x[3]-myField.vectPot('X',t))+(x[4]-myField.vectPot('Y',t))*(x[4]-myField.vectPot('Y',t))+(x[5]-myField.vectPot('Z',t))*(x[5]-myField.vectPot('Z',t));

  asympEnergy=Vsq/2+myPotential->potentialEnergy(x);

  return asympEnergy;


}



#endif
