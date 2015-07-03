#ifndef MOLECULE_H
#define MOLECULE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>

#include"electrostaticpotential.h"

enum enumOrientations{W,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3};

//Structure which can both represente a molecule orientation as both an enum type (myEnum) and a string type (myString)
struct moleculeOrientation
{

moleculeOrientation(enumOrientations myEnum) : myEnum(myEnum)
{
switch(myEnum)
 {
 case W :
  myString="W";
  break;
 case X1 :
  myString="X1";
  break;
 case X2 :
  myString="X2";
  break;
 case X3 :
  myString="X3";
  break;
 case Y1 :
  myString="Y1";
  break;
 case Y2 :
  myString="Y2";
  break;
 case Y3 :
  myString="Y3";
  break;
 case Z1 :
  myString="Z1";
  break;
 case Z2 :
  myString="Z2";
  break;
 case Z3 :
  myString="Z3";
  break;
 }  
}

enumOrientations operator()()
{
return myEnum;
}

enumOrientations myEnum;
std::string myString;

};

//Define characteristics and properties of an electrostatic potential for a given molecule
template<typename state_type> class Molecule : public ElectrostaticPotential<state_type>
{

 public:

  //We declare an orientation for the chiral molecule
  moleculeOrientation myOrientation;

  //Nuclei charges and covalent bond length
  //For orientation 'W'
  double c[4]; 
  double l[4]; 
  //For others orientations (1,2,3 means the atoms respectively along axis x,y,z; 0 is the one at the origin)
  double charge[4];
  double bondLength[4]; 
  //Some variables
  double inverseRadialDistanceCube[4];
  double sumInverseRadialDistanceCube;

  //Constructor
  Molecule(moleculeOrientation myOrientation=moleculeOrientation(W));
  
  //Necessary method for copy procedure when using polymophism (abstract parent class)
  Molecule* Clone()
  {
  return new Molecule(*this);
  }
  
  //Method which sets the molecule orientation
  void setOrientation(moleculeOrientation myOrientation);

  //Method which builds the molecule and sets the molecule orientation
  void setMoleculeOrientation(moleculeOrientation myOrientation);
  
  //Method which computes some quantities in advance
  void preparePotential(const state_type &x);

  //Methods which return the coulomb force on each components 
  double operator()(char component, const state_type &x);

  //Method which returns the electrostatic potential energy
  double potentialEnergy(const state_type &x);

};

template<typename state_type>
Molecule<state_type>::Molecule(moleculeOrientation myOrientation) : myOrientation(myOrientation)
{
  //Nuclei charges
  c[0]=0.4;
  c[1]=0.2;
  c[2]=0.2;
  c[3]=0.2;
  //Covalent bond length
  l[1]=2.0;
  l[2]=1.0;
  l[3]=3.0;
  
  setMoleculeOrientation(myOrientation.myEnum);

  //Ionization potential
  this->IP=0.32772257;
}

//Method which builds the molecule
template<typename state_type>
void Molecule<state_type>::setMoleculeOrientation(moleculeOrientation _myOrientation_)
{
myOrientation=_myOrientation_;
switch(myOrientation.myEnum)
 {
 case W :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[2];
  charge[3]=c[3];
  bondLength[1]=l[1];
  bondLength[2]=l[2];
  bondLength[3]=l[3];
  break;
 case X1 :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[3];
  charge[3]=c[2];
  bondLength[1]=l[1];
  bondLength[2]=l[3];
  bondLength[3]=-l[2];
  break;
 case X2 :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[2];
  charge[3]=c[3];
  bondLength[1]=l[1];
  bondLength[2]=-l[2];
  bondLength[3]=-l[3];
  break;
 case X3 :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[3];
  charge[3]=c[2];
  bondLength[1]=l[1];
  bondLength[2]=-l[3];
  bondLength[3]=l[2];
  break;
 case Y1 :
  charge[0]=c[0];
  charge[1]=c[3];
  charge[2]=c[2];
  charge[3]=c[1];
  bondLength[1]=l[3];
  bondLength[2]=l[2];
  bondLength[3]=-l[1];
  break;
 case Y2 :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[2];
  charge[3]=c[3];
  bondLength[1]=-l[1];
  bondLength[2]=l[2];
  bondLength[3]=-l[3];
  break;
 case Y3 :
  charge[0]=c[0];
  charge[1]=c[3];
  charge[2]=c[2];
  charge[3]=c[1];
  bondLength[1]=-l[3];
  bondLength[2]=l[2];
  bondLength[3]=l[1];
  break;
 case Z1 :
  charge[0]=c[0];
  charge[1]=c[2];
  charge[2]=c[1];
  charge[3]=c[3];
  bondLength[1]=-l[2];
  bondLength[2]=l[1];
  bondLength[3]=l[3];
  break;
 case Z2 :
  charge[0]=c[0];
  charge[1]=c[1];
  charge[2]=c[2];
  charge[3]=c[3];
  bondLength[1]=-l[1];
  bondLength[2]=-l[2];
  bondLength[3]=l[3];
  break;
 case Z3 :
  charge[0]=c[0];
  charge[1]=c[2];
  charge[2]=c[1];
  charge[3]=c[3];
  bondLength[1]=l[2];
  bondLength[2]=-l[1];
  bondLength[3]=l[3];
  break;
 }  
}


template<typename state_type>
void Molecule<state_type>::preparePotential(const state_type &x)
{

  inverseRadialDistanceCube[0] = charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,3./2.);
  inverseRadialDistanceCube[1] = charge[1]/pow((x[0]-bondLength[1])*(x[0]-bondLength[1])+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,3./2.) ;
  inverseRadialDistanceCube[2] = charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[2])*(x[1]-bondLength[2])+x[2]*x[2]+this->softParameter*this->softParameter,3./2.) ;
  inverseRadialDistanceCube[3] = charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[3])*(x[2]-bondLength[3])+this->softParameter*this->softParameter,3./2.) ;
  sumInverseRadialDistanceCube = inverseRadialDistanceCube[0]+inverseRadialDistanceCube[1]+inverseRadialDistanceCube[2]+inverseRadialDistanceCube[3];

}



template<typename state_type>
double Molecule<state_type>::operator()(char component, const state_type &x)
{

  switch(component)
    {
    case 'X' :
      return -x[0]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[1]*bondLength[1];

    case 'Y' :
      return -x[1]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[2]*bondLength[2];

    case 'Z' :
      return -x[2]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[3]*bondLength[3];
    }


}


template<typename state_type>
double Molecule<state_type>::potentialEnergy(const state_type& x)
{

  return -charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[1]/pow((x[0]-bondLength[1])*(x[0]-bondLength[1])+x[1]*x[1]+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[2]/pow(x[0]*x[0]+(x[1]-bondLength[2])*(x[1]-bondLength[2])+x[2]*x[2]+this->softParameter*this->softParameter,1./2.)-charge[3]/pow(x[0]*x[0]+x[1]*x[1]+(x[2]-bondLength[3])*(x[2]-bondLength[3])+this->softParameter*this->softParameter,1./2.);

}



#endif
