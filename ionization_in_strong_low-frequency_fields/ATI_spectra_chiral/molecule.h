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

enum enumOrientations{W,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,X3Z2,X3Z3};

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
 case X3Z2 :
  myString="X3Z2";
  break;
 case X3Z3 :
  myString="X3Z3";
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
  moleculeOrientation toBeRemoved;
  int myOrientation;

  //Nuclei charges and covalent bond length
  double charge[4];
  double bondLength[4]; 
  //Some variables
  double inverseRadialDistanceCube[4];
  double sumInverseRadialDistanceCube;
  double X[4], Y[4], Z[4];
  double xi,yi,zi;
  double lebedevWeight;

  //Constructor
  Molecule(moleculeOrientation toBeRemoved=moleculeOrientation(W));
  
  //Necessary method for copy procedure when using polymophism (abstract parent class)
  Molecule* Clone()
  {
  return new Molecule(*this);
  }
  
  void setLebedevOrientation(int myOrientation);

  //Method which computes some quantities in advance
  void preparePotential(const state_type &x);

  //Methods which return the coulomb force on each components 
  double operator()(char component, const state_type &x);

  //Method which returns the electrostatic potential energy
  double potentialEnergy(const state_type &x);

};

template<typename state_type>
Molecule<state_type>::Molecule(moleculeOrientation toBeRemoved): toBeRemoved(toBeRemoved)
{
  //Nuclei charges
  charge[0]=0.2;
  charge[1]=0.3;
  charge[2]=0.3;
  charge[3]=0.2;
  //Covalent bond length
  bondLength[1]=2.0;
  bondLength[2]=1.0;
  bondLength[3]=3.0;
  setLebedevOrientation(3);

  //Ionization potential
  //this->IP=0.32772257;
  this->IP=0.31631047;
}

//

template<typename state_type>
void Molecule<state_type>::setLebedevOrientation(int _myOrientation_)
{
myOrientation=_myOrientation_;
std::ifstream lebedev_table("lebedev_table", std::ios::in);

if(lebedev_table)  
{      
for(int i=1; i<=myOrientation; i++)
lebedev_table >> xi >> yi >> zi >> lebedevWeight;
 
lebedev_table.close();
}

else 
std::cerr << "Error when trying to open lebedev_table" << std::endl;
std::cout<<xi<<" "<<yi<<" "<<zi<<std::endl;

if(xi==0 && yi==0 && zi==1)
{
X[1]=bondLength[1];
Y[1]=0.;
Z[1]=0.;

X[2]=0.;
Y[2]=0.;
Z[2]=-bondLength[2];

X[3]=0.;
Y[3]=bondLength[3];
Z[3]=0.;
}
if(xi==0 && yi==0 && zi==-1)
{
X[1]=bondLength[1];
Y[1]=0.;
Z[1]=0.;

X[2]=0.;
Y[2]=0.;
Z[2]=bondLength[2];

X[3]=0.;
Y[3]=-bondLength[3];
Z[3]=0.;
}
if(xi!=0 || yi!=0)
{
X[1]=bondLength[1]*yi/sqrt(xi*xi+yi*yi);
Y[1]=bondLength[1]*xi;
Z[1]=-bondLength[1]*xi*zi/sqrt(xi*xi+yi*yi);

X[2]=-bondLength[2]*xi/sqrt(xi*xi+yi*yi);
Y[2]=bondLength[2]*yi;
Z[2]=-bondLength[2]*yi*zi/sqrt(xi*xi+yi*yi);

X[3]=0.;
Y[3]=bondLength[3]*zi;
Z[3]=bondLength[3]*sqrt(xi*xi+yi*yi);
}
}

template<typename state_type>
void Molecule<state_type>::preparePotential(const state_type &x)
{

  inverseRadialDistanceCube[0] = charge[0]/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],3./2.);
  inverseRadialDistanceCube[1] = charge[1]/pow((x[0]-X[1])*(x[0]-X[1])+(x[1]-Y[1])*(x[1]-Y[1])+(x[2]-Z[1])*(x[2]-Z[1]),3./2.) ;
  inverseRadialDistanceCube[2] = charge[2]/pow((x[0]-X[2])*(x[0]-X[2])+(x[1]-Y[2])*(x[1]-Y[2])+(x[2]-Z[2])*(x[2]-Z[2]),3./2.) ;
  inverseRadialDistanceCube[3] = charge[3]/pow((x[0]-X[3])*(x[0]-X[3])+(x[1]-Y[3])*(x[1]-Y[3])+(x[2]-Z[3])*(x[2]-Z[3]),3./2.) ;
  sumInverseRadialDistanceCube = inverseRadialDistanceCube[0]+inverseRadialDistanceCube[1]+inverseRadialDistanceCube[2]+inverseRadialDistanceCube[3];

}



template<typename state_type>
double Molecule<state_type>::operator()(char component, const state_type &x)
{

      switch(component)
    {
    case 'X' :
      return -x[0]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[1]*X[1]+inverseRadialDistanceCube[2]*X[2]+inverseRadialDistanceCube[3]*X[3];

    case 'Y' :
      return -x[1]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[1]*Y[1]+inverseRadialDistanceCube[2]*Y[2]+inverseRadialDistanceCube[3]*Y[3];

    case 'Z' :
      return -x[2]*sumInverseRadialDistanceCube+inverseRadialDistanceCube[1]*Z[1]+inverseRadialDistanceCube[2]*Z[2]+inverseRadialDistanceCube[3]*Z[3];
    }



}


template<typename state_type>
double Molecule<state_type>::potentialEnergy(const state_type& x)
{

    return -charge[0]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-charge[1]/sqrt((x[0]-X[1])*(x[0]-X[1])+(x[1]-Y[1])*(x[1]-Y[1])+(x[2]-Z[1])*(x[2]-Z[1]))-charge[2]/sqrt((x[0]-X[2])*(x[0]-X[2])+(x[1]-Y[2])*(x[1]-Y[2])+(x[2]-Z[2])*(x[2]-Z[2]))-charge[3]/sqrt((x[0]-X[3])*(x[0]-X[3])+(x[1]-Y[3])*(x[1]-Y[3])+(x[2]-Z[3])*(x[2]-Z[3]));

}



#endif
