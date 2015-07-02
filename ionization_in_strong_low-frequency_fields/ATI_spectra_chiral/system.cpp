#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"


using namespace std;

template<typename state_type, typename potential_type>
System<state_type, potential_type>::System(/*ElectrostaticPotential<state_type>*/ potential_type& _myPotential_, ElectricField &myField) : myPotential(_myPotential_), myField(myField)
{
  //ElectrostaticPotential<state_type> *myPotential=_myPotential_.Clone();
}
/*
//Copy constructor (in order to use System as firstprivate in OPEN MP)
template<typename state_type, typename potential_type>
System<state_type, potential_type>::System(System const & source)
{
  //Copy constructor in ElectrostaticPotential is the default copy constructor since no data member is a pointer; however since it uses polymorphism, it needs the method "Clone" (See in ElectrostaticPotential and derived classes)
  ElectrostaticPotential<state_type> *myPotential=(source.myPotential)->Clone();
  
  //Copy constructor in ElectricField is the default copy constructor since no data member is a pointer
  ElectricField myField;
  myField=source.myField;
}

//Copy assignement operator
template<typename state_type, typename potential_type>
System<state_type, potential_type> & System<state_type, potential_type>::operator = (const System & source)
{
  ElectrostaticPotential<state_type> *myPotential=(source.myPotential)->Clone();
  ElectricField myField;
  myField=source.myField;
  return *this;
}  
  */
//Destructor
template<typename state_type, typename potential_type>
System<state_type, potential_type>::~System()
{
//delete myPotential;
}

template<typename state_type, typename potential_type>
void System<state_type, potential_type>::operator() (const state_type &x , state_type &dxdt , const double&  t)
{    

  {
    //myPotential->preparePotential(x);
    myPotential.preparePotential(x);

    dxdt[0]=x[3];
    dxdt[1]=x[4];
    dxdt[2]=x[5];
 /*  
    dxdt[3]=(*myPotential)('X',x)-myField('X',t);
    dxdt[4]=(*myPotential)('Y',x)-myField('Y',t);
    dxdt[5]=(*myPotential)('Z',x)-myField('Z',t);
  */
    dxdt[3]=(myPotential)('X',x)-myField('X',t);
    dxdt[4]=(myPotential)('Y',x)-myField('Y',t);
    dxdt[5]=(myPotential)('Z',x)-myField('Z',t);
  
  
  }
}

  

template class System <vector<double>, Hydrogen< vector<double> > >;
template class System <vector<double>, Molecule< vector<double> > >;
template class System <double[6], Hydrogen<double[6]> >;
template class System <double[6], Molecule<double[6]> >;
//For more details about templates, see https://www.cs.umd.edu/class/fall2002/cmsc214/Projects/P2/proj2.temp.html
