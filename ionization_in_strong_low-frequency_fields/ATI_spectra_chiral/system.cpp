#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"


using namespace std;

template<typename state_type>
System<state_type>::System(ElectrostaticPotential<state_type> *myPotential, ElectricField &myField) : myPotential(myPotential), myField(myField) {}

template<typename state_type>
void System<state_type>::operator() (const state_type &x , state_type &dxdt , const double&  t)
{      
#pragma omp critical
  {
    myPotential->preparePotential(x);

    dxdt[0]=x[3];
    dxdt[1]=x[4];
    dxdt[2]=x[5];
    dxdt[3]=(*myPotential)('X',x)-myField('X',t);
    dxdt[4]=(*myPotential)('Y',x)-myField('Y',t);
    dxdt[5]=(*myPotential)('Z',x)-myField('Z',t);
  }
}



template class System <vector<double> >;
template class System <double[6]>;
//For more details about templates, see https://www.cs.umd.edu/class/fall2002/cmsc214/Projects/P2/proj2.temp.html
