#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"
#include"parameters.h"

using namespace std;

template<typename state_type>
System<state_type>::System(string whichField="withoutField") : whichField(whichField)
  { 
    K=Z;
    double field[3]={0.,0.,0.};
    double softParameter=0.;
  }

template<typename state_type>
    void System<state_type>::operator() (const state_type &x , state_type &dxdt , const double  t)
    {      

    //If it has been requested, we switch on the electric field
 if(whichField=="withField")
      {
	field[0]=0.;
	field[1]=0.;
	field[2]=-fieldAmpl*cos(pulsation*t+phase);
      }

    double a=pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter,3./2.); 
    dxdt[0]=x[3];
    dxdt[1]=x[4];
    dxdt[2]=x[5];
    dxdt[3]=-K*x[0]/a+field[0];
    dxdt[4]=-K*x[1]/a+field[1];
    dxdt[5]=-K*x[2]/a+field[2];
    }



template class System <vector<double> >;
template class System <double* >;
//For more details about templates, see https://www.cs.umd.edu/class/fall2002/cmsc214/Projects/P2/proj2.temp.html
