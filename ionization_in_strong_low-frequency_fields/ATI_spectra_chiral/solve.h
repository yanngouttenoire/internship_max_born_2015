#ifndef SOLVE_H
#define SOLVE_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"

//Forward declaration
template<typename state_type>
class System;

template<typename state_type>
class Solve
{

public:

//We implement a method which solves an ode system using Runge-Kutta of order 5 with adaptive step-size algorithm 
//The adaptive step-size algorithm uses Runge-Kutta of order 4 embedded in RK5
//For more details: http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c16-2.pdf

void controlledRK5(System<state_type> &system, state_type &x, double &t, double &dt, double &error, double desiredErrorMin, double desiredErrorMax);


};



//We implement the rk5 iteration for the integration of ode
//We also implement a rk4 iteration embedded in the rk5 one. Then, we can compute the error and adapt the step size
//For more details: http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c16-2.pdf

template<typename state_type>
void Solve<state_type>::controlledRK5(System<state_type> &system, state_type &x, double &t, double &dt, double &error, double desiredErrorMin, double desiredErrorMax)
{

  //The iteration of yn is
  //yn+1=yn+c0*k0+c1*k1+c2*k2+c3*k3+c4*k4+c5*k5
  //where
  //k0=dt*f(xn,yn)
  //k1=dt*f(xn+a1*dt,yn+b10*k0)
  //k2=dt*f(xn+a2*dt,yn+b20*k0+b21*k1)
  //k3=dt*f(xn+a3*dt,yn+b30*k0+b31*k1+b32*k2)
  //k4=dt*f(xn+a4*dt,yn+b40*k0+b41*k1+b42*k2+b43*k3)
  //k5=dt*f(xn+a5*dt,yn+b50*k0+b51*k1+b52*k2+b53*k3+b54*k4)
  //an, bn and cn are called Butcher's coefficients

  //Temporary variables
  int i,j,p;

  double k[6][6];
  double xm[7][6];

  //First derivative of the ODE
  double dxdt[6][6];
 
  //Butcher's coefficients

  double a[6]={0.,1/5.,3/10.,3/5.,1.,7/8.};
  double b[6][5]={{0.,0.,0.,0.,0.},{1/5.,0.,0.,0.,0.},{3/40.,9/40.,0.,0.,0.},{3/10.,-9/10.,6/5.,0.,0.},{-11/54.,5/2.,-70/27.,35/27.,0.},{1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.}};
  
  //Fifth-order Runge-Kuta method
  double c5[6]={37/378.,0.,250/621.,125/594.,0.,512/1771.};
  //Embedded four-order Runge kutta method  
  double c4[6]={2825/27648.,0.,18575/48384.,13525/55296.,277/14336.,1/4.};


  //RK4 solution variable
  double x4[6];
  //RK5 solution variable
  double x5[6];

int kk=0;

 do
    {
     
  for(i=0; i<6; i++)
    {
      x4[i]=x[i];
      x5[i]=x[i];
    }

       for(i=0; i<6; i++)
	{
	  xm[0][i]=x[i];
	}
   
      for(j=1 ; j<=6; j++)
	{
	  system(xm[j-1],dxdt[j-1],t+dt*a[j-1]);
      
	  for(i=0; i<6; i++)
	    {
	      xm[j][i]=x[i];
	      k[j-1][i]=dt*dxdt[j-1][i];

	      for(p=0; p<j && j<6 ; p++)
		{
		  xm[j][i]=xm[j][i]+b[j][p]*k[p][i];
		}
	    }
	}

      for(i=0; i<6; i++)
	{
	  for(j=0; j<6; j++)
	    {
	      x4[i]=x4[i]+c4[j]*k[j][i];
	    }
	}

      for(i=0; i<6; i++)
	{
	  for(j=0; j<6; j++)
	    {
	      x5[i]=x5[i]+c5[j]*k[j][i];
	    }
	}
   
     //We set the "error" as the result of the difference between rk5 solution and rk4 solution
      error=fabs(pow(x5[0]*x5[0]+x5[1]*x5[1]+x5[2]*x5[2],1./2.)-pow(x4[0]*x4[0]+x4[1]*x4[1]+x4[2]*x4[2],1./2.))/pow(x5[0]*x5[0]+x5[1]*x5[1]+x5[2]*x5[2],1./2.);
   
    //Since error terms in rk4 and rk5 are respectly of order 5 and 6 according the step dt, the "error" we just computed scales as the step dt to the power of 5
    //Then, the new step as to be estimated as the following: newStep=oldStep*pow(desiredError/oldError,1/5.)

  if(error>desiredErrorMax)
	dt=dt*pow(desiredErrorMax/error,0.2);
    
  if(error<desiredErrorMin && error!=0.)
	dt=dt*pow(desiredErrorMin/error,0.2);

 kk=kk+1;
if(kk>100)
error=desiredErrorMax;
     }

  //We do not exit from the loop while error is not the one expected
  while(error>desiredErrorMax || error<desiredErrorMin && error!=0.);
     

  //We keep the rk5 solution
  for(i=0; i<6; i++)
    {
      x[i]=x5[i];
    }

  //We increment time
  t=t+dt;

}



#endif
