#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"solve.h" 

Solve::Solve(double* a_q, double *a_t, Ode* a_ode)
{
  q=a_q;
  m_ode=a_ode;
  t=a_t;
}

void Solve::setOde(Ode* a_ode)
{
  m_ode=a_ode;
}

//We implement the rk4 iteration for the integration of ode
void Solve::rk4(double dt)
{

  //yn+1=yn+c0*k0+c1*k1+c2*k2+c3*k3
  //k0=dt*f(xn,yn)
  //k1=dt*f(xn+a1*dt,yn+b10*k0)
  //k2=dt*f(xn+a2*dt,yn+b20*k0+b21*k1)
  //k3=dt*f(xn+a3*dt,yn+b30*k0+b31*k1+b32*k2)
 
  int i,j,p;
  double qp[4][6];
  double qm[4][6];
  double k[4][6];

  double a[4]={0.,1/2.,1/2.,1.};
  double b[4][4]={{0.,0.,0.,0.},{1/2.,0.,0.,0.},{0.,1/2.,0.,0.},{0.,0.,1.,0.}};
  double c[4]={1/6.,1/3.,1/3.,1/6.};
  
  for(i=0; i<6; i++)
    {
      qm[0][i]=q[i];
    }
   
  for(j=1 ; j<=4; j++)
    {
      m_ode->diff(qm[j-1],*t+dt*a[j-1]);
      
      for(i=0; i<6; i++)
	{
	  k[j-1][i]=dt*m_ode->qp[i];
	  for(p=j-1; p<j && j<4 ; p++)
	    {
	      qm[j][i]=q[i];
	      qm[j][i]=qm[j][i]+b[j][p]*k[p][i];
	    }
	}
	
    }


  for(i=0; i<6; i++)
    {
      for(j=0; j<4; j++)
	{
	  q[i]=q[i]+c[j]*k[j][i];
	}
    }

 *t=*t+dt;

}
