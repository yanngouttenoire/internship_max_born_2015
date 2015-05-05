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

void Solve::rk4(double dt)
{
  int i,p;
  double qm[6];
  double k[4][6];
  
  m_ode->diff(q,*t);
  for(i=0; i<6; i++)
    {
      k[0][i] = dt*m_ode->qp[i];
      qm[i]=q[i]+k[0][i]/2;
      *t=*t+dt/2.;
    }

  m_ode->diff(qm,*t);
  for(i=0; i<6; i++)
    {
      k[1][i] = dt*m_ode->qp[i];
      qm[i]=q[i]+k[1][i]/2;
    }

  m_ode->diff(qm,*t);
  for(i=0; i<6; i++)
    {
      k[2][i] = dt*m_ode->qp[i];
      qm[i]=q[i]+k[2][i];
      *t=*t+dt/2.;
    }

  m_ode->diff(qm,*t);
  for(i=0; i<6; i++)
    {
      k[3][i] = dt*m_ode->qp[i];
    }

  for(i=0; i<6; i++)
    {
      q[i]= q[i] + (k[0][i]+2.*k[1][i]+2.*k[2][i]+k[3][i])/6.;
    }
 
}

