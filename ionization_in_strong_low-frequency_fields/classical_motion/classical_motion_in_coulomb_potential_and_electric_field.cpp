#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>


using namespace std;

double uaIntensity=3.5094451E16;
double uaTime=24.1888421562712E-18;
double uaEnergy=27.211608;
double lightSpeed=2.99792458E8;
double waveLenght=2E-6;
double fieldAmpl=sqrt(30*(1E13)/uaIntensity);
double pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
double IP=24.587/uaEnergy;
double opticalCycle=2.*M_PI/pulsation;
double pulseDuration=3.*opticalCycle;
double phase=0.;
double ellipticity=0.;
double Z=1.; //Positive charge of the nucleus


//Ordinary differential equations of the dynamic of an electron in coulomb potential 
void coulomb(double* q, double t, double* qp)
{
  double K=Z;

 //We introduce a spatial threshold below which we switch off the coulomb potential
  double thresholdLenght=1.;
  if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])<thresholdLenght)
    {
      qp[0]=q[3];
      qp[1]=q[4];
      qp[2]=q[5];
      qp[3]=0.;
      qp[4]=0.;
      qp[5]=0.;  
      return;
    }

  double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],3./2.);  //Do not forget the dots in 3./2. !
  qp[0]=q[3];
  qp[1]=q[4];
  qp[2]=q[5];
  qp[3]=-K*q[0]/a;
  qp[4]=-K*q[1]/a;
  qp[5]=-K*q[2]/a;

}

//We add the electric field
void coulomb_field(double* q, double t, double* qp)
{
  double K=Z;

    //We introduce a spatial threshold below which we switch off the coulomb potential
  double thresholdLenght=1.;
  if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])<thresholdLenght)
    {
      qp[0]=q[3];
      qp[1]=q[4];
      qp[2]=q[5];
      qp[3]=fieldAmpl*pow(sin(M_PI*t/2./pulseDuration),2)*cos(ellipticity)*sin(pulsation*t+phase);
      qp[4]=-fieldAmpl*pow(sin(M_PI*t/2./pulseDuration),2)*sin(ellipticity)*cos(pulsation*t+phase);
      qp[5]=0.;  
      return;
    }

  double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],3./2.);
  qp[0]=q[3];
  qp[1]=q[4];
  qp[2]=q[5];
  qp[3]=-K*q[0]/a+fieldAmpl*pow(sin(M_PI*t/2./pulseDuration),2)*cos(ellipticity)*sin(pulsation*t+phase);
  qp[4]=-K*q[1]/a-fieldAmpl*pow(sin(M_PI*t/2./pulseDuration),2)*sin(ellipticity)*cos(pulsation*t+phase);
  qp[5]=-K*q[2]/a;

}

//We implement the rk4 iteration for the integration of ode
void rk4(void (&diff)(double*,double,double*), double* q, double &t, double dt)
{
  int i;
  double qp[6];
  double qm[6];
  double k[4][6];

  diff(q,t,qp);
  for(i=0; i<6; i++)
    {
      k[0][i] = dt*qp[i];
      qm[i]=q[i]+k[0][i]/2.;
      t=t+dt/2.;
    }

  diff(qm,t,qp);
  for(i=0; i<6; i++)
    {
      k[1][i] = dt*qp[i];
      qm[i]=q[i]+k[1][i]/2.;
    }

  diff(qm,t,qp);
  for(i=0; i<6; i++)
    {
      k[2][i] = dt*qp[i];
      qm[i]=q[i]+k[2][i];
      t=t+dt/2.;
    }

  diff(qm,t,qp);	
  for(i=0; i<6; i++)
    {
      k[3][i] = dt*qp[i];
    }

  for(i=0; i<6; i++)
    {
      q[i]= q[i] + (k[0][i]+2.*k[1][i]+2.*k[2][i]+k[3][i])/6.;
    }
}

//We build the function which sets the initial conditions
void IC(double* q)
{

  //Position of the electron after tunneling
  double x_birth=sqrt(IP*IP-4*Z*abs(fieldAmpl*cos(phase)));
  x_birth=(x_birth+IP)/(fieldAmpl*cos(phase))/2.;         

  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldAmpl*cos(phase)/sqrt(2.*IP));  

  //Perpendicular velocity after tunneling
  double vPerp=4.*sigma_V*(2.*(double(rand())/double(RAND_MAX))-1.);

  //We display some informations
  cout<<"sigma_V= "<<sigma_V<<endl;
  cout<<"xbirth= "<<x_birth<<endl;
  cout<<"vPerp= "<<vPerp<<endl;
  cout<<" "<<endl;

  //We set the IC
  q[0]=x_birth;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
  q[4]=0.;
  q[5]=0.;

}

//Keplerian energy
double energy(double* q)
{
  return 0.5*(q[3]*q[3]+q[4]*q[4]+q[5]*q[5])-Z/pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],1./2.);
}

int main()
{

  //We display some informations
  cout<<" "<<endl;
  cout<<"fieldIntensity= "<<fieldAmpl<<endl;
  cout<<"opticalCycle= "<<2*M_PI/pulsation<<endl;
  cout<<"cycleNumber= "<<pulseDuration/opticalCycle<<endl;
  cout<<"ellipticity= "<<ellipticity<<endl;
  cout<<" "<<endl;

  //We open a file with a view to writing in it
  fstream dat("data.dat",ios::out);

  int n, i,j;

  //We fix the step
  double dt=0.001;

  //We fix the number of points
  int N=100000;

  //We create the array in which we will store the orbit
  double* q=new double[6];

  //We create the time variable
  double t;

  //We call a function which will fix the IC
  IC(q);
  t=0.;
  
  for(n=1; n<=N; n++)
    { 

      //We call the rk4 function which solve eq of the motion
      rk4(coulomb_field,q,t,dt);
    
      //We write the position of the electron in the phase space
      for (i=0; i<6; i++)
	{
	  dat<<q[i]<<" ";
	}
      dat<<" "<<endl;

      //We display some informations
      if(n%int(N/3)==0)
	{
	  cout<<"Energy= "<<energy(q)<<endl;
	  cout<<"Distance= "<<pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],1./2.)<<endl;
	  cout<<" "<<endl;
	}
    }


  //We plot a circle which somehow corresponds to the ionic core

 double u;
 
  dat<<" "<<endl;

  for(i=0; i<=100; i++)
    {
      u=double(i)*2.*M_PI/100.;
      dat<<cos(u)<<" "<<sin(u)<<" "<<0<<endl;
    }

  //We display the duration of the simulation in ua
  cout<<"t= "<<t<<endl;
 
    
  //We plot the orbit (and the "ionic circle") with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"splot 'data.dat' using 1:2:3 with line title 'Coulomb Orbit'"<<endl;

  gnu<<"pause -1"<<endl;
  gnu<<"set term postscript"<<endl;
  gnu<<"set output 'courbe.ps'"<<endl;
  gnu<<"replot"<<endl;
  gnu.close();

  system("gnuplot data.gnu");

  return 0;
}

