#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<iomanip>


using namespace std;

double uaIntensity=3.5094451E16;
double uaTime=24.1888421562712E-18;
double uaEnergy=27.211608;
double lightSpeed=2.99792458E8;
double waveLenght=2E-6;
double fieldAmpl=sqrt(30*(1E13)/uaIntensity);
double pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
double IP=0.5;
double opticalCycle=2.*M_PI/pulsation;
double pulseDuration=3.*opticalCycle;
double phase=0.;
double ellipticity=0.;
double Z=1.; //Positive charge of the nucleus


//Ordinary differential equations of the dynamic of an electron in coulomb potential 
void coulomb(double* q, double t, double* qp)
{
  double K=Z;

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

  //yn+1=yn+c0*k0+c1*k1+c2*k2+c3*k3
  //k0=dt*f(xn,yn)
  //k1=dt*f(xn+a1*dt,yn+b10*k0)
  //k2=dt*f(xn+a2*dt,yn+b20*k0+b21*k1)
  //k3=dt*f(xn+a3*dt,yn+b30*k0+b31*k1+b32*k2)

  int i,j,p;
  int order=4;
  double qp[order][6];
  double qm[order][6];
  double k[order][6];

  double a[4]={0.,1/2.,1/2.,1.};
  double b[4][4]={{0.,0.,0.,0.},{1/2.,0,0,0},{0.,1/2.,0.,0.},{0.,0.,1.,0}};
  double c[4]={1/6.,1/3.,1/3.,1/6.};


  for(i=0; i<6; i++)
    {
    qm[0][i]=q[i];
    }

  for(j=1 ; j<=order; j++)
    {
cout<<"qm= "<<qm[j][i]<<" ";
      diff(qm[j-1],t+dt*a[j-1],qp[j-1]);
      
      if(j!=order)
	{

	  for(i=0; i<6; i++)
	    {
	      cout<<"qp= "<<qp[j][i]<<" ";
	      qm[j][i]=q[i];
	      k[j][i]=dt*qp[j][i];

	      for(p=0; p<j; p++)
		{
		  qm[j][i]=qm[j][i]+b[j][p]*k[p][i];
		}
	    }
	  cout<<" "<<endl;
 
	}

    }
cout<<" "<<endl;

  for(i=0; i<6; i++)
    {

      for(j=0; j<order; j++)
	{
	  cout<<k[j][i]<<" ";
	  q[i]= q[i] + c[j]*k[j][i];;

	}
      cout<<" "<<endl;
    }
 cout<<" "<<endl;

  t=t+dt;
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
  cout<<" "<<endl;

  //We set the IC
  q[0]=2.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
  q[4]=0.1;
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
  double dt=0.0001;

  //We fix the number of points
  int N=100000;

  //We create the array in which we will store the orbit
  double* q=new double[6];

  //We create the time variable
  double t;

  //We call a function which will fix the IC
  IC(q);
  t=0.;
  
  time_t start=time(NULL);

  for(n=1; n<=N; n++)
    { 

      //We call the rk4 function which solve eq of the motion
      rk4(coulomb,q,t,dt);
    
      //We write the position of the electron in the phase space
      for (i=0; i<6; i++)
	{
	  dat<<q[i]<<" ";
	}
      dat<<" "<<endl;

      //We display some informations
      if(n%100==0)
	{
	  cout<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F";
	  cout<<"Energy= "<<energy(q)<<endl;
	  cout<<"Distance= "<<pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],1./2.)<<endl;
	  time_t end=time(NULL);
	  //We display the elapsed time and the load bar
	  cout<<"Elapsed time= "<<difftime(end,start)<<endl;
	  cout<<"\r"<<"progression= "<<setprecision(3)<<(double)(n)/(double)(N)*100.<<"%"<<setprecision(cout.precision())<<"                "<<endl;
	}
    }

  //We plot a circle which somehow corresponds to the ionic core

 double u;
 
  dat<<" "<<endl;

  for(i=0; i<=100; i++)
    {
      u=double(i)*2.*M_PI/100.;
      dat<<0.01*cos(u)<<" "<<0.01*sin(u)<<" "<<0<<endl;
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

