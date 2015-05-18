#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>

using namespace std;

//Plateau in above-threashold-ionization spectra and chaotic behavior in rescattering processes
//Bambi Hu, Jie Liu, Shi-gang Chen

//VARIABLES DECLARATION

double uaIntensity=3.5094451E16;
double uaTime=24.1888421562712E-18;
double uaEnergy=27.211608;
double lightSpeed=2.99792458E8;
double waveLenght=1064E-9;
double fieldAmpl=sqrt(13*(1E13)/uaIntensity);
double pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
double IP=13.605804/uaEnergy;
double opticalCycle=2.*M_PI/pulsation;
double pulseDuration=3.*opticalCycle;
double phase=0;
double ellipticity=0.1;
double Z=1.; //Positive charge of the nucleus

//We declare a softening parameter for the soft coulomb potential
double softParameter=0.1;

//We consider the rotating frame of reference (xPrime, yPrime, zPrime) where the electric field holds the xPrime direction
double vPerp, weight;
int iField, iVPerp=0;
int nField=100, nVPerp=100;

//We fix the step
double dt=0.00001; 

//We declare runge kutta error, its max, its min and the desired error (desired accuracy) 
double error;
double errorMax=1E-25;
double desiredErrorMax=1E-10;
double desiredErrorMin=desiredErrorMax/10.;

//We declare the number of points for trajectories
int nBirth;
int nEnd;

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare a variable for the electron position after tunneling
double rhoBirth;
double rhoBirthOld;

//We create the time variable
double t;

//We create the array in which we will store the orbit
double* q=new double[6];

//We declare containers which will contains electron spectra
map<int,double> asymptVel;

//Bins width
double binsWidth;

//We open files with a view to writing in them 

fstream dat("data.dat",ios::out);


//FUNCTIONS DECLARATION


//Ordinary differential equations of the dynamic of an electron in coulomb potential 
void coulomb(double* q, double t, double* qp)
{
  double K=Z;

  double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter*softParameter,3./2.); 
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
 
  double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter*softParameter,3./2.);
  qp[0]=q[3];
  qp[1]=q[4];
  qp[2]=q[5];
  qp[3]=-K*q[0]/a-fieldAmpl*cos(pulsation*t+phase);
  qp[4]=-K*q[1]/a;
  qp[5]=-K*q[2]/a;

}


//We implement the rk5 iteration for the integration of ode
//We also imprement a embedded rk4 with a view to computing the error and then to adapting the step size
//http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c16-2.pdf
void rk5(void (&diff)(double*,double,double*), double* q, double &t, double &dt, double &error, double desiredErrorMin, double desiredErrorMax)
{

  //yn+1=yn+c0*k0+c1*k1+c2*k2+c3*k3+c4*k4+c5*k5
  //k0=dt*f(xn,yn)
  //k1=dt*f(xn+a1*dt,yn+b10*k0)
  //k2=dt*f(xn+a2*dt,yn+b20*k0+b21*k1)
  //k3=dt*f(xn+a3*dt,yn+b30*k0+b31*k1+b32*k2)
  //k4=dt*f(xn+a4*dt,yn+b40*k0+b41*k1+b42*k2+b43*k3)
  //k5=dt*f(xn+a5*dt,yn+b50*k0+b51*k1+b52*k2+b53*k3+b54*k4)

  int i,j,p;
  double qp[6][6];
  double qm[7][6];
  double k[6][6];

  double a[6]={0.,1/5.,3/10.,3/5.,1.,7/8.};
  double b[6][5]={{0.,0.,0.,0.,0.},{1/5.,0.,0.,0.,0.},{3/40.,9/40.,0.,0.,0.},{3/10.,-9/10.,6/5.,0.,0.},{-11/54.,5/2.,-70/27.,35/27.,0.},{1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.}};
  
  //Fifth-order Runge-Kuta method
  double c5[6]={37/378.,0.,250/621.,125/594.,0.,512/1771.};
  //Embedded four-order Runge kutta method  
  double c4[6]={2825/27648.,0.,18575/48384.,13525/55296.,277/14336.,1/4.};


  double q4[6];
  double q5[6];

 do
    {
     
  for(i=0; i<6; i++)
    {
      q4[i]=q[i];
      q5[i]=q[i];
    }

       for(i=0; i<6; i++)
	{
	  qm[0][i]=q[i];
	}
   
      for(j=1 ; j<=6; j++)
	{
	  diff(qm[j-1],t+dt*a[j-1],qp[j-1]);
      
	  for(i=0; i<6; i++)
	    {
	      qm[j][i]=q[i];
	      k[j-1][i]=dt*qp[j-1][i];

	      for(p=0; p<j && j<6 ; p++)
		{
		  qm[j][i]=qm[j][i]+b[j][p]*k[p][i];
		}
	    }
	}

      for(i=0; i<6; i++)
	{
	  for(j=0; j<6; j++)
	    {
	      q4[i]=q4[i]+c4[j]*k[j][i];
	    }
	}

      for(i=0; i<6; i++)
	{
	  for(j=0; j<6; j++)
	    {
	      q5[i]=q5[i]+c5[j]*k[j][i];
	    }
	}
   
      error=fabs(pow(q5[0]*q5[0]+q5[1]*q5[1]+q5[2]*q5[2],1./2.)-pow(q4[0]*q4[0]+q4[1]*q4[1]+q4[2]*q4[2],1./2.))/pow(q5[0]*q5[0]+q5[1]*q5[1]+q5[2]*q5[2],1./2.);
   
    
  if(error>desiredErrorMax)
	dt=dt*pow(desiredErrorMax/error,0.2);
    
  if(error<desiredErrorMin && error!=0.)
	dt=dt*pow(desiredErrorMin/error,0.2);

     }

  while(error>desiredErrorMax || error<desiredErrorMin && error!=0.);
     

  for(i=0; i<6; i++)
    {
      q[i]=q5[i];
    }

  t=t+dt;

}


//We compute the asymptoticVelocity
double asymptoticVelocityOld(double* q)
{

  double Vsq=q[3]*q[3]+q[4]*q[4]+q[5]*q[5];
  double dist=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter*softParameter);
  double asymptVsq=Vsq-2.*Z/dist;

  if(asymptVsq<0)
    return 0.;
else
  return sqrt(asymptVsq);

}


//We compute the asymptoticVelocity according the article
double asymptoticVelocity(double* q)
{
  double vectPot=fieldAmpl/pulsation*sin(pulsation*t+phase);
  double Vsq=(q[3]+vectPot)*(q[3]+vectPot)+q[4]*q[4]+q[5]*q[5];
  double dist=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter*softParameter);
  double asymptVsq=Vsq-2.*Z/dist;

  if(asymptVsq<0)
    return 0.;
else
  return sqrt(asymptVsq);

}


//We set the initial perpendicular velocity
void setInitialVPerp()
{

  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldBirth/sqrt(2.*IP));

  //Perpendicular velocity after tunneling
  vPerp=2.*double(iVPerp)/double(nVPerp)*sigma_V;
  vPerp=2.*vPerp/4.;

  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
  weight=4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerp)/fieldBirth*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/fieldBirth);

}

//We set the initial ionization time
void setTBirth()
{

  //We set the initial time tBirth with a view to scanning the phase values from -pi/4 to pi/4
  tBirth=-M_PI/4./pulsation+double(iField)/double(nField)*M_PI/2./pulsation;

}

//We set the initial field value
void setFieldBirth()
{

  //We set the initial field value
  fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
  fieldBirth=fabs(fieldBirth);

  /*if(fieldBirth<1E-7)
    fieldBirth=1E-7;*/

}

//We set the electron position after tunneling
void setRhoBirth()
{
 
  double rhoBirthSq=IP*IP-4*Z*fabs(fieldBirth);
  rhoBirth=sqrt(IP*IP-4*Z*fabs(fieldBirth));
  rhoBirth=(rhoBirth+IP)/(fieldBirth)/2.; 

}

//We set the electron position after tunneling according the article
void setRhoBirthArt()
{
  /*
 double q,p,s,theta,eta,e;
 e=fieldBirth;
  
 p=-1./3./e/e+2./e;
 q=-2./9./e/e/e+2./3/e/e+1/e;
 s=sqrt(-pow(p/3.,3));
 theta=1/3.*acos(-q/2./s);
 eta=1/3./e+2.*pow(s,-3./2.)*cos(theta);

 rhoBirth=fabs(eta/2.);
 cout<<" e= "<<e<<" p= "<<p<<" q= "<<q<<" s= "<<s<<" theta= "<<theta<<" eta= "<<eta<<" rhoBirth= "<<rhoBirth<<endl;
  */


  // we find the root of the function which gives us the position of the electron after tunneling
 double x=0.5;
 for(int k=0; k<=50; k++)                              
	{
	  x=x-(-1./4./x-1./8./x/x-1./8.*fabs(fieldBirth)*x+1./8.)/(+1./4./x/x+1./4./x/x/x-fabs(fieldBirth)/8.);
	}

 rhoBirth=-x/2.*fabs(fieldBirth)/fieldBirth;


}

//We build the function which sets the initial conditions
void IC(double* q) 
{

  //We set the initial time tBirth
  setTBirth();
  t=tBirth; 

  //We set the initial field value
  setFieldBirth();
 
  //We set the electron position after tunneling
  setRhoBirthArt();

  double rhoBirthSq=IP*IP-4*Z*fabs(fieldBirth);
  rhoBirthOld=sqrt(IP*IP-4*Z*fabs(fieldBirth));
  rhoBirthOld=(rhoBirth+IP)/(fieldBirth)/2.; 

  //We set vPerp and the weight associated for the two directions perpendicular to the field       
  setInitialVPerp();
 
  //We set the IC
  q[0]=rhoBirth;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
  q[4]=vPerp;
  q[5]=0.;


}

//We set the number of point for the trajectory
void setTrajParameters()
{
  //We make the initial loop control variable nBirth to correspond with the ionization time tBirth
  nBirth=int(tBirth/dt); 
}

//We set the histogram intervals width
void setBinsWidth(int binsNumber=1000)
{

  double ponderomotiveEnergy=fieldAmpl*fieldAmpl/4./pulsation/pulsation;
  binsWidth=10./double(binsNumber);

}

//We store asymptotic velocities in containers of map type
void storeDataBinning()
{

  setBinsWidth();

  int range;
  
      //We compute the x value of the new point in the histogram
      range=int(asymptoticVelocity(q)/binsWidth);
      
      if(range==0) return;
      if(asymptoticVelocity(q)>10.) return;

      //We declare a container of map type
      //It can contains pairs with are couple of objects
      //First the range, second the asymptic velocity value

      map<int,double>::iterator findRange= asymptVel.find(range);
   
      if(findRange==asymptVel.end())
	{
	  asymptVel[range]=weight;
	}
      else
	{
	  asymptVel[range]=asymptVel[range]+weight;
	}
}

//Finally we write all the data binning in a file
void writeDataBinning()
{
  
  map<int,double>::iterator it= asymptVel.begin();

  for(1; it!=asymptVel.end(); it++)
    {
      dat<<(it->first)*binsWidth<<" "<<log10(it->second)<<endl;
    }

}

       
//We implement a load bar with a view to displaying the remaining time

//First top
time_t start = time (NULL);  

static inline void loadbar(int i, int np)
{

  //If it the first time, we execute it we leave three line break
  static int k=0;
  if(k==0)
    {
  cout<<" "<<endl;
  cout<<" "<<endl;
  cout<<" "<<endl;
    }
  k++;

  //Load bar
  double ratio =(double)(i)/(double)(np);
  double p=ratio*(100); 
 
  int s=cout.precision();

  cout<<"\r"<<"progression= "<<setprecision(3)<<p<<"%"<<setprecision(s)<<"                "<<endl;
 
  //second top
  time_t end = time (NULL);         

  //Elapsed time since the beginning of the simulation
  double T0=difftime(end,start);

  int T1=(int)(T0/3600);
  int T2=(int)(T0/60)-T1*60;
  int T3=(int)(T0)-T1*3600-T2*60;

  cout<<"elapsed time= "<<T1<<" h "<<T2<<" min "<<T3<<" s "<<"             "<<endl;


  //Remaining time since the beginning of the simulation
  double Z0=(100-p)*T0/p;
  int Z1=(int)(Z0/3600);
  int Z2=(int)(Z0/60)-Z1*60;
  int Z3=(int)(Z0)-Z1*3600-Z2*60; 
  
  if(fabs(Z1)>1000 || fabs(Z2)>60 || fabs(Z3)>60)
    {
      Z1=0;
      Z2=0;
      Z3=0;
    }

  cout<<"\r"<<"estimated remaining time= "<<Z1<<" h "<<Z2<<" min "<<Z3<<" s "<<"    "<<endl;

  cout<<"rhoBirthArt= "<<rhoBirth<<" rhoBirthOld= "<<rhoBirthOld<<" diff= "<<rhoBirth-rhoBirthOld<<"      "<<endl;
  cout<<"ponderomotriveMin= "<<fieldBirth*fieldBirth/4./pulsation/pulsation<<"      "<<endl;
  cout<<"ponderomotriveMax= "<<fieldAmpl*fieldAmpl/4./pulsation/pulsation<<"       "<<endl;      
  cout<<"asymptoticVelocityArt= "<<asymptoticVelocity(q)<<" asymptoticVelocityOld= "<<asymptoticVelocityOld(q)<<" diff= "<<asymptoticVelocity(q)-asymptoticVelocityOld(q)<<"          "<<endl;
  cout<<"weight= "<<weight<<"        "<<endl;
  cout<<"phaseBirth= "<<pulsation*tBirth*180./M_PI<<"    "<<endl;
  cout<<"vPerp= "<<vPerp<<endl;

  //We ask the cursor to go three lines up
  if(int(i)!=int(np))
    cout<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F";

}


//FUNCTION MAIN

int main()
{

  //We perform three loops
  //first, for each ionization time
  //second, for each perpendicular velocity along yPrime
  //third, for each perpendicular velocity along zPrime

  for(iField=0; iField<=nField; iField++)
    {
      for(iVPerp=0; iVPerp<=nVPerp; iVPerp++)
	{
	  
	  //We call the function which will fix the initial condition
	  IC(q);

	  //We set the number of points for the trajectory
	  setTrajParameters();
	   
	  //We compute the trajectory
	  //No final loop control variable: the loop stops when the electron is completely ionized
	  //We declare a boolean variable	  
	  bool stop=true;
	  for(int n=nBirth; stop  ; n++)
	    { 
		  
	      //We call the function which solve eq of the motion
	      rk5(coulomb_field,q,t,dt,error,desiredErrorMin,desiredErrorMax);

	      //if(nBirth%100==0) {cout<<"+1";}

	      //We stop when the electron is fully ionized
	      if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])>40.)
		stop=false;

	      if(n==nBirth+50000000) stop=false;
	    }

	  //We store the asymptotic velocity in a container of map type with a view to make a data binning
	  storeDataBinning();
	      
	  //We update the load bar
	  loadbar(iVPerp+(iField-1)*nVPerp, nField*nVPerp);

	}
    }

  //Finally we write the data binning in a file
  writeDataBinning();

  //We display the spectra with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"set xtics rotate out"<<endl;
  gnu<<"set multiplot  layout 2, 1"<<endl;
  gnu<<"set key on outside left bmargin box title 'nField="<<nField<<", nVPerp="<<nVPerp<<", dt="<<dt<<"'"<<endl;
  gnu<<"plot 'data.dat' using 1:2 w l title 'Spectra'"<<endl;
  gnu<<"unset multiplot"<<endl;
  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
  gnu<<"set output 'spectra.eps'"<<endl;
 
  gnu<<"replot"<<endl;
  
  gnu.close();
  system("gnuplot data.gnu");


  return 0;
}

