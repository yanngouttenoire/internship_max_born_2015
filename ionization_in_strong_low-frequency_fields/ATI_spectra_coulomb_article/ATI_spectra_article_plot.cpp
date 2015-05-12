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
double fieldAmpl=sqrt(12*(1E13)/uaIntensity);

double pulsation=2.*M_PI*lightSpeed/waveLenght*uaTime;
double IP=13.605804/uaEnergy;
double opticalCycle=2.*M_PI/pulsation;
double pulseDuration=3.*opticalCycle;
double phase=0;
double ellipticity=0.1;
double Z=1.; //Positive charge of the nucleus

//We declare a softening parameter for the soft coulomb potential
double softParameter=0.001;

//We consider the rotating frame of reference (xPrime, yPrime, zPrime) where the electric field holds the xPrime direction
double vPerp, weight;
int iField, iVPerp=0;
int nField=0, nVPerp=0;
double tBirth0=0.;
double vPerp0=0.123;
double rhoBirth0=11.;

//We fix the step
double dt=0.00001; 

//We declare the number of points for trajectories
int nBirth;
int nEnd;

//We create the time variable
double t;

//We create the array in which we will store the orbit
double* q=new double[6];

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare a variable for the electron position after tunneling
double rhoBirth;
double rhoBirthOld;

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
 
  double a=pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter,3./2.);
  qp[0]=q[3];
  qp[1]=q[4];
  qp[2]=q[5];
  qp[3]=-K*q[0]/a;
  qp[4]=-K*q[1]/a;
  qp[5]=-K*q[2]/a-fieldAmpl*cos(pulsation*t+phase);

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
    }
  t=t+dt/2.; //Out of the loop !

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
    }
  t=t+dt/2.;

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


//We compute the asymptoticVelocity
double asymptoticVelocityOld(double* q)
{

  double Vsq=q[3]*q[3]+q[4]*q[4]+q[5]*q[5];
  double dist=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
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
  double dist=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter);
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

  vPerp=vPerp0;

  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
  weight=4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerp)/fieldBirth/M_PI*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/fieldBirth);

}

//We set the initial ionization time
void setTBirth()
{

  //We set the initial time tBirth with a view to scanning the phase values from -pi/4 to pi/4
  // tBirth=-M_PI/4./pulsation+double(iField)/double(nField)*M_PI/2./pulsation;
  tBirth=tBirth0;
}

//We set the initial field value
void setFieldBirth()
{

  //We set the initial field value
  fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
  fieldBirth=fabs(fieldBirth);

  cout<<"fieldBirth= "<<fieldBirth<<endl;
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
  //Not working, there is a factor 2
 double q,p,s,theta,eta,e;
 e=fieldBirth;
  
 p=-1./3./e/e+2./e;
 q=-2./9./e/e/e+2./3./e/e+1./e;
 s=sqrt(-pow(p/3.,3));
 theta=1/3.*acosh(-q/2./s);
 eta=1/3./e+2.*pow(s,-3./2.)*cosh(theta);

 rhoBirth=fabs(eta/2.);

 //Factor 2
 rhoBirth=-rhoBirth/2.;

 cout<<"rhoBirthAnalytic= "<<rhoBirth<<endl;
  
  cout<<"Newton-Raphson-AnalyticTest= "<<-1./4./eta-1./8./eta/eta-1./8.*fieldBirth*eta+1./8.<<endl;

  // we find the root of the function which gives us the position of the electron after tunneling
  double x=0.5;
 for(int k=0; k<=50; k++)                              
	{
	  x=x-(-1./4./x-1./8./x/x-1./8.*fieldBirth*x+1./8.)/(+1./4./x/x+1./4./x/x/x-fieldBirth/8.);
	}

 rhoBirth=-x/2.;

 cout<<"rhoBirthNewton= "<<rhoBirth<<endl;
 cout<<"Newton-Raphson-NewtonTest= "<<-1./4./x-1./8./x/x-1./8.*fieldBirth*x+1./8.<<endl;

 /*for(double x=-20.; x<20; x=x+0.1)
   {
     if(!(x>-0.5 && x<0.5)) 
       {
   dat<<x<<" "<<-1./4./x-1./8./x/x-1./8.*fieldBirth*x+1./8.<<endl;
   }
   }*/
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
  q[0]=0.;
  q[1]=0.;
  q[2]=rhoBirth;
  q[3]=vPerp;
  q[4]=0.;
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
      //   dat<<(it->first)*binsWidth<<" "<<log10(it->second)<<endl;
    }

}

//Keplerian energy
double energy(double* q)
{
  return 0.5*(q[3]*q[3]+q[4]*q[4]+q[5]*q[5])-Z/pow(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+softParameter,1./2.);
}
       
//We implement a load bar with a view to displaying the remaining time

//First top
time_t start = time (NULL);  

static inline void loadbar(int i, int np)
{

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
  cout<<" "<<endl;

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
	  bool stop=true; //stop
	  cout<<" "<<endl;
	    for(int n=nBirth; n<=nBirth+5000000  ; n++)
	    { 
		  
	      //We call the function which solve eq of the motion
	      rk4(coulomb,q,t,dt);

	      //cout<<"\033[F"<<energy(q)<<endl;

	      dat<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<t<<" "<<energy(q)<<" "<<fabs(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])-softParameter)/(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])-softParameter)<<endl;
	      // dat<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;

	      //dat<<t<<" "<<-fieldAmpl*cos(pulsation*t+phase)<<endl;

	      // if(n%100==0) cout<<"asymptoticVelocityArt= "<<asymptoticVelocity(q)<<endl;
	      // dat<<t<<" "<<asymptoticVelocity(q)<<" "<<sqrt(q[3]*q[3]+q[4]*q[4]+q[5]*q[5])<<" "<<asymptoticVelocity(q)-sqrt(q[3]*q[3]+q[4]*q[4]+q[5]*q[5])<<endl;

	      //We stop when the electron is fully ionized
	      if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])>600.)
		stop=false;

	    }

	  //We store the asymptotic velocity in a container of map type with a view to make a data binning
	  storeDataBinning();
	      
	  //We update the load bar
	  loadbar(iVPerp+(iField-1)*nVPerp, nField*nVPerp);

	}
    }

  //Finally we write the data binning in a file
  writeDataBinning();


  //We plot a circle which somehow corresponds to the ionic core

  double u;
 
  dat<<" "<<endl;
 dat<<" "<<endl;

  for(int i=0; i<=100; i++)
    {
      u=double(i)*2.*M_PI/100.;
         dat<<softParameter*1000*cos(u)<<" "<<0<<" "<<softParameter*1000*sin(u)<<endl;
    }
 
  //We display the spectra with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"set xtics rotate out"<<endl;
 

  //gnu<<"set xlabel 'X'"<<endl;
 gnu<<"set zlabel 'Y'"<<endl;
 gnu<<"set zlabel 'Z'"<<endl;
 
 // gnu<<"set tmargin 0"<<endl;
 //gnu<<"set bmargin 5"<<endl;
 // gnu<<"set style data boxes"<<endl;
    gnu<<"set multiplot  layout 2, 2"<<endl;
  gnu<<"set key on outside left bmargin box title 'nField="<<nField<<", nVPerp="<<nVPerp<<", dt="<<dt<<"'"<<endl;
  //  gnu<<"plot 'data.dat' using 1:2 w l title 'Spectra'"<<endl;

 gnu<<"plot 'data.dat' using 1:3 w l title 'Spectra'"<<endl;
 //gnu<<"plot 'data.dat' using 1:2 w l title 'Spectra'"<<endl;
   gnu<<"plot 'data.dat' using 1:3 w l title 'Spectra'"<<endl;
   gnu<<"plot 'data.dat' index 0 using 4:5 w l title 'Spectra'"<<endl;
 
   gnu<<"plot 'data.dat' index 0 using 4:6 w l title 'Spectra'"<<endl;
     gnu<<"unset multiplot"<<endl;
  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,20'"<<endl;
  gnu<<"set output 'spectra.eps'"<<endl;
 
  gnu<<"replot"<<endl;
  
  gnu.close();
  system("gnuplot data.gnu");


  return 0;
}

