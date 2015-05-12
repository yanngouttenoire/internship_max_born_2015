#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>

using namespace std;

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

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare a variable for the electron position after tunneling
double rhoBirth;

//We consider the rotating frame of reference (xPrime, yPrime, zPrime) where the electric field holds the xPrime direction
double vPerpYPrime, vPerpZPrime, weightYPrime, weightZPrime;
int iField, iVyPrime=0, iVzPrime=0;
int nField=100, nVyPrime=10, nVzPrime=10;

//We create the time variable
double t;

//We create the array in which we will store the orbit
double* q=new double[6];

//We fix the step
double dt=1.; 

//We declare the number of points for trajectories
int nBirth;
int nEnd;

//We declare containers which will contains electron spectra
map<int,double> asymptVelZUp;
map<int,double> asymptVelZDown;

//Bins width
double binsWidth;

//We open files with a view to writing in them
fstream dat("data.dat",ios::out);


//FUNCTIONS DECLARATION


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

//We set the initial perpendicular velocity
void setInitialVPerp(double &vPerp, double &weight, int iVPerp, int nVPerp)
{

  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldBirth/sqrt(2.*IP));

  //Perpendicular velocity after tunneling
  vPerp=-2.*sigma_V+4.*double(iVPerp)/double(nVPerp)*sigma_V;
  vPerp=2.*vPerp/4.;

  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields
  weight=4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerp)/fieldBirth*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/fieldBirth);
}

//We set the initial ionization time
void setTBirth(int iField, int nField)
{
 
  //We set the initial time tBirth
  tBirth=double(iField)/double(nField)*2.*pulseDuration; 

}

//We set the initial field value
void setFieldBirth()
{
 
  //We set the initial field value
  fieldBirth=fieldAmpl*pow(sin(M_PI*tBirth/2./pulseDuration),2)*sqrt(pow(cos(ellipticity)*sin(pulsation*tBirth+phase),2)+pow(sin(ellipticity)*cos(pulsation*tBirth+phase),2));
  fieldBirth=fabs(fieldBirth);

  if(fieldBirth==0)
    fieldBirth=1E-6;
}

//We set the electron position after tunneling
void setRhoBirth()
{
 
  double rhoBirthSq=IP*IP-4*Z*fabs(fieldBirth);
  rhoBirth=sqrt(IP*IP-4*Z*fabs(fieldBirth));
  rhoBirth=(rhoBirth+IP)/(fieldBirth)/2.; 

}

//We build the function which sets the initial conditions
void IC(double* q) 
{
  //We set the initial time tBirth
  setTBirth(iField, nField);
  t=tBirth; 

  //We set the initial field value
  setFieldBirth();

  //We set the electron position after tunneling
  setRhoBirth();

  //We set vPerp and the weight associated for the two directions perpendicular to the field       
  setInitialVPerp(vPerpYPrime, weightYPrime, iVyPrime, nVyPrime);
  setInitialVPerp(vPerpZPrime, weightZPrime, iVzPrime, nVzPrime);

  //We set the IC
  q[0]=rhoBirth*cos(ellipticity)*sin(pulsation*tBirth+phase);
  q[1]=-rhoBirth*sin(ellipticity)*cos(pulsation*tBirth+phase);
  q[2]=0.;
  q[3]=vPerpYPrime*sin(ellipticity)*cos(pulsation*tBirth+phase);
  q[4]=vPerpYPrime*cos(ellipticity)*sin(pulsation*tBirth+phase);
  q[5]=vPerpZPrime;

}

//We set the number of point for the trajectory
void setTrajParameters()
{

  //The trajectory will stop at the end of the impulsion
  nEnd=int(2.*pulseDuration/dt); 

  //We compute the correponding step dt
  dt=2.*pulseDuration/double(nEnd);

  //We make the initial loop control variable nBirth to correspond with the ionization time tBirth
  nBirth=int(tBirth/dt); 

}

//We compute the asymptoticVelocity
double asymptoticVelocity(double* q)
{

  double Vsq=q[3]*q[3]+q[4]*q[4]+q[5]*q[5];
  double dist=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
  double asymptVsq=Vsq-2.*Z/dist;

  if(asymptVsq<0)
    return 0.;
else
  return sqrt(asymptVsq);

}

//We set the histogram intervals width
void setBinsWidth(int binsNumber=10000)
{

  double ponderomotiveEnergy=fieldAmpl*fieldAmpl/4./pulsation/pulsation;
  binsWidth=20./double(binsNumber);

}

//We store asymptotic velocities in containers of map type
void storeDataBinning()
{

  setBinsWidth();

  double weight=weightYPrime*weightZPrime;
  double range;
  
  if(vPerpZPrime>=0.)
    {

      //We compute the x value of the new point in the histogram
      range=int(asymptoticVelocity(q)/binsWidth);

      //We declare a container of map type
      //It can contains pairs with are couple of objects
      //First the range, second the asymptic velocity value
      map<int,double>::iterator findRange= asymptVelZUp.find(range);
   
      if(findRange==asymptVelZUp.end())
	{
	  asymptVelZUp[range]=weight;
	}
      else
	{
	  asymptVelZUp[range]=asymptVelZUp[range]+weight;
	}

    }

  if(vPerpZPrime<0.)
    {
     
      range=int(asymptoticVelocity(q)/binsWidth);
      map<int,double>::iterator findRange= asymptVelZDown.find(range);

      if(findRange==asymptVelZDown.end())
	{
	  asymptVelZDown[range]=weight;
	}
      else
	{
	  asymptVelZDown[range]=asymptVelZDown[range]+weight;
	}

    }
}

//Finally we write all the data binning in a file
void writeDataBinning()
{
  
  map<int,double>::iterator itUp= asymptVelZUp.begin();
  map<int,double>::iterator itDown= asymptVelZDown.begin();

  for(1; itUp!=asymptVelZUp.end(); itUp++)
    {
      dat<<(itUp->first)*binsWidth<<" "<<itUp->second<<endl;
    }

  dat<<" "<<endl;
  dat<<" "<<endl;

  for(1; itDown!=asymptVelZDown.end(); itDown++)
    {
      dat<<(itDown->first)*binsWidth<<" "<<itDown->second<<endl;
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

  cout<<"ponderomotriveMin= "<<fieldBirth*fieldBirth/4./pulsation/pulsation<<"      "<<endl;
  cout<<"ponderomotriveMax= "<<fieldAmpl*fieldAmpl/4./pulsation/pulsation<<"       "<<endl;      cout<<"i= "<<i<<" np= "<<np<<endl;
  cout<<"asymptoticVelocity= "<<asymptoticVelocity(q)<<"           "<<endl;

  //We ask the cursor to go three lines up
  if(int(i)!=int(np))
    cout<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F";
}


//FUNCTION MAIN

int main()
{

  //We perform three loops
  //first, for each ionization time
  //second, for each perpendicular velocity along yPrime
  //third, for each perpendicular velocity along zPrime

  for(iField=1; iField<=nField; iField++)
    {
      for(iVyPrime=1; iVyPrime<=nVyPrime; iVyPrime++)
	{
	  for(iVzPrime=1; iVzPrime<=nVzPrime; iVzPrime++)
	    {
	  
	      //We call the function which will fix the initial condition
	      IC(q);

	      //We set the number of points for the trajectory
	      setTrajParameters();
	   
	      //We compute the trajectory
	      for(int n=nBirth; n<=nEnd; n++)
		{ 
		  
		  //We call the function which solve eq of the motion
		  rk4(coulomb_field,q,t,dt);

		}

	      //We store the asymptotic velocity in a container of map type with a view to make a data binning
	      storeDataBinning();
	      
	      //We update the load bar
	      loadbar(iVzPrime+nVzPrime*((iVyPrime-1)+(iField-1)*nVyPrime), nField*nVyPrime*nVzPrime);

	    }
	}
    }

  //Finally we write the data binning in a file
  writeDataBinning();

  //We display the spectra with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);


  gnu<<"set xtics rotate out"<<endl;

  gnu<<"set multiplot  layout 2, 1"<<endl;
  gnu<<"set key on outside left bmargin box title 'nField="<<nField<<", nVyPrime="<<nVyPrime<<", nVzPrime="<<nVzPrime<<"'"<<endl;

  gnu<<"plot 'data.dat' index 0 using 1:2 w l title 'Spectra Z up'"<<endl;
  gnu<<"plot 'data.dat' index 1 using 1:2 w l title 'Spectra Z down'"<<endl;

  gnu<<"unset multiplot"<<endl;

  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
  gnu<<"set output 'spectra.eps'"<<endl;
  gnu<<"replot"<<endl;
  gnu.close();

  system("gnuplot data.gnu");

  return 0;
}

