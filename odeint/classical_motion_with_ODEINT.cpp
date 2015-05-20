#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<iomanip>
#include<map>
#include <vector>

using namespace std;

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

//Plateau in above-threashold-ionization spectra and chaotic behavior in rescattering processes
//Bambi Hu, Jie Liu, Shi-gang Chen

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

typedef vector<double> state_type;

//We consider the rotating frame of reference (xPrime, yPrime, zPrime) where the electric field holds the xPrime direction
double vPerp, weight;
double sigmaV;
int iField, iVPerp=0;
int nField=1, nVPerp=1;

//We create the array in which we will store the orbit
state_type x(6);

//We create the time variable
double t;

//We fix the step
double dt=0.000001;

//We declare a variable for the initial time
double tBirth;

//We declare a variable for the initial field value
double fieldBirth;

//We declare a variable for the electron position after tunneling
double rhoBirth;

//We declare containers which will contains electron spectra
map<int,double> asymptVel;

//Bins width
double binsWidth;
 
//We open a file with a view to writing in it
fstream dat("data.dat",ios::out);

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//Ordinary differential equations of the dynamic of an electron in coulomb potential 
class coulomb
{
  double K;
  double field[3];
  string whichField;

public:

  coulomb(string whichField="withoutField") : whichField(whichField)
  { 
    K=Z;
    double field[3]={0.,0.,0.};
  }

    void operator() ( const state_type &x , state_type &dxdt , const double  t  )
    {      

    //If it has been requested, we switch on the electric field
 if(whichField=="withField")
      {
	field[0]=0.;
	field[1]=0.;
	field[2]=-fieldAmpl*cos(pulsation*t+phase);
      }

    double a=pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],3./2.); 
    dxdt[0]=x[3];
    dxdt[1]=x[4];
    dxdt[2]=x[5];
    dxdt[3]=-K*x[0]/a+field[0];
    dxdt[4]=-K*x[1]/a+field[1];
    dxdt[5]=-K*x[2]/a+field[2];
    }
};

//integrate_observer
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

//We compute the asymptoticVelocity according the article
double asymptoticEnergy(state_type x)
{
  double vectPot=fieldAmpl/pulsation*sin(pulsation*t+phase);
  double Vsq=(x[5]+vectPot)*(x[5]+vectPot)+x[4]*x[4]+x[3]*x[3];
  double dist=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  double asymptVsq=Vsq-2.*Z/dist;

  if(asymptVsq<0)
    return 0.;
else
  return asymptVsq/2.;

}

//We set the electron position after tunneling according the article
void setRhoBirth()
{

  //We set the initial field value
  double  fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
  fieldBirth=fabs(fieldBirth);

  // we find the root of the function which gives us the position of the electron after tunneling
  double x=100.;

  for(int k=0; k<50; k++)                       
    {
       x=x-(-1./4./x-1./8./x/x-1./8.*fabs(fieldBirth)*x+1./8.)/(+1./4./x/x+1./4./x/x/x-fabs(fieldBirth)/8.);
    }
 
  rhoBirth=-x/2.*fabs(fieldBirth)/fieldBirth;

}


//We set the initial perpendicular velocity
void setVPerpBirth()
{
  //REF: J. Liu, Classical Trajectory Perspective of Atomic Ionization in Strong Laser Fields

  //Width of the velocity distributions after tunneling
  sigmaV=sqrt(fieldBirth/sqrt(2.*IP));

  //Perpendicular velocity after tunneling
  vPerp=2.*double(iVPerp)/double(nVPerp)*sigmaV;
  vPerp=2.*vPerp/4.;
  
  //Ionization rate with a given field and a given transverse velocity according ADK distribution
  weight=fabs(4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerp)/fieldBirth*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/fieldBirth));

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

}

  //We build the function which sets the initial conditions
void IC(state_type& x)
{

  //We set the initial time tBirth
  setTBirth();

  tBirth=0.;

  t=tBirth; 

  //We set the initial field value
  setFieldBirth();

  fieldBirth=0.06;

  //Position of the electron after tunneling
  setRhoBirth(); 

  //We set vPerp and the weight associated for the two directions perpendicular to the field       
  setVPerpBirth();

  //We set the IC

   vPerp=6E-5;
   //0.000304043 0 -0.000524244 -14.5082 0 55.61  

   rhoBirth=0.1;
   vPerp=0.123;
   
      x[0]=0.;
      x[1]=0.;
      x[2]=rhoBirth;
      x[3]=vPerp;
      x[4]=0.;
      x[5]=0.;
      /*
    x[0]=0.000304043;
   x[1]=0.;
   x[2]=-0.000524244;
   x[3]=-14.5082;
   x[4]=0.;
   x[5]=55.61;
      */

}


//We set the histogram intervals width
void setBinsWidth(int binsNumber=100)
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
  range=int(asymptoticEnergy(x)/binsWidth);
  
  if(range<=0) return;
  //  if(asymptoticVelocity(q)>10.) return;

  //We declare a container of map type
      //It can contains pairs with are couple of objects
      //First the range, second the asymptic velocity value

  dat<<"asymptoticEnergy= "<<asymptoticEnergy(x)<<" range= "<<range<<endl;

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

  //We normalize the distribution law
  double sum=0.;
  for(1; it!=asymptVel.end(); it++)
    {
          sum=sum+it->second;
	  dat<<"sum= "<<sum<<endl;
    }

  for(it=asymptVel.begin(); it!=asymptVel.end(); it++)
    {
      dat<<(it->first)*binsWidth<<" "<<log10( (it->second)/sum )<<endl;
    }

}


//Keplerian energy
double energy(state_type x)
{
  return 0.5*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])-Z/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],1./2.);
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

  cout<<"rhoBirth="<<"\t"<<rhoBirth<<endl;
  cout<<"ponderomotriveMin="<<"\t"<<fieldBirth*fieldBirth/4./pulsation/pulsation<<"      "<<endl;
  cout<<"ponderomotriveMax="<<"\t"<<fieldAmpl*fieldAmpl/4./pulsation/pulsation<<"       "<<endl;      cout<<"asymptoticEnergy="<<"\t"<<asymptoticEnergy(x)<<endl;
 cout<<"KeplerianEnergy="<<"\t"<<energy(x)<<endl;
  cout<<"weight="<<"\t"<<weight<<"        "<<endl;
  cout<<"phaseBirth="<<"\t"<<pulsation*tBirth*180./M_PI<<"    "<<endl;
  cout<<"vPerp= "<<vPerp<<endl;
  cout<<"sigmaV= "<<sigmaV<<endl;
  cout<<"vPerp="<<"\t"<<vPerp<<endl;
  cout<<"weight="<<"\t"<<weight<<endl;

  //We ask the cursor to go three lines up
  if(int(i)!=int(np))
    cout<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F";

}

//FUNCTION MAIN


int main()
{

  //ODEINT objects

  // define_adapt_stepper
  // typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
  // typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
   typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

  // integrate_adapt
  typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
  controlled_stepper_type controlled_stepper(1.0E-6);

  coulomb m_coulomb("withoutField");
   
  vector<state_type> x_vec;
  vector<double> times;
  size_t steps;
   

  //We perform three loops
  //first, for each ionization time
  //second, for each perpendicular velocity along yPrime
  //third, for each perpendicular velocity along zPrime

  for(iField=1; iField<=nField; iField++)
    {
      for(iVPerp=1; iVPerp<=nVPerp; iVPerp++)
	{
	  
	   //We call a function which will fix the IC
	  IC(x);

	  //We compute the trajectory

	  //We call the rk4 function which solve eq of the motion
	  steps=integrate_adaptive( controlled_stepper , m_coulomb , x , 0.0 , 500.0 , dt,push_back_state_and_time( x_vec , times ));

	  // steps = integrate(m_coulomb ,  x , 0.0 , 600.0 , 0.0001, push_back_state_and_time( x_vec , times )  );

  	  	 
	   //We write the position of the electron in the phase space in a file
	  for(int j=0; j<=steps; j++)
	  {
	  for (int i=0; i<6; i++)
	  {
	    dat<<x_vec[j][i]<<"\t";
	  }
	  /*	  dat<<times[j]<<"\t"<<-fieldAmpl*cos(pulsation*times[j]+phase)<<"\t"<<energy(x_vec[j])<<"\t"<<asymptoticEnergy(x_vec[j])<<"\t"<<fabs(energy(x_vec[j])-asymptoticEnergy(x_vec[j]));*/
	  dat<<" "<<endl;
	  }
	  
	  //We store the asymptotic velocity in a container of map type with a view to make a data binning

	  storeDataBinning();
	  

	  //We update the load bar
	  loadbar(iVPerp+(iField-1)*nVPerp, nField*nVPerp);
	  
	  
	}
    }




  //Finally we write the data binning in a file
  dat<<" "<<endl;
  dat<<" "<<endl;
  writeDataBinning();
 
    
 //We display the spectra with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"set xtics rotate out"<<endl;
  gnu<<"set multiplot  layout 1, 1"<<endl;
  gnu<<"set bmargin  1.5"<<endl;
  gnu<<"set tmargin  0.5"<<endl;
  gnu<<"set lmargin  8"<<endl;
  gnu<<"set rmargin  2"<<endl;
  gnu<<"set key on outside left bmargin box title 'nField="<<nField<<", nVPerp="<<nVPerp<<", dt="<<dt<<"'"<<endl;
  gnu<<"set key off"<<endl;
  gnu<<"plot 'data.dat' index 0 using 1:3 w l title 'Spectra'"<<endl;
  //  gnu<<"plot 'data.dat' index 0 using 7:8 w l title 'Spectra'"<<endl;
  // gnu<<"plot 'data.dat' index 0 using 7:9 w l title 'Spectra'"<<endl;
  gnu<<"unset multiplot"<<endl;
  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
  gnu<<"set output 'spectra.eps'"<<endl;
 
  gnu<<"replot"<<endl;
  
  gnu.close();
  system("gnuplot data.gnu");

  return 0;
}

