#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<iomanip>
#include <vector>

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;


using namespace std;
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

  //We create the array in which we will store the orbit
  state_type x(6);

  //We create the time variable
  double t;
double tBirth=-6.356;

double rhoBirth;

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;


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



//[ integrate_observer
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
//]

//We compute the asymptoticVelocity according the article
double asymptoticEnergy(state_type x)
{
  double vectPot=fieldAmpl/pulsation*sin(pulsation*t+phase);
  double Vsq=(x[3]+vectPot)*(x[3]+vectPot)+x[4]*x[4]+x[5]*x[5];
  double dist=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  double asymptVsq=Vsq-2.*Z/dist;

  if(asymptVsq<0)
    return 0.;
else
  return asymptVsq/2.;

}

//We set the electron position after tunneling according the article
void setRhoBirthArt()
{


  //We set the initial field value
 double  fieldBirth=fieldAmpl*cos(pulsation*tBirth+phase);
  fieldBirth=fabs(fieldBirth);

   // we find the root of the function which gives us the position of the electron after tunneling
 double x=100.;
 for(int k=0; k<=50; k++)                              
	{
	  x=x-(-1./4./x-1./8./x/x-1./8.*fabs(fieldBirth)*x+1./8.)/(+1./4./x/x+1./4./x/x/x-fabs(fieldBirth)/8.);
	}

 rhoBirth=-x/2.*fabs(fieldBirth)/fieldBirth;

}

//We build the function which sets the initial conditions
void IC(state_type& x)
{

  //Position of the electron after tunneling
  setRhoBirthArt();      

  //Width of the velocity distributions after tunneling
  double sigma_V=sqrt(fieldAmpl*cos(phase)/sqrt(2.*IP));  

  //Perpendicular velocity after tunneling
  double vPerp=4.*sigma_V*(2.*(double(rand())/double(RAND_MAX))-1.);

  //We display some informations
  cout<<"sigma_V= "<<sigma_V<<endl;
  cout<<"vPerp= "<<vPerp<<endl;
  cout<<" "<<endl;
  cout<<" "<<endl;

  //We set the IC
  x[0]=0.;
  x[1]=0.;
  x[2]=rhoBirth;
  x[3]=0.130;
  x[4]=0.;
  x[5]=0.;

  t=tBirth;
}

//Keplerian energy
double energy(double* x)
{
  return 0.5*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])-Z/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],1./2.);
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

   //We fix the step
  double dt=0.0001;

   //We call a function which will fix the IC
  IC(x);
  
  time_t start=time(NULL);

   // define_adapt_stepper
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
  
    // integrate_adapt
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper;

    coulomb m_coulomb("withoutField");
   
    vector<state_type> x_vec;
    vector<double> times;
    size_t steps;
   
    dt=0.01;

	//We call the rk4 function which solve eq of the motion
  steps=integrate_adaptive( controlled_stepper , m_coulomb , x , 0.0 , 500.0 , dt,  push_back_state_and_time( x_vec , times ) );

    // steps = integrate(m_coulomb ,  x , 0.0 , 600.0 , 0.0001, push_back_state_and_time( x_vec , times )  );

    

    
      /*  for( size_t i=0; i<=steps; i++ )
      {
        cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] <<endl;
      }*/

       //We write the position of the electron in the phase space
    for(int j=0; j<=steps; j++)
      {
	for (int i=0; i<6; i++)
	  {
	    dat<<x_vec[j][i]<<" ";
	  }
	dat<<" "<<endl;
	}

    cout<<"asymptoticEnergy= "<<asymptoticEnergy(x)<<endl;

      //We display some informations
      /* if(n%100==0)
	{
	  cout<<"\033[F"<<"\033[F"<<"\033[F"<<"\033[F";
	  cout<<"Energy= "<<energy(q)<<endl;
	  cout<<"Distance= "<<pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],1./2.)<<endl;
	  time_t end=time(NULL);
	  //We display the elapsed time and the load bar
	  cout<<"Elapsed time= "<<difftime(end,start)<<endl;
	  cout<<"\r"<<"progression= "<<setprecision(3)<<(double)(n)/(double)(N)*100.<<"%"<<setprecision(cout.precision())<<"                "<<endl;
	}
      */


  //We display the duration of the simulation in ua
  cout<<"t= "<<t<<endl;
 
    
  //We plot the orbit (and the "ionic circle") with gnuplot

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"plot 'data.dat' using 3:1 with line title 'Coulomb Orbit'"<<endl;

  gnu<<"pause -1"<<endl;
  gnu<<"set term postscript"<<endl;
  gnu<<"set output 'courbe.ps'"<<endl;
  gnu<<"replot"<<endl;
  gnu.close();

  system("gnuplot data.gnu");

  return 0;
}

