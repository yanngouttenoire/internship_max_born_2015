#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

#include"system.h"
#include"solve.h"
#include"spectra.h"
#include"ic.h"
#include"display.h"
#include"plot.h"
#include"parameters.h"
#include"observer.h"


using namespace std;

//Plateau in above-threashold-ionization spectra and chaotic behavior in rescattering processes
//Bambi Hu, Jie Liu, Shi-gang Chen

//VARIABLES DECLARATION


//We consider the rotating frame of reference (xPrime, yPrime, zPrime) where the electric field holds the xPrime direction
int iFieldBirth, iVPerpBirth=0;

//We declare a variable for the step in controlledRK5
double dt=0.001; 

//We declare runge kutta error, its max, its min and the desired error (desired accuracy) 
double error;
double errorMax=1E-25;

//We declare a variable which contains the current nbr of trajectories
int nTraj;

//We declare the duration during which we want to let the electron propagates 
double trajDuration=500.;

//We declare boolean control 
bool stopStepper;
bool isdtMin;
double distMin;

//We declare the time variable
double t;

//We declare the array in which we will store the orbit
typedef vector<double> state_type;
state_type x(6);

//Bins width
double binsWidth;

//We open files with a view to writing in them 

fstream dataFile("data.dat",ios::out);



//FUNCTION MAIN

int main()
{


//We instantiate the differente object implemented in the others .cpp files
//Each object will execute a specific task
IC myIC;
System<state_type> mySystem("withField");	
Solve<state_type> mySolve;
Display myDisplay;
Spectra<state_type> mySpectra(0.1);
Plot myPlot;

  //We perform two loops
  //first, for each ionization time (initial field value)
  //second, for each perpendicular velocity

  for(iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(iVPerpBirth=1; iVPerpBirth<=nVPerpBirth; iVPerpBirth++)
	{
	  
	  //INITIAL CONDITIONS
          //We set the ionization time
	  myIC.setTBirth(iFieldBirth, nFieldBirth);
	  myIC.setFieldBirth();
	  myIC.setVPerpBirth(iVPerpBirth, nVPerpBirth);
	  myIC.setRhoBirth();
	  myIC.setIC<state_type>(x,t);


	  //We compute the trajectory

	
	  //We initialise the boolean controls	  
	  stopStepper=true; 
	  isdtMin=false;

	  distMin=100.;


	    for(nTraj=0; stopStepper; nTraj++)
	    { 
		  
	      //We call the function which solve eq of the motion
	      mySolve.controlledRK5(mySystem,x,t,dt,error,desiredErrorMin,desiredErrorMax);


	      //The trajectory stops when the electron has sufficiently propagated 
	      if((t-myIC.tBirth)>trajDuration)
		stopStepper=false;

              //We look for the best approach (not important)
	      if(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<distMin) 
		distMin=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

              //We check if the step is no too small (otherwise the simulation will take too much time)
	      if(dt<dtMin)
		{
		  stopStepper=false;
		  isdtMin=true;
		}

	    }

	  //We store the asymptotic velocity in a container of map type with a view to make a data binning
	  mySpectra.storeDataBinning(x, myIC.weightIonization, isdtMin);
	         
	    if(isdtMin==true)
	      dtMinNbr+=1;

	  //We update the load bar
	  myDisplay.loadbar(iVPerpBirth+(iFieldBirth-1)*nVPerpBirth, nFieldBirth*nVPerpBirth);

	}
    }

  //Finally we write the data binning in the file "dataFile"
   mySpectra.writeDataBinning(dataFile);

   fstream gnuFile("data.gnu", ios::out);

   myPlot.gnuplot(gnuFile);
 
   gnuFile.close();
   system("gnuplot data.gnu");

  return 0;
}

