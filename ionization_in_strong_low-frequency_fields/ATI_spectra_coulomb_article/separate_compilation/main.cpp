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
#include"electrostaticpotential.h"
#include"electricfield.h"


using namespace std;

//Plateau in above-threashold-ionization spectra and chaotic behavior in rescattering processes
//Bambi Hu, Jie Liu, Shi-gang Chen

//VARIABLES DECLARATION

//Numbers of computed points
int nFieldBirth=100, nVPerpBirth=1000;
int iFieldBirth, iVPerpBirth;

//We declare the time variable
double t;

//We declare the array in which we will store the orbit
typedef vector<double> state_type;
state_type x(6);

//We declare a variable which contains the current nbr of trajectories
int nTraj;

//We declare the duration during which we want to let the electron propagates 
double trajDuration=500.;

//We declare a variable which will contain the minimum approach distance
double distMin;

//We declare a variable for the step in controlledRK5, the min allowed value and a counter which count how many events have not been accepted because dt was smaller than dtMin
double dt=0.001; 
double dtMin=1E-20;
int dtMinReachedNbr=0;

//We declare runge kutta error, its max allowed value, and the desired error min and max
double error;
double errorMax=1E-25;
double desiredErrorMax=1E-13;
double desiredErrorMin=desiredErrorMax/10.;

//We declare boolean controls
bool stopStepper;
bool isdtMin;

//Bins width
double binsWidth;

//We open files with a view to writing in them 
fstream dataFile("data.dat",ios::out);


//FUNCTION MAIN

int main()
{

//We leave few lines break
cout<<" " <<endl;
cout<<" "<<endl;

//We instantiate the differente object implemented in the others .cpp files
//Each object will execute a specific task

//Contains the electrostatic potential properties
ElectrostaticPotential<state_type> myPotential;

//Contains the electric field properties
ElectricField myField;

//Sets the initial condition for the ionization probability, perpendicular velocity, field at birth, electron position at birth
IC<state_type> myIC(myPotential, myField);

//Contains the ordinary differential system of equations of motion of the electron in the electrostatic potential and the electric field
System<state_type> mySystem(myPotential, myField);	

//Contains the method which solves ODE: runge kutta 5 with controlled step-size algorithm
Solve<state_type> mySolve;

//Contains methods for outputting data in terminal
Display myDisplay;

//Contains methods for doing a binning procedure and build a spectrum
Spectra<state_type> mySpectra(myPotential, myField);

//Contains methods for drawing curves
Plot myPlot;


  //We perform two loops
  //first, for each ionization time (initial field value)
  //second, for each perpendicular velocity

  for(iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(iVPerpBirth=1; iVPerpBirth<=nVPerpBirth; iVPerpBirth++)
	{
	  
          //We move the cursor back up with a view to rewriting on previous script and displaying a stable output
          myDisplay.moveCursorBackUp();	

	  //INITIAL CONDITIONS
          //We set the ionization time
	  myIC.setTBirth(iFieldBirth, nFieldBirth);
	  myIC.setFieldBirth();
	  myIC.setVPerpBirth(iVPerpBirth, nVPerpBirth);
	  myIC.setRhoBirth();
	  myIC.setIC(x,t);

	  //We initialise the boolean controls	  
	  stopStepper=true; 
	  isdtMin=false;
	  distMin=100.;

	  //We compute the trajectory

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

	  //We store the asymptotic velocity in a container of map type with a view to making a data binning
	  mySpectra.storeDataBinning(x, t, myIC.weightIonization, isdtMin);
	         
	    if(isdtMin==true)
	      dtMinReachedNbr+=1;

	  //We update the load bar and display some informations
          if(iVPerpBirth%10==0)
          {
	  myDisplay.loadbar(iVPerpBirth+(iFieldBirth-1)*nVPerpBirth, nFieldBirth*nVPerpBirth);
	  myDisplay("distMin",distMin);
          myDisplay("rhoBirth",myIC.rhoBirth);
          myDisplay("phaseBirth",myField.pulsation*myIC.tBirth*180./M_PI);
          myDisplay("vPerpBirth", myIC.vPerpBirth);
          myDisplay("step", dt);
          myDisplay("stepMin", dtMin);
          myDisplay("error", error);
          myDisplay("errorMin", desiredErrorMin);
          myDisplay("asymptoticEnergy",mySpectra.asymptoticEnergy(x,t));
          myDisplay("weightIonization",myIC.weightIonization);
          myDisplay("spectraPointNbr", mySpectra.spectraPointsNbr);
          myDisplay("dtMinReachedNbr", dtMinReachedNbr/mySpectra.spectraPointsNbr*100., "%");

          }
	}
    }

  //Finally we write the data binning in the file "dataFile"
   mySpectra.writeDataBinning(dataFile);

  //We build the legend of the plot
   myPlot.addKey("nField",nFieldBirth);
   myPlot.addKey("nVPerp",nVPerpBirth);
   myPlot.addKey("ErrorMax",desiredErrorMax);
   myPlot.addKey("dtMin",dtMin);
   myPlot.addKey("dtMinReachedNbr",dtMinReachedNbr);
   myPlot.addKey("linear field");
   myPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
   myPlot.addKey("waveLenght",myField.waveLenght*1.E9, "nm");
   myPlot.addKey("duration",myDisplay.elapsedTime);
 
   myPlot.gnuplot();


  return 0;
}

