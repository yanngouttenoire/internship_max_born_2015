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
#include"molecule.h"
#include"hydrogen.h"
#include"electricfield.h"


using namespace std;

//Plateau in above-threashold-ionization spectra and chaotic behavior in rescattering processes
//Bambi Hu, Jie Liu, Shi-gang Chen

//VARIABLES DECLARATION

//Numbers of computed points
int nFieldBirth=1000, nVYPerpBirth=1, nVZPrimPerpBirth=1000;
int iFieldBirth, iVYPerpBirth, iVZPrimPerpBirth;

//We declare the time variable
double t;

//We declare the array in which we will store the orbit
typedef double state_type[6];
state_type x;

//We declare a variable for the step in controlledRK5, the min allowed value
double dt=0.1; 
double dtMin=1E-20;

//We declare a counter which counts how many events have not been accepted because dt was smaller than dtMin
int unexpectedStopNbr=0;
//And a counter which counts how many initial conditions have not been accepted because the probality of ionization was too small
int weightTooSmallNbr=0;

//We declare runge kutta error, its max allowed value, and the desired error min and max
double error;
double desiredErrorMax=1E-10;
double desiredErrorMin=desiredErrorMax/10.;

//We declare a minimum threshold value for the probability of ionization
double weightMinThreshold=1E-6;

//We declare boolean controls
bool stopStepper;
bool unexpectedStop;
bool isWeightTooSmall;

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
//ElectrostaticPotential<state_type> myPotential;
Hydrogen<state_type> myPotential;
//Molecule<state_type> myPotential;

//Contains the electric field properties
ElectricField myField(0.);

//Sets the initial condition for the ionization probability, perpendicular velocity, field at birth, electron position at birth
IC<state_type> myIC(myPotential, myField);

//Contains the ordinary differential system of equations of motion of the electron in the electrostatic potential and the electric field
System<state_type> mySystem(myPotential, myField);	

//Contains the method which solves ODE: runge kutta 5 with controlled step-size algorithm
Solve<state_type> mySolve;

//Contains methods for outputting data in terminal
Display myDisplay;

//Contains methods for doing a binning procedure and build a spectrum
Spectra<state_type> mySpectra(myPotential, myField,0.005);

//Contains methods for drawing curves
Plot myPlot;


  //We perform two loops
  //first, for each ionization time (initial field value)
  //second, for each perpendicular velocity

  for(iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(iVZPrimPerpBirth=0; iVZPrimPerpBirth<nVZPrimPerpBirth; iVZPrimPerpBirth++)
       {
         for(iVYPerpBirth=0; iVYPerpBirth<nVYPerpBirth; iVYPerpBirth++)
	   {

	    //We move the cursor back up with a view to rewriting on previous script and displaying a stable output
          myDisplay.moveCursorBackUp();         	

          //We initialise the boolean controls	  
 	  stopStepper=false; 
	  unexpectedStop=false;
          isWeightTooSmall=false;

          //We update the step dt
          dt=0.001; 

	  //INITIAL CONDITIONS
          //We set the ionization time
	  myIC.setTBirth(iFieldBirth, nFieldBirth);
	  myIC.setFieldBirth();
	  myIC.setVYPerpBirth(iVYPerpBirth, nVYPerpBirth);
	  myIC.setVXZPerpBirth(iVZPrimPerpBirth, nVZPrimPerpBirth);
          myIC.setWeightIonization();

          //We check if weightIonization is big enough
          if(myIC.weightIonization < weightMinThreshold)
          { 
           isWeightTooSmall=true;
           stopStepper=true;
          }
         
	  myIC.setRhoBirth();
          myIC.setPolarCoordBirth();
	  myIC.setIC(x,t);

	  //We compute the trajectory

	    for(int nTraj=0; !stopStepper ; nTraj++)
	    { 

	      //We call the function which solve eq of the motion
	      mySolve.controlledRK5(mySystem,x,t,dt,error,desiredErrorMin,desiredErrorMax);

             //If the electron is always bonded to the attractor, we do not consider the event 
	      if((t-myIC.tBirth)>10.*myField.opticalCycle)
                {
                unexpectedStop=true;
		stopStepper=true;
                }
                 
             //We stop when the electron is fully ionized
	      if(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])>300.)
		stopStepper=true;

              //We check if the step is no too small (otherwise the simulation will take too much time)
	      if(dt<dtMin)
		{
		  stopStepper=true;
		  unexpectedStop=true;
		}
	    }

	  //We store the asymptotic velocity in a container of map type with a view to making a data binning
	  mySpectra.storeDataBinning(x, t, myIC.weightIonization, unexpectedStop || isWeightTooSmall);
	   
            
	    if(unexpectedStop==true)
	      unexpectedStopNbr+=1;
            if(isWeightTooSmall==true)
	      weightTooSmallNbr+=1;

		    //We update the load bar and display some informations
          if(iVYPerpBirth+nVYPerpBirth*((iVZPrimPerpBirth-1)+(iFieldBirth-1)*nVZPrimPerpBirth)%1000==0)
          {
	  myDisplay.loadbar(iVYPerpBirth+nVYPerpBirth*((iVZPrimPerpBirth-1)+(iFieldBirth-1)*nVZPrimPerpBirth),nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
          myDisplay("rhoBirth",myIC.rhoBirth);
          myDisplay("phaseBirth",myField.pulsation*myIC.tBirth*180./M_PI);
          myDisplay("vPerpBirth", myIC.vPerpBirth);
          myDisplay("vYPerpBirth", myIC.vYPerpBirth);
          myDisplay("vZPrimPerpBirth", myIC.vZPrimPerpBirth);
          myDisplay("step", dt);
          myDisplay("stepMin", dtMin);
          myDisplay("error", error);
          myDisplay("errorMin", desiredErrorMin);
          myDisplay("asymptoticEnergy",mySpectra.asymptoticEnergy(x,t));
          myDisplay("weightIonization",myIC.weightIonization);
          myDisplay("binsWidth",mySpectra.binsWidth);
          myDisplay("spectraPointNbr", mySpectra.spectraPointsNbr);
          myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
          myDisplay("unexpectedStopNbr", double(unexpectedStopNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          myDisplay("weightMinThreshold",weightMinThreshold);
          myDisplay("weightTooSmallNbr", double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          }

	 }         
      }
    }

  //Finally we write the data binning in the file "dataFile"
   mySpectra.writeDataBinning(dataFile);

  //We build the legend of the plot
   myPlot.addKey("nField",nFieldBirth);
   myPlot.addKey("nVYPerp",nVYPerpBirth);
   myPlot.addKey("nVZPrimPerp",nVZPrimPerpBirth);
   myPlot.addKey("weightMinThreshold",weightMinThreshold);
   myPlot.addKey("weightTooSmallNbr", double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
   myPlot.addKey("spectraPointNbr", mySpectra.spectraPointsNbr);
   myPlot.addKey("unexpectedStopNbr",unexpectedStopNbr);
   myPlot.addKey("binsWidth",mySpectra.binsWidth);
   myPlot.addKey("ErrorMax",desiredErrorMax);
   myPlot.addKey("dtMin",dtMin);
   myPlot.addKey("linear field");
   myPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
   myPlot.addKey("waveLenght",myField.waveLenght*1.E9, "nm");
   myPlot.addKey("duration",myDisplay.elapsedTime);
     
   myPlot.addInstruction("set xlabel 'Asymptotic energy (au)'");
   myPlot.addInstruction("set ylabel 'Probability (log)'");
   myPlot.addInstruction("set xrange [0:1]");
 
   myPlot.setPlotType("plot");
   myPlot.addPlot("'data.dat' index 0 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY positive'");
   myPlot.addPlot("'data.dat' index 1 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY negative'");

   myPlot.gnuplot();

  return 0;
}

