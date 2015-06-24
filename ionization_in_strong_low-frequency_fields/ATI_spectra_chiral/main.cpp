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
int nFieldBirth=100, nVYPerpBirth=1, nVZPrimPerpBirth=100;
int iFieldBirth, iVYPerpBirth, iVZPrimPerpBirth;

//Closest approach
double closestApproach;

//We declare the time variable
double t;

//We declare the array in which we will store the orbit
typedef double state_type[6];
state_type x;


//We declare a variable for the step in controlledRK5, the min allowed value
double step; 
double stepMin=1E-20;

//We declare a counter which counts how many events have not been accepted because step was smaller than stepMin
int stepTooSmallNbr=0;
//And a counter which counts how many initial conditions have not been accepted because the probality of ionization was too small
int weightTooSmallNbr=0;

//We declare runge kutta error, its max allowed value, and the desired error min and max
double error;
double desiredErrorMax=1E-14;
double desiredErrorMin=desiredErrorMax/10.;

//We declare a minimum threshold value for the probability of ionization
double weightThreshold=0.;
//We declare boolean controls
bool stopStepper;
bool isStepTooSmall;
bool isWeightTooSmall;

//Bins width
double binsWidth;

//We open files with a view to writing in them 
fstream plotFile("plot.dat",ios::out);
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
Hydrogen<state_type> *myPotential=new Hydrogen<state_type>;
myPotential->setIP(0.5792);
//Molecule<state_type> *myPotential=new Molecule<state_type>();

//Contains the electric field properties
ElectricField myField(0.0);

//Sets the initial condition for the ionization probability, perpendicular velocity, field at birth, electron position at birth
IC<state_type> myIC(myPotential, myField);

//Contains the ordinary differential system of equations of motion of the electron in the electrostatic potential and the electric field
System<state_type> mySystem(myPotential, &myField);	

//Contains the method which solves ODE: runge kutta 5 with controlled step-size algorithm
Solve<state_type> mySolve;

//Contains methods for outputting data in terminal
Display myDisplay;

//Contains methods for doing a binning procedure and build a spectrum
Spectra<state_type> mySpectra(myPotential, myField,0.002);

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
	  isStepTooSmall=false;
          isWeightTooSmall=false;
          
          //We prepare closest approach variable
          closestApproach=100.;

          //We update the step step
          step=0.001; 

	  //INITIAL CONDITIONS
          //We set the ionization time
	  myIC.setTBirth(iFieldBirth,nFieldBirth);
	  myIC.setFieldBirth();
	  myIC.setVYPerpBirth(iVYPerpBirth,nVYPerpBirth);
	  myIC.setVXZPerpBirth(iVZPrimPerpBirth, nVZPrimPerpBirth);
          myIC.setWeightIonization();

          //We check if weightIonization is big enough
          if(myIC.weightIonization < weightThreshold)
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
	     mySolve.controlledRK5(mySystem,x,t,step,error,desiredErrorMin,desiredErrorMax);

             //plotFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<x[4]<<" "<<x[5]<<" "<<t<<" "<<myField('Z',t)<<endl;
	
             //If the electron is always bonded to the attractor, we do not consider the event 
	     //if(t>10.*myField.cyclesNbr*myField.opticalCycle)
	     if(t>8.*myField.cyclesNbr*myField.opticalCycle)                
	     stopStepper=true;

	     if(t>myField.opticalCycle/2. && sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<closestApproach)
	     closestApproach=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
                 
             //We check if the step is no too small (otherwise the simulation will take too much time)
	     if(step<stepMin)
		{
		  stopStepper=true;
		  isStepTooSmall=true;
		}
	    }


	  //We write data in a file
	    mySpectra.writeData(x, t, myIC.tBirth, myIC.fieldBirth, myIC.vPerpBirth, closestApproach, myIC.weightIonization, dataFile, isStepTooSmall || isWeightTooSmall);
	   
            
	    if(isStepTooSmall==true)
	      stepTooSmallNbr+=1;
            if(isWeightTooSmall==true)
	      weightTooSmallNbr+=1;

		    //We update the load bar and display some informations
          //if(iVYPerpBirth+nVYPerpBirth*((iVZPrimPerpBirth-1)+(iFieldBirth-1)*nVZPrimPerpBirth)%1000==0)
          {
	  myDisplay.loadbar(iVYPerpBirth+nVYPerpBirth*((iVZPrimPerpBirth-1)+(iFieldBirth-1)*nVZPrimPerpBirth),nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
          myDisplay("ellipticity", myField.ellipticity);
          myDisplay("rhoBirth",myIC.rhoBirth);
          myDisplay("phaseBirth",myField.pulsation*myIC.tBirth*180./M_PI);
          myDisplay("vPerpBirth", myIC.vPerpBirth);
          myDisplay("fieldpBirth", myIC.fieldBirth);
          myDisplay("vYPerpBirth", myIC.vYPerpBirth);
          myDisplay("vZPrimPerpBirth", myIC.vZPrimPerpBirth);
         // myDisplay.variableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
          //myDisplay.variableArg<double>("bondLength", 3, myPotential->bondLength[0],  myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);
          myDisplay("step", step);
          myDisplay("stepMin", stepMin);
          myDisplay("error", error);
          myDisplay("errorMin", desiredErrorMin);
          myDisplay("asymptoticEnergy",mySpectra.asymptoticEnergy(x,t));
          myDisplay("weightIonization",myIC.weightIonization);
          myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
          myDisplay("spectraPointNbr", double(mySpectra.spectraPointsNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          myDisplay("trappedElectronNbr", double(mySpectra.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          myDisplay("stepTooSmallNbr", double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          myDisplay("weightTooSmallNbr", double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
          myDisplay("weightThreshold",weightThreshold);
          myDisplay("binsWidth",mySpectra.binsWidth);
          }


plotFile<<" "<<endl;
plotFile<<" "<<endl;



	 }         
      }
    }


  //We build the legend of the plot
  /* myPlot.addKey("nField",nFieldBirth);
   myPlot.addKey("nVYPerp",nVYPerpBirth);
   myPlot.addKey("nVZPrimPerp",nVZPrimPerpBirth);*/
  // myPlot.addKeyVariableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
  // myPlot.addKeyVariableArg<double>("bondLength", 3, myPotential->bondLength[0],  myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);

 /*  myPlot.addKey("spectraPointNbr", double(mySpectra.spectraPointsNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
   myPlot.addKey("stepTooSmallNbr",int(double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
   myPlot.addKey("weightTooSmallNbr", int(double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
   myPlot.addKey("weightThreshold",weightThreshold);
   myPlot.addKey("binsWidth",mySpectra.binsWidth);*/
   
   /*
   myPlot.addKey("tBirth",myIC.tBirth);
   myPlot.addKey("vPerpBirth",myIC.vPerpBirth);
   myPlot.addKey("weightIonization",myIC.weightIonization);
   myPlot.addKey("asymptoticEnergy",mySpectra.asymptoticEnergy(x,t));
   //myPlot.addKey("ErrorMax",desiredErrorMax);
   myPlot.addKey("stepMin",stepMin);
   myPlot.addKey("ellipticity", myField.ellipticity);
   myPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
   myPlot.addKey("waveLenght",myField.waveLenght*1.E9, "nm");
   myPlot.addKey("duration",myDisplay.elapsedTime);
*/

//Specific to the article 8
   //myPlot.addInstruction("set xrange [0:8]");


   myPlot.addInstruction("set multiplot  layout 3, 2 title 'Searching for soft-recollision trajectory'");

   myPlot.addInstruction("set style line 1 lc rgb '#db0000' pt 6 ps 1 lt 1 lw 1 "); //red
   myPlot.addInstruction("set style line 2 lc rgb '#062be5' pt 6 ps 1 lt 1 lw 1 "); //blue '#0060ad'
 
   //myPlot.setPlotType("plot");
   //myPlot.addPlot("'data.dat' index 0 using 3:2 w l ls 1 title 'Photo-electron trajectories, errorMax=1E-6'");

/*
std::ostringstream asymptEnergy0;
asymptEnergy0<<"plot 'data.dat' index 0 using 3:2 w l ls 1 title '(a) tBirth=8.050, vPerpBirth=0.123, asymptEnergy="<<asymptEnergy[0]*37.3<<" (eV)'";

myPlot.addInstruction(asymptEnergy0.str());
*/

   //IMPORTANT: in order to align plots in column we have to fix the xrange, the lmargin and rmargin
   std::ostringstream xrange;
   xrange<<"set xrange [0:"<<mySpectra.spectraPointsNbr<<"]";
   myPlot.addInstruction(xrange.str());

   myPlot.addInstruction("unset xtics");
   myPlot.addInstruction("unset ytics");

   myPlot.addInstruction("set bmargin 1");
   myPlot.addInstruction("set lmargin 12");
   myPlot.addInstruction("set rmargin 3");
   myPlot.addInstruction("set ytics");
   myPlot.addInstruction("set ylabel 'fieldBirth (au)'");
   myPlot.addInstruction("plot 'data.dat' using 1:2 w l ls 1 notitle");

   myPlot.addInstruction("set lmargin 0");
   myPlot.addInstruction("set rmargin 13");
   myPlot.addInstruction("unset ytics");
   myPlot.addInstruction("unset ylabel");
   myPlot.addInstruction("set y2tics");
   myPlot.addInstruction("set y2label 'velocityBirth (au)'");
   myPlot.addInstruction("plot 'data.dat' using 1:3 w l ls 1 notitle");

   myPlot.addInstruction("set tmargin 0");
   myPlot.addInstruction("set lmargin 12");
   myPlot.addInstruction("set rmargin 3");
   myPlot.addInstruction("unset y2tics");
   myPlot.addInstruction("unset y2label");
   myPlot.addInstruction("set ytics");
   myPlot.addInstruction("set ylabel 'closestApproach (au)'");
   myPlot.addInstruction("set yrange [0:10]");
   myPlot.addInstruction("plot 'data.dat' using 1:4 w l ls 1 notitle");

   myPlot.addInstruction("set lmargin 0");
   myPlot.addInstruction("set rmargin 13");
   myPlot.addInstruction("unset ytics");
   myPlot.addInstruction("unset ylabel");
   myPlot.addInstruction("unset yrange");
   myPlot.addInstruction("set y2tics");
   myPlot.addInstruction("set y2label 'weightIonization'");
   myPlot.addInstruction("plot 'data.dat' using 1:5 w l ls 1 notitle");

   myPlot.addInstruction("set bmargin 4");
   myPlot.addInstruction("set lmargin 12");
   myPlot.addInstruction("set rmargin 3");
   myPlot.addInstruction("unset y2tics");
   myPlot.addInstruction("unset y2label");
   myPlot.addInstruction("set ytics");
   myPlot.addInstruction("set ylabel 'asymptoticEnergy (eV)' ");
   myPlot.addInstruction("set yrange [0:8]");
   myPlot.addInstruction("set xlabel 'Trajectory ID' ");
   myPlot.addInstruction("set xtics offset 0,0.3");
   myPlot.addInstruction("plot 'data.dat' using 1:6 w l ls 1 notitle");

   myPlot.addInstruction("set lmargin 0");
   myPlot.addInstruction("set rmargin 13");
   myPlot.addInstruction("unset ytics");
   myPlot.addInstruction("unset ylabel");
   //myPlot.addInstruction("set key on outside bmargin");
  // myPlot.addInstruction("plot 'data.dat' using 1:6 w l ls 1  title 'asymptoticEnergy (eV)',\\");
  // myPlot.addInstruction("'data.dat' using 1:4 w l ls 2 title 'closestApproach (au)'");
   myPlot.addInstruction("set xrange [0:10]");
   myPlot.addInstruction("set xlabel 'closestApproach (au)' ");
   myPlot.addInstruction("set y2range [0:20]");
   myPlot.addInstruction("set y2tics");
   myPlot.addInstruction("set y2label 'asymptEnergy (eV)' ");
   myPlot.addInstruction("plot 'data.dat' using 4:6 w l ls 1 notitle");

   myPlot.addInstruction("unset multiplot");


   myPlot.gnuplot();

  return 0;
}

