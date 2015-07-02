#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#ifdef _OPENMP
#include<omp.h>
#endif

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

//VARIABLES DECLARATION

//Numbers of computed points
int nFieldBirth=100, nVYPerpBirth=1, nVZPrimPerpBirth=100;
int iFieldBirth, iVYPerpBirth, iVZPrimPerpBirth;

//We declare some variables for OPENMP information
int threadID;
int threadsNbr;
int threadsNbrMax;

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
//a counter which counts how many initial conditions have not been accepted because the probality of ionization was too small
int weightTooSmallNbr=0;

//We declare runge kutta error, its max allowed value, and the desired error min and max
double error;
double desiredErrorMax=1E-12;
double desiredErrorMin=desiredErrorMax/10.;

//We declare a minimum threshold value for the probability of ionization
double weightThresholdRatio=5.;
double weightThreshold;

//We declare boolean controls
bool stopStepper;
bool isStepTooSmall;
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

  /********We instantiate the differente object implemented in the others .cpp files***********/
  /********Each object will execute a specific task********************************************/

  //Contains the electrostatic potential properties
  Hydrogen<state_type> *myPotential=new Hydrogen<state_type>;
  myPotential->setIP(0.5792);
  //Molecule<state_type> *myPotential=new Molecule<state_type>();

  //Contains the electric field properties
  ElectricField myField(0.5);

  //Sets the initial condition for the ionization probability, perpendicular velocity, field at birth, electron position at birth
  IC<state_type> myIC(myPotential, myField);

  //Contains the ordinary differential system of equations of motion of the electron in the electrostatic potential and the electric field
  System<state_type> mySystem(myPotential, myField);	

  //Contains the method which solves ODE: runge kutta 5 with controlled step-size algorithm
  Solve<state_type> mySolve;

  //Contains methods for outputting data in terminal
  Display myDisplay;

  //Contains methods for doing a binning procedure and build a spectrum
  Spectra<state_type> mySpectra(myPotential, myField, &myIC, 180., 0.01);

  //Contains methods for drawing curves
  Plot myPlot;

  //We set the ionization rate threshold
  weightThreshold=myIC.getMaxWeightIonization(2)/weightThresholdRatio;
  
  //We declare a 2-dimensional array in which we will store the data
  //The first dimension accounts for the number of the trajectory
  //The second dimension accounts for the asymptotic energy (index=0), the tunneling ionization probability (index=1), a bool which tell if the trajectory is correct (index=2) and the profile of the trajectory (index=3), for instance we can distinguish electrons detected upward than those detected downward according the laser field polarization
double **dataBinning = new double*[nFieldBirth*nVYPerpBirth*nVZPrimPerpBirth];
for(int i = 0; i < nFieldBirth*nVYPerpBirth*nVZPrimPerpBirth; ++i) 
{
    dataBinning[i] = new double[4];
}


  /************************************We perform 3 loops**********************************/
  /***************First, for each ionization time (initial field value)********************/
  /***************Second and third, for each perpendicular velocity************************/

#ifdef _OPENMP
  threadsNbrMax=omp_get_max_threads();
  omp_set_num_threads(threadsNbrMax);
#endif

#pragma omp parallel for schedule(dynamic) collapse(3) private(x,t,error,step,stopStepper,isStepTooSmall,isWeightTooSmall,myDisplay) firstprivate(myIC,desiredErrorMax,desiredErrorMin) 

  for(iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(iVYPerpBirth=0; iVYPerpBirth<nVYPerpBirth; iVYPerpBirth++)
	{

	  for(iVZPrimPerpBirth=0; iVZPrimPerpBirth<nVZPrimPerpBirth; iVZPrimPerpBirth++)
	    {

	      /**************************We set the initial conditions********************************/
	      myIC.setTBirth(iFieldBirth, nFieldBirth);
	      myIC.setFieldBirth();
	      myIC.setVYPerpBirth(iVYPerpBirth, nVYPerpBirth);
	      myIC.setVXZPerpBirth(iVZPrimPerpBirth, nVZPrimPerpBirth);
	      myIC.setWeightIonization();
         
	      myIC.setRhoBirth();
	      myIC.setPolarCoordBirth();
	      myIC.setIC(x,t);

	      /**************************We compute the trajectory***********************************/

	      //We initialise the boolean controls	  
	      stopStepper=false; 
	      isStepTooSmall=false;
	      isWeightTooSmall=false;

	      //We check if weightIonization is big enough
	      if(myIC.weightIonization < weightThreshold)
		{ 
		  isWeightTooSmall=true;
		  stopStepper=true;
		}

	      //We update the step step
	      step=0.0001; 

	      for(int nTraj=0; !stopStepper ; nTraj++)
		{ 
		  //We call the function which solve eq of the motion
		  mySolve.controlledRK5(mySystem,x,t,step,error,desiredErrorMin,desiredErrorMax);

		  //We wait long enough for the end of the pulse
		  if(t>4.*myField.cyclesNbr*myField.opticalCycle)
		    stopStepper=true;

		  //We check if the step is no too small (otherwise the simulation will take too much time)
		  if(step<stepMin)
		    {
		      stopStepper=true;
		      isStepTooSmall=true;
		    }
		}
            
	      if(isStepTooSmall==true)
		stepTooSmallNbr+=1;
	      if(isWeightTooSmall==true)
		weightTooSmallNbr+=1;

	      //We store the asymptotic energy, the tunneling rate and bool which tell us if the trajectory has a good profile and which profile
	      int index=iVYPerpBirth+nVYPerpBirth*(iVZPrimPerpBirth+(iFieldBirth-1)*nVZPrimPerpBirth);
	      dataBinning[index][0]=mySpectra.asymptoticEnergy(x,t);
	      dataBinning[index][1]=myIC.weightIonization;
	      dataBinning[index][2]=mySpectra.hasTrajectoryGoodProfile(x,t,isStepTooSmall || isWeightTooSmall);
	      dataBinning[index][3]=mySpectra.whichProfile(x,t);
	      	      	      
/*#pragma omp critical
	      mySpectra.storeDataBinning(x,t, myIC.weightIonization, isStepTooSmall || isWeightTooSmall);*/
	   

	      if(iVYPerpBirth+nVYPerpBirth*(iVZPrimPerpBirth+(iFieldBirth-1)*nVZPrimPerpBirth)%500==0)
		{
		#pragma omp critical
		 {
		  //We move the cursor back up with a view to rewriting on previous script and displaying a stable output
	     	  myDisplay.moveCursorBackUp();  
	     	  
	          //We update the load bar and display some informations
	    	  myDisplay.loadbar(iVYPerpBirth+nVYPerpBirth*(iVZPrimPerpBirth+(iFieldBirth-1)*nVZPrimPerpBirth),nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
		  myDisplay("rhoBirth",myIC.rhoBirth);
		  myDisplay("phaseBirth",myField.pulsation*myIC.tBirth*180./M_PI);
		  myDisplay("vPerpBirth", myIC.vPerpBirth);
		  myDisplay("fieldpBirth", myIC.fieldBirth);
		  myDisplay("vYPerpBirth", myIC.vYPerpBirth);
		  myDisplay("vZPrimPerpBirth", myIC.vZPrimPerpBirth);
		  myDisplay("ellipticity", myField.ellipticity);

		  // myDisplay.variableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
		  //myDisplay.variableArg<double>("bondLength", 3, myPotential->bondLength[0],  myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);

		  myDisplay("step", step);
		  myDisplay("stepMin", stepMin);
		  myDisplay("error", error);
		  myDisplay("errorMin", desiredErrorMin);

		  myDisplay("asymptoticEnergy",mySpectra.asymptoticEnergy(x,t));
		  myDisplay("angleDetection",mySpectra.angleDetection, "degree");
		  myDisplay("binsWidth",mySpectra.binsWidth);

		  myDisplay("weightIonization",myIC.weightIonization);
		  myDisplay.variableArg<double>("weightThreshold",2,weightThreshold, weightThresholdRatio);
		  myDisplay("weightTooSmallNbr", double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");

		  myDisplay("trappedElectronNbr", double(mySpectra.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("stepTooSmallNbr", double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("spectraPointNbr", double(mySpectra.spectraPointsNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
#ifdef _OPENMP
		  myDisplay("omp_get_thread_num", omp_get_thread_num());
		  myDisplay("omp_get_num_threads", omp_get_num_threads());
		  myDisplay("omp_get_max_threads", threadsNbrMax);
		  myDisplay("omp_in_parallel", omp_in_parallel());
		  myDisplay("omp_get_dynamic", omp_get_dynamic());
		  myDisplay("omp_get_nested", omp_get_nested());
#endif
		 }
		}
	    }         
	}
    }

  //We perform the data binning process
  mySpectra.storeDataBinning(dataBinning,nFieldBirth*nVYPerpBirth*nVZPrimPerpBirth); 
  //Finally we write the data binning in the file "dataFile"
  mySpectra.writeDataBinning(dataFile);

  //We build the legend of the plot
  myPlot.addKey("nField",nFieldBirth);
  myPlot.addKey("nVYPerp",nVYPerpBirth);
  myPlot.addKey("nVZPrimPerp",nVZPrimPerpBirth);
  // myPlot.addKeyVariableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
  // myPlot.addKeyVariableArg<double>("bondLength", 3, myPotential->bondLength[0],  myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);
  myPlot.addKey("weightTooSmallNbr", int(double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPlot.addKey("stepTooSmallNbr",int(double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPlot.addKey("trappedElectronNbr", int(double(mySpectra.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPlot.addKey("angleTooLargeNbr", int(double(mySpectra.angleTooLargeNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPlot.addKey("spectraPointsNbr", double(mySpectra.spectraPointsNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
  myPlot.addKey("weightThresholdRatio",weightThresholdRatio);
  #ifdef _OPENMP
  myPlot.addKey("threadsNbrMax", threadsNbrMax);
  #endif
  myPlot.addKey("binsWidth",mySpectra.binsWidth);
  myPlot.addKey("ErrorMax",desiredErrorMax);
  myPlot.addKey("stepMin",stepMin);
  myPlot.addKey("angleDetection",mySpectra.angleDetection, "deg");
  myPlot.addKey("ellipticity", myField.ellipticity);
  myPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
  myPlot.addKey("waveLenght",myField.waveLenght*1.E6, "micro-m");
  myPlot.addKey("duration",myDisplay.elapsedTime());
     
  myPlot.addInstruction("set xlabel 'Asymptotic energy (eV)' offset 0,4");
  myPlot.addInstruction("set ylabel 'Probability (linear scale)'");

  //myPlot.addInstruction("set xrange [0:8]");
  myPlot.addInstruction("set xtics offset 0,0.3");

  myPlot.addInstruction("set style line 1 lc rgb '#db0000' pt 6 ps 1 lt 1 lw 2 "); //red
  myPlot.addInstruction("set style line 2 lc rgb '#062be5' pt 6 ps 1 lt 1 lw 2 "); //blue '#0060ad'
 
  myPlot.setPlotType("plot");

  //We consider electrons differently depending if they are detected along the polarization of the field or not
  myPlot.addPlot("'data.dat' index 0 using 1:2 w l ls 1 title 'Photo-electron detected upward according the initial field'");
  myPlot.addPlot("'data.dat' index 1 using 1:2 w l ls 2 title 'Photo-electron detected downward according the initial field'");

  myPlot.gnuplot();

  return 0;
}

