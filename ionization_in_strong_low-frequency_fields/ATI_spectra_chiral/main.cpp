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
int nFieldBirth=100, nVYPerpBirth=1000, nVZPrimPerpBirth=500;
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

//We select which electrostatic potential we want to use
#define MOLECULE
//#define HYDROGEN

#ifdef HYDROGEN
typedef Hydrogen<state_type> potential_type;
double IP=0.5792;
#endif

#ifdef MOLECULE
typedef Molecule<state_type> potential_type;
int myOrientation=38;
#endif

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
double binsWidth=0.05;

//Ellipticity
double ellipticity=0.1;

//Angle between velocity vector and field polarization within which we detect electrons
double angleDetection=180.;


//FUNCTION MAIN

int main()
{

  //We leave few lines break
  cout<<" " <<endl;
  cout<<" "<<endl;

  /********We instantiate the differente object implemented in the others .cpp files****************/
  /********Each object will execute a specific task*************************************************/

  //Contains the electrostatic potential properties
  potential_type *myPotential=new potential_type;
  
#ifdef HYDROGEN
   myPotential->setIP(IP);
#endif

#ifdef MOLECULE
    myPotential->setLebedevOrientation(myOrientation);
#endif
  
  //Contains the electric field properties
  ElectricField myField(ellipticity);

  //Sets the initial condition for the ionization probability, perpendicular velocity, field at birth, electron position at birth
  IC<state_type> myIC(myPotential, myField);

  //Contains the ordinary differential system of equations of motion of the electron in the electrostatic potential and the electric field
  System<state_type,potential_type > mySystem(*myPotential, myField);	

  //Contains the method which solves ODE: runge kutta 5 with controlled step-size algorithm
  Solve<state_type,potential_type > mySolve;

  //Contains methods for outputting data in terminal
  Display myDisplay;

  //Contains methods for drawing curves
  Plot myPlot;

  //We set the ionization rate threshold
  weightThreshold=myIC.getMaxWeightIonization(2)/weightThresholdRatio;
  

  /************************************We perform 3 loops**********************************/
  /***************First, for each ionization time (initial field value)********************/
  /***************Second and third, for each perpendicular velocity************************/

#ifdef _OPENMP
  threadsNbrMax=omp_get_max_threads();
  omp_set_num_threads(threadsNbrMax);
#endif

#ifndef _OPENMP
  //Contains methods for doing a binning procedure and build a spectrum
  Spectra<state_type> mySpectrum(myPotential, myField, &myIC, angleDetection, binsWidth);
#endif

#ifdef _OPENMP
  vector<Spectra<state_type> > mySpectra(threadsNbrMax,Spectra<state_type>(myPotential, myField, &myIC, angleDetection, binsWidth));
#endif

#pragma omp parallel for schedule(dynamic) collapse(3) private(x,t,error,step,stopStepper,isStepTooSmall,isWeightTooSmall) firstprivate(myIC,desiredErrorMax,desiredErrorMin, mySystem) 

  for(iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(iVYPerpBirth=0; iVYPerpBirth<=nVYPerpBirth; iVYPerpBirth++)
	{

	  for(iVZPrimPerpBirth=0; iVZPrimPerpBirth<=nVZPrimPerpBirth; iVZPrimPerpBirth++)
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

 		 //We perform the data binning process and store data in map container
#ifndef _OPENMP
   mySpectrum.storeDataBinning(x,t, myIC.weightIonization, isStepTooSmall || isWeightTooSmall);
#endif

#ifdef _OPENMP
   mySpectra[omp_get_thread_num()].storeDataBinning(x,t, myIC.weightIonization, isStepTooSmall || isWeightTooSmall);
#endif
	     
	      if(iVYPerpBirth+nVYPerpBirth*(iVZPrimPerpBirth+(iFieldBirth-1)*nVZPrimPerpBirth)%50000==0)
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
#ifdef MOLECULE
                  myDisplay("Molecule orientation", myPotential->myOrientation);
		  myDisplay.variableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
		  myDisplay.variableArg<double>("bondLength", 3, myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);
#endif		  

		  myDisplay("step", step);
		  myDisplay("stepMin", stepMin);
		  myDisplay("error", error);
		  myDisplay("errorMin", desiredErrorMin);

		  myDisplay("weightIonization",myIC.weightIonization);
		  myDisplay.variableArg<double>("weightThreshold",2,weightThreshold, weightThresholdRatio);
		  myDisplay("weightTooSmallNbr", double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");

#ifndef _OPENMP
	          myDisplay("asymptoticEnergy",mySpectrum.asymptoticEnergy(x,t));
		  myDisplay("angleDetection",mySpectrum.angleDetection, "degree");
		  myDisplay("binsWidth",mySpectrum.binsWidth);
		  myDisplay("trappedElectronNbr", double(mySpectrum.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("stepTooSmallNbr", double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("spectraPointNbr", double(mySpectrum.electronsDetectedNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
#endif
		  
#ifdef _OPENMP
		  myDisplay("asymptoticEnergy", mySpectra[omp_get_thread_num()].asymptoticEnergy(x,t));
		  myDisplay("angleDetection", mySpectra[omp_get_thread_num()].angleDetection, "degree");
		  myDisplay("binsWidth", mySpectra[omp_get_thread_num()].binsWidth);
		  myDisplay("trappedElectronNbr", double( mySpectra[omp_get_thread_num()].trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("stepTooSmallNbr", double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("spectraPointNbr", double( mySpectra[omp_get_thread_num()].electronsDetectedNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);

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

/*******We merge all the data binning computed by each thread and write it in a file***********************************/

#ifdef _OPENMP
    Spectra<state_type> mySpectrum(myPotential, myField, &myIC, angleDetection, binsWidth);
    //We merge all the data binning computed by each thread
    mySpectrum.mergeSpectra(mySpectra);
#endif

   //We open a file with a view to writing in it
   std::ostringstream dataFileStream;
   dataFileStream<<"leb_orient_"<<myPotential->myOrientation<<".dat";
   fstream dataFile(dataFileStream.str().c_str(),ios::out);

   //We write the data binning in the file "dataFile"
   mySpectrum.writeDataBinning(dataFile);

/*****************************************************PLOT*****************************************************/
/****************************************We build the key of the plot******************************************/
  myPlot.addKey("nField",nFieldBirth);
  myPlot.addKey("nVYPerp",nVYPerpBirth);
  myPlot.addKey("nVZPrimPerp",nVZPrimPerpBirth);
                
#ifdef MOLECULE
  myPlot.addKey("Molecule orientation",myPotential->myOrientation);
  myPlot.addKeyVariableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
  myPlot.addKeyVariableArg<double>("bondLength", 3, myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);
#endif

#ifdef HYDROGEN
  myPlot.addKey("hydrogen potential");
  myPlot.addKey("IP", myPotential->IP);
#endif
  
  myPlot.addKey("weightTooSmallNbr", int(double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPlot.addKey("stepTooSmallNbr",int(double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  
    myPlot.addKey("trappedElectronNbr", int(double(mySpectrum.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
    
  if(mySpectrum.angleDetection<180)
    myPlot.addKey("angleTooLargeNbr", int(double(mySpectrum.angleTooLargeNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
    myPlot.addKey("electronsDetectedNbr", double(mySpectrum.electronsDetectedNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
    myPlot.addKey("binsWidth",mySpectrum.binsWidth,"eV");
    
  if(mySpectrum.angleDetection<180)
  myPlot.addKey("angleDetection",mySpectrum.angleDetection, "deg");
  
  myPlot.addKey("weightThresholdRatio",weightThresholdRatio);
  
  #ifdef _OPENMP
  myPlot.addKey("threadsNbrMax", threadsNbrMax);
  #endif

  myPlot.addKey("ErrorMax",desiredErrorMax);
  myPlot.addKey("stepMin",stepMin);
  
  myPlot.addKey("ellipticity", myField.ellipticity);
  myPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
  myPlot.addKey("waveLenght",myField.waveLenght*1.E6, "micro-m");
  myPlot.addKey("duration",myDisplay.elapsedTime());
  if(myField.ellipticity>0)
  myPlot.addKey("laser polarization rotates clockwise when we look along y");  
  if(myField.ellipticity<0)
  myPlot.addKey("laser polarization rotates counterclockwise when we look along y"); 
  
  //We put the keys on
  myPlot.setKeysOn();
    
  /****************************************We build the build the different parts of plot******************************************/
  
  myPlot.addInstruction("set terminal epscairo font 'Gill Sans,9' rounded fontscale 0.4");
  myPlot.addInstruction("set multiplot layout 1, 1");
  //Remove border on bottom and right, these borders are useless and make it harder to see plotted lines near the border
  //Also, put it in grey, no need for so much emphasis on a border 
  myPlot.addInstruction("set border 6 back linestyle 100");
  myPlot.addInstruction("set grid linestyle 101");
  //We also put the key box in grey
  myPlot.addInstruction("set key box ls 100");

  myPlot.addInstruction("set xtics scale 0");
  myPlot.addInstruction("set format x '' ");
  myPlot.addInstruction("set xrange [0:20]");
  myPlot.addInstruction("set x2tics 1");
  myPlot.addInstruction("set x2label 'Asymptotic energy (eV)'");
  myPlot.addInstruction("set ylabel 'Probability (linear scale)'");
   
  /**************************************We indicate which curves we want to plot************************************/   
   
  
  //We consider electrons differently depending if they are detected in y>0 or y<0
  std::ostringstream plot1;
  plot1<<"plot 'data.dat' index 0 using 1:2 w l ls 1 title 'Photo-electrons detected in the upper half-space (y>0), number="<< double(mySpectrum.electronsDetectedUPNbr)/(mySpectrum.electronsDetectedNbr)*100.<<" %', \\";
  std::ostringstream plot2;
  plot2<<"'data.dat' index 1 using 1:2 w l ls 2 title 'Photo-electrons detected in the lower half-space (y<0), number="<< double(mySpectrum.electronsDetectedDOWNNbr)/(mySpectrum.electronsDetectedNbr)*100.<<" %'";

  myPlot.addInstruction(plot1.str());
  myPlot.addInstruction(plot2.str());
  
  myPlot.addInstruction("unset multiplot");

  myPlot.gnuplot();

  return 0;
}

