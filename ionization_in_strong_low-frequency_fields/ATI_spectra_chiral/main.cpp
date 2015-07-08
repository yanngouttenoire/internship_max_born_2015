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
int nFieldBirth=10000, nVYPerpBirth=1, nVZPrimPerpBirth=10000;

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
//#define MOLECULE
#define HYDROGEN

#ifdef HYDROGEN
typedef Hydrogen<state_type> potential_type;
double IP=0.5; //0.5792;
#endif

#ifdef MOLECULE
typedef Molecule<state_type> potential_type;
moleculeOrientation myOrientation(X3);
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
double weightThresholdRatio=10000.;
double weightThreshold;

//We declare boolean controls
bool stopStepper;
bool isStepTooSmall;
bool isWeightTooSmall;

//Bins width
double binsWidth=2.;

//Ellipticity
double ellipticity=0.;

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
    myPotential->setMoleculeOrientation(myOrientation);
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

  for(int iFieldBirth=1; iFieldBirth<=nFieldBirth; iFieldBirth++)
    {
      for(int iVYPerpBirth=0; iVYPerpBirth<nVYPerpBirth; iVYPerpBirth++)
	{

	  for(int iVZPrimPerpBirth=0; iVZPrimPerpBirth<nVZPrimPerpBirth; iVZPrimPerpBirth++)
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
		 // if(t>4.*myField.cyclesNbr*myField.opticalCycle)
		 if((t-myIC.tBirth)>2*myField.opticalCycle*myField.cyclesNbr)
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
#ifdef MOLECULE
                  myDisplay("Molecule orientation", myPotential->myOrientation.myString);
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
		  myDisplay("binsWidthEnergy",mySpectrum.binsWidthEnergy);
		  myDisplay("binsWidthAngle",mySpectrum.binsWidthAngle);
		  myDisplay("trappedElectronNbr", double(mySpectrum.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("stepTooSmallNbr", double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("spectraPointNbr", double(mySpectrum.electronsDetectedNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
		  myDisplay("ptsNumber", nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth);
#endif
		  
#ifdef _OPENMP
		  myDisplay("asymptoticEnergy", mySpectra[omp_get_thread_num()].asymptoticEnergy(x,t));
		  myDisplay("angleDetection", mySpectra[omp_get_thread_num()].angleDetection, "degree");
		  myDisplay("binsWidthEnergy", mySpectra[omp_get_thread_num()].binsWidthEnergy);
		  myDisplay("binsWidthAngle", mySpectra[omp_get_thread_num()].binsWidthAngle);
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
   fstream PESFile("pes.dat",ios::out);
   fstream ARPESFile("arpes.dat",ios::out);
   
   //We write the data binning in the file "PESFile"
   mySpectrum.writeDataBinning(PESFile,ARPESFile);

/*****************************************************PLOT*****************************************************/

  //Contains methods for drawing curves
  Plot myPESPlot;
  Plot myARPESPlot;

/****************************************We build the key of the plot******************************************/
  myPESPlot.addKey("nField",nFieldBirth);
  myPESPlot.addKey("nVYPerp",nVYPerpBirth);
  myPESPlot.addKey("nVZPrimPerp",nVZPrimPerpBirth);
                
#ifdef MOLECULE
  myPESPlot.addKey("Molecule orientation",myPotential->myOrientation.myString);
  myPESPlot.addKeyVariableArg<double>("charges", 4, myPotential->charge[0],  myPotential->charge[1],  myPotential->charge[2],  myPotential->charge[3]);
  myPESPlot.addKeyVariableArg<double>("bondLength", 3, myPotential->bondLength[1],  myPotential->bondLength[2],  myPotential->bondLength[3]);
#endif

#ifdef HYDROGEN
  myPESPlot.addKey("hydrogen potential");
  myPESPlot.addKey("IP", myPotential->IP);
#endif
  
  myPESPlot.addKey("weightTooSmallNbr", int(double(weightTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  myPESPlot.addKey("stepTooSmallNbr",int(double(stepTooSmallNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
  
    myPESPlot.addKey("trappedElectronNbr", int(double(mySpectrum.trappedElectronNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
    
  if(mySpectrum.angleDetection<180)
    myPESPlot.addKey("angleTooLargeNbr", int(double(mySpectrum.angleTooLargeNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*1000.)/10., "%");
    myPESPlot.addKey("electronsDetectedNbr", double(mySpectrum.electronsDetectedNbr)/(nFieldBirth*nVZPrimPerpBirth*nVYPerpBirth)*100., "%");
    myPESPlot.addKey("binsWidthEnergy",mySpectrum.binsWidthEnergy,"eV");
    myPESPlot.addKey("binsWidthAngle",mySpectrum.binsWidthAngle,"°");
    
  if(mySpectrum.angleDetection<180)
  myPESPlot.addKey("angleDetection",mySpectrum.angleDetection, "deg");
  
  myPESPlot.addKey("weightThresholdRatio",weightThresholdRatio);
  
  #ifdef _OPENMP
  myPESPlot.addKey("threadsNbrMax", threadsNbrMax);
  #endif

  myPESPlot.addKey("ErrorMax",desiredErrorMax);
  myPESPlot.addKey("stepMin",stepMin);
  
  myPESPlot.addKey("ellipticity", myField.ellipticity);
  myPESPlot.addKey("fieldAmplMax",myField.fieldAmpl, "au");
  myPESPlot.addKey("waveLenght",myField.waveLenght*1.E6, "micro-m");
  myPESPlot.addKey("duration",myDisplay.elapsedTime());
  if(myField.ellipticity>0)
  myPESPlot.addKey("laser polarization rotates clockwise when we look along y");  
  if(myField.ellipticity<0)
  myPESPlot.addKey("laser polarization rotates counterclockwise when we look along y"); 
  
  //We put the keys on
  myPESPlot.setKeysOn();
   
  /***********************We add some instructions and indicate which curves we want to plot*******************************/   
   
  /***************************************for the PES*****************************************/ 
  
      
  //We load a file with gnuplot instructions already written
  myPESPlot.setLoadFile("epscairo2d.gnu");
  
  myPESPlot.addInstruction("set multiplot layout 1, 1");
  //myPESPlot.addInstruction("set xrange [0:20]");
  
  //We consider electrons differently depending if they are detected in y>0 or y<0
  std::ostringstream plot1;
  plot1<<"plot 'pes.dat' index 0 using ($1/13.5):(log($2)) w l ls 3 title 'Photo-electrons detected in the upper half-space (y>0), number="<< double(mySpectrum.electronsDetectedUPNbr)/(mySpectrum.electronsDetectedNbr)*100.<<" %', \\";
  std::ostringstream plot2;
  plot2<<"'pes.dat' index 1 using ($1/13.5):2 w l ls 4 title 'Photo-electrons detected in the lower half-space (y<0), number="<< double(mySpectrum.electronsDetectedDOWNNbr)/(mySpectrum.electronsDetectedNbr)*100.<<" %'";

  myPESPlot.addInstruction(plot1.str());
  myPESPlot.addInstruction(plot2.str());
  
  myPESPlot.addInstruction("unset multiplot");

  myPESPlot.gnuplot("pes.gnu","pes.eps");
  
   /***************************************for the ARPES*****************************************/ 
    
    //We load a file with gnuplot instructions already written
    myARPESPlot.setLoadFile("epscairopolar.gnu");
    
    myARPESPlot.addInstruction("set multiplot layout 2, 2");
  
//Margins for each row resp. column
  myARPESPlot.addInstruction("TMARGIN = 'set tmargin at screen 0.90; set bmargin at screen 0.55'");
  myARPESPlot.addInstruction("BMARGIN = 'set tmargin at screen 0.45; set bmargin at screen 0.10'");
  myARPESPlot.addInstruction("LMARGIN = 'set lmargin at screen 0.25; set rmargin at screen 0.65'");
  myARPESPlot.addInstruction("RMARGIN = 'set lmargin at screen 0.55; set rmargin at screen 0.95'");

  myARPESPlot.addInstruction("set multiplot layout 2,2");

  myARPESPlot.addInstruction("set_key(x,text) = sprintf(\"set label %f '%s' at graph 0.5,-0.1 center\", x,text)");

  //Pay attention that no space caractere follows \\ in instruction

  myARPESPlot.addInstruction("eval set_key(1,'(a) PAD for region [0, 2Up]')");
  myARPESPlot.addInstruction("@TMARGIN; @LMARGIN");
  myARPESPlot.addInstruction("plot 'arpes.dat' index 0 using 1:2 w l ls 3 notitle, \\"); 
  myARPESPlot.addInstruction("'arpes.dat' index 0 using (-$1):2 w l ls 3 notitle");
  myARPESPlot.addInstruction("unset label 1");

  myARPESPlot.addInstruction("eval set_key(2,'(b) PAD for region [2Up, 8Up]')");
  myARPESPlot.addInstruction("@TMARGIN; @RMARGIN");
  myARPESPlot.addInstruction("plot 'arpes.dat' index 1 using 1:2 w l ls 3  notitle, \\");
  myARPESPlot.addInstruction("'arpes.dat' index 1 using (-$1):2 w l ls 3  notitle");
  myARPESPlot.addInstruction("unset label 2");

  myARPESPlot.addInstruction("eval set_key(3,'(c) PAD for region [8Up, 10Up]')");
  myARPESPlot.addInstruction("@BMARGIN; @LMARGIN");
  myARPESPlot.addInstruction("plot 'arpes.dat' index 2 using 1:2 w l ls 3  notitle, \\");
  myARPESPlot.addInstruction("'arpes.dat' index 2 using (-$1):2 w l ls 3  notitle");
  myARPESPlot.addInstruction("unset label 3");
  
  myARPESPlot.addInstruction("eval set_key(4,'(b) PAD near 4Up')");
  myARPESPlot.addInstruction("@BMARGIN; @RMARGIN");
  myARPESPlot.addInstruction("plot 'arpes.dat' index 3 using 1:2 w l ls 3  notitle, \\");
  myARPESPlot.addInstruction("'arpes.dat' index 3 using (-$1):2 w l ls 3  notitle");
   
  myARPESPlot.addInstruction("unset multiplot");

  myARPESPlot.gnuplot("arpes.gnu","arpes.eps");
  
  return 0;
}

