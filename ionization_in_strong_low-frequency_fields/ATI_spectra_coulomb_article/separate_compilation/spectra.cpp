#include"spectra.h"
#include"parameters.h"

using namespace std;

//We set the histogram intervals width
template<typename state_type>
Spectra<state_type>::Spectra(double binsWidth) : binsWidth(binsWidth) {};


//We store asymptotic energies in containers of map type with a view to make a data binning
template<typename state_type>
void Spectra<state_type>::storeDataBinning(const state_type& x, const double& weightIonization, const bool& inexpectedStop)
{

  int range;
  
  //We compute the x value of the new point in the histogram
  range=int(asymptoticEnergy(x)/binsWidth);
  
  if(inexpectedStop) return;

  //We declare a container of map type
  //It can contains pairs which are couple of objects
  //First the range, second the asymptic velocity value

  map<int,double>::iterator findRange= asymptEnergy.find(range);
     
  spectraPointsNbr++;
    
  if(findRange==asymptEnergy.end())
    {
      asymptEnergy[range]=weightIonization;
    }
  else
    {
      asymptEnergy[range]=asymptEnergy[range]+weightIonization;
    }
}



  //Finally we write all the data binning in a file
template<typename state_type>
void Spectra<state_type>::writeDataBinning(fstream& dataFile)
{
  
  map<int,double>::iterator it= asymptEnergy.begin();

  //We normalize the distribution law
  double sum=0.;
  for(1; it!=asymptEnergy.end(); it++)
    {
          sum=sum+it->second;
    }

  //We write the histogram in a file	
  for(it=asymptEnergy.begin(); it!=asymptEnergy.end(); it++)
    {
      dataFile<<(it->first)*binsWidth<<" "<<log10( (it->second)/sum )<<endl;
    }

}


//We compute asymptotic energy
template<typename state_type>
double Spectra<state_type>::asymptoticEnergy(const state_type &x)
{
  double vectPot=fieldAmpl/pulsation*sin(pulsation*t+phase);
  double Vsq=x[3]*x[3]+x[4]*x[4]+(x[5]+vectPot)*(x[5]+vectPot);
  double dist=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+softParameter*softParameter);
  double E=Vsq/2-Z/dist;

  return E;

}


template class Spectra <vector<double> >;
template class Spectra <double* >;
//For more details about templates, see https://www.cs.umd.edu/class/fall2002/cmsc214/Projects/P2/proj2.temp.html
