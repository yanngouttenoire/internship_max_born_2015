#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>

using namespace std;

int main()

{
int ni=100;
int nj=100;
double fieldBirthMax=0.06;
double IP=0.5;
double vPerpMax=0.3;
double weight;
double fieldBirth;
double vPerp;

  //We open a file with a view to writing in it
  fstream dat("data.dat", ios::out);

for(int i=0; i<ni; i++)
{
for(int j=0; j<nj; j++)
{
 fieldBirth=i*fieldBirthMax/ni;
 vPerp=j*vPerpMax/nj;
 weight=4./fieldBirth*exp(-2.*pow(2.*IP,3./2.)/3./fieldBirth)*fabs(vPerp)/fieldBirth*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/fieldBirth);
 dat<<fieldBirth<<" "<<vPerp<<" "<<weight<<endl;
}
}

  //We open a file with a view to writing in it
  fstream gnu("data.gnu", ios::out);

  gnu<<"set xlabel 'fieldBirth'"<<endl;
  gnu<<"set ylabel 'vperp'"<<endl;
  gnu<<"splot 'data.dat' using 1:2:3 w l title 'weightIonization'"<<endl;
  gnu<<"unset multiplot"<<endl;

  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
  gnu<<"set output 'spectra.eps'"<<endl;
  gnu<<"replot"<<endl;
  gnu.close();

  system("gnuplot data.gnu");

  return 0;
}

