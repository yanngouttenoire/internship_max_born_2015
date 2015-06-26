#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>

using namespace std;

double fieldAmpl=0.0534;
double IP=0.5792;
double vPerpMax=0.3;
double cyclesNbr=2.;
double pulsation=0.0228;
double opticalCycle=2.*M_PI/pulsation;
double phase=0.;

double field(double t)
{
return fabs(exp(-2.*log(2.)*pow(t/cyclesNbr/opticalCycle,2))*fieldAmpl*cos(pulsation*t+phase));
}

double weight(double t, double vPerp, int i)
{
if(i==1)
return 4./field(t)*exp(-2.*pow(2.*IP,3./2.)/3./field(t))*fabs(vPerp)/field(t)*exp(-pow(2.*IP,1./2.)*vPerp*vPerp/field(t));
if(i==2)
return 1./field(t)/field(t)*exp(-2.*pow((2.*IP+vPerp*vPerp),3./2.)/3/field(t))*vPerp/sqrt(1+vPerp*vPerp/2./IP);
}

int main()

{
int ni=100;
int nj=100;
double fieldBirth;
double vPerp;
double max1=0.;
double max2=0.;
double ti;

  //We open a file with a view to writing in it
  fstream dat("data.dat", ios::out);

for(int i=0; i<ni; i++)
{
for(int j=0; j<nj; j++)
{
 ti=-M_PI/pulsation/90.*25+i/double(ni)*2.*M_PI/pulsation/90.*25;
 vPerp=0.01+j*vPerpMax/nj;

 dat<<field(ti)<<" "<<vPerp<<" "<<weight(ti, vPerp,1)<<endl;

if(weight(ti, vPerp,1)>max1)
max1=weight(ti, vPerp,1);

if(weight(ti, vPerp,2)>max2)
max2=weight(ti, vPerp,2);
}
}

cout<<"max1= "<<max1<<" "<<weight(0.,sqrt(field(0)/2.),1)<<endl;
cout<<"max2= "<<max2<<" "<<weight(0.,sqrt(field(0)/2.),2)<<endl;

cout<<"t=0 "<<scientific<<weight(0.,sqrt(field(0)/2.),1)<<" "<<weight(0.,sqrt(field(0)/2.),2)<<endl;
cout<<"t=20° "<<scientific<<weight(-M_PI/pulsation/90.*20.,sqrt(field(-M_PI/pulsation/90.*20.)/2.),1)<<" "<<weight(-M_PI/pulsation/90.*20.,sqrt(field(-M_PI/pulsation/90.*20.)/2.),2)<<endl;
cout<<"t=25° "<<scientific<<weight(-M_PI/pulsation/90.*25.,sqrt(field(-M_PI/pulsation/90.*25.)/2.),1)<<" "<<weight(-M_PI/pulsation/90.*25.,sqrt(field(-M_PI/pulsation/90.*25.)/2.),2)<<endl;
cout<<"t=-T/2 "<<weight(-M_PI/pulsation,sqrt(field(-M_PI/pulsation)/2.),1)<<" "<<weight(-M_PI/pulsation,sqrt(field(-M_PI/pulsation)/2.),2)<<endl;


cout<<"ratio_20° "<<fixed<<weight(0.,sqrt(field(0)/2.),1)/weight(-M_PI/pulsation/90.*20.,sqrt(field(-M_PI/pulsation/90.*20.)/2.),1)<<" "<<weight(0.,sqrt(field(0)/2.),2)/weight(-M_PI/pulsation/90.*20.,sqrt(field(-M_PI/pulsation/90.*20.)/2.),2)<<endl;


cout<<"ratio_T/2° "<<fixed<<weight(0.,sqrt(field(0)/2.),1)/weight(-M_PI/pulsation,sqrt(field(-M_PI/pulsation)/2.),1)<<" "<<weight(0.,sqrt(field(0)/2.),2)/weight(-M_PI/pulsation,sqrt(field(-M_PI/pulsation)/2.),2)<<endl;

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

 // system("gnuplot data.gnu");

  return 0;
}

