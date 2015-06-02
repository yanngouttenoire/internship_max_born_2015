#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
using namespace std;

int main()
{

  double dt, N;

  cout<<"ptNumber ?"<<endl;
  cin>>N;
  cout<<"dt ?"<<endl;
  cin>>dt;
  
  fstream gnu("data.gnu", ios::out);

  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,15'"<<endl;
  gnu<<"set output 'spectra_"<<"ptsNumber="<<N<<"_dt="<<dt<<"t="<<N*dt<<"ua"<<".eps'"<<endl;
  gnu<<"set xtics rotate out"<<endl;
  gnu<<"set key on outside left bmargin box title 'ptsNumber="<<N<<"_dt="<<dt<<"t="<<N*dt<<"ua'"<<endl;
  gnu<<"plot 'data.dat' index 0 using 1:2 w l "<<endl;
  gnu<<"unset output"<<endl;
 
  
  gnu.close();
  system("gnuplot data.gnu");

}
