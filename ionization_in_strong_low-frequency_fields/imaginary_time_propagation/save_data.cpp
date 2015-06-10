#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>


using namespace std;

int main()
{

  fstream gnu("data.gnu", ios::out);

  gnu<<"set xrange [-2:2]"<<endl;
  gnu<<"set yrange [-2:2]"<<endl;
  gnu<<"set multiplot layout 2, 1 "<<endl;
  gnu<<"set key on outside left bmargin box "<<endl;
  gnu<<"splot 'ground_state_wf.dat' index 0 using 1:2:4 w l title 'wave_function'"<<endl;
  gnu<<"unset multiplot"<<endl;
  gnu<<"pause -1"<<endl;
  gnu<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
  gnu<<"set output 'courbe.eps'"<<endl;
  gnu<<"replot"<<endl;
  gnu.close();
  system("gnuplot data.gnu");

}
