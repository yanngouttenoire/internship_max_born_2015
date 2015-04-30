#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string>
#include"plot.h" 
using namespace std;

void Plot::plotIt(string file, string title)
{
  //system("mkdir -p data"); //Option -p: no error if directory already exists
fstream gnu("data.gnu", ios::out);

gnu<<"splot '"<<file<<"'  u 1:2:3 w l title '"<<title<<"'"<<endl;

gnu<<"pause -1"<<endl;
gnu<<"set term postscript"<<endl;
gnu<<"set output 'courbe.ps'"<<endl;
gnu<<"replot"<<endl;
gnu.close();

system("gnuplot data.gnu");
//system("rm data.dat data.gnu");
}
