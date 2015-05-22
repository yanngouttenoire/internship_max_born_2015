#include"plot.h"
#include"parameters.h"
using namespace std;

 //We display the spectra with gnuplot

  //We open a file with a view to writing in it
void Plot::gnuplot(fstream& gnuFile)
{

gnuFile<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
gnuFile<<"set output 'spectrum.eps'"<<endl;

gnuFile<<"set xlabel 'Asymptotic energy (au)'"<<endl;
gnuFile<<"set ylabel 'Probability (log)'"<<endl;

gnuFile<<"set multiplot  layout 1, 1 
gnuFile<<"set key on outside left bmargin box title 'nFieldBirth="<<nFieldBirth<<", nVPerpBirth="<<nVPerpBirth<<", errorMax="<<desiredErrorMax<<", errorMin="<<desiredErrorMin<<", dtMin="<<dtMin<<", dtMinNbr="<<dtMinNbr<<" \n linear field, fieldAmpl="<<fieldAmpl<<" au, waveLenght="<<waveLenght<<", vPerp=0->sigmaV, phase=-45->45 \n duration="<<elapsedTime<<"'"<<endl;
gnuFile<<"plot 'data.dat' using 1:2 lc rgb 'violet' title 'Photo-electron spectrum with rk5 controlled stepper (rk4 embedded)'"<<endl;
gnuFile<<"unset multiplot"<<endl;

gnuFile<<"pause -1"<<endl;
gnuFile<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
gnuFile<<"set output 'courbe.eps'"<<endl;
gnuFile<<"replot"<<endl;


}
