#include"plot.h"
using namespace std;



 //We set the plot type
void Plot::setPlotType(std::string m_plotType)
{
plotType=m_plotType;
}

 //We add a plot 
void Plot::addPlot(std::string instruction)
{
plot.push_back(instruction);
}

//We display the spectra with gnuplot
void Plot::gnuplot()
{

fstream gnuFile("data.gnu", ios::out);

gnuFile<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
gnuFile<<"set output 'spectrum.eps'"<<endl;

gnuFile<<"set xlabel 'Asymptotic energy (au)'"<<endl;
gnuFile<<"set ylabel 'Probability (log)'"<<endl;

gnuFile<<"set multiplot  layout 1, 1"<<endl;


gnuFile<<"set key on outside left bmargin box title sprintf(\"";
vector<string>::iterator it=legend.begin();
for(int k=1; it!=legend.end(); it++, k++)
{
gnuFile<<*it<<", ";
if(k==5) 
gnuFile<<"\\n ";
}
gnuFile<<"\")"<<endl;


it = plot.begin();
gnuFile<<plotType<<" ";
for(int k=0; it!=plot.end(); it++, k++)
{
if(k!=0)
gnuFile<<", \\"<<endl;
gnuFile<<*it;
}
gnuFile<<" "<<endl;


gnuFile<<"unset multiplot"<<endl;

gnuFile<<"pause -1"<<endl;
gnuFile<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<endl;
gnuFile<<"set output 'courbe.eps'"<<endl;
gnuFile<<"replot"<<endl;

    
gnuFile.close();
system("gnuplot data.gnu");

}
