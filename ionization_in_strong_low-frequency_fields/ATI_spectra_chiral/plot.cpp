#include"plot.h"
using namespace std;


//We add a plot instruction
void Plot::addInstruction(std::string instruction)
{
  instructions.push_back(instruction);
}


//We display the spectra with gnuplot
void Plot::gnuplot()
{

  fstream gnuFile("data.gnu", ios::out);

  gnuFile<<"set output 'spectrum.eps'"<<endl;

//Line style for axes, grey color
  gnuFile<<"set style line 100 linecolor rgb '#808080'"<<endl;

//Line style for grid, dashed, grey color
  gnuFile<<"set style line 101 linecolor rgb '#808080' linetype 0 "<<endl;

//Line styles for curves
//RED
  gnuFile<<"set style line 1 lc rgb '#A00000' pt 6 ps 1 lt 1 lw 2"<<endl;
//GREEN-BLUE (cyan) (Complementary Color)
  gnuFile<<"set style line 2 lc rgb '#00A0A0' pt 6 ps 1 lt 1 lw 2"<<endl;

//GREEN
  gnuFile<<"set style line 3 lc rgb '#00A000' pt 6 ps 1 lt 1 lw 2"<<endl;
//RED-BLUE (purple) (Complementary Color)
  gnuFile<<"set style line 6 lc rgb '#A000A0' pt 6 ps 1 lt 1 lw 2"<<endl;

//BLUE
  gnuFile<<"set style line 4 lc rgb '#0000A0' pt 6 ps 1 lt 1 lw 2"<<endl; 
//RED-GREEN (gold) (Complementary Color)
  gnuFile<<"set style line 5 lc rgb '#A0A000' pt 6 ps 1 lt 1 lw 2"<<endl;

  vector<string>::iterator it;
  if(isKeysOn==true)
  {
   int size=0;
   it=keys.begin();
  gnuFile<<"set key on outside center bmargin Left reverse box title sprintf(\"";
  for(int k=1; it!=keys.end(); it++, k++)
    {
    size+=it->size();
      if(k==6 || size>120) 
        {
	 gnuFile<<"\\n ";
	 k=1;
	 size=0;
	} 
      if(it==--keys.end())
        gnuFile<<*it;
      else
      gnuFile<<*it<<", ";
    }

  gnuFile<<"\")"<<endl;
  }
  
  
  for(it=instructions.begin(); it!=instructions.end(); it++)
    {
      gnuFile<<*it<<endl;
    }

  gnuFile.close();
  system("gnuplot data.gnu && evince spectrum.eps &");

}
