#include"plot.h"
using namespace std;


//We add a plot instruction
void Plot::addInstruction(std::string instruction)
{
  instructions.push_back(instruction);
}


//We display the spectra with gnuplot
void Plot::gnuplot(string gnuFileName, string output)
{

  fstream gnuFile(gnuFileName, ios::out);

  gnuFile<<"set output '"<<output<<"'"<<endl;
  //We load a file with some gnuplot intructions already written
  gnuFile<<"load '"<<loadFile<<"'"<<endl;
  

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
  
  std::ostringstream shell;
  shell<<"gnuplot "<<gnuFileName<<" && evince "<<output<<" &";
  system((shell.str()).c_str());

}
