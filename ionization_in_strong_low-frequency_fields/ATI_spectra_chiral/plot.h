#ifndef PLOT_H
#define PLOT_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<string>
#include<sstream>


//Preprocessor function
//Display x (not the content of x)
#define quote(x) #x

class Plot
{

public:
std::vector<std::string> keys;
std::vector<std::string> plot;
std::vector<std::string> instructions;
std::string plotType;

template<typename T>
void addKey(std::string keyName, const T &key, std::string unit="")
{
std::ostringstream ostringKey;
ostringKey<<keyName<<" = "<<key<<unit;
keys.push_back(ostringKey.str());
}

void addKey(std::string key)
{
keys.push_back(key);
}
  
  //We choose the type of plot
void setPlotType(std::string plotType);

  //We add a plot instruction
void addInstruction(std::string instruction);

 //We add a plot 
void addPlot(std::string instruction);

 //We display the spectra with gnuplot
 void gnuplot();


};

#endif
