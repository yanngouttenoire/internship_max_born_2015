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
std::vector<std::string> legend;

template<typename T>
void addKey(std::string keyName, const T &key, std::string unit="")
{
std::ostringstream ostringKey;
ostringKey<<keyName<<" = "<<key<<unit;
legend.push_back(ostringKey.str());
}

void addKey(std::string key)
{
legend.push_back(key);
}

 //We display the spectra with gnuplot
 void gnuplot();


};

#endif
