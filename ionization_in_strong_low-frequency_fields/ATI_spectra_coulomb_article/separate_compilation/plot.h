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

class Plot
{

public:

 //We display the spectra with gnuplot
 void gnuplot(std::fstream& gnuFile);


};

#endif
