#ifndef DISPLAY_H
#define DISPLAY_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>

//Preprocessor function
//Display x (not the content of x)
#define quote(x) #x

class Display
{

public:

int displayCounter;
std::string elapsedTime;

Display() {displayCounter=0;};

//We implement a method which put the cursor back up in a view to displaying a stationnary output
void moveCursorBackUp();
   
//We implement a load bar with a view to displaying the remaining time
void loadbar(int currentPoint, int finalPoint);

//We implement a method which outputs a data in the terminal
template<typename T>
void operator()(std::string dataName, const T& data, std::string unit="")
{
  std::cout<<dataName<<"= "<<data<<" "<<unit<<"                  "<<std::endl;
  displayCounter++;	
}

void operator()(std::string data)
{
  std::cout<<data<<std::endl;
  displayCounter++;	
}

};

#endif
