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

int displayCounter;

public:

Display() {displayCounter=0;};

//We implement a method which put the cursor back up in a view to displaying a stationnary output
void moveCursorBackUp();
   
//We implement a load bar with a view to displaying the remaining time
void loadbar(int currentPoint, int finalPoint);

//We implement a method which outputs a data in the terminal
template<typename T>
void display(T &data)
{
  std::cout<<quote(data)<<"= "<<data<<std::endl;
  displayCounter++;	
}


};

#endif
