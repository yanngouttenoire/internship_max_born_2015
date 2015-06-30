#ifndef DISPLAY_H
#define DISPLAY_H
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<iomanip>
#include<vector>
#include<stdarg.h>

//Preprocessor function
//Display x (not the content of x)
#define quote(x) #x

class Display
{

 public:

  int displayCounter;
  time_t start;
  time_t end;

  Display() {displayCounter=0; start=time(NULL);};

  //We implement a method which put the cursor back up in a view to displaying a stationnary output
  void moveCursorBackUp();
   
  //We implement a load bar with a view to displaying the remaining time
  void loadbar(int currentPoint, int finalPoint);

  //We return the elasped time since the creation of the instance of Display
  std::string elapsedTime();


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

  template<typename T>
    void variableArg(std::string dataName, int n,...)
    {
      std::cout<<dataName<<"= ";
      va_list ap;
      va_start(ap,n);
      for(int i=0; i<n; i++) std::cout<<va_arg(ap,T)<<" ";
      va_end(ap);
      std::cout<<" "<<std::endl;
      displayCounter++;
    }


};

#endif
