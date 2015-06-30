#include<sstream>

#include"display.h"


using namespace std;

//We implement a method which put the cursor back up in a view to displaying a stationnary output
void Display::moveCursorBackUp()
{

  while(displayCounter!=0)
    {
      cout<<"\033[F";
      --displayCounter;
    }

}
       
//We implement a load bar with a view to displaying the remaining time

void Display::loadbar(int i, int np)
{

  //1//We build the progress bar in pourcent
  double ratio =(double)(i)/(double)(np);
  double progress=ratio*(100); 
 
  int s=cout.precision();

  cout<<"progression= "<<setprecision(3)<<progress<<"%"<<setprecision(s)<<"                "<<endl;
  displayCounter++;        

  //2//Elapsed time since the beginning of the simulation
  cout<<"elapsed time= "<<elapsedTime()<<"             "<<endl;
  displayCounter++;

  //3//Remaining time since the beginning of the simulation
  time_t end = time (NULL); 
  double T0=difftime(end,start);
  double Z0=(100-progress)*T0/progress;
  int Z1=(int)(Z0/3600);
  int Z2=(int)(Z0/60)-Z1*60;
  int Z3=(int)(Z0)-Z1*3600-Z2*60; 
  
  if(fabs(Z1)>1000 || fabs(Z2)>60 || fabs(Z3)>60)
    {
      Z1=0;
      Z2=0;
      Z3=0;
    }

  cout<<"estimated remaining time= "<<Z1<<" h "<<Z2<<" min "<<Z3<<" s "<<"    "<<endl;
  displayCounter++;


}

//We return a string variable for the elapsed Time
std::string Display::elapsedTime()
{
  time_t end = time (NULL); 
  double T0=difftime(end,start);

  int T1=(int)(T0/3600);
  int T2=(int)(T0/60)-T1*60;
  int T3=(int)(T0)-T1*3600-T2*60;

  ostringstream ostringElapsedTime;
  ostringElapsedTime <<T1<<" h "<<T2<<" min "<<T3<<" s ";
  return ostringElapsedTime.str();
}



