#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string>
#include"simulation.h"

using namespace std;

Simulation::Simulation(Ode* a_ode)
{

//We create the array in which we will store the orbit
q=new double[6];
t=new double; 

//Starting from the argument of the constructor, we initialise the object which contains the dynamic equations of the motion
m_ode=a_ode;

//We allocate an instance of IC: handles the initial conditions
m_ic=new Ic(q);

//We allocate an instance of solve: contains methods which can solve the differential system of equation "m_ode"
m_solve=new Solve(q, t, m_ode);

//We allocate an instance of Plot for ploting action
m_plot=new Plot();

}

Simulation::~Simulation()
{

//We delete all the pointers which we allocated in the constructor
delete[] q;
delete t;
delete m_solve;
delete m_plot;

}

void Simulation::setOde(Ode* a_ode)
{
m_ode=a_ode;
}

void Simulation::run(string algorithm, double step, unsigned int ptsNbr)
{

//We open a file in a view to performing output operation
fstream dat("data.dat",ios::out);
 
//We call the method of IC which fixes the initial conditions
m_ic->setIc(0.,0.,1.,0.,1.,0., t);

for(unsigned int n=1; n<=ptsNbr; n++)
{ 

//We call the method rk4 of Solve with a step dt
if(algorithm=="rk4")
m_solve->rk4(step);
else 
cout<<"The algorithm you want to use is not available"<<endl;

//We write coordinates of the orbit in a file
for (unsigned int i=0; i<6; i++)
{
dat<<q[i]<<" ";
}
dat<<" "<<endl;
}

//We call the method of Plot which displays the orbit with gnuplot
m_plot->plotIt("data.dat", "Coulomb orbit");
}

