#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string>
#include"simulation.h"

using namespace std;


int main()
{
//We declare an instance of Coulomb which contains ODE of Coulomb potential
//Argument:1 -> atomics units
Coulomb m_coulomb(1);
//We declare an instance of Simulation which contains handles the numerics integration of the electron in a Coulomb potential 
Simulation m_simulation(&m_coulomb);
m_simulation.run("rk4",0.1, 100);

return 0;
}

