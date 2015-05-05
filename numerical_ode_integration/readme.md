%Numerical ODE integration using OOP
%Yann Gouttenoire

# Motivated interest

  - Implementation of a code for integration of ordinary differentials equations.
  - Using **object-oriented programming** which has the advantages to be modular, 
extensible and reusable. 

# Explanation
## Classes

The main goal of **OOP** is to first define the different elements of your code as __objects__ and second, to make these objects interact together.  
In OOP, each object in an __instances__ of a program-code-template named __classes__.  Each class is designed and programmed to accomplish one particular thing.  
A class contains first, others objects (called __member variables__ or __fields__) and second, procedures called __methods__.  
Member variables are the data of the object and methods are the interface the object presents to the outside world.  
Each class is associated with two files:

  - the header file (\*.h) where we define member variables and methods,
  - the source file (\*.cpp) where we implement the methods.

An advantage of OOP is that a simple reading of the header files (\*.h) can be enough to understand the code.

## Contents of the code

Here are implemented:

  * A class __Ode__ which can contains ordinary differential equations through the method __diff__. Since __diff__ is not implemented (it is a __virtual__ method), __Ode__ is an abstract class, and then, we cannot create any object with it.
  * A class __Coulomb__ which inherits __Ode__, it implements the virtual method __diff__ and defines an ode, so that, the dynamic ode for an electron in a coulomb potential.
  * A class  __Solve__ which contains an __Ode__ object and a method called __rk4__ for the integration of the ode using rk4 iteration.
  * A class __Ic__ which set the initial condition for the variable of ode.  
  * A class __Plot__ which contains a method called __plotIt__ which can plot a data file whith gnuplot.
  * A class __Simulation__ which contains a method called __run__ which start the numerical simulation. 

## Class diagram

The whole structure of the system is described in the __class diagram__.  
It contains one box for each class. A box contains at the top the name of the class, at the middle the members variables and at the bottom the methods the class executes.  
Lines graphically represents the relationships between classes:

  - A line ended with an arrow head means a inheritance relationship.
  - A line ended with a diamond head means a composition relationship. It happens when a class contains on instance of an other class as a member variable.
