#Description
This repository gathers together my work as intern at the Max-Born Institute. I am representated by M.I (Max-Born Institute, Berlin) and E.C (ISMO, Orsay).
The subject deals with the recollion process. It is a process in 3 step:
  
  -The ionization of the electron with an intense electric field. In function of the field parameters (Intensity, pulsation, initial phase), the electron reaches the continuum of energy after a multiphoton absorption or after he tunnels through the classically forbidden region caused by the coulomb potential. 
  -Then, the electron propagates and oscillates in the field. We can do the Strong Field Approximation which consists to consider that if the field is sufficiently intense, we can ignore the ionic potential and consider the electron behaves as a plane waves.
  -In function of the field parameters (Intensity, pulsation, initial phase),
the electron can scatter the ionic core and then can ionize a second electron, recombinate, of just passes through.

#Accomplished work
 Implementation of a general code for computing solutions of ordinary differential equations system with Runge Kutta 4 iteration in C++ and in Matlab.


#Objective
 To compute classical trajectories of an electron in a Coulomb potential and in an intense electric field. Then, to study the influence of a weak ellipticity compared to an linear one. I already dispose of a code which computes the trajectories of the electron and determines which are the conditions on the initial phase of the field for having the closest approach to the nucleus. The orbit starts just after the ionization process. The initial position is the position of the electron just after he finishes tunneling. The initial parallel (to the field) velocity is taken null. The initial perpendicular velocity follows a velocity distribution (ADK distribution ?) centered on zero. 
I will try to replace the linear field by a weakly elliptic one. 
Then, we could try to replace the Coulomb hydrogenoid potential by the potential of a chiral molecule. 