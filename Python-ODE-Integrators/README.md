# Python ODE Integrators
This assignment introduce python ODE solvers. Using the new tool, simulate two planets orbiting a sun. 

## Objectives
Accomplish the following by Fri. Feb. 14, 2014. 

  1. Write the equations of motion for a two planet system
  2. Transform those equations of motion into a function that can be handled by the ODE solving machinery.
  3. Set the masses and initial conditions such that the system falls into stable orbits. 
  4. Have facilities that allow you to reason about the output.
  5. Write the equation for the energy of the system.
  6. This is the basis of the n-body problem . The simple loops over two planets here can be generalized to loops over many bodies. 

## Equations of Motion

[Equations](http://wiki.cs.umt.edu/classes/cs477/index.php/Python_ODE_Integrators "Equations")

## Assignment Tasks
To receive full credit for the assignment complete the following tasks. 
  1. Begin by simulating the system with the Runge-Kutta 4th order method that you wrote. 
  Note the time step, and the total time required for simulation. Compare this to what is used for subsequent steps. 
  2. As a group, modify the templates on the [Integrating ODEs](http://wiki.cs.umt.edu/classes/cs477/index.php/Integrating_ODEs "Integrating ODEs") to simulate this two-planet system and create a wiki page to post your source code, images of the output, and commentary on these questions. 
  3. Address the conservation of energy and momentum as specifically as possible with detailed plots. 
  4. Investigate the relation between the atol and rtol parameters, and the conservation. 
  5. Another interesting dynamical system is a planet orbiting about two fixed stars of equal mass. 
  In this case there are no closed orbits, but the orbits can be classified as either stable or unstable. 
  Stable orbits may be open loops that encircle both stars, figure eights, or Kepler-like orbits the encircle only one star. 
  Unstable orbits will eventually collide with one of the stars. 
  Modify the program written for two planets to simulate the double star system, with the first star located at the origin and the second star of equal mass located at (2,0). 
  Place the planet at (1.1,1) and systematically vary the x and y components of the velocity to obtain different types of orbits. 
  Try other initial conditions. 
  Provide plots of each unique behavior that you discover. 

