# Homework 1
Refer to [Homework 1](http://wiki.cs.umt.edu/classes/cs477/index.php/Homework_1 "homework 1") 

## Problem 1
Find the value of the error function using euler, euler-richardson, and rk4. Using the generated ODE error function. Simulate it and then compare it to the python built in error function. Use the sum square error to compare the functions using the different solvers. 

    \textrm{erf}(x)=\frac{2}{\sqrt \pi} \int_0^{x} e^{-x'^2} dx',

## Problem 2
Simulate early Early Modern Man and Neanderthal populations given limited resources. Construct an ODE and solve ti to determine the best parameters that results in a 10,000 year extinction period for the Neanderthals.

    \frac{dN}{dt} = N[A-D(N+E)-B]

    \frac{dE}{dt} = E[A-D(N+E)-sB]

## Problem 3 Extra Credit
Alter the formulas such that they make the assumption that the Neanderthals were being eaten by the humans.

Assume that interaction is weak (1 x 10^4) and Neanderthal consumption conversion to humans is 10%. 


