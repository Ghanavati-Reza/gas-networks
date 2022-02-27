# Dynamic Gas Network Model 

In this repository all the files are writen in MATLAB

Gas Network is defined as ODEs and solved using ode solver

We can define any pipeline network that contains pipes, nodes, supplies, and deliveries.

Also, it is possible to define multiple events such as increase or decrease the pressure in each supply or mass flow in each delivery 

The definition of files is as follows:

- examle.m
used to define the structure of the gas network and the characteristics of the fluid

- ode-solver.m
used to find number of state variables, initialize the states, and define the ode solver

- ode-fun.me
used to define ODEs that will be solved by the ode solver

- plot_results.m
used to show pressure and mass flow profiles of each pipeline in the network
