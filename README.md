# PXY_MP_model_solver
This was written by Dr Natasha Savage and edited by Kristine Bagdassaria.

Solver for the model 'A MP-PXY negative feedback loop may stabilize the auxin maxima in the cambium during secondary growth'

This is a MATLAB code used to solve a system of equations to steady state. The system of equations describes the dynamics of auxin, cytokinin, PIN active transport, MP, PXY and TDIF. It is a reaction-diffusion system, thus the reaction and diffusion parts are solved separately. The reaction part is solved using the Euler method, while the diffusion part is solved using a modification to Euler which consists of a diffusion matrix. The output is a table with parameters on the lefthandside, three columns of 999 (for visual separation) and the obtained concentrations at steady state.

The parameters are read from a file. The programme then runs these parameters through the solver and outputs a table of parameters, three columns of 999 for visual separation, and the final steady state concentrations of the relevant components. The programme uses a function that ensures the concentrations reach a steady state. At every step dt, Euler_Function_GitHub.m calculates the reaction part of the equations using the Euler method. This is then
updated with a diffusion matrix Hop. The above calclations are repeated until a steady state is reached. 


Main programme: Numerical_Solver_GitHub.m
Steady state solver function: SS_Simulation_GitHub.m
Euler function: Euler_Function_GitHub.m
Reaction function: Reaction_Part_GitHub.m
