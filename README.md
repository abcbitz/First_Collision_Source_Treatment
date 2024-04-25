
Surface Method Raytracer (SMR) , v. 1.0
Author: Alex Christensen
Date: April 2024
LLNL-CODE-862911

This work was completed as part of a PhD project to examine a new method for solving the uncollided flux. 

Credit to Greg von Winckel for lglnodes and lgwt functions and Manuel A. Diaz for his GenerateCentralDiffsCoefs function.

-------------------------------------------------

This code is a MATLAB project used to solve for the uncollided neutron flux using ray tracing techniques. It contains two different approaches - the interior point method and the surface method. 

Contained are two folders for two different example problems - one a square pin of highly absorbing material surrounded by a low absorbing material, while the other uses a functionally defined cross section distribution. 

To run the interior point method, run the main_int.m file. 

To run the surface method, run the main_surf.m file. 

To run a series of tests and compare the results of the two different methods, use the repeated.m file. 

To create different test problems, modify the Functions/Case_Functions/xs_generator.m file.

If using a functionally defined cross section and want to directly integrate the function rather than trace through the mesh, set use_func = 1 in the main


