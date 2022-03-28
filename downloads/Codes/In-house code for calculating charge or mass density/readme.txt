This code (mem_density.c) is initially written by Xubo Lin using the template "xdrfile-1.1.4.tar.gz" downloaded from www.gromacs.org, which is modified by Xiu Li for the current purpose. If you use this code for calculating the charge/mass density along the lipid membrane, please cite at least one of the following papers:
1. Xubo Lin*, Alemayehu A. Gorfe*. Transmembrane Potential of Physiologically Relevant Model Membranes: Effects of Membrane Asymmetry. The Journal of Chemical Physics, 2020, 153, 105103. 
2. Xubo Lin, Vinay Nair, Yong Zhou, Alemayehu A. Gorfe*. Membrane Potential and Dynamics in a Ternary Lipid Mixture: Insights from Molecular Dynamics Simulations. Physical Chemistry Chemical Physics, 2018, 20, 15841-15851.

Usage£º
1. gcc mem_density.c xdrfile.c xdrfile_xtc.c -o mem_density -lm
2. ./mem_density num1 num2 char      (num1 and num2 are integer numbers for time with the unit of "ps"; char can be "charge" or "mass")

Code description£º
In order to use this code, a txt file containing atom index (1st column), charge (2nd column) and mass (3rd column) derived from the molecule.itp file needs to be prepared. When running this code over the trajectory (*.xtc), it will fisrtly read the information from these txt files and then calculate the center-of-mass of the lipid membrane and density distribution with the membrane center as the origin for each frame, which is averaged over selected time region.