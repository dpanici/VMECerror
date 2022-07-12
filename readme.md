


Set of functions and script used to calculate VMEC fixed-boundary equilibrium force error for [Panici 2022](https://arxiv.org/abs/2203.17173)

Installation
____________

- clone this repository

- Add this folder (`VMECerror`) and the subfolders `VMECerror/plotting` and `VMECerror/fxns` to your [MATLAB path](https://www.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html) in order to run the scripts here and in the Zenodo files for the Panici 2022 paper.

**Must also have [matlabVMEC](https://github.com/lazersos/matlabVMEC) cloned locally and added on the MATLAB path**

-------------------------------------------------------------------------
force_error.m is the main script which will read in the data from the VMEC
output files (wout) and get the R,Z,lambda fourier series coefficients,
and calculate the necessary B,J and derivatives for findinf the force
error.


debug_plot_quants.m plots quantities and compares them to what 
matlabVMEC calculates. Uncomment lines to compare different quantities.
In these plots, "VMEC" labelled quantities are quantitites found using the
fourier series coefficients from the VMEC wout file for that quantity (like
bsupumnc, jcurrumnc for B^u and J^u), while the quantities labelled "my XX" are
found from starting with just R,Z,and lamba calculated from the VMEC wout file
and then multiplying them together in real space.

example.m shows using the scripts to calculate force error FSA for an example
VMEC equilibrium
