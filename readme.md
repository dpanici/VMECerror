Set of functions and script used to calculate VMEC fixed-boundary equilibrium force error for [Panici 2022](https://arxiv.org/abs/2203.17173)

Add this folder (`VMECerror`) and the subfolders `VMECerror/plotting` and `VMECerror/fxns` to your MATLAB path in order to run the scripts here and in the Zenodo files for the Panici 2022 paper.

**Must also have [matlabVMEC](https://github.com/lazersos/matlabVMEC) cloned locally and added on the MATLAB path**

-------------------------------------------------------------------------
force_error.m is the main script which will read in the data from the VMEC
output files (wout) and get the R,Z,lambda fourier series coefficients,
and calculate the necessary B,J and derivatives for findinf the force
error.


debug_plot_quants.m plots a bunch of quantities and compares them to what 
matlabVMEC calculates
	Can use VMECplot command to plot quantities that readVMEC does not
	give, but matlabVMEC has internally (such as g)
