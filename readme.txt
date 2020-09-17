VMEC Error Analysis

force_error.m is the main script which will read in the data from the VMEC
output files (wout) and get the R,Z,lambda fourier series coefficients,
and calculate the necessary B,J and derivatives for findinf the force
error.

plot_force_error.m plots the force error in the volume at the v=0 toroidal plane
(i.e. plots it in R-Z at phi=0)

debug_plot_quants.m plots a bunch of quantities and compares them to what 
matlabVMEC calculates
	Can use VMECplot command to plot quantities that readVMEC does not
	give, but matlabVMEC has internally (such as g)
