%% runs force_error script and plots force error
close all
clearvars

u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

% file = 'VMECfiles/wout_DSHAPE_Vac_s512_M15.nc';
file = 'VMECfiles/wout_HELIOTRON_s256_M12_N3.nc';
% file = 'VMECfiles/wout_HELIOTRON_s512_M12_N3.nc'; % This HELIOTRON case agrees
% very well everywhere except near axis

% file = 'VMECfiles/wout_HELIOTRON_Vac.nc';
% file = 'VMECfiles/wout_SOLOVEV.nc';
% file = 'VMECfiles/wout_W7X_standard.nc';
% file = 'VMECfiles/wout_W7X_standard_M15_N12_s256.nc'

% load in VMEC data then run force error to calculate everything
data = read_vmec(file);
force_error
% plot force_error will do as its name implies
plot_force_error

% debug_plot_quants script has a lot of plotting scripts for checking
% intermediate quantities against VMEC
debug_plot_quants

