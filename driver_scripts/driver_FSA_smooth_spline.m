%% runs force_error script and plots force error
% uses smoothing spline
close all
clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set calculation parameters, such as which radial derivative method to use,can leave as the defaults
use_lambda_full_mesh=true; % interpolate lambda from half mesh to full mesh
%% deriv_method determines which method to use for the radial derivatives of R,Z and Lambda%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% These work fine and should be used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deriv_method='finite difference'; % central differences
% deriv_method='finite difference 4th'; % 4th order central differences
% deriv_method = 'spline' % MATLAB cubic splines
%
% deriv_method = 'spapi' % MATLAB spline of arbitrary order
global SPLINE_ORDER_SPAPI % order of spline for spapi (degree = order -1)
SPLINE_ORDER_SPAPI=3;
%
 deriv_method = 'smooth_spline'
global SMOOTH_FACTOR % controls how smooth spline is, let it be a negative number to let MATLAB decide how smooth the spline should be (should let MATLAB decide)
SMOOTH_FACTOR=0.99999999; % set to -1 to allow matlab to select
%
use_real_space_radial_derivs = false; % if true, use radial derivatives obtained from R,Z after transforming to real space, not radial derivs of the RMNC, ZMNS fourier coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% these do not work well/have implementation issues and should not be used %%%
% deriv_method='poly'

% deriv_method='tension_spline'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if one of these use_XX variables is true, will use that for the radial derivatives
% and will ignore deriv_method
%% Least Squares Fit Piecewise Polynomial
use_piecewise_lsq = false; % least squares fit a piecewise polynomial to the data (will smooth out discontinuities/roughness)
global POLY_LSQ_WINDOW_SIZE
POLY_LSQ_WINDOW_SIZE=16; % number of data points to consider at a time in the polynomial fit
global POLY_LSQ_ORDER
POLY_LSQ_ORDER = 5; % polynomial order to fit the data window with
%% One-sided Cubic Spline Basis
use_my_cubic_spline=false; % use a one-sided spline basis to fit the data
global MY_SPLINE_END_CONDITION
MY_SPLINE_END_CONDITION = 'natural'; %% either 'natural' (2nd deriv=0 at endpoints) or 'not-a-knot' (3rd deriv = 0 at endpoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select VMEC output filename to calculate force error of



file='../example_files/wout_W7X_s256_M12_N12_f12_cpu1_32GB.nc';

SPLINE_ORDER_SPAPI = 3;
  try
    data = read_vmec(file);
    force_error
[val,s_ind] = min(abs(data.phi-0.1));
    plot_force_error

    p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;
   % F_avgg = F_V_avg / p_V_avg
    E = VMEC_W_B + W_p;
    F_rhos = trapz(u,trapz(v,F.*abs_g_vmec,3),2);
    p_rhos = trapz(u,trapz(v,abs(presr).*sqrt(gSS).*abs_g_vmec,3),2);
    F_fsa = F_rhos ./ p_V_avg; % flux surface avg
 catch e
     warning('Error, likely no wout for %s',file);
     fprintf(1,'The identifier was:\n%s',e.identifier);
     fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
 end


figure
plot(s,F_fsa)
set(gca, 'YScale', 'log')