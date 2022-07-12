%% Example use of VMECerror matlab scripts
% calculate and plot the force error residuals of a W7-X VMEC
% fixed-boundary equilibrium

close all
clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set calculation parameters, such as which radial derivative method to use,can leave as the defaults
use_lambda_full_mesh=true; % interpolate lambda from half mesh to full mesh
%% deriv_method determines which method to use for the radial derivatives of R,Z and Lambda%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% These work fine and should be used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv_method='finite difference'; % central differences
% deriv_method='finite difference 4th'; % 4th order central differences
% deriv_method = 'spline' % MATLAB cubic splines
%
% deriv_method = 'spapi' % MATLAB spline of arbitrary order
global SPLINE_ORDER_SPAPI % order of spline for spapi (degree = order -1)
SPLINE_ORDER_SPAPI=6;
%
% deriv_method = 'smooth_spline'
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


file='./example_files/wout_W7X_s256_M12_N12_f12_cpu1_32GB.nc';


%% read VMEC file using matlabVMEC read_VMEC function
data = read_vmec(file);
%% run force_error script to calculate force error
force_error

%% calculate energy and Flux Surface Avg of Force
E = W_B + W_p;
F_rhos = trapz(u,trapz(v,F.*abs_g_vmec,3),2) ./ trapz(u,trapz(v,abs_g_vmec,3),2);
vol = trapz(v,trapz(u,trapz(s,abs_g_vmec)));

% choose normalization for FSA based on if pressure or vacuum
p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ vol;
if max(data.presf) > 1e-3 % normalize by pressure gradient
    normalize = p_V_avg
else % normalize by magnetic pressure gradient
    grad_B_pres_s = (BU.*Bu_s + Bu.*BU_s + Bv.*BV_s + BV.*Bv_s)*0.5; % contravariant s component of the magnetic pressure gradient
    grad_B_pres_u = (BU.*Bu_u + Bu.*BU_u + Bv.*BV_u + BV.*Bv_u).*0.5; % contravariant u component of the magnetic pressure gradient
    grad_B_pres_v = (BU.*Bu_v + Bu.*BU_v + Bv.*BV_v + BV.*Bv_v).*0.5; % contravariant v component of the magnetic pressure gradient
    grad_B_pres = (grad_B_pres_s .* eS + grad_B_pres_u .* eU + grad_B_pres_v .* eV);

    mag_grad_B_pres = sqrt((grad_B_pres_s .* dot(grad_B_pres,eS,4) + grad_B_pres_u .* dot(grad_B_pres,eU,4) + grad_B_pres_v .* dot(grad_B_pres,eV,4) ) );
    normalize = mag_grad_B_pres
end
F_fsa = F_rhos ./ p_V_avg; % flux surface avg normalized by volume average of pressure gradient


% plot force error and flux surfaces
plot_force_error

figure
plot(s(1:data.ns),F_fsa(1:data.ns))
set(gca, 'YScale', 'log')

ylabel('Normalized Force Error FSA')
xlabel('s')
