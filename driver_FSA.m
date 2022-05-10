%% runs force_error script and plots force error
close all
clearvars -except F_fsa2
% F_fsa2 = F_fsa;

%% deriv_method determines which method to use for the radial derivatives of R,Z and Lambda%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% These work fine and should be used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv_method='finite difference'; % central differences
% deriv_method='finite difference 4th'; % 4th order central differences
% deriv_method = 'spline' % MATLAB cubic splines
% deriv_method = 'spapi' % MATLAB spline of arbitrary order
global SPLINE_ORDER_SPAPI % order of spline for spapi (degree = order -1)
SPLINE_ORDER_SPAPI=6
% deriv_method = 'smooth_spline'
global SMOOTH_FACTOR % controls how smooth spline is, let it be a negative number to let MATLAB decide how smooth the spline should be (should let MATLAB decide)
SMOOTH_FACTOR=0.99999999; % set to -1 to allow matlab to select

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

% file = 'VMECfiles/wout_W7X_s512_M16_N16_f12_cpu1.nc'
% file='VMECfiles/wout_W7X_s500_M16_N16.nc'
% file='VMECfiles/wout_W7X_s729_M16_N16_f12.nc'
% file='VMECfiles/wout_W7X_s10_1024_M16_N16.nc';
% file='wout_W7X_s500_M16_N16_f12_cpu1_32GB_not_pow2.nc';
% file='VMECfiles/wout_W7X_s1000_M16_N16_f12_cpu1_32GB_not_pow2.nc';
% file = 'VMECfiles/wout_VMEC_initial_12.nc';
% file='VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu8_diff_ftol.nc';
figure
for MN=16:16
file=sprintf('VMECfiles/wout_W7X_s512_M%d_N%d_f12_cpu1.nc',MN,MN);
%% read VMEC file 
data = read_vmec(file);
%% run force_error script to calculate force error
force_error

% plot_force_error %uncomment if want to plot Force cross-sectional plot

%% calculate energy and Flux Surface Avg of Force
E = W_B + W_p;
F_rhos = trapz(u,trapz(v,F.*abs_g_vmec,3),2);
p_rhos = trapz(u,trapz(v,abs(presr).*sqrt(gSS).*abs_g_vmec,3),2);
F_fsa = F_rhos ./ p_rhos; % flux surface avg
p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;

% dlmwrite('VMECfiles/W7X_s1000_higher_ftol_F_fsa.txt',F_fsa);

% dlmwrite('VMECfiles/wout_initial_F.txt',F(:,:,nfp_v_index));
% dlmwrite('VMECfiles/wout_initial_R.txt',R(:,:,nfp_v_index));
% dlmwrite('VMECfiles/wout_initial_Z.txt',Z(:,:,nfp_v_index));
% 
% dlmwrite('VMECfiles/wout_initial_F_norm.txt',F(:,:,nfp_v_index)./p_V_avg);


plot_label=sprintf('MN=%d',MN);
% figure
plot(s(1:data.ns),F_fsa(1:data.ns),'DisplayName',plot_label)
set(gca, 'YScale', 'log')

ylabel('Normalized Force Error FSA')
xlabel('s')
hold on
end
legend
% plot(s,F_fsa2(1:end),'DisplayName','higher ftol')
% a1 = gca
% for ns = [3,9,27,81]
% 
%     f2 = figure
%     a2 = copyobj(a1,f2)
%    for loc = linspace(0,1,ns)
%    if loc == 0
%     xline(loc,'--','DisplayName',sprintf('ns=%d Surfaces',ns))
%     hold on
%     else
%     xline(loc,'--','HandleVisibility','off')
%     hold on
%     end
%    end
% legend
% end


% figure
% title('Iota Profile')
% plot(data.phi,data.iotaf)
% ylabel('Iota')
% xlabel('s')