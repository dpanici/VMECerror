%% runs force_error script and plots force error
close all
clearvars

u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at
% deriv_method='finite difference';
deriv_method='spapi'
F_avgs = zeros([1,length(3:15)])
global SPLINE_ORDER_SPAPI
for i = 3:15
SPLINE_ORDER_SPAPI = i
% file = 'VMECfiles/wout_DSHAPE_Vac_s512_M15.nc';
% file = 'VMECfiles/wout_HELIOTRON_s256_M12_N3.nc';
% file = 'VMECfiles/wout_HELIOTRON_s512_M12_N3.nc'; % This HELIOTRON case agrees
% very well everywhere except near axis

% file = 'VMECfiles/wout_HELIOTRON_Vac.nc';
% file = 'VMECfiles/wout_SOLOVEV.nc';
% file = 'VMECfiles/wout_W7X_standard.nc';
% file = 'VMECfiles/wout_W7X_standard_M15_N12_s256.nc'

file = 'VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu1.nc';
% load in VMEC data then run force error to calculate everything
data = read_vmec(file);
force_error
% plot force_error will do as its name implies
plot_force_error

[val,s_ind] = min(abs(data.phi-0.1));


F_V_avg = trapz(v,trapz(u,trapz(s(s_ind:end),F(s_ind:end,:,:).*abs_g_vmec(s_ind:end,:,:)))) ./ data.Volume;
p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;
F_avgg = F_V_avg / p_V_avg
E = VMEC_W_B + W_p;
F_avgs(SPLINE_ORDER_SPAPI-2) = F_avgg
% debug_plot_quants script has a lot of plotting scripts for checking
% intermediate quantities against VMEC
% debug_plot_quants
end
figure
plot(3:15,F_avgs)
