%% runs force_error script and plots force error
close all
clearvars
% F_FSA_fac = F_fsa
u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

path_to_W7X = '/p/desc-usr/dpanici/vmec/W7X_pressure/';
% load in VMEC data then run force error to calculate everything

% plot force_error will do as its name implies
% plot_force_error
% deriv_method='factor difference';
deriv_method='finite difference';

% deriv_method = 'spline'

% deriv_method = 'smooth_spline'
% deriv_method='tension_spline'
Ns = [1 2 3];
Ms = [8 10 12 14 16];
%Ms = [14 16]
%ftols = [4 8 12];
ftols = [4 8 12];

s_nums = [128 256 512 1024 2048]; %no s1024 needed for convergence sweep over M,N

cpus = [1 2 4 8];
iiii = 1;
force_filename = 'F_spline_FSA_norm_by_p_rho_avg.txt';
global SPLINE_ORDER_SPAPI
global SPLINE_TENSION
global FACTOR_S
global SMOOTH_FACTOR

FACTOR_S = 0.3
SPLINE_TENSION = 1
SMOOTH_FACTOR=1 % set to -1 to allow matlab to select
for curr_s=[2048]
    cpu = 1;
    if curr_s==2048
        cpu=8;
    end
    file_pre = sprintf('VMECfiles/wout_W7X_s%d_M16_N16_f12_cpu%d' ,curr_s,cpu)  ;
    file = sprintf('%s.nc',file_pre)
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
        F_fsa = F_rhos ./ p_rhos; % flux surface avg
%         dlmwrite(sprintf('VMECfiles/s%d_%s',curr_s,force_filename),F_fsa);
     catch e
         warning('Error, likely no wout for %s',file);
         fprintf(1,'The identifier was:\n%s',e.identifier);
         fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
     end
end

figure
plot(s,F_fsa)
hold on
% plot(s,F_FSA_fac,'DisplayName','dim_ang=120')

set(gca, 'YScale', 'log')
figure
JS_rhos = trapz(u,trapz(v,JS.*abs_g_vmec,3),2);
plot(s,JS_rhos)
xlabel('s')
ylabel('J^s')
title('J^s FSA')