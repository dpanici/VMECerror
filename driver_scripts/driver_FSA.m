%% runs force_error script and plots force error
close all
clearvars

s_index=6; % index of s to start plotting quantities at

path_to_W7X = '/p/desc-usr/dpanici/vmec/W7X_pressure/';
% load in VMEC data then run force error to calculate everything

% plot force_error will do as its name implies
% plot_force_error
% deriv_method='factor difference';
deriv_method='finite difference';
% deriv_method = 'spapi'
% deriv_method = 'spline'
% deriv_method='poly'
% deriv_method = 'smooth_spline'
% deriv_method='tension_spline'
use_piecewise_lsq = false; % 
use_my_cubic_spline=false;
Ns = [1 2 3];
Ms = [8 10 12 14 16];
%Ms = [14 16]
%ftols = [4 8 12];
ftols = [4 8 12];

s_nums = [128 256 512 1024 2048]; %no s1024 needed for convergence sweep over M,N

cpus = [1 2 4 8];
iiii = 1;
force_filename = 'spline_6_F_FSA_norm_by_p_rho_avg.txt';
global SPLINE_ORDER_SPAPI
global SPLINE_TENSION
global FACTOR_S
global SMOOTH_FACTOR
global POLYFIT_DEGREE
global MY_SPLINE_END_CONDITION
MY_SPLINE_END_CONDITION = 'natural';
%% LSQ params
global POLY_LSQ_WINDOW_SIZE
POLY_LSQ_WINDOW_SIZE=16; % 16 was good and 6 poly order
global POLY_LSQ_ORDER
POLY_LSQ_ORDER = 5; % polynomial order
%%

POLYFIT_DEGREE=10

FACTOR_S = 0.3
SPLINE_ORDER_SPAPI=6
SPLINE_TENSION = 1
SMOOTH_FACTOR=0.99999999 % set to -1 to allow matlab to select
% SMOOTH_FACTOR=-1
figure
for deriv_method=["finite difference","spline"]
    
for POLY_LSQ_WINDOW_SIZE=16

for curr_s=[2048]
    cpu = 1;
    if curr_s==2048
        cpu=8;
    end
    file_pre = sprintf('VMECfiles/wout_W7X_s%d_M16_N16_f12_cpu%d' ,curr_s,cpu)  ;
    file = sprintf('%s.nc',file_pre);
%     file = 'VMECfiles/wout_DSHAPE_s1024_M14_N0.nc';
%     file = 'VMECfiles/wout_DSHAPE_s256_M14_N1_LBSUBS_T.nc'
%     file = 'VMECfiles/wout_HELIOTRON_s1024_M12_N2.nc'
%     file = 'VMECfiles/wout_HELIOTRON_s1024_M6_N3.nc'
%     file = 'VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu8_const_iota.nc' %
%     const noble irrational
%     file = 'VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu8_zero_iota.nc' %

%     file = 'wout_DSHAPE_s1024_M14_N1_iota_1.nc'
% file = 'VMECfiles/wout_DSHAPE_s256_M14_N0_zero_iota_initial.nc'
% file = 'VMECfiles/wout_W7X_vac_s1024_M16_N16_f12_cpu8.nc'
      try
        data = read_vmec(file);
        force_error
	[val,s_ind] = min(abs(data.phi-0.1));
%         plot_force_error
        
        p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;
       % F_avgg = F_V_avg / p_V_avg
        E = VMEC_W_B + W_p;
        F_rhos = trapz(u,trapz(v,F.*abs_g_vmec,3),2);
        p_rhos = trapz(u,trapz(v,abs(presr).*sqrt(gSS).*abs_g_vmec,3),2);
        F_fsa = F_rhos ./ p_rhos; % flux surface avg
%         force_filename = 'DSHAPE_F_FSA.txt';
%         dlmwrite(sprintf('VMECfiles/s%d_%s',curr_s,force_filename),F_fsa);
     catch e
         warning('Error, likely no wout for %s',file);
         fprintf(1,'The identifier was:\n%s',e.identifier);
         fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
     end
end
% plot_rationals
plot_label=sprintf('%s',deriv_method)
figure
plot(s(80:data.ns-20),F_fsa(80:data.ns-20),'DisplayName',plot_label)
hold on
set(gca, 'YScale', 'log')
% legend
ylabel('Normalized Force Error FSA')
xlabel('s')
% if exist('F_FSA_fac','var')
% 
end
end
% plot_label2='winsize=9, k=5'
% % 
% plot(s,F_FSA_fac,'--','DisplayName',plot_label2)
% % % 
legend
% % end
% title('W7X F FSA for Two Iota Profiles')
% for i=unique(rat_s)
%     xline(i)
%     hold on
% end

%% plot F error (not normalized)
% figure
% plot(s,F_FSA_fac,'DisplayName','Vacuum')
% hold on
% plot(s,F_rhos,'DisplayName','Pressure')
% legend
% ylabel('Force Error FSA (N/m)')
% xlabel('s')
% set(gca, 'YScale', 'log')
% 


figure
title('Iota Profile')
plot(data.phi,data.iotaf)
ylabel('Iota')
xlabel('s')

F_s_fsa = trapz(u,trapz(v,F_s.*abs_g_vmec,3),2);
F_beta_fsa = trapz(u,trapz(v,F_beta.*abs_g_vmec,3),2);
gSS_fsa = trapz(u,trapz(v,gSS.*abs_g_vmec,3),2);
mag_beta_fsa = trapz(u,trapz(v,new_mag_beta.*abs_g_vmec,3),2);



figure
subplot(5,1,1)
ax1 = gca;
plot(s,F_rhos)
title('Force Error in (N/m^3)')
set(gca, 'YScale', 'log')
subplot(5,1,2)
ax2 = gca;
plot(s,abs(F_s_fsa))
title('F_s')
set(gca, 'YScale', 'log')
subplot(5,1,3)
ax3 = gca;
plot(s,abs(F_beta_fsa))
title('F_beta')
set(gca, 'YScale', 'log')
subplot(5,1,4)
ax4 = gca;
plot(s,gSS_fsa)
title('gSS')
set(gca, 'YScale', 'log')
subplot(5,1,5)
ax5 = gca;
plot(s,mag_beta_fsa)
title('|beta|')
set(gca, 'YScale', 'log')
linkaxes([ax1 ax2 ax3 ax4 ax5],'x')





figure
JS_rhos = trapz(u,trapz(v,JS.*abs_g_vmec,3),2);
plot(s,JS_rhos,'DisplayName','JS FSA')
legend
xlabel('s')
ylabel('J^s')
title('J^s FSA')