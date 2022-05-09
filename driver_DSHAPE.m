%% runs force_error script and plots force error
close all
clearvars

u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

path_to_W7X = '/p/desc-usr/dpanici/vmec/DSHAPE/ftol_sweep/ns_sweep/';
path_to_HELIOTRON = '../HELIOTRON/HELIO_1CPU/';

% load in VMEC data then run force error to calculate everything

% plot force_error will do as its name implies
% plot_force_error
%deriv_method='finite difference';

global SPLINE_ORDER_SPAPI
global SPLINE_TENSION
global FACTOR_S
global SMOOTH_FACTOR
global POLYFIT_DEGREE
global WINDOW_SIZE

use_lambda_full_mesh=true;

deriv_method='finite difference';
use_piecewise_lsq=false;
Ns = [1 2 3];
Ms = [8 10 12 14 16 18 20 22];
%Ms = [14 16]
%ftols = [4 8 12];
ftols = [4 8 12];
ftols = 13:18;
ftols = [4 8 12 13 14 15 16];
s_nums = [128 256 512 1024 2048]; %no s1024 needed for convergence sweep over M,N

s_nums=[128 256];
Ms=[18 20 22];
ftols=[12]

Ms=[16]
Ns=[0]
s_nums=[16 32 64 128 256 512 1024]
s_nums=[1024]
Ms=[6 10 14 16 20 24 28 32 36]
cpus = [1];
iiii = 1;
force_filename = 'F_findif_full_vol_normalized_vol_weighted_correctly_with_erho_01s_099s_lam_half_mesh.txt';
force_filename = 'F_findif_normalized_01_099s_lam_full_mesh_all.txt';
% for N = Ns
for s_num = s_nums
for ftol = ftols
iiii=1
for M = Ms

F_avg = zeros([1 size(s_nums)]);
Es = zeros([1 size(s_nums)]);

for cpu=cpus
    file_id = sprintf('s%d_M%d_N%d',s_num,M,0);
    foldername = sprintf('%s%s',path_to_W7X,file_id);
    filename = sprintf('%s/wout_DSHAPE_%s.nc',foldername,file_id);
    file = filename
    % probably skipped ftol = 8 and M<12 of cpu4?     
%    if isfile(sprintf('%s/%s',foldername,force_filename)) % dont run if F.txt already exists
%break
%    end 
      try
        data = read_vmec(file);
        force_error
        [val,s_ind] = min(abs(data.phi-0.1));
	[val,s_ind_end] = min(abs(data.phi-0.99));
        plot_force_error
        
%         F_avgg = nanmean(F(:,:,:).*abs(abs_g_vmec),'all')/nanmean(mag_grad_B_pres(:,:,:).*abs(abs_g_vmec),'all')*mu0;
        %F_avgg = nanmean(F(:,:,:).*abs(abs_g_vmec)./data.Volume,'all')/nanmean(abs(presr(:,:,:)),'all');
        % must multiply F by g, then integrate over volume, then divide by
        % volume, do same for presr, then divide F_V_avg by p_V_avg
        
%        F_V_avg = trapz(v,trapz(u,trapz(s(s_2:end),F(2:end,:,:).*abs_g_vmec(2:end,:,:)))) ./ data.Volume;
%        p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;
        vol = trapz(v,trapz(u,trapz(s(s_ind:s_ind_end),abs_g_vmec(s_ind:s_ind_end,:,:))))*data.nfp
        F_V_avg = trapz(v,trapz(u,trapz(s(s_ind:s_ind_end),F(s_ind:s_ind_end,:,:).*abs_g_vmec(s_ind:s_ind_end,:,:)))) ./ vol;
        p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume
        p_V_avg2 = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ vol
        F_avgg = F_V_avg / p_V_avg
        F_avgg2 = F_V_avg / p_V_avg2
        E = VMEC_W_B + W_p;
        
        dlmwrite(sprintf('%s/%s',foldername,force_filename),F_avgg);
%         dlmwrite(sprintf('%s/E.txt',foldername),E);
     catch e
         warning('Error, likely no wout for %s',file);
         fprintf(1,'The identifier was:\n%s',e.identifier);
         fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
     end
    end

end
end
% debug_plot_quants script has a lot of plotting scripts for checking
% intermediate quantities against VMEC
% debug_plot_quants



end

% end
% end
