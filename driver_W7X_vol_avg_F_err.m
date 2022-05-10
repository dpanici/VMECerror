%% runs force_error script and plots force error
close all
clearvars

u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

path_to_W7X = '/p/desc-usr/dpanici/vmec/W7X_pressure/';
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

Ms = [8 10 12 14 16 18 20 22];
ftols = [4 8 12 13 14];
s_nums = [128 256 512 1024 2048]; 

cpus = [1];
iiii = 1;
force_filename = 'F_findif_normalized_01_099s_lam_full_mesh_all.txt';
% for N = Ns
for s_num = s_nums
for ftol = ftols
iiii=1
for M = Ms


for cpu=cpus
    file_id = sprintf('s%d_M%d_N%d_f%d_cpu%d',s_num,M,M,ftol,cpu);
    if ftol < 13 && M < 17
    
    file_id = sprintf('s%d_M%d_N%d_f%d_cpu%d_32GB',s_num,M,M,ftol,cpu);
    end
    foldername = sprintf('%s%s',path_to_W7X,file_id);
    filename = sprintf('%s/wout_W7X_%s.nc',foldername,file_id);
    file = filename
      try
        data = read_vmec(file);
        force_error
        [val,s_ind] = min(abs(data.phi-0.1));
	[val,s_ind_end] = min(abs(data.phi-0.99));
        plot_force_error
        
        % must multiply F by g, then integrate over volume, then divide by
        % volume, do same for presr, then divide F_V_avg by p_V_avg
        
        vol = trapz(v,trapz(u,trapz(s(s_ind:s_ind_end),abs_g_vmec(s_ind:s_ind_end,:,:))))*data.nfp
        F_V_avg = trapz(v,trapz(u,trapz(s(s_ind:s_ind_end),F(s_ind:s_ind_end,:,:).*abs_g_vmec(s_ind:s_ind_end,:,:)))) ./ vol;
        p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;
        F_avgg = F_V_avg / p_V_avg
        E = VMEC_W_B + W_p;
        
         dlmwrite(sprintf('%s/%s',foldername,force_filename),F_avgg);
     catch e
         warning('Error, likely no wout for %s',file);
         fprintf(1,'The identifier was:\n%s',e.identifier);
         fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
     end
    end

end
end



end
