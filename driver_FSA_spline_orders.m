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
deriv_method='finite difference';
Ns = [1 2 3];
Ms = [8 10 12 14 16];
%Ms = [14 16]
%ftols = [4 8 12];
ftols = [4 8 12];

s_nums = [128 256 512 1024]; %no s1024 needed for convergence sweep over M,N

cpus = [1 2 4 8];
iiii = 1;
force_filename = 'F_findif_FSA_norm_by_p_V_avg.txt';
global SPLINE_ORDER_SPAPI
% for N = Ns
for s_num = s_nums
for ftol = ftols
iiii=1
for M = Ms

F_avg = zeros([1 size(s_nums)]);
Es = zeros([1 size(s_nums)]);

for cpu=cpus
    file_id = sprintf('s%d_M%d_N%d_f%d_cpu%d',s_num,M,M,ftol,cpu);
    foldername = sprintf('%s%s',path_to_W7X,file_id);
    filename = sprintf('%s/wout_W7X_%s.nc',foldername,file_id);
    file = 'VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu1.nc'   
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
        dlmwrite(sprintf('%s/order%d_%s',foldername,SPLINE_ORDER_SPAPI,force_filename),F_fsa);
%         dlmwrite(sprintf('%s/E.txt',foldername),E);
     catch e
         warning('Error, likely no wout for %s',file);
         fprintf(1,'The identifier was:\n%s',e.identifier);
         fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
     end
    end

end
end




end

% end
% end
