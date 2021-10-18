%% runs force_error script and plots force error
% specifically for the HELIO_1CPU folder
close all
clearvars

u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

% load in VMEC data then run force error to calculate everything

deriv_method='finite difference';

force_filename = 'F_findif_full_vol_normalized_vol_weighted_correctly_with_erho_01s.txt';



filename = 'VMECfiles/wout_W7X_s128_M16_N16_f12_cpu1.nc';
file = filename
  try
    data = read_vmec(file);
    force_error
[val,s_ind] = min(abs(data.phi-0.03));
    plot_force_error

%         F_avgg = nanmean(F(:,:,:).*abs(abs_g_vmec),'all')/nanmean(mag_grad_B_pres(:,:,:).*abs(abs_g_vmec),'all')*mu0;
    %F_avgg = nanmean(F(:,:,:).*abs(abs_g_vmec)./data.Volume,'all')/nanmean(abs(presr(:,:,:)),'all');
    % must multiply F by g, then integrate over volume, then divide by
    % volume, do same for presr, then divide F_V_avg by p_V_avg

%        F_V_avg = trapz(v,trapz(u,trapz(s(s_2:end),F(2:end,:,:).*abs_g_vmec(2:end,:,:)))) ./ data.Volume;
%        p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume;

    F_V_avg = trapz(v,trapz(u,trapz(s(s_ind:end),F(s_ind:end,:,:).*abs_g_vmec(s_ind:end,:,:)))) ./ data.Volume;
    p_V_avg = trapz(v,trapz(u,trapz(s,abs(presr).*sqrt(gSS).*abs_g_vmec))) ./ data.Volume
    % this is too low by a factor of data.nfp!! must divide result in F_Avgg by
    % data.nfp to get actual vol avg of p
    F_avgg = F_V_avg / p_V_avg
    E = VMEC_W_B + W_p

    %dlmwrite(sprintf('%s',force_filename),F_avgg);
    %dlmwrite(sprintf('E.txt'),E);
 catch e
     warning('Error, likely no wout for %s',file);
     fprintf(1,'The identifier was:\n%s',e.identifier);
     fprintf(1,'\nThere was an error! The message was:\n%s',e.message);
 end


