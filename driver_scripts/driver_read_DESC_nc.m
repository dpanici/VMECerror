%% runs force_error script and plots force error
close all
clearvars
n = 21;
Es = zeros([1 n]);
Es_mine = zeros([1 n]);

F_means = zeros([1 n]);

dxs = linspace(-5e-2,5e-2,n);


u_index= 11; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=6; % index of s to plot quantities at

% method used to calculate radial derivatives
% deriv_method='makima';
deriv_method='finite difference';
% deriv_method= 'spline'
% deriv_method='factor_spline'
% deriv_method = 'pchip';
% deriv_method = 'makima';
% deriv_method='poly' % currently takes forever, try smaller poly? 

% file = 'VMECfiles/wout_DSHAPE_s512_M13_N0.nc';
% file = '../VMECerror/VMECfiles/wout_DSHAPE.nc';
% file = 'VMECfiles/wout_HELIOTRON_s16_M6_N3.nc';
% file = 'VMECfiles/wout_HELIOTRON_s512_M12_N3.nc'; % This HELIOTRON case agrees
% very well everywhere except near axis

% file = 'VMECfiles/wout_HELIOTRON_Vac.nc';
% file = 'VMECfiles/wout_SOLOVEV.nc';
% file = 'VMECfiles/wout_W7X_standard.nc';
% file = 'VMECfiles/wout_W7X_s758_M15_N15.nc'
% file = 'VMECfiles/wout_W7X_standard_M15_N12_s256.nc'
% file = 'VMECfiles/wout_W7X_s256_M15_N15_LBSUBS_F.nc'

% file = 'VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu1.nc' % from VMEC runs
% file = 'DESC_wout_W7X_M16_N16_ansi.nc';

file = 'DESC_wout_W7X_s2048_M16_N16_ansi.nc';
%% CURRENTLY compare_g.m has things for converting g from VMEC -> DESC coordinate system

% load in VMEC data then run force error to calculate everything
iii = 1
% for dx=dxs
data = read_vmec(file);
% data.rmnc(1,:) = data.rmnc(1,:) + dx;
force_error

F_rhos = trapz(u,trapz(v,F.*abs_g_vmec,3),2);
p_rhos = abs(trapz(u,trapz(v,presr.*sqrt(gSS).*abs_g_vmec,3),2));
plot_force_error




if max(presf)<1e-3
p_B_rhos = abs(trapz(u,trapz(v,mag_grad_B_pres.*abs_g_vmec,3),2))/mu0;
figure
plot(s(2:end),F_rhos(2:end)./p_B_rhos(2:end))
ylabel('Normalized <F> by gradient of B^2')
xlabel('s')
else
figure
plot(s(5:end),F_rhos(5:end)./p_rhos(5:end))
ylabel('Normalized <F> by |grad(p)|')
xlabel('s')  
end
% % plot force_error will do as its name implies
% plot_force_error
% 
% Es(iii) = VMEC_W_B + W_p;
% Es_mine(iii) = W_B + W_p;
% F_means(iii) = nanmean(F(8:end,:,:),'all');
% iii = iii + 1
% end
% debug_plot_quants script has a lot of plotting scripts for checking
% intermediate quantities against VMEC
% debug_plot_quants

% figure()
% plot(dxs,F_means)
% title('average |F| over volume versus dr')
% xlabel('dr')
% ylabel('F [N/m^3]')
% ylim([0 1e4])
% 
% figure()
% plot(dxs,Es)
% title('E versus dr')
% xlabel('dr')
% ylabel('E [J]')
% 
% figure()
% plot(dxs,Es_mine)
% title('My calc E versus dr')
% xlabel('dr')
% ylabel('E [J]')
%    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index) ./ abs(presr(i_s,:,nfp_v_index)) ./ sqrt(gSS(i_s,:,nfp_v_index))))

for v_ind=[1,16,31,46]
    R_sec = R(:,:,v_ind);
    Z_sec = Z(:,:,v_ind);
    F_sec = F(:,:,v_ind);
    presr_sec = presr(:,:,v_ind);
    gSS_sec = gSS(:,:,v_ind);
    
%     writematrix(R_sec,sprintf('W7X_calc_data/R/R_phi_%1.3f.csv',v(v_ind)))
%     writematrix(Z_sec,sprintf('W7X_calc_data/Z/Z_phi_%1.3f.csv',v(v_ind)))
%     writematrix(F_sec,sprintf('W7X_calc_data/F/F_phi_%1.3f.csv',v(v_ind)))
%     writematrix(presr_sec,sprintf('W7X_calc_data/presr/presr_phi_%1.3f.csv',v(v_ind)))
%     writematrix(gSS_sec,sprintf('W7X_calc_data/gSS/gSS_phi_%1.3f.csv',v(v_ind)))
    
end

% writematrix(F_rhos./p_rhos,sprintf('F_FSA_DESC_VMEC_%s.csv',deriv_method));
figure
JS_rhos = trapz(u,trapz(v,JS.*abs_g_vmec,3),2);
plot(s,JS_rhos)
xlabel('s')
ylabel('J^s')
title('J^s FSA')