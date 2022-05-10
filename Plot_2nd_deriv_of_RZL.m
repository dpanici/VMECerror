%% This script will calculate 1st and 2nd derivatives of R,Z,L fourier coefficients from VMEC
% and plot them versus s. Currently has the spike locations for W7X M=N=16
% hardcoded in but that can be easily removed
close all
clearvars
%% read in the vmec file using matlabVMEC
vmecfile = 'VMECfiles/wout_W7X_s2048_M16_N16_f12_cpu8.nc';
% vmecfile = 'VMECfiles/wout_HELIOTRON_s1024_M6_N3.nc';
% vmecfile = 'VMECfiles/wout_ESTELL.nc';
% vmecfile='VMECfiles/wout_W7X_s512_M16_N16_f12_cpu8_no_s32.nc'
% vmecfile='VMECfiles/wout_W7X_s2187_M16_N16_f12.nc'
% vmecfile='VMECfiles/wout_W7X_s625_M16_N16_f12.nc'

data = read_vmec(vmecfile);

data.phi = data.phi/data.phi(end);% normalize data.phi by phi_b to get s
% need to do this becauuse the s_deriv and s2_deriv function use data.phi
% but assume it is the normalized flux s

N = data.ns;
s = linspace(0,1,N);

%% these are the peaks in F that occur for W7X M=N=16 (hardcoded here)
% they are equispaced apart which is fishy, and occur at same spots for
% ns=512,1024,2048 and for pressure/vacuum cases and for different iota
% profiles
force_peak_s_locations=[0.7498,0.781,0.8123,0.8436,0.8749,0.9062,0.9374];

spacing = force_peak_s_locations(2)-force_peak_s_locations(1);
data.ns*spacing


start_mode_ind = 1;
%% range of s to plot on
start_s_val = 0.3;
end_s_val = 0.98;
[val,sind_1] = min(abs(data.phi-start_s_val));
[val,sind_2] = min(abs(data.phi-end_s_val));



%% rmnc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate 1st and 2nd derivatives of RMNC
% can change the deriv method by just changing the 3rd argument, possible
% deriv methods are in the s_deriv and s2_deriv functions, but
% main ones are:
% 'finite difference', 'finite difference 4th', 'spline'
y_2deriv_spline = s2_deriv(data.rmnc,data,'spline');
y_2deriv_findif = s2_deriv(data.rmnc,data,'finite difference'); 
% y_2deriv_findif=y_2deriv_spline;
y_deriv_findif = s_deriv(data.rmnc,data,'finite difference');

%% plot the value and radial derivatives of the coefficient with the largest 2nd derivative
[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_findif(:,sind_1:sind_2)),[],2));% find the mode with the largest 2nd derivative to plot the value of vs s

figure
plot(data.phi(sind_1:sind_2),data.rmnc(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of R_{mnc} vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))

for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of dR_{mnc}/ds vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_2deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of d^2R_{mnc}/ds^2 vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
% for loc=force_peak_s_locations
%     [val,sind_peak] = min(abs(data.phi-loc));
%     if loc == force_peak_s_locations(1)
%     xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
%     hold on
%     else
%     xline(data.phi(sind_peak),'--','HandleVisibility','off')
%     hold on
%     end
% end


legend
%% plot the second derivative of all the rmnc modes

%%%%%%%%%%%%%%%%% with m labelled %%%%%%%%%%%%%%%%%%%
figure

cmap = jet(max(data.xm) +1);
colormap(cmap)

for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline(i,sind_1:sind_2);
    plot(s(sind_1:sind_2),data_window,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC R_{mnc}',length(data.xm)))
c = colorbar; caxis([0, max(data.xm)]), c.Label.String = '|m|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend
%%%%%%%%%%%%%%%%% with n labelled %%%%%%%%%%%%%%%%%%%
figure

cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline(i,sind_1:sind_2);
    plot(s(sind_1:sind_2),data_window,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC R_{mnc}',length(data.xm)))
c = colorbar; caxis([0, max(data.xm)]), c.Label.String = '|n|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend

%% zmns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate 1st and 2nd derivatives of ZMNS
% can change the deriv method by just changing the 3rd argument, possible
% deriv methods are in the s_deriv and s2_deriv functions, but
% main ones are:
% 'finite difference', 'finite difference 4th', 'spline'
y_2deriv_spline_z = s2_deriv(data.zmns,data,'spline');
y_2deriv_findif_z = s2_deriv(data.zmns,data,'finite difference'); 
y_deriv_findif_z = s_deriv(data.rmnc,data,'finite difference');
%% plot the value and radial derivatives of the coefficient with the largest 2nd derivative
[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_findif_z(:,sind_1:sind_2)),[],2));% find the mode with the largest 2nd derivative to plot the value of vs s

figure
plot(data.phi(sind_1:sind_2),data.zmns(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
ylabel('Value')
xlabel('s')
title(sprintf('Value of Z_{mns} vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
ylabel('Value')
xlabel('s')
title(sprintf('Value of dZ_{mns}/ds vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_2deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of d^2Z_{mns}/ds^2 vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend

%% plot the second derivative of all the zmns modes
%%%%%%%%%%%%%%%%% with m labelled %%%%%%%%%%%%%%%%%%%
figure
cmap = jet(max(data.xm) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_z(i,sind_1:sind_2);
    plot(s(sind_1:sind_2),data_window,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC Z_{mns}',length(data.xm)))
c = colorbar; caxis([0, max(data.xn/data.nfp)]), c.Label.String = '|m|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend
%%%%%%%%%%%%%%%%% with n labelled %%%%%%%%%%%%%%%%%%%
figure

cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_z(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(s(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC Z_{mns}',length(data.xm)))
c = colorbar; caxis([0, max(data.xn/data.nfp)]), c.Label.String = '|n|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend

%% lmns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate 1st and 2nd derivatives of LMNS
% can change the deriv method by just changing the 3rd argument, possible
% deriv methods are in the s_deriv and s2_deriv functions, but
% main ones are:
% 'finite difference', 'finite difference 4th', 'spline'
y_2deriv_spline_L = s2_deriv(data.lmns,data,'spline');
y_2deriv_findif_L = s2_deriv(data.lmns,data,'finite difference'); 
y_deriv_findif_L = s_deriv(data.lmns,data,'finite difference');

%% plot the value and radial derivatives of the coefficient with the largest 2nd derivative
[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_findif_L(:,sind_1:sind_2)),[],2));% find the mode with the largest 2nd derivative to plot the value of vs s

figure
plot(data.phi(sind_1:sind_2),data.lmns(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
ylabel('Value')
xlabel('s')
title(sprintf('Value of L_{mns} vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
ylabel('Value')
xlabel('s')
title(sprintf('Value of dL_{mns}/ds vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend
figure
plot(data.phi(sind_1:sind_2),y_2deriv_findif(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of d^2L_{mns}/ds^2 vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc));
    if loc == force_peak_s_locations(1)
    xline(data.phi(sind_peak),'--','DisplayName','Location of Spike')
    hold on
    else
    xline(data.phi(sind_peak),'--','HandleVisibility','off')
    hold on
    end
end
legend

%% plot the second derivative of all the lmns modes
%%%%%%%%%%%%%%%%% with m labelled %%%%%%%%%%%%%%%%%%%
figure
cmap = jet(max(data.xm) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_L(i,sind_1:sind_2);
    plot(s(sind_1:sind_2),data_window,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC L_{mns}',length(data.xm)))
c = colorbar; caxis([0, max(data.xm)]), c.Label.String = '|m|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend
%%%%%%%%%%%%%%%%% with n labelled %%%%%%%%%%%%%%%%%%%
figure

cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_L(i,sind_1:sind_2);
    plot(s(sind_1:sind_2),data_window,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
end
title(sprintf('Second radial deriv all %d modes of VMEC L_{mns}',length(data.xm)))
c = colorbar; caxis([0, max(data.xn/data.nfp)]), c.Label.String = '|n|';
xlabel('s')
ylabel('Value')
for loc=force_peak_s_locations
    if loc == force_peak_s_locations(1)
    xline(loc,'--','DisplayName','Location of Spike')
    hold on
    else
    xline(loc,'--','HandleVisibility','off')
    hold on
    end
end
legend