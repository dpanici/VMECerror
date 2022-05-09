close all
% data = read_vmec('../VMECerror/VMECfiles/wout_W7X_s2048_M16_N16_f12_cpu8.nc');
data = read_vmec('../VMECerror/VMECfiles/wout_W7X_s1024_M16_N16_f12_cpu8_shift_iota_up_02.nc')
% data = read_vmec('VMECfiles/wout_W7X_vac_s1024_M16_N16_f12_cpu8.nc');

data.phi = data.phi/data.phi(end);

global POLY_LSQ_WINDOW_SIZE
POLY_LSQ_WINDOW_SIZE=16; % 16 was good and 6 poly order
global POLY_LSQ_ORDER
POLY_LSQ_ORDER = 6; % polynomial order

y_2deriv_spline = s2_deriv(data.rmnc,data,'spline');
y_2deriv_findif = s2_deriv(data.rmnc,data,'finite difference 4th');
[y_deriv_lsq,y2_deriv_lsq]=least_squares_fit_coeffs(data.rmnc,data,'none');

y_deriv_findif_1st = s_deriv(data.rmnc,data,'finite difference 1st');

N = data.ns;
x = linspace(0,1,N);

%% aside: plot s2 deriv of all R modes and just look at that window
% maybe normalize them to the max value in that window (so avoiding edges)
force_peak_s_locations=[0.7498,0.781,0.8123,0.8436,0.8749,0.9062,0.9374];
start_mode_ind = 1;

%% rmnc
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_findif(:,sind_1:sind_2)),[],2));
figure
plot(data.phi(sind_1:sind_2),data.rmnc(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
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
plot(data.phi(sind_1:sind_2),y_deriv_findif_1st(mode_ind_large_deriv,sind_1:sind_2),'.-','HandleVisibility','off')
ylabel('Value')
xlabel('s')
title(sprintf('Value of dR_{mnc}/ds vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc))
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
% plot(data.phi(sind_1:sind_2),y_2deriv_spline(mode_ind_large_deriv,sind_1:sind_2),'.-','DisplayName','Spline')
% hold on
ylabel('Value')
xlabel('s')
title(sprintf('Value of d^2R_{mnc}/ds^2 vs s near edge for m=%d n=%d',data.xm(mode_ind_large_deriv),data.xn(mode_ind_large_deriv)/data.nfp))
for loc=force_peak_s_locations
    [val,sind_peak] = min(abs(data.phi-loc))
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

cmap = jet(max(data.xm) +1);
colormap(cmap)

for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'--','Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
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
%%%%%%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%%%
figure
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
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

%% zmns
figure
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
[y_deriv_lsq,y2_deriv_lsq] = least_squares_fit_coeffs(data.zmns,data,'none');
y_2deriv_spline_z = s2_deriv(data.zmns,data,'spline');

[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_spline_z(:,sind_1:sind_2)),[],2));
figure
plot(data.phi(sind_1:sind_2),data.zmns(mode_ind_large_deriv,sind_1:sind_2),'-')
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



cmap = jet(max(data.xm) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_z(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'--','Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
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
%%%%%%%%%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%%%
figure
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
y_2deriv_spline_z = s2_deriv(data.zmns,data,'spline');


cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_z(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'--','Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
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

%% lmns
figure
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
y_2deriv_spline_L = s2_deriv(data.lmns,data,'spline');
[y_deriv_lsq,y2_deriv_lsq]= least_squares_fit_coeffs(data.lmns,data,'none');

[val,mode_ind_large_deriv] = max(max(abs(y_2deriv_spline_L(:,sind_1:sind_2)),[],2));
figure
plot(data.phi(sind_1:sind_2),data.lmns(mode_ind_large_deriv,sind_1:sind_2),'.')
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

cmap = jet(max(data.xm) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_L(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'--','Color',cmap(abs(data.xm(i)) + 1, :),'HandleVisibility','off')
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
%%%%%%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%
figure
[val,sind_1] = min(abs(data.phi-0.7));
[val,sind_2] = min(abs(data.phi-0.98));
y_2deriv_spline_L = s2_deriv(data.lmns,data,'spline');

cmap = jet(max(data.xn/data.nfp) +1);
colormap(cmap)
for i=start_mode_ind:length(data.xm)
    data_window = y_2deriv_spline_L(i,sind_1:sind_2);
    data_window_normed = data_window;%./ max(abs(data_window));
    plot(x(sind_1:sind_2),data_window_normed,'Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
    hold on
    plot(x(sind_1:sind_2),y2_deriv_lsq(i,sind_1:sind_2),'--','Color',cmap(abs(data.xn(i)/data.nfp) + 1, :),'HandleVisibility','off')
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