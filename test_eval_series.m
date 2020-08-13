% test eval_series
close all, clear all
data = read_vmec('wout_HELIOTRON_16x4x4.nc');

dimS = data.ns;
dimU = 100;
dimV = 100;

s = linspace(0,1,dimS);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

suvgrid = ndgrid(s,u,v);
tic
R = eval_series2(suvgrid,data.rmnc,data,'c');
Z = eval_series2(suvgrid,data.zmns,data,'s');
toc

figure()
for i=1:dimS
    plot(R(i,:,1),Z(i,:,1))
    hold on
end
title('phi = 0')
xlabel('R')
ylabel('Z')
ylim([-1.25,1.25])
xlim([8.25,11.5])
