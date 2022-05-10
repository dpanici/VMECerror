% test eval_series
close all, clear all
% data = read_vmec('wout_HELIOTRON_16x4x4.nc');
data = read_vmec('wout_HELIOTRON.nc');

dimS = data.ns;
dimU = 100;
dimV = 100;

s = linspace(0,1,dimS);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

suvgrid = ndgrid(s,u,v);
tic
rr = data.rmnc;
pad_len = length(data.xn_nyq)-length(data.xn);
rr_pad = padarray(rr,[pad_len 0],0);
rs = rr_pad(1 + pad_len:end,:);
zz = data.zmns;
zz_pad = padarray(zz,[pad_len 0],0);
zs = zz_pad(1 + pad_len:end,:);

data.xn_nyq = [data.xn zeros(1,pad_len)];

data.xm_nyq = [data.xm zeros(1,pad_len)];

R = eval_series_nyq(suvgrid,rs,data,'c');
Z = eval_series_nyq(suvgrid,zs,data,'s');
toc

figure()
for i=1:dimS
    plot(R(i,:,10),Z(i,:,10))
    hold on
end
title('phi = 0')
xlabel('R')
ylabel('Z')
ylim([-1.25,1.25])
xlim([8.25,11.5])
