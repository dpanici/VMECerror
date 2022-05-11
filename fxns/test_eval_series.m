% test eval_series
close all, clear all
file='../example_files/wout_W7X_s256_M12_N12_f12_cpu1_32GB.nc';
data=read_vmec(file);
dimS = data.ns;
dimU = 100;
dimV = 100;

s = linspace(0,1,dimS);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

suvgrid = ndgrid(s,u,v);
tic
R = eval_series(suvgrid,data.rmnc,data,'c');
Z = eval_series(suvgrid,data.zmns,data,'s');
toc

figure()
for i=1:dimS
    plot(R(i,:,1),Z(i,:,1))
    hold on
end
title('phi = 0')
xlabel('R')
ylabel('Z')

