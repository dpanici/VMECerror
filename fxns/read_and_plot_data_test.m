clear all, clearvars
% Read tutorial VMEC results and do some sort of analysis of outputs
data = read_vmec('wout_HELIOTRON_16x4x4.nc');


% VMECplot(data) % GUI to plot various quantities such as flux surfaces, p
% or iota profiles, etc, does not need an input, will let you select from
% wout files in current directory

% Fourier coefficients are returned as a 2D array in the "NESCOIL" (nu+nv)
% format

%rmnc, zmns store the R and Z coefficients, and these arrays are 2D
% first index is over the radial coordinate s (discrete radial surfaces),
% the second index is vectorized over the poloidal/toroidal modes
% the xn, xm variables contain the toroidal and poloidal mode numbers for a
% given second index of rmnc/zmns.

% INDEXING is actually reversed, the rows are for the modes and columns are
% the radial s coordinate so (:,s_ind) gives you the modes for the given
% s_ind flux surface

% these coefficients are for fourier series of the form R = sum_0,-N^ M,N of rmnc(s) * cos(mu -
% nvNfp), where Nfp is the number of field periods in the device (Nfp=0 for
%axisymmetric)

nfp = data.nfp;

for i = 1:length(data.xm)
    fprintf('index = %d, m=%d, n=%d\n',i,data.xm(i),data.xn(i))
end

%to get R,Z for a given flux surface with label s1, you take the sum of
%rmnc(:,1) * cos(mu - nvNfp)

%as an example let's fix v=0 to get just a single toroidal XS of the field.
u = linspace(0,2*pi,100);
v = linspace(0,2*pi,6*19+1);
figure()
for k=1:
%figure() 
for j=1:data.ns
R = zeros([1,length(u)]);
Z = zeros([1,length(u)]);

for i = 1:length(data.xm)
    R = R + data.rmnc(i,j) .* cos(data.xm(i).*u - data.xn(i).*v(k));
    Z = Z + data.zmns(i,j) .* sin(data.xm(i).*u - data.xn(i).*v(k));
end

plot(R,Z)
hold on
end
if k==1
    %scatter(sum(data.raxiscc),data.zaxiscs(1))
end
title(sprintf('phi = %f',v(k)))
ylabel('Z(m)')
ylim([-1.25,1.25])
xlim([8.25,11.5])
xlabel('R(m)')
end


