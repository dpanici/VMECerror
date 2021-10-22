close all
FACTOR_S = 0.1;

file = 'VMECfiles/wout_W7X_s128_M16_N16_f12_cpu1.nc';
data = read_vmec(file);
phia = data.phi(end)
data.phi = data.phi/data.phi(end);
rsmnc_fac = s_deriv(data.rmnc,data,'factor difference');
rsmnc_fin = s_deriv(data.rmnc,data,'finite difference');
rsmnc_smooth = s_deriv(data.rmnc,data,'smooth_spline');
rsmnc_spline = s_deriv(data.rmnc,data,'spline');
rsmnc_fin_4th = s_deriv(data.rmnc,data,'finite difference 4th');


%% test 4th order findif with test fxn
fake_dat = zeros(size(data.rmnc));
fake_dat(1,:) = data.phi.^2 + data.phi;
fake_dat_deriv_fin_4th = s_deriv(fake_dat,data,'finite difference 4th');
fake_dat_deriv_fin = s_deriv(fake_dat,data,'finite difference');
fake_dat_deriv2_fin_4th = s2_deriv(fake_dat,data,'finite difference 4th');
fake_dat_deriv2_fin = s2_deriv(fake_dat,data,'finite difference');

real_deriv = 2.*data.phi + 1;
real_deriv2 = 2.*ones(size(data.phi));

figure
title('Testing if 4th order fin dif is correct')
plot(data.phi,real_deriv,'DisplayName','real deriv')
hold on
plot(data.phi,fake_dat_deriv_fin(1,:),'--','DisplayName','findif')
hold on
plot(data.phi,fake_dat_deriv_fin_4th(1,:),'DisplayName','findif 4th')
ylabel('value')
xlabel('s')
legend
figure
title('Testing if 4th order 2nd deriv fin dif is correct')
plot(data.phi,real_deriv2,'DisplayName','real deriv')
hold on
plot(data.phi,fake_dat_deriv2_fin(1,:),'--','DisplayName','findif')
hold on
plot(data.phi,fake_dat_deriv2_fin_4th(1,:),'DisplayName','findif 4th')
ylabel('value')
xlabel('s')
legend




i=200
figure
plot(data.phi,rsmnc_fac(i,:),'DisplayName','Factored')
hold on
plot(data.phi,data.rsmnc(i,:),'--','DisplayName','VMEC output')
hold on
plot(data.phi,data.rmnc(i,:),'--','DisplayName','VMEC RMNC')
hold on
plot(data.phi,rsmnc_fin(i,:),'b-','DisplayName','findif')
hold on
plot(data.phi,rsmnc_smooth(i,:),'--','DisplayName','smooth')
hold on
plot(data.phi,rsmnc_spline(i,:),'--','DisplayName','spline')
hold on
plot(data.phi,rsmnc_fin_4th(i,:),'--','DisplayName','findif 4th')
xlabel('s')
title(sprintf(' coeffs of rmnc(%d) vs s',i))


legend
i=2
figure
plot(rsmnc_fac(:,i),'DisplayName','Factored')
hold on
xlabel('MN index')
plot(data.rsmnc(:,i),'--','DisplayName','VMEC')
hold on
plot(rsmnc_fin(:,i),'--','DisplayName','findif')
hold on
plot(rsmnc_smooth(:,i),'--','DisplayName','smooth')
hold on
plot(rsmnc_spline(:,i),'--','DisplayName','spline')

title(sprintf('fourier coeff s deriv amplitudes rsmnc at s=%f',data.phi(i)))
legend

%% take coeffs, div by rho m and mult by rho m and see if same

im = 500
coeffs_im =data.rmnc(im,:);
coeffs_im_fac = coeffs_im./sqrt(data.phi).^data.xm(im);
coeffs_im_fac_unfac = coeffs_im_fac.*sqrt(data.phi).^data.xm(im);
figure
xlabel('s')
plot(data.phi,coeffs_im,'DisplayName','original')
hold on
plot(data.phi,coeffs_im_fac_unfac,'--','DisplayName','factored then unfactored')
ylabel('fourier Amplitude')
title(sprintf('Original, then factored by rhom and re-multiplied coeffs of rmnc(%d)',im))
legend

%% compare second derivs

figure
im=20

rssmnc_fac = s2_deriv(data.rmnc,data,'factor difference');
% rssmnc_fac_2sd = s_deriv(rsmnc_smooth,data,'factor difference');
rssmnc_fin = s2_deriv(data.rmnc,data,'finite difference');
rssmnc_smooth = s2_deriv(data.rmnc,data,'smooth_spline');
rssmnc_spline = s2_deriv(data.rmnc,data,'spline');
plot(data.phi,rssmnc_fac(im,:),'DisplayName','Factored')
hold on
plot(data.phi,data.rmnc(im,:),'.--','DisplayName','VMEC RMNC')
hold on
plot(data.phi,rssmnc_fin(im,:),'--','DisplayName','findif')
hold on
plot(data.phi,rssmnc_smooth(im,:),'--','DisplayName','smooth')
hold on
plot(data.phi,rssmnc_spline(im,:),'--','DisplayName','spline')
% hold on
% plot(data.phi,rssmnc_fac_2sd(im,:),'DisplayName','Factored sderiv twice')
legend
xlabel('s')
title(sprintf('second deriv of rmnc(%d) with m = %d',im,data.xm(im)))
