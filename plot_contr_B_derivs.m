close all
s_index = 30
u_index = 1
v_nfp_index = 1

if exist('BU_vmec','var')== false
BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
end
if exist('BV_vmec','var')==false
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
end

%% Plot s derivs
figure()

plot(data.phi(s_index:end),BU_sa(s_index:end,u_index,v_nfp_index))

title(sprintf('d(BU)/ds vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
ylabel('Deriv')
xlabel('s')

hold on
plot(data.phi(s_index:end),BU_s(s_index:end,u_index,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(data.phi(s_index:end),BU_vmec(s_index:end,u_index,v_nfp_index),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC BU')

figure()
plot(data.phi(s_index:end),BV_sa(s_index:end,u_index,v_nfp_index))
title(sprintf('d(BV)/ds vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')

hold on
plot(data.phi(s_index:end),BV_s(s_index:end,u_index,v_nfp_index))
ylabel('Deriv')
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(data.phi(s_index:end),BV_vmec(s_index:end,u_index,v_nfp_index),'k')
plot(data.phi(s_index:end),BV(s_index:end,u_index,v_nfp_index),'r--')
legend('Analytic deriv', 'Numerical deriv','Zero', 'VMEC B_v','My Bv')

%% Plot u derivs
figure()

plot(u,BU_ua(s_index,:,v_nfp_index))
title(sprintf('d(BU)/du vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))
ylabel('Deriv')
xlabel('u')

hold on
plot(u,BU_u(s_index,:,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(u,BU_vmec(s_index,:,v_nfp_index),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC BU')

figure()

plot(u,BV_ua(s_index,:,v_nfp_index))
title(sprintf('d(BV)/du vs u at s =%f,  v=%f',data.phi(s_index),v(v_nfp_index)))
xlabel('u')
ylabel('Deriv')

hold on
plot(u,BV_u(s_index,:,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(u,BV_vmec(s_index,:,v_nfp_index),'k')
plot(u,BV(s_index,:,v_nfp_index),'r--')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC BV','My BV')


%% Plot v derivs
figure()

plot(v,reshape(BU_va(s_index,u_index,:),size(v)))
title(sprintf('d(BU)/dv vs v at s=%f, u=%f',data.phi(s_index),u(u_index)))
xlabel('v')
ylabel('Deriv')
hold on
plot(v,reshape(BU_v(s_index,u_index,:),size(v)))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(v,reshape(BU_vmec(s_index,u_index,:),size(v)),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC BU')

figure()

plot(v,reshape(BV_va(s_index,u_index,:),size(v)))
title(sprintf('d(BV)/dv vs v at s=%f, u=%f',data.phi(s_index),u(u_index)))
xlabel('v')
ylabel('Deriv')

hold on
plot(v,reshape(BV_v(s_index,u_index,:),size(v)))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(v,reshape(BV_vmec(s_index,u_index,:),size(v)),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC BV')

