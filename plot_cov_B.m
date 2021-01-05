% close all
s_index = 30
u_index = 1
v_nfp_index = 1
if exist('Bs_vmec','var')== false
Bs_vmec = eval_series_nyq(suvgrid,data.bsubsmns,data,'s');
end
if exist('Bu_vmec','var')== false
Bu_vmec = eval_series_nyq(suvgrid,data.bsubumnc,data,'c');
end
if exist('Bv_vmec','var')==false
Bv_vmec = eval_series_nyq(suvgrid,data.bsubvmnc,data,'c');
end

% Plot s derivs
figure()

plot(data.phi(s_index:end),Bu_sa(s_index:end,u_index,v_nfp_index))

title(sprintf('d(Bu)/ds vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
ylabel('Deriv')
xlabel('s')

hold on
plot(data.phi(s_index:end),Bu_s(s_index:end,u_index,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(data.phi(s_index:end),Bu_vmec(s_index:end,u_index,v_nfp_index),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC B_u')

figure()
plot(data.phi(s_index:end),Bv_sa(s_index:end,u_index,v_nfp_index))
title(sprintf('d(Bv)/ds vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')

hold on
plot(data.phi(s_index:end),Bv_s(s_index:end,u_index,v_nfp_index))
ylabel('Deriv')
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(data.phi(s_index:end),Bv_vmec(s_index:end,u_index,v_nfp_index),'k')
plot(data.phi(s_index:end),Bv(s_index:end,u_index,v_nfp_index),'r--')
legend('Analytic deriv', 'Numerical deriv','Zero', 'VMEC B_v','My Bv')

% Plot u derivs
figure()

plot(u,Bs_ua(s_index,:,v_nfp_index))
title(sprintf('d(Bs)/du vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))
ylabel('Deriv')
xlabel('u')

hold on
plot(u,Bs_u(s_index,:,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(u,Bs_vmec(s_index,:,v_nfp_index),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC B_s')

figure()

plot(u,Bv_ua(s_index,:,v_nfp_index))
title(sprintf('d(Bv)/du vs u at s =%f,  v=%f',data.phi(s_index),v(v_nfp_index)))
xlabel('u')
ylabel('Deriv')

hold on
plot(u,Bv_u(s_index,:,v_nfp_index))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(u,Bv_vmec(s_index,:,v_nfp_index),'k')
plot(u,Bv(s_index,:,v_nfp_index),'r--')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC B_v','My Bv')


% Plot v derivs
figure()

plot(v,reshape(Bs_va(s_index,u_index,:),size(v)))
title(sprintf('d(Bs)/dv vs v at s=%f, u=%f',data.phi(s_index),u(u_index)))
xlabel('v')
ylabel('Deriv')
hold on
plot(v,reshape(Bs_v(s_index,u_index,:),size(v)))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(v,reshape(Bs_vmec(s_index,u_index,:),size(v)),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC B_s')

figure()

plot(v,reshape(Bu_va(s_index,u_index,:),size(v)))
title(sprintf('d(Bu)/dv vs v at s=%f, u=%f',data.phi(s_index),u(u_index)))
xlabel('v')
ylabel('Deriv')

hold on
plot(v,reshape(Bu_v(s_index,u_index,:),size(v)))
hold on
yline(0,'--')
hold on
yyaxis right
ylabel('Value')
plot(v,reshape(Bu_vmec(s_index,u_index,:),size(v)),'k')
legend('Analytic deriv','Numerical deriv','Zero', 'VMEC B_u')

