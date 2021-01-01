% run after force error script has ran
if exist('g_vmec','var') == false
g_vmec = eval_series_nyq(suvgrid,data.gmnc,data,'c');
end
figure()

plot(data.phi(10:end),g(10:end,u_index,v_nfp_index))
title(sprintf('g vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')

hold on
plot(data.phi(10:end),g_vmec(10:end,u_index,v_nfp_index))
title(sprintf('g vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')
legend('My Calc','VMEC')

figure()
subplot(1,2,1)
plot(v,reshape(g(s_index,u_index,:),size(v)))
title(sprintf('g vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('g')

subplot(1,2,2)
plot(v,reshape(g_vmec(s_index,u_index,:),size(v)))
title(sprintf('g vmec vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('g')