%% guu
% vs s
figure()

plot(data.phi(s_index:end),g_uu(s_index:end,u_index,v_nfp_index))
title(sprintf('g_{uu} vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g_uu')

% vs u
figure()

plot(u,g_uu(s_index,:,v_nfp_index))
title(sprintf('g_{uu} vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('g_uu')


% vs v
figure()

plot(v,reshape(g_uu(s_index,u_index,:),size(v)))
title(sprintf('g_uu vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('g_uu')

%% BU*guu
% vs s
figure()

plot(data.phi(s_index:end),BU(s_index:end,u_index,v_nfp_index).*g_uu(s_index:end,u_index,v_nfp_index))
title(sprintf('BU*g_uu vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('BU*g_uu')

% vs u
figure()

plot(u,BU(s_index,:,v_nfp_index).*g_uu(s_index,:,v_nfp_index))
title(sprintf('BU*g_uu vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('BU*g_uu')


% vs v
figure()

plot(v,reshape(BU(s_index,u_index,:).*g_uu(s_index,u_index,:),size(v)))
title(sprintf('BU*g_uu vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('BU*g_uu')


%% gvu
figure()

plot(data.phi(s_index:end),g_vu(s_index:end,u_index,v_nfp_index))
title(sprintf('g_vu vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g_vu')

% vs u
figure()

plot(u,g_vu(s_index,:,v_nfp_index))
title(sprintf('g_vu vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('g_vu')


% vs v
figure()

plot(v,reshape(g_vu(s_index,u_index,:),size(v)))
title(sprintf('g_vu vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('g_vu')

%% BV*gvu
% vs s
figure()

plot(data.phi(s_index:end),BV(s_index:end,u_index,v_nfp_index).*g_vu(s_index:end,u_index,v_nfp_index))
title(sprintf('BV*g_vu vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('BV*g_vu')

% vs u
figure()

plot(u,BV(s_index,:,v_nfp_index).*g_vu(s_index,:,v_nfp_index))
title(sprintf('BV*g_vu vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('BV*g_uu')


% vs v
figure()

plot(v,reshape(BU(s_index,u_index,:).*g_uu(s_index,u_index,:),size(v)))
title(sprintf('BV*g_vu vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('BV*g_vu')

%% BU*g_uu + BV*gvu
% vs s
figure()

plot(data.phi(s_index:end),BU(s_index:end,u_index,v_nfp_index).*g_uu(s_index:end,u_index,v_nfp_index) + BV(s_index:end,u_index,v_nfp_index).*g_vu(s_index:end,u_index,v_nfp_index))
title(sprintf('B_u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_u = BU * g_uu + BV*g_vu')
if exist('Bu_vmec','var') == false
Bu_vmec = eval_series_nyq(suvgrid,data.bsubumnc,data,'c');
end
hold on
plot(data.phi(s_index:end), Bu_vmec(s_index:end,u_index,v_nfp_index))
legend('my B_u','VMEC B_u')