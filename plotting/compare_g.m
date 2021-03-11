% run after force error script has ran
if exist('g_vmec','var') == false
g_vmec = eval_series_nyq(suvgrid,data.gmnc,data,'c');
end

% vs s
figure()

plot(data.phi(s_index:end),g(s_index:end,u_index,v_nfp_index))
title(sprintf('g vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')

hold on
plot(data.phi(s_index:end),g_vmec(s_index:end,u_index,v_nfp_index))
title(sprintf('g vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')
legend('My Calc','VMEC')

% sweep over u and plot diff in my g and VMEC vs s at various u
figure()
for i = 1:size(u)
    diff =abs(g_vmec(s_index:end,i,v_nfp_index)-g(s_index:end,i,v_nfp_index));
    if min(ismembertol(g_vmec(s_index:end,i,v_nfp_index),g(s_index:end,i,v_nfp_index))) == 0
    plot(data.phi(s_index:end),diff,'DisplayName',sprintf('u=%f',u(i)))
    hold on
    end
end
title(sprintf('diff in g and VMEC g vs s at various u for v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('diff in g')
legend

% vs u
figure()

plot(u,g(s_index,:,v_nfp_index))
title(sprintf('g vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))

hold on
plot(u,g_vmec(s_index,:,v_nfp_index),'k--')
xlabel('u')
ylabel('B_s')
legend('My Calc','VMEC')

% sweep over s and plot diff in my g and VMEC vs u at various s
figure()
for i = 1:data.ns
    diff =abs(g_vmec(i,:,v_nfp_index)-g(i,:,v_nfp_index));
    if max(diff) > 0.01
    plot(u,diff,'DisplayName',sprintf('s=%f',data.phi(i)))
    hold on
    end
end
title(sprintf('diff in g and VMEC g vs u at various s for v = %f',v(v_nfp_index)))
xlabel('u')
ylabel('diff in g')
legend

% vs v
figure()
plot(v,reshape(g(s_index,u_index,:),size(v)))
title(sprintf('g vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('g')
hold on
plot(v,reshape(g_vmec(s_index,u_index,:),size(v)),'k--')
legend('My Calc','VMEC')

%% plot 2d
% g
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index), abs((g(:,:,nfp_v_index) - g_vmec(:,:,nfp_v_index))),10)
colormap jet
caxis([0,0.01])
colorbar;
title(sprintf('Difference in my g and VMEC g at v = %f',v(v_nfp_index)))
xlabel('R')
ylabel('Z')
axis equal

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),g(:,:,nfp_v_index),10)
colormap jet
% caxis([0,0.01])
colorbar;
title(sprintf('Difference in my g and VMEC g at v = %f',v(v_nfp_index)))
xlabel('R')
ylabel('Z')
axis equal



% s,u
figure()
pcolor(u,s,abs(g(:,:,nfp_v_index) - g_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.01])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my g and VMEC g at v = %f',v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(g(:,u_index,:) - g_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my g and VMEC g at u = %f',u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(g(s_index,:,:) - g_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my g and VMEC g at s = %f',data.phi(s_index)))

