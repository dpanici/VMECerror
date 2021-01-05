if exist('Bs_vmec','var') == false
Bs_vmec = eval_series_nyq(suvgrid,data.bsubsmns,data,'s');
end
if exist('Bu_vmec','var') == false
Bu_vmec = eval_series_nyq(suvgrid,data.bsubumnc,data,'c');
end
if exist('Bv_vmec','var') == false
Bv_vmec = eval_series_nyq(suvgrid,data.bsubvmnc,data,'c');
end
%% Plot Bs
% vs s
figure()

plot(data.phi(s_index:end),Bs(s_index:end,u_index,v_nfp_index))
title(sprintf('B_s vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_s')


hold on
plot(data.phi(s_index:end),Bs_vmec(s_index:end,u_index,v_nfp_index),'k--')
title(sprintf('B_s vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_s')
legend('My Calc','VMEC')

% sweep over u and plot diff in my B_s and VMEC vs s at various u
figure()
for i = 1:size(u)
    diff =abs(Bs_vmec(s_index:end,i,v_nfp_index)-Bs(s_index:end,i,v_nfp_index));
    if min(ismembertol(Bs_vmec(s_index:end,i,v_nfp_index),Bs(s_index:end,i,v_nfp_index))) == 0
    plot(data.phi(s_index:end),diff,'DisplayName',sprintf('u=%f',u(i)))
    hold on
    end
end
title(sprintf('diff in B_s and VMEC B_s vs s at various u for v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('diff in B_s')
legend

% vs u
figure()

plot(u,Bs(s_index,:,v_nfp_index))
title(sprintf('B_s vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))

hold on
plot(u,Bs_vmec(s_index,:,v_nfp_index),'k--')
xlabel('u')
ylabel('B_s')
legend('My Calc','VMEC')

% sweep over s and plot diff in my B_s and VMEC vs u at various s
figure()
for i = 1:data.ns
    diff =abs(Bs_vmec(i,:,v_nfp_index)-Bs(i,:,v_nfp_index));
    if max(diff) > 0.01
    plot(u,diff,'DisplayName',sprintf('s=%f',data.phi(i)))
    hold on
    end
end
title(sprintf('diff in B_s and VMEC B_s vs u at various s for v = %f',v(v_nfp_index)))
xlabel('u')
ylabel('diff in B_s')
legend


% vs v
figure()

plot(v,reshape(Bs(s_index,u_index,:),size(v)))
title(sprintf('B_s vs v  at s=%f, u==%f',data.phi(s_index),u(u_index)))

hold on
plot(v,reshape(Bs_vmec(s_index,u_index,:),size(v)),'k--')
xlabel('v')
ylabel('B_s')
legend('My Calc','VMEC')

%% Plot Bu
% vs s
figure()

plot(data.phi(s_index:end),Bu(s_index:end,u_index,v_nfp_index))
title(sprintf('B_u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_u')


hold on
plot(data.phi(s_index:end),Bu_vmec(s_index:end,u_index,v_nfp_index),'k--')
title(sprintf('B_u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_u')
legend('My Calc','VMEC')


figure()
for i = 1:size(u)
    diff =abs(Bu_vmec(s_index:end,i,v_nfp_index)-Bu(s_index:end,i,v_nfp_index));
    if max(diff) > 0.005
    plot(data.phi(s_index:end),diff,'DisplayName',sprintf('u=%f',u(i)))
    hold on
    end
end
title(sprintf('diff in B_u and VMEC B_u vs s at various u for v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('diff in B_u')
legend

% vs u
figure()

plot(u,Bu(s_index,:,v_nfp_index))
title(sprintf('B_u vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))

hold on
plot(u,Bu_vmec(s_index,:,v_nfp_index),'k--')
xlabel('u')
ylabel('B_u')
legend('My Calc','VMEC')

figure()
for i = 1:data.ns
    diff =abs(Bu_vmec(i,:,v_nfp_index)-Bu(i,:,v_nfp_index));
    if max(diff) > 0.01
    plot(u,diff,'DisplayName',sprintf('s=%f',data.phi(i)))
    hold on
    end
end

title(sprintf('diff in B_u and VMEC B_u vs u at various s for v = %f',v(v_nfp_index)))
xlabel('u')
ylabel('B_u')
legend



% vs v
figure()

plot(v,reshape(Bu(s_index,u_index,:),size(v)))
title(sprintf('B_u vs v  at s=%f, u==%f',data.phi(s_index),u(u_index)))

hold on
plot(v,reshape(Bu_vmec(s_index,u_index,:),size(v)),'k--')
xlabel('v')
ylabel('B_u')
legend('My Calc','VMEC')

% Plot Bv
% vs s
figure()

plot(data.phi(s_index:end),Bv(s_index:end,u_index,v_nfp_index))
title(sprintf('B_v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_v')


hold on
plot(data.phi(s_index:end),Bv_vmec(s_index:end,u_index,v_nfp_index),'k--')
title(sprintf('B_v vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_v')
legend('My Calc','VMEC')

% vs u
figure()

plot(u,Bv(s_index,:,v_nfp_index))
title(sprintf('B_v vs u at s=%f, nfp*phi=%f',data.phi(s_index),v(v_nfp_index)))

hold on
plot(u,Bv_vmec(s_index,:,v_nfp_index),'k--')
xlabel('u')
ylabel('B_v')
legend('My Calc','VMEC')

figure()
for i = 1:data.ns
    diff =abs(Bv_vmec(i,:,v_nfp_index)-Bv(i,:,v_nfp_index));
    if max(diff) > 0.01
    plot(u,diff,'DisplayName',sprintf('s=%f',data.phi(i)))
    hold on
    end
end

title(sprintf('diff in B_v and VMEC B_v vs u at various s for v = %f',v(v_nfp_index)))
xlabel('u')
ylabel('B_v')
legend

% vs v
figure()

plot(v,reshape(Bv(s_index,u_index,:),size(v)))
title(sprintf('B_v vs v  at s=%f, u==%f',data.phi(s_index),u(u_index)))

hold on
plot(v,reshape(Bv_vmec(s_index,u_index,:),size(v)),'k--')
xlabel('v')
ylabel('B_v')
legend('My Calc','VMEC')

%% plot 2d
% Bs
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index), abs((Bs(:,:,nfp_v_index) - Bs_vmec(:,:,nfp_v_index))))%./Bu(:,:,v_nfp_index) .* 100))
colormap jet
caxis([0,0.005])
colorbar;
title(sprintf('Difference in my B_s and VMEC B_s at v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('u')

axis equal


% Bu
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index), abs((Bu(:,:,nfp_v_index) - Bu_vmec(:,:,nfp_v_index))))%./Bu(:,:,v_nfp_index) .* 100))
colormap jet
caxis([0,0.005])
colorbar;
title(sprintf('Difference in my B_u and VMEC B_u at v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('u')

axis equal

% Bv
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index), abs((Bv(:,:,nfp_v_index) - Bv_vmec(:,:,nfp_v_index))))%./Bu(:,:,v_nfp_index) .* 100))
colormap jet
caxis([0,0.005])
colorbar;
title(sprintf('Difference in my B_v and VMEC B_v at v = %f',v(v_nfp_index)))
xlabel('s')
ylabel('u')

axis equal
