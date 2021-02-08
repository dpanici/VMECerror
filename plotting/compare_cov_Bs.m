if exist('Bs_vmec','var') == false
Bs_vmec = eval_series_nyq(suvgrid,data.bsubsmns,data,'s');
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
quant = Bs;
quant_str = 'Bs';
quant_vmec = Bs_vmec;

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.01])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% ratios 
clims_ratio = [0,0.25];

%s,u 
figure()
pcolor(u,s,abs((quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)) ./ quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('u')
ylabel('s')
title(sprintf(' Abs Pct Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs((quant(:,u_index,:) - quant_vmec(:,u_index,:)) ./ quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs pct Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs((quant(s_index,:,:) - quant_vmec(s_index,:,:)) ./ quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Pct Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

