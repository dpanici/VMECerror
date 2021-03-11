if exist('Bv_vmec','var') == false
Bv_vmec = eval_series_nyq(suvgrid,data.bsubvmnc,data,'c');
end

%% Plot Bv
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
quant = Bv;
quant_str = 'Bv';
quant_vmec = Bv_vmec;

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
caxis([0,0.01])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.01])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% ratios 

clims_ratio = [0,0.1];
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
