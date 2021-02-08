if exist('BV_vmec,','var') == false
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
end
% vs s
figure()

plot(data.phi(s_index:end),BV(s_index:end,u_index,v_nfp_index))
title(sprintf('B^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^v')


hold on
plot(data.phi(s_index:end),BV_vmec(s_index:end,u_index,v_nfp_index))
title(sprintf('B^v vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^v')
legend('My Calc','VMEC')

% vs u
figure()

plot(u,BV(s_index,:,v_nfp_index))
title(sprintf('B^v vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('u')
ylabel('B^u')

hold on
plot(u,BV_vmec(s_index,:,v_nfp_index))
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('u')
ylabel('B^v')
legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(BV(s_index,u_index,:),size(v)))
title(sprintf('B^v vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^v')

hold on
plot(v,reshape(BV_vmec(s_index,u_index,:),size(v)))
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^v')
legend('My Calc','VMEC')


% figure()
% subplot(1,2,1)
% plot(v,reshape(BV(s_index,u_index,:),size(v)))
% title(sprintf('B^v vs v at s=%f, u=%f',s(s_index),u(u_index)))
% xlabel('v')
% ylabel('B^v')

% subplot(1,2,2)
% plot(v,reshape(BV_vmec(s_index,u_index,:),size(v)))
% title(sprintf('B^v vmec vs v at s=%f, u=%f',s(s_index),u(u_index)))
% xlabel('v')
% ylabel('B^v')
% 
% figure()
% plot(u,BV(2,:,v_nfp_index))
% title(sprintf('B^v vs u at s=%f, nfp*phi=%f',s(2),v(v_nfp_index)))
% xlabel('u')
% ylabel('B^v')

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BV(:,:,nfp_v_index)-BV_vmec(:,:,nfp_v_index))))
colorbar; 
colormap jet
 caxis([-5 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('log scale difference btwn my B^v and VMEC at nfp*phi=%f',v(nfp_v_index)))


%% plot surfaces
quant = BV;
quant_str = 'BV';
quant_vmec = BV_vmec;
clims_ratio = [0,0.1];
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
