if exist('Z_s_vmec,','var') == false
Z_s_vmec = eval_series(suvgrid,data.zsmns,data,'s');
end

quant = Z_s;
quant_str = 'Z_s';
quant_vmec = Z_s_vmec;

% vs s
figure()
plot(data.phi(s_index:end),quant(s_index:end,u_index,v_nfp_index))
title(sprintf('%s vs s at u=%f, nfp*phi=%f',quant_str,u(u_index),v(v_nfp_index)))



hold on
plot(data.phi(s_index:end),quant_vmec(s_index:end,u_index,v_nfp_index))
xlabel('s')
ylabel(sprintf('%s',quant_str))
legend('My Calc','VMEC')

% vs u
figure()

plot(u,quant(s_index,:,v_nfp_index))
title(sprintf('%s vs v at u=%f, s=%f',quant_str,u(u_index),s(s_index)))
xlabel('u')
ylabel(sprintf('%s',quant_str))

hold on
plot(u,quant_vmec(s_index,:,v_nfp_index))
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(quant(s_index,u_index,:),size(v)))
title(sprintf('%s vs v at u=%f, s=%f',quant_str,u(u_index),s(s_index)))
xlabel('v')
ylabel(sprintf('%s',quant_str))

hold on
plot(v,reshape(quant_vmec(s_index,u_index,:),size(v)))

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
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(quant(:,:,nfp_v_index)-quant_vmec(:,:,nfp_v_index))))
colorbar; 
colormap jet
 caxis([-5 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('log scale difference btwn my %s and VMEC at nfp*phi=%f',quant_str,v(nfp_v_index)))


%% plot surfaces


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


%s,v
figure()
pcolor(v,s,reshape(quant(:,u_index,:),[dimS,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('My %s at u = %f',quant_str,u(u_index)))
%s,v
figure()
pcolor(v,s,reshape(quant_vmec(:,u_index,:),[dimS,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('VMEC %s at u = %f',quant_str,u(u_index)))



