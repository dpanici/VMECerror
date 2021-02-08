if exist('magB_vmec','var') == false
% magB_vmec = sqrt((BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4));
magB_vmec = eval_series_nyq(suvgrid,data.bmnc,data,'c');
end
if exist('magB','var') == false   
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4));% + (BU.*BV).*dot(eu,ev,4) + (BU.*BV).*dot(ev,eu,4));
end

figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB(:,:,v_nfp_index))
% caxis([0.15 0.3])
colormap jet
colorbar
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('||B|| I calculate at nfp*phi=%f',v(v_nfp_index)))

%vs u
figure()

plot(u,magB(s_index,:,v_nfp_index))
title(sprintf('||B|| vs u at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('u')
ylabel('||B||')
% ylim([0,1.1*max(magB(s_index,:,v_nfp_index))])

hold on
plot(u,magB_vmec(s_index,:,v_nfp_index),'--')
%title(sprintf('||B|| vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))

legend('My Calc','VMEC')

%vs s

figure()

plot(data.phi,magB(:,u_index,v_nfp_index))
title(sprintf('||B|| vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('||B||')
ylim([0,1.1*max(magB(:,u_index,v_nfp_index))])

hold on
plot(data.phi,magB_vmec(:,u_index,v_nfp_index),'--')
%title(sprintf('||B|| vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))

legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(magB(s_index,u_index,:),size(v)))
title(sprintf('||B|| vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('||B||')
ylim([0,1.1*max(magB(s_index,u_index,:))])

hold on
plot(v,reshape(magB_vmec(s_index,u_index,:),size(v)))
% title(sprintf('||B|| vmec vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('||B||')
ylim([0,1.1*max(magB(s_index,u_index,:))])
legend('My Calc', 'VMEC')
%% plot surfaces

quant = magB;
quant_str = '|B|';
quant_vmec = magB_vmec;

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

%u,v
figure()
pcolor(v,u,reshape(quant(s_index,:,:),[dimU,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('My %s at s = %f',quant_str,data.phi(s_index)))

%u,v
figure()
pcolor(v,u,reshape(quant_vmec(s_index,:,:),[dimU,dimV]))
colormap jet
% caxis([0,0.01])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('VMEC %s at s = %f',quant_str,data.phi(s_index)))

quant = magB;
quant_str = '|B|';
quant_vmec = magB_vmec;

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


