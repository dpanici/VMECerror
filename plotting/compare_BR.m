if exist('BU_vmec,','var') == false
BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
end
if exist('BV_vmec,','var') == false
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
end
if exist('R_u_vmec,','var') == false
R_u_vmec = eval_series(suvgrid,data.rumns,data,'s');
end
if exist('R_v_vmec,','var') == false
R_v_vmec = eval_series(suvgrid,data.rvmns,data,'s');
end
if exist('B_R_vmec,','var') == false
B_R_vmec = BU_vmec .* R_u_vmec + BV_vmec .* R_v_vmec;
end
if exist('B_R,','var') == false
B_R = BU .* R_u + BV .* R_v;
end



figure()

plot(data.phi(s_index:end),B_R(s_index:end,u_index,v_nfp_index))
title(sprintf('BR vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('BR')

hold on
plot(data.phi(s_index:end),B_R_vmec(s_index:end,u_index,v_nfp_index),'k--')
%title(sprintf('B^u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
legend('My Calc','VMEC')

% vs u
figure()

plot(u,B_R(s_index,:,v_nfp_index))
title(sprintf('B_R vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('B_R')

hold on
plot(u,B_R_vmec(s_index,:,v_nfp_index),'k--')
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(B_R(s_index,u_index,:),size(v)))
title(sprintf('B_R vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B_R')

hold on
plot(v,reshape(B_R_vmec(s_index,u_index,:),size(v)),'k--')
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
legend('My Calc','VMEC')

% 2D plots
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),B_R_vmec(:,:,nfp_v_index))
colorbar;
colormap jet
% caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B_R from matlabVMEC at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),B_R(:,:,nfp_v_index))
colorbar; 
colormap jet
% caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B_R I calculate at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(B_R(:,:,nfp_v_index)-B_R_vmec(:,:,nfp_v_index)))
colorbar; 
colormap jet
%  caxis([-5 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('Abs Val Difference btwn my B_R and VMEC at nfp*phi=%f',v(nfp_v_index)))

%% plot surfaces
quant = B_R;
quant_str = 'B_R';
quant_vmec = B_R_vmec;

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
