if exist('BV_vmec,','var') == false
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
end
if exist('R_vmec,','var') == false
R_vmec = eval_series(suvgrid,data.rmnc,data,'c');
end
if exist('B_Phi_vmec,','var') == false
B_Phi_vmec = BV_vmec .* R_vmec;
end
if exist('B_Phi,','var') == false
B_Phi = BV .* R;
end



figure()

plot(data.phi(s_index:end),B_Phi(s_index:end,u_index,v_nfp_index))
title(sprintf('Bphi vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('Bphi')

hold on
plot(data.phi(s_index:end),B_Phi_vmec(s_index:end,u_index,v_nfp_index),'k--')
%title(sprintf('B^u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
legend('My Calc','VMEC')

% vs u
figure()

plot(u,B_Phi(s_index,:,v_nfp_index))
title(sprintf('B_Phi vs u at v=%f, s=%f',v(v_nfp_index),s(s_index)))
xlabel('u')
ylabel('B_Phi')

hold on
plot(u,B_Phi_vmec(s_index,:,v_nfp_index),'k--')
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(B_Phi(s_index,u_index,:),size(v)))
title(sprintf('B_Phi vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B_Phi')

hold on
plot(v,reshape(B_Phi_vmec(s_index,u_index,:),size(v)),'k--')
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
legend('My Calc','VMEC')

% 2D plots
figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(B_Phi_vmec(:,:,nfp_v_index)))
colorbar;
colormap jet
% caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B_Phi from matlabVMEC at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),B_Phi(:,:,nfp_v_index))
colorbar; 
colormap jet
% caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B_Phi I calculate at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(B_Phi(:,:,nfp_v_index)-B_Phi_vmec(:,:,nfp_v_index)))
colorbar; 
colormap jet
%  caxis([-5 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('Abs Val Difference btwn my B_Phi and VMEC at nfp*phi=%f',v(nfp_v_index)))

%% plot surfaces
quant = B_Phi;
quant_str = 'B_Phi';
quant_vmec = B_Phi_vmec;

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
