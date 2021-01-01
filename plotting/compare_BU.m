if exist('BU_vmec,','var') == false
BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
end
figure()

plot(data.phi(s_index:end),BU(s_index:end,u_index,v_nfp_index))
title(sprintf('B^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')

hold on
plot(data.phi(s_index:end),BU_vmec(s_index:end,u_index,v_nfp_index))
%title(sprintf('B^u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')
legend('My Calc','VMEC')
% vs u
figure()

plot(u,BU(s_index,:,v_nfp_index))
title(sprintf('B^u vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('u')
ylabel('B^u')

hold on
plot(u,BU_vmec(s_index,:,v_nfp_index))
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('u')
ylabel('B^u')
legend('My Calc','VMEC')

% vs v
figure()

plot(v,reshape(BU(s_index,u_index,:),size(v)))
title(sprintf('B^u vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^u')

hold on
plot(v,reshape(BU_vmec(s_index,u_index,:),size(v)))
%title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^u')
legend('My Calc','VMEC')

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BU_vmec(:,:,nfp_v_index))))
colorbar;
colormap jet
caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B^u from matlabVMEC at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BU(:,:,nfp_v_index))))
colorbar; 
colormap jet
caxis([-2 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('B^u I calculate at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BU(:,:,nfp_v_index)-BU_vmec(:,:,nfp_v_index))))
colorbar; 
colormap jet
 caxis([-5 0])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('Difference btwn my B^u and VMEC at nfp*phi=%f',v(nfp_v_index)))
