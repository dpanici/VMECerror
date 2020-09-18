%% Plot various quantities for comparison with matlabVMEC
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m

close all


% nfp_index=0;
u_index= 1; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at

BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
magB_vmec = sqrt((BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4));
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4));
% figure()
% contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),magB(:,:,nfp_v_index))
% colormap jet
% colorbar
% caxis([0 0.5])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('||B|| I calculate at nfp*phi=%f',v(nfp_v_index)))


figure()
subplot(1,2,1)
plot(data.phi,magB(:,u_index,v_nfp_index))
title(sprintf('||B|| vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('||B||')
ylim([0,max(magB(:,u_index,v_nfp_index))])

subplot(1,2,2)
plot(data.phi,magB_vmec(:,u_index,v_nfp_index))
title(sprintf('||B|| vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('||B||')
ylim([0,max(magB(:,u_index,v_nfp_index))])

ev_ev = dot(ev,ev,4);

figure()
plot(data.phi,ev_ev(:,u_index,v_nfp_index))
title(sprintf('dot(ev,ev) vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('ev dot ev')


% figure()
% contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BU_vmec(:,:,nfp_v_index))))
% colorbar;
% colormap jet
% caxis([-2 0])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('B^u from matlabVMEC at nfp*phi=%f',v(nfp_v_index)))

% figure()
% contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(abs(BU(:,:,nfp_v_index))))
% colorbar; 
% colormap jet
% caxis([-2 0])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('B^u I calculate at nfp*phi=%f',v(nfp_v_index)))

figure()
subplot(1,2,1)
plot(data.phi(3:end),BU(3:end,u_index,v_nfp_index))
title(sprintf('B^u vs s at u=%f, nfp*phi=%f\n mine is almost exactly 2*pi larger than matlabVMEC',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')

subplot(1,2,2)
plot(data.phi(3:end),BU_vmec(3:end,u_index,v_nfp_index))
title(sprintf('B^u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')


figure()
plot(data.phi(2:end),g(2:end,u_index,v_nfp_index))
title(sprintf('g vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')



figure()
subplot(1,2,1)
plot(data.phi(3:end),BV(3:end,u_index,v_nfp_index))
title(sprintf('B^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^v')


subplot(1,2,2)
plot(data.phi(3:end),BV_vmec(3:end,u_index,v_nfp_index))
title(sprintf('B^v vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^v')
% figure()
% plot(u,BV(2,:,v_nfp_index))
% title(sprintf('B^v vs u at s=%f, nfp*phi=%f',s(2),v(v_nfp_index)))
% xlabel('u')
% ylabel('B^v')


% figure()
% plot(data.phi,L(:,1,v_nfp_index))
% title(sprintf('L vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
% xlabel('s')
% ylabel('L')


% figure()
% contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),log10(abs(BV_vmec(:,:,v_nfp_index))))
% c=colorbar; 
% caxis([-1 2])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('B^v from matlabVMEC at nfp*phi=%f',v(v_nfp_index)))

% figure()
% contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),log10(abs(BV(:,:,v_nfp_index))))
% c=colorbar; 
% caxis([-1,2])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('B^v I calculate at nfp*phi=%f',v(v_nfp_index)))



% 
% figure()
% contourf(R(:,:,nfp_v_index),Z(:,:,v_nfp_index),log10(magB_vmec(:,:,v_nfp_index)))
% c=colorbar; 
% colormap jet
% caxis([-1 2])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('||B|| from matlabvmec at nfp*phi=%f',v(v_nfp_index)))
