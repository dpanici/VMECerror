%% Plot various quantities for comparison with matlabVMEC
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m

close all


% nfp_index=0;
u_index= 1; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=3; % index of s to plot quantities at (that arent plotted versus s)

g_vmec = eval_series_nyq(suvgrid,data.gmnc,data,'c');
BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
% JU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
magB_vmec = sqrt((BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4));
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4));

%% Surface plot ||B||
figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB(:,:,v_nfp_index))
colormap jet
colorbar
caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('||B|| I calculate at nfp*phi=%f',v(v_nfp_index)))

%% Plot ||B||
figure()
subplot(1,2,1)
plot(data.phi,magB(:,u_index,v_nfp_index))
title(sprintf('||B|| vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('||B||')
ylim([0,1.1*max(magB(:,u_index,v_nfp_index))])

subplot(1,2,2)
plot(data.phi,magB_vmec(:,u_index,v_nfp_index))
title(sprintf('||B|| vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('||B||')
ylim([0,1.1*max(magB(:,u_index,v_nfp_index))])

figure()
subplot(1,2,1)
plot(v,reshape(magB(s_index,u_index,:),size(v)))
title(sprintf('||B|| vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('||B||')
ylim([0,1.1*max(magB(s_index,u_index,:))])

subplot(1,2,2)
plot(v,reshape(magB_vmec(s_index,u_index,:),size(v)))
title(sprintf('||B|| vmec vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('||B||')
ylim([0,1.1*max(magB(s_index,u_index,:))])


%% Plot B^u
figure()
subplot(1,2,1)
plot(data.phi(3:end),BU(3:end,u_index,v_nfp_index))
title(sprintf('B^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')

subplot(1,2,2)
plot(data.phi(3:end),BU_vmec(3:end,u_index,v_nfp_index))
title(sprintf('B^u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B^u')

figure()
subplot(1,2,1)
plot(v,reshape(BU(s_index,u_index,:),size(v)))
title(sprintf('B^u vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^u')

subplot(1,2,2)
plot(v,reshape(BU_vmec(s_index,u_index,:),size(v)))
title(sprintf('B^u vmec vs v at u=%f, s=%f',u(u_index),s(s_index)))
xlabel('v')
ylabel('B^u')

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


%% Plot B^v
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
title(sprintf('Difference btwn my B^u and VMEC at nfp*phi=%f',v(nfp_v_index)))


% Plot J^u, J^v (Units?)
JU_vmec = cfunct(u,v,data.currumnc,data.xm_nyq,data.xn_nyq);

figure()
subplot(1,2,1)
plot(data.phi(3:end),JU(3:end,u_index,v_nfp_index))
title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u')

subplot(1,2,2)
plot(data.phi(3:end),JU_vmec(3:end,u_index,v_nfp_index))
title(sprintf('J^u VMEC vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u VMEC')

figure()
plot(data.phi(3:end),JU(3:end,u_index,v_nfp_index))
title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u')
ylim([0,1e6])
hold on
plot(data.phi(3:end),JU_vmec(3:end,u_index,v_nfp_index))
legend('mine','VMEC')


figure()
plot(data.phi(3:end),JV(3:end,u_index,v_nfp_index))
title(sprintf('J^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^v')

% Plot g

figure()
subplot(1,2,1)
plot(data.phi(10:end),g(10:end,u_index,v_nfp_index))
title(sprintf('g vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')

subplot(1,2,2)
plot(data.phi(10:end),g_vmec(10:end,u_index,v_nfp_index))
title(sprintf('g vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('g')

figure()
subplot(1,2,1)
plot(v,reshape(g(s_index,u_index,:),size(v)))
title(sprintf('g vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('g')

subplot(1,2,2)
plot(v,reshape(g_vmec(s_index,u_index,:),size(v)))
title(sprintf('g vmec vs v at s=%f, u=%f',s(s_index),u(u_index)))
xlabel('v')
ylabel('g')



figure()
plot(data.phi,L(:,1,v_nfp_index))
title(sprintf('L vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('L')


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
