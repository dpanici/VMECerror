if exist('JU_vmec','var') == false
JU_vmec = eval_series_nyq(suvgrid,data.currumnc,data,'c');
end
if exist('JV_vmec','var') == false
JV_vmec = eval_series_nyq(suvgrid,data.currvmnc,data,'c');
end
if exist('g_vmec','var') == false
g_vmec = eval_series_nyq(suvgrid,data.gmnc,data,'c');
end
plot_rational_iotas=false;

% plot J * g
figure()

plot(data.phi(s_index:end),JU(s_index:end,u_index,v_nfp_index).*g(s_index:end,u_index,v_nfp_index))
title(sprintf('sqrt(g)*J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('sqrt(g)*J^u')

hold on
plot(data.phi(s_index:end),JU_vmec(s_index:end,u_index,v_nfp_index))
title(sprintf('sqrt(g)*J^u VMEC vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('sqrt(g)*J^u')
legend('My Calc','VMEC')
% ylim([min(JU_vmec(s_index:end,u_index,v_nfp_index)),max(JU_vmec(s_index:end,u_index,v_nfp_index))]);
if exist('rat_s','var') == false
plot_rationals
end
% 
% for i=unique(rat_s)
%     xline(i)
%     hold on
% end
% % plot J
figure()

plot(data.phi(s_index:end),JU(s_index:end,u_index,v_nfp_index))
title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))


hold on
plot(data.phi(s_index:end),JU_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index))
% title(sprintf('J^u VMEC vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u')
legend('My Calc','VMEC')
ylim([min(JU_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index)),max(JU_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index))]);
% ylim([-1.5e4,3e5])
if plot_rational_iotas
for i=unique(rat_s)
    xline(i)
    hold on
end
end
% plot J * g
figure()
plot(data.phi(s_index:end),JV(s_index:end,u_index,v_nfp_index).*g(s_index:end,u_index,v_nfp_index))

title(sprintf('sqrt(g)*J^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('sqrt(g)*J^v')
hold on
plot(data.phi(s_index:end),JV_vmec(s_index:end,u_index,v_nfp_index))
legend('mine','VMEC')
ylim([min(JV_vmec(s_index:end,u_index,v_nfp_index)),max(JV_vmec(s_index:end,u_index,v_nfp_index))]);

legend('My Calc','VMEC')
% for i=unique(rat_s)
%     xline(i,'HandleVisibility','off')
%     hold on
% end
%plot J
figure()
% plot(data.phi(s_index:end),JV(s_index:end,u_index,v_nfp_index).*g(s_index:end,u_index,v_nfp_index))
plot(data.phi(s_index:end),JV(s_index:end,u_index,v_nfp_index))
title(sprintf('J^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^v')
hold on
plot(data.phi(s_index:end),JV_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index))
legend('My Calc','matlabVMEC')
% hold on
% plot(data.phi(s_index:end),JV_smoothed(s_index:end,u_index,v_nfp_index),'DisplayName','Smooth Spline')
% legend('finite difference','VMEC output','Smooth Spline')
ylim([min(JV_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index)),max(JV_vmec(s_index:end,u_index,v_nfp_index)./g_vmec(s_index:end,u_index,v_nfp_index))]);
% ylim([2.45e4,3.5e4])
if plot_rational_iotas
for i=unique(rat_s)
    xline(i)
    hold on
end
end
% plot our JS
%plot J
figure()
% plot(data.phi(s_index:end),JV(s_index:end,u_index,v_nfp_index).*g(s_index:end,u_index,v_nfp_index))
plot(sqrt(data.phi(s_index:end)),JS(s_index:end,u_index,v_nfp_index))
title(sprintf('J^s vs rho at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('rho')
ylabel('J^s')
v_ind=1;
%% saving so can read in with python to compare with DESC
% writematrix(R(:,:,1),sprintf('../VMEC_W7X_s2048_R_JS_phi_%1.3f.csv',v(v_ind)))
% writematrix(Z(:,:,1),sprintf('../VMEC_W7X_s2048_Z_JS_phi_%1.3f.csv',v(v_ind)))
% writematrix(JS(:,:,1),sprintf('../VMEC_W7X_s2048_JS_phi_%1.3f.csv',v(v_ind)))




%% 2D
clims = [0,1e4];
clims_ratio = [0,0.3];
%% JS
figure()
h=pcolor(u,s,JS(:,:,nfp_v_index))
set(h, 'EdgeColor', 'none');
colormap jet
caxis([-4e3,4e3])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('My JS at v = %f',v(v_nfp_index)))

figure()
hh=pcolor(v,u,reshape(JS(200,:,:),[dimU,dimV]))
set(hh, 'EdgeColor', 'none');
colormap jet
% caxis(clims)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('My calc JS at s = %f',data.phi(200)))


%% JU
% s,u
figure()
pcolor(u,s,abs(g(:,:,nfp_v_index).*JU(:,:,nfp_v_index) - JU_vmec(:,:,nfp_v_index)))
colormap jet
caxis(clims)
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my g*JU and VMEC g*JU at v = %f',v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(g(:,u_index,:).*JU(:,u_index,:) - JU_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis(clims)
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my g*JU and VMEC g*JU at u = %f',u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(g(s_index,:,:).*JU(s_index,:,:) - JU_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis(clims)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my g*JU and VMEC g*JU at s = %f',data.phi(s_index)))

%ratios
% s,u
figure()
pcolor(u,s,abs((g(:,:,nfp_v_index).*JU(:,:,nfp_v_index) - JU_vmec(:,:,nfp_v_index)) ./ JU_vmec(:,:,nfp_v_index)))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('u')
ylabel('s')
title(sprintf('pct diff my g*JU and VMEC g*JU at v = %f',v(v_nfp_index)))

%s,v
figure()
pcolor(v,s,reshape(abs((g(:,u_index,:).*JU(:,u_index,:) - JU_vmec(:,u_index,:)) ./ JU_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('s')
title(sprintf('pct diff my g*JU and VMEC g*JU at u = %f',u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs((g(s_index,:,:).*JU(s_index,:,:) - JU_vmec(s_index,:,:)) ./ JU_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('pct diff my g*JU and VMEC g*JU at s = %f',data.phi(s_index)))



%% JV
% s,u
figure()
pcolor(u,s,abs(g(:,:,nfp_v_index).*JV(:,:,nfp_v_index) - JV_vmec(:,:,nfp_v_index)))
colormap jet
caxis(clims)
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my g*JV and VMEC g*JV at v = %f',v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(g(:,u_index,:).*JV(:,u_index,:) - JV_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis(clims)
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my g*JV and VMEC g*JV at u = %f',u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(g(s_index,:,:).*JV(s_index,:,:) - JV_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis(clims)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my g*JV and VMEC g*JV at s = %f',data.phi(s_index)))
% ratios
% s,u
figure()
pcolor(u,s,abs((g(:,:,nfp_v_index).*JV(:,:,nfp_v_index) - JV_vmec(:,:,nfp_v_index)) ./ JV_vmec(:,:,nfp_v_index)))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('u')
ylabel('s')
title(sprintf('pct diff my g*JV and VMEC g*JV at v = %f',v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs((g(:,u_index,:).*JV(:,u_index,:) - JV_vmec(:,u_index,:)) ./ JV_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('s')
title(sprintf('pct diff my g*JV and VMEC g*JV at u = %f',u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs((g(s_index,:,:).*JV(s_index,:,:) - JV_vmec(s_index,:,:)) ./ JV_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis(clims_ratio)
colorbar
xlabel('v')
ylabel('u')
title(sprintf('pct diff my g*JV and VMEC g*JV at s = %f',data.phi(s_index)))
