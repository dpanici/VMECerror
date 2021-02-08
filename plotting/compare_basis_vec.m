eseS = dot(es,eS,4); 
eueU = dot(eu,eU,4);
eveV = dot(ev,eV,4);

eseU = dot(es,eU,4); 
eseV = dot(es,eV,4);

eueS = dot(eu,eS,4);
eueV = dot(eu,eV,4);

eveS = dot(ev,eS,4);
eveU = dot(ev,eU,4);


eseS_is_one = ismembertol(eseS,ones(size(eseS)));
eueU_is_one = ismembertol(eueU,ones(size(eueU)));
eveV_is_one = ismembertol(eveV,ones(size(eveV)));

%% plot eseS


%% plot surfaces
% quant = eseS;
% quant_str = 'es*eS';
% quant_vmec = ones(size(quant));
% 
% % s,u
% figure()
% pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('u')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
% %s,v
% figure()
% pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))
% 
% %u,v
% figure()
% pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('u')
% title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))
% 
% %% plot surfaces eueU
% quant = eueU;
% quant_str = 'eu*eU';
% quant_vmec = ones(size(quant));
% 
% % s,u
% figure()
% pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('u')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
% %s,v
% figure()
% pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))
% 
% %u,v
% figure()
% pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('u')
% title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))
% 
% %% plot surfaces eveV
% quant = eveV;
% quant_str = 'ev*eV';
% quant_vmec = ones(size(quant));
% 
% % s,u
% figure()
% pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('u')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
% %s,v
% figure()
% pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('s')
% title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))
% 
% %u,v
% figure()
% pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
% colormap jet
% caxis([0,0.005])
% colorbar
% xlabel('v')
% ylabel('u')
% title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))


%% check cross terms

%% plot eseU


% plot surfaces
quant = eseU;
quant_str = 'es*eU';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% plot surfaces eseV
quant = eseV;
quant_str = 'es*eV';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% plot surfaces eueS
quant = eueS;
quant_str = 'eu*eS';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% plot surfaces eueV
quant = eueV;
quant_str = 'eu*eV';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% plot surfaces eveS
quant = eveS;
quant_str = 'ev*eS';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

%% plot surfaces eveU
quant = eveU;
quant_str = 'ev*eU';
quant_vmec = zeros(size(quant));

% s,u
figure()
pcolor(u,s,abs(quant(:,:,nfp_v_index) - quant_vmec(:,:,nfp_v_index)))
colormap jet
caxis([0,0.005])
colorbar
xlabel('u')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at v = %f',quant_str,v(v_nfp_index)))
%s,v
figure()
pcolor(v,s,reshape(abs(quant(:,u_index,:) - quant_vmec(:,u_index,:)),[dimS,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('s')
title(sprintf('Abs Difference in my %s and VMEC at u = %f',quant_str,u(u_index)))

%u,v
figure()
pcolor(v,u,reshape(abs(quant(s_index,:,:) - quant_vmec(s_index,:,:)),[dimU,dimV]))
colormap jet
caxis([0,0.005])
colorbar
xlabel('v')
ylabel('u')
title(sprintf('Abs Difference in my %s and VMEC at s = %f',quant_str,data.phi(s_index)))

