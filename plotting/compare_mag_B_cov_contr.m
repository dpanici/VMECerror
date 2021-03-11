% run after force error script has ran
if exist('Bs_vmec','var')== false
Bs_vmec = eval_series_nyq(suvgrid,data.bsubsmns,data,'s');
end
if exist('Bu_vmec','var')== false
Bu_vmec = eval_series_nyq(suvgrid,data.bsubumnc,data,'c');
end
if exist('Bv_vmec','var')==false
Bv_vmec = eval_series_nyq(suvgrid,data.bsubvmnc,data,'c');
end
if exist('magB_cov','var') == false
% magB_cov = sqrt((Bs.^2).*dot(eS,eS,4) + (Bu.^2).*dot(eU,eU,4) + (Bv.^2).*dot(eV,eV,4));
magB_cov = sqrt((Bs.^2).*dot(eS,eS,4) + (Bu.^2).*dot(eU,eU,4) + (Bv.^2).*dot(eV,eV,4));
end
if exist('magB','var') == false
% magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4))    
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4) + (BU.*BV).*dot(eu,ev,4) + (BU.*BV).*dot(ev,eu,4));
end
if exist('magB_vmec','var') == false
% magB_vmec = sqrt((BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4));
% magB_vmec = sqrt((BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4) + (BU_vmec.*BV_vmec).*dot(eu,ev,4) + (BU_vmec.*BV_vmec).*dot(ev,eu,4));
magB_vmec = eval_series_nyq(suvgrid,data.bmnc,data,'c');
end
if exist('magB_cov_vmec','var') == false
magB_cov_vmec = sqrt((Bs_vmec.^2).*dot(eS,eS,4) + (Bu_vmec.^2).*dot(eU,eU,4) + (Bv_vmec.^2).*dot(eV,eV,4));
end


% compare my cov B and contr B mags
figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB(:,:,v_nfp_index) - magB_cov(:,:,v_nfp_index) )
colormap jet
colorbar
% caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('diff in ||B|| from cov B and contrav B  at nfp*phi=%f',v(v_nfp_index)))

% compare my cov B mag to VMEC cov B mag

figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB_cov_vmec(:,:,v_nfp_index) - magB_cov(:,:,v_nfp_index) )
colormap jet
colorbar
% caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('diff in ||B|| from VMEC cov B and my cov B  at nfp*phi=%f',v(v_nfp_index)))

% compare VMEC cov B mag to VMEC contr B mag
figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB_cov_vmec(:,:,v_nfp_index) - magB_vmec(:,:,v_nfp_index) )
colormap jet
colorbar
% caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('diff in ||B|| from VMEC cov B and VMEC contrav B  at nfp*phi=%f',v(v_nfp_index)))

% mag B with VMEC cov B
figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB_cov_vmec(:,:,v_nfp_index) )
colormap jet
colorbar
% caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('||B|| from VMEC cov B at nfp*phi=%f',v(v_nfp_index)))

% mag B with VMEC contr B
figure()
contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),magB_vmec(:,:,v_nfp_index) )
colormap jet
colorbar
% caxis([0.2 0.5])
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('||B|| from VMEC contr B at nfp*phi=%f',v(v_nfp_index)))


quant = magB_cov;
quant_str = '|B| from covB';
quant_vmec = magB_cov_vmec;

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

clims_ratio = [0,0.25];
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

