%% Plot the force error
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m
close all
% Plot just the v=0 plane
ds = s(2)-s(1);
du = u(2)-u(1);
dv = v(2)-v(1);

nfp_v_index = 1;
u_index = 2;

figure()
for i=1:size(R,1)
    plot(R(i,:,nfp_v_index),Z(i,:,nfp_v_index))
    hold on
    if i == 1
        plot(R(i,:,nfp_v_index),Z(i,:,nfp_v_index),'r.')
    end
end
xlabel('R (m)')
ylabel('Z (m)')
% xlim([8.7,11.3])
% ylim([-1.3,1.3])
axis equal


title(sprintf('Flux Surfaces at nfp*phi=%f',v(nfp_v_index)))

%should I divide by (presr * gss) or just presr? 
for nfp_v_index=1:1
figure()
% contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index) ./ abs(presr(floor(data.ns/4),:,nfp_v_index) .* sqrt(gss(floor(data.ns/4),:,nfp_v_index)))))
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
c=colorbar; 
% caxis([0 1]);

xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('Log Force error [N/m^3] at nfp*phi=%f',v(nfp_v_index)))
end

i_s = floor(data.ns/4);
if max(data.presf) > 1e-3
    

    figure()
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index) ./ abs(presr(i_s,:,nfp_v_index)))) %./ sqrt(gss(floor(data.ns/2),:,nfp_v_index)))))
    % contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
    c=colorbar; 
    caxis([-3 1.5]);
    hold on
    plot(R(8,:,nfp_v_index),Z(8,:,nfp_v_index),'r')
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    title(sprintf('Force error normalized by grad(p(s=%f)) at nfp*phi=%f, on log scale',s(floor(data.ns/4)),v(nfp_v_index)))

else % vacuum case, use ratio of gradient of magnetic pressure to magnetic field line tension as error metric
    B_tension = (BU.*BU_u + BV.*BU_v).*eu + (BU.*BV_u + BV.*BV_v).*ev;
    B_tension_U = (BU.*BU_u + BV.*BU_v);
    B_tension_V = (BU.*BV_u + BV.*BV_v);
    mag_B_tension = sqrt( (BU.*BU_u + BV.*BU_v).* dot(B_tension,eu,4) + (BU.*BV_u + BV.*BV_v).*dot(B_tension,ev,4)) ;
%     mag_B_tension = sqrt( ((BU.*BU_u + BV.*BU_v).^2) .* dot(eu,eu,4) + ((BU.*BV_u + BV.*BV_v).^2) .*dot(ev,ev,4) + ...
%                     2 .* (BU.*BU_u + BV.*BU_v).*(BU.*BV_u + BV.*BV_v) .*dot(eu,ev,4) ) ;
    

    grad_B_pres_s = (BU.*Bu_s + Bu.*BU_s + Bv.*BV_s + BV.*Bv_s)*0.5; % contravariant s component of the magnetic pressure gradient
    grad_B_pres_u = (BU.*Bu_u + Bu.*BU_u + Bv.*BV_u + BV.*Bv_u).*0.5; % contravariant u component of the magnetic pressure gradient
    grad_B_pres_v = (BU.*Bu_v + Bu.*BU_v + Bv.*BV_v + BV.*Bv_v).*0.5; % contravariant v component of the magnetic pressure gradient
    grad_B_pres = (grad_B_pres_s .* eS + grad_B_pres_u .* eU + grad_B_pres_v .* eV);

    mag_grad_B_pres = sqrt((grad_B_pres_s .* dot(grad_B_pres,eS,4) + grad_B_pres_u .* dot(grad_B_pres,eU,4) + grad_B_pres_v .* dot(grad_B_pres,eV,4) ) );

%     vac_F = (grad_B_pres - B_tension)./mu0;
% magnitude of vac_F is mag(grad_B_pres) + mag(B_tension) - 2 * B_pressure
% dot B_tension
    mag_vac_F = sqrt(mag_grad_B_pres.^2 + mag_B_tension.^2 - 2.*dot(grad_B_pres,B_tension,4))./mu0;
    
    figure()
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),(abs(mag_grad_B_pres(:,:,nfp_v_index) - mag_B_tension(:,:,nfp_v_index))./mag_grad_B_pres(:,:,nfp_v_index)),20)
    % contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
    hold on
    c=colorbar; 
    caxis([0 1]);
    hold on
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    title(sprintf('(|\\nabla B^2 /2| - |(B \\cdot \\nabla)B|) / |\\nabla B^2 /2| ','Interpreter','latex'))

    figure()
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(mu0.*F(:,:,nfp_v_index))./nanmean(mag_grad_B_pres(:,1,nfp_v_index)),20)
    % contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
    hold on
    c=colorbar; 
    caxis([0 1]);
    hold on
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    title(sprintf('(|F|) / |\\nabla B^2 /2| ','Interpreter','latex'))

    figure()
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(mu0.*F(:,:,nfp_v_index))./nanmean(mag_B_tension(:,1,nfp_v_index)),20)
    % contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
    hold on
    c=colorbar; 
    caxis([0 1]);
    hold on
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    title(sprintf('(|F|) / |(B \\cdot \\nabla)B| ','Interpreter','latex'))
    
    
    figure()
    pcolor(u,s,mag_grad_B_pres(:,:,1)./mu0)
    colormap jet
%     caxis([min(abs(F(:,:,nfp_v_index)),[],'all') 1e-3*max(abs(F(:,:,nfp_v_index)),[],'all')])
    caxis([0,1e5])
    colorbar
    xlabel('u')
    ylabel('s')
    title(sprintf(' |\\nabla B^2 /2mu0| ','Interpreter','latex'))
    
    figure()
%     pcolor(u,s,repmat(avg_B_tension./mu0,[1 dimU]))
    pcolor(u,s,mag_B_tension(:,:,1)./mu0)
    colormap jet
%     caxis([min(abs(F(:,:,nfp_v_index)),[],'all') 1e-3*max(abs(F(:,:,nfp_v_index)),[],'all')])
    caxis([0,1e5])
    colorbar
    xlabel('u')
    ylabel('s')
    title(sprintf(' |(B \\cdot \\nabla)B| ','Interpreter','latex'))
    
    ratio = mag_B_tension ./ mag_grad_B_pres;
    
    figure()
%     pcolor(u,s,repmat(avg_B_tension./mu0,[1 dimU]))
    pcolor(u,s(10:end), abs(ratio(10:end,:,1)))
    colormap jet
%     caxis([min(abs(F(:,:,nfp_v_index)),[],'all') 1e-3*max(abs(F(:,:,nfp_v_index)),[],'all')])
%     caxis([0,1e2])
    colorbar
    xlabel('u')
    ylabel('s')
    title(sprintf('|(B \\cdot \\nabla)B| / |\\nabla B^2 /2| ','Interpreter','latex'))
    
    figure()

    plot(data.phi(1:end),mag_B_tension(1:end,u_index,1))

    hold on
    plot(data.phi(1:end),mag_grad_B_pres(1:end,u_index,1))
    hold on
    plot(data.phi(1:end),ratio(1:end,u_index,1))
    xlabel('s')
    ylabel('Magnitude')
    legend('B tension','B Pressure','B ten / B pres')
    
    
    figure()

    plot(data.phi(1:end),F(1:end,u_index,1))

    hold on
    plot(data.phi(1:end),mag_grad_B_pres(1:end,u_index,1)./mu0)
    xlabel('s')
    ylabel('Magnitude')
    legend('F','B Pressure')
    
    figure()

    plot(data.phi(1:end),F(1:end,u_index,1))

    hold on
    plot(data.phi(1:end),mag_vac_F(1:end,u_index,1))
    hold on
    plot(data.phi(1:end),(mag_grad_B_pres(1:end,u_index,1) - mag_B_tension(1:end,u_index,1))./mu0 )
    xlabel('s')
    ylabel('Magnitude')
    legend('F','Vac F','B Pressure - B tension')

end


