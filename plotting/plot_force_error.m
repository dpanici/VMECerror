%% Plot the force error
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m
close all
% Plot just the v=0 plane
% ds = s(2)-s(1);
% du = u(2)-u(1);
% dv = v(2)-v(1);

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

axis equal


title(sprintf('Flux Surfaces at nfp*phi=%f',v(nfp_v_index)))

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
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index) ./ abs(presr(i_s,:,nfp_v_index))))
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

    

    grad_B_pres_s = (BU.*Bu_s + Bu.*BU_s + Bv.*BV_s + BV.*Bv_s)*0.5; % contravariant s component of the magnetic pressure gradient
    grad_B_pres_u = (BU.*Bu_u + Bu.*BU_u + Bv.*BV_u + BV.*Bv_u).*0.5; % contravariant u component of the magnetic pressure gradient
    grad_B_pres_v = (BU.*Bu_v + Bu.*BU_v + Bv.*BV_v + BV.*Bv_v).*0.5; % contravariant v component of the magnetic pressure gradient
    grad_B_pres = (grad_B_pres_s .* eS + grad_B_pres_u .* eU + grad_B_pres_v .* eV);

    mag_grad_B_pres = mu0.*sqrt((grad_B_pres_s .* dot(grad_B_pres,eS,4) + grad_B_pres_u .* dot(grad_B_pres,eU,4) + grad_B_pres_v .* dot(grad_B_pres,eV,4) ) );


    mag_vac_F = sqrt(mag_grad_B_pres.^2 + mag_B_tension.^2 - 2.*dot(grad_B_pres,B_tension,4))./mu0;
    

    figure()
    contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),abs(F(:,:,nfp_v_index))./nanmean(mag_grad_B_pres(:,1,nfp_v_index)),20)
    % contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index)))
    hold on
    c=colorbar; 
    hold on
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    title(sprintf('(|F|) / |\\nabla B^2 /2| ','Interpreter','latex'))



end


