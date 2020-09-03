%% Plot the force error
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid

% Plot just the v=0 plane
nfp_v_index = 1;

figure()
for i=1:size(R,1)
    plot(R(i,:,nfp_v_index),Z(i,:,nfp_v_index))
    hold on
end
xlabel('R (m)')
ylabel('Z (m)')
xlim([8.7,11.3])
ylim([-1.3,1.3])

title(sprintf('Flux Surfaces at nfp*phi=%f',v(nfp_v_index)))

figure()
contourf(R(:,:,nfp_v_index),Z(:,:,nfp_v_index),log10(F(:,:,nfp_v_index) ./ abs(presr(data.ns/2,:,nfp_v_index))))
c=colorbar; 
caxis([0 1]);
xlabel('R (m)')
ylabel('Z (m)')
axis equal
title(sprintf('Force error at nfp*phi=%f',v(nfp_v_index)))

