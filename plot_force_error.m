%% Plot the force error
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid
close all
% Plot just the v=0 plane

figure()
for i=1:data.ns
    plot(R(i,:,1),Z(i,:,1))
    hold on
end
xlabel('R (m)')
ylabel('Z (m)')

title('Force error at phi=0')

figure()
contourf(R(:,:,1),Z(:,:,1),log10(F(:,:,1)))
c=colorbar; 
xlabel('R (m)')
ylabel('Z (m)')

title('Force error at phi=0')

