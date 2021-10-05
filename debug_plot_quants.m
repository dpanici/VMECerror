%% Plot various quantities for comparison with matlabVMEC
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m

close all

u_index= 9; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=1; % index of s to plot quantities at (that arent plotted versus s)


%% R,Z derivs
% run('plotting/compare_R_Z_v_derivs.m')

%% basis vector dot products

% run('plotting/compare_basis_vec')

%% Plot R derivs
% run('plotting/compare_R_s')
% run('plotting/compare_Z_s')


%% surface plot check that my cov B and contrav B give same result for |B|

% run('plotting/compare_mag_B_cov_contr')

%% plot ||B||

% run('plotting/compare_mag_B')

%% Plot B^2 derivatives
% run('plotting/check_grad_Bsq')

%% Plot B^u

% run('plotting/compare_BU')

%% Plot B^v

% run('plotting/compare_BV')

%% Plot J^u, J^v (Units?)
run('plotting/compare_J')
% 
%% Plot Cyl B

% run('plotting/compare_BR')
% run('plotting/compare_BPhi')
% run('plotting/compare_BZ')

%% Plot g

% run('plotting/compare_g')

%% Plot Lambda

% run('plotting/compare_L')

%% plot covariant B components 

% run('plotting/compare_cov_Bs')
% run('plotting/compare_cov_Bu')
% run('plotting/compare_cov_Bv')

%% plot metric tensor
% run('plotting/plot_metric_tensor')

%% compare covariant B derivs I calculate here to ones I get with contravariant B

% derivplot(Bs,Bs_u,u)
% derivplot(Bs,Bs_v,v)
% 
% derivplot(Bu,Bu_s,s)
% derivplot(Bu,Bu_v,v)
% 
% derivplot(Bv,Bv_u,u)
% derivplot(Bv,Bv_s,s)


%% plot current density

% J = JS.*es + JU.*eu + JV.*ev;
% figure()
% contourf(R(:,:,v_nfp_index),Z(:,:,v_nfp_index),abs(J(:,:,v_nfp_index)))
% colormap jet
% colorbar
% caxis([1e5 6e5])
% xlabel('R (m)')
% ylabel('Z (m)')
% axis equal
% title(sprintf('Current Density I calculate at nfp*phi=%f',v(v_nfp_index)))
