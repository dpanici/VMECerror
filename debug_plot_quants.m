%% Plot various quantities for comparison with matlabVMEC
% this script assumes R,Z, and F have been defined and calculated already
% on the suvgrid after having run force_error.m

close all

u_index= 5; % index of u to plot quantities at
v_nfp_index=1; % index of v to plot quantities at
nfp_v_index = v_nfp_index;
s_index=15; % index of s to plot quantities at (that arent plotted versus s)


%% scatter plot for basis vector dot products

% run('plotting/compare_basis_vec')

%% Plot R derivs
% run('plotting/compare_R_s')
% run('plotting/compare_Z_s')


%% surface plot check thatZmy cov B and contrav B give same result for |B|

% run('plotting/compare_mag_B_cov_contr')

%% plot ||B||

% run('plotting/compare_mag_B')

%% Plot B^u

% run('plotting/compare_BU')

%% Plot B^v

% run('plotting/compare_BV')

%% Plot J^u, J^v (Units?)
% run('plotting/compare_J')
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
run('plotting/compare_cov_Bu')
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


%% try plotting my covariant J versus vmec
% VMEC USES J^X * g, which is why I was off before
% Ju = JS .* g_su + JU.*g_uu + JV .* g_vu;
% Jv = JS .* g_sv + JU.*g_uv + JV .* g_vv;
% 
% figure()
% 
% plot(data.phi(s_index:end),Ju(s_index:end,u_index,v_nfp_index))
% title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
% xlabel('s')
% ylabel('J^u')
% 
% hold on
% plot(data.phi(s_index:end),JU_vmec(s_index:end,u_index,v_nfp_index))
% title(sprintf('J^u VMEC vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
% xlabel('s')
% ylabel('J^u VMEC')
% legend('My Calc J_u','VMEC J^u supposedly')
% 
% figure()
% plot(data.phi(s_index:end),Jv(s_index:end,u_index,v_nfp_index))
% title(sprintf('J^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
% xlabel('s')
% ylabel('J^v')
% hold on
% plot(data.phi(s_index:end),JV_vmec(s_index:end,u_index,v_nfp_index))
% legend('My Calc J_v','VMEC J^v supposedly')
