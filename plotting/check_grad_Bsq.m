%magB_sq already exists from get_energy.m
% derivmeth = 'finite difference';
derivmeth='spline';
Bsq_s = real_space_deriv(magB_sq,s,derivmeth);
Bsq_u = real_space_deriv(magB_sq,u,derivmeth);
Bsq_v = real_space_deriv(magB_sq,v,derivmeth);

figure
plot(s,2*grad_B_pres_s(:,u_index,v_nfp_index),'DisplayName','Analytic')
hold on
plot(s,Bsq_s(:,u_index,v_nfp_index),'Linestyle','--','DisplayName','Numerical')
xlabel('s')
ylabel('Grad |B|^2 wrt s')
legend

figure
plot(u,2*grad_B_pres_u(s_index,:,v_nfp_index),'DisplayName','Analytic')
hold on
plot(u,Bsq_u(s_index,:,v_nfp_index),'Linestyle','--','DisplayName','Numerical')
xlabel('u')
ylabel('Grad |B|^2 wrt u')
legend

figure
plot(v,reshape(2*grad_B_pres_v(s_index,u_index,:),size(v)),'DisplayName','Analytic')
hold on
plot(v,reshape(Bsq_v(s_index,u_index,:),size(v)),'Linestyle','--','DisplayName','Numerical')
xlabel('v')
ylabel('Grad |B|^2 wrt v')
legend

