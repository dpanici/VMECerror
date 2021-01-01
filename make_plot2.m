% Plot Bs
figure(11)

plot(data.phi(s_index:end),Bs(s_index:end,u_index,v_nfp_index))
title(sprintf('B_s vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_s')


hold on

% Plot Bu
figure(12)

plot(data.phi(s_index:end),Bu(s_index:end,u_index,v_nfp_index))
title(sprintf('B_u vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_u')
legend('Analytic','Numerical')
% Plot Bv
figure(13)

plot(data.phi(s_index:end),Bv(s_index:end,u_index,v_nfp_index))
title(sprintf('B_v vmec vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('B_v')
legend('Analytic','Numerical')