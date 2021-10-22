close all
F_s_fsa = trapz(u,trapz(v,F_s.*abs_g_vmec,3),2);
F_beta_fsa = trapz(u,trapz(v,F_beta.*abs_g_vmec,3),2);
gSS_fsa = trapz(u,trapz(v,gSS.*abs_g_vmec,3),2);
mag_beta_fsa = trapz(u,trapz(v,new_mag_beta.*abs_g_vmec,3),2);
s_plot = s(s_ind:end-4);

figure
subplot(5,1,1)
ax1 = gca;
plot(s_plot,F_rhos(s_ind:end-4))
title('FSA of Force Error in (N/m^3)')
set(gca, 'YScale', 'log')
subplot(5,1,2)
ax2 = gca;
plot(s_plot,abs(F_s_fsa(s_ind:end-4)))
title('F_s')
set(gca, 'YScale', 'log')
subplot(5,1,3)
ax3 = gca;
plot(s_plot,abs(F_beta_fsa(s_ind:end-4)))
title('F_beta')
set(gca, 'YScale', 'log')
subplot(5,1,4)
ax4 = gca;
plot(s_plot,gSS_fsa(s_ind:end-4))
title('gSS')
set(gca, 'YScale', 'log')
subplot(5,1,5)
ax5 = gca;
plot(s_plot,mag_beta_fsa(s_ind:end-4))
title('|beta|')
set(gca, 'YScale', 'log')
linkaxes([ax1 ax2 ax3 ax4 ax5],'x')


u_index=1;
v_index=1;
s_ind=15;
s_plot = s(s_ind:end-4);

fsize=16

figure
subplot(5,1,1)
ax1 = gca;
plot(s_plot,abs(F(s_ind:end-4,u_index,v_index)))
title('$$ |F| = \sqrt{F_{s}^{2}gSS + F_{\beta}^2 |\beta|}$$ at u=0, v=0','interpreter','latex','FontSize',fsize)
set(gca, 'YScale', 'log')
subplot(5,1,2)
ax2 = gca;
plot(s_plot,abs(F_s(s_ind:end-4,u_index,v_index)))
title('$$|F_s|$$','interpreter','latex','FontSize',fsize)
set(gca, 'YScale', 'log')
subplot(5,1,3)
ax3 = gca;
plot(s_plot,abs(F_beta(s_ind:end-4,u_index,v_index)))
title('$$|F_{\beta}|$$','interpreter','latex','FontSize',fsize)
% set(gca, 'YScale', 'log')
subplot(5,1,4)
ax4 = gca;
plot(s_plot,gSS(s_ind:end-4,u_index,v_index))
title('$$g^{ss}$$','interpreter','latex','FontSize',fsize)
set(gca, 'YScale', 'log')
subplot(5,1,5)
ax5 = gca;
plot(s_plot,new_mag_beta(s_ind:end-4,u_index,v_index))
title('$$|\beta|$$','interpreter','latex','FontSize',fsize)
set(gca, 'YScale', 'log')
linkaxes([ax1 ax2 ax3 ax4 ax5],'x')


% F_s components
s_plot = s(s_ind:end-4);
figure
subplot(6,1,1)
ax1 = gca;
plot(s_plot,abs(F_s(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title("$$ |F_s| = |\sqrt{g}(J^vB^u - J^u B^v) + p'|$$",'interpreter','latex')
subplot(6,1,2)
ax2 = gca;
plot(s_plot,abs(JV(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title('|JV|')
subplot(6,1,3)
ax3 = gca;
plot(s_plot,abs(BU(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title('|BU|')
subplot(6,1,4)
ax4 = gca;
plot(s_plot,abs(JU(s_ind:end-4,u_index,v_index)))
title('|JU|')
set(gca, 'YScale', 'log')

subplot(6,1,5)
ax5 = gca;
plot(s_plot,abs(BV(s_ind:end-4,u_index,v_index)))
title('|BV|')
set(gca, 'YScale', 'log')

subplot(6,1,6)
ax6 = gca;
plot(s_plot,abs(g(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')
title('$$|\sqrt(g)|$$','interpreter','latex')
linkaxes([ax1 ax2 ax3 ax4 ax5,ax6],'x')

%%%
figure
subplot(4,1,1)
ax1 = gca;
plot(s_plot,abs(F_s(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title("$$ |F_s| = |\sqrt{g}(J^vB^u - J^u B^v) + p'|$$ at u=0, v=0",'interpreter','latex')
subplot(4,1,2)
ax2 = gca;
gJVBU = g.*JV.*BU;
plot(s_plot,gJVBU(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

gJUBV = g.*JU.*BV;

title('$$\sqrt{g}(J^vB^u)$$','interpreter','latex')
subplot(4,1,3)
ax3 = gca;

plot(s_plot,-gJUBV(s_ind:end-4,u_index,v_index),'r','DisplayName','Negative')
set(gca, 'YScale', 'log')
legend
title('$$\sqrt{g}(J^uB^v)$$','interpreter','latex')

subplot(4,1,4)
ax4 = gca;

plot(s_plot,abs(presr(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')
legend
title('$$|\nabla p|$$','interpreter','latex')


linkaxes([ax1 ax2 ax3 ax4],'x')

%% check currents
figure
subplot(4,1,1)
ax1 = gca;
gJVBU = g.*JV.*BU;
plot(s_plot,gJVBU(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$\sqrt{g}(J^vB^u)$$','interpreter','latex')
subplot(4,1,2)
ax2 = gca;

plot(s_plot,JV(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$J^v$$','interpreter','latex')

subplot(4,1,3)
ax3 = gca;

plot(s_plot,BU(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$B^u$$','interpreter','latex')

subplot(4,1,4)
ax4 = gca;

plot(s_plot,abs(g(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title('$$|\sqrt{g}|$$','interpreter','latex')


linkaxes([ax1 ax2 ax3 ax4],'x')

%% check JV
%JV = (Bu_s - Bs_u) ./ mu0 ./ g;
figure
subplot(4,1,1)
ax1 = gca;
plot(s_plot,JV(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$J^v = (\frac{\partial B_u}{\partial s} - \frac{\partial B_s}{\partial u})/\sqrt{g} / \mu_0$$','interpreter','latex')
subplot(4,1,2)
ax2 = gca;

plot(s_plot,-Bu_s(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$-\frac{\partial B_u}{\partial s}$$','interpreter','latex')

subplot(4,1,3)
ax3 = gca;

plot(s_plot,-Bs_u(s_ind:end-4,u_index,v_index))
set(gca, 'YScale', 'log')

title('$$-\frac{\partial B_s}{\partial u}$$','interpreter','latex')

subplot(4,1,4)
ax4 = gca;

plot(s_plot,abs(1./g(s_ind:end-4,u_index,v_index)))
set(gca, 'YScale', 'log')

title('$$-\frac{1}{\sqrt{g}}$$','interpreter','latex')


linkaxes([ax1 ax2 ax3 ax4],'x')

%% cov B radial derivative
%Bu_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4) + dot(BU.*eu + BV.*ev,eus,4);
term1 = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4);
term1_str = '$$[\left(\frac{\partial B^u}{\partial s}\mathbf{e}_u + B^u \mathbf{e}_{us} + \frac{\partial B^v}{\partial s}\mathbf{e}_v + B^v \mathbf{e}_{vs}\right) \cdot \mathbf{e}_u]$$';
term2 = dot(BU.*eu + BV.*ev,eus,4);
term2_str = '$$(B^u \mathbf{e}_u + B^v \mathbf{e}_v) \cdot \mathbf{e}_{us}$$';
figure
subplot(3,1,1)
ax1 = gca;
plot(s_plot,Bu_s(s_ind:end-4,u_index,v_index))

title('$$ \frac{\partial B_u}{\partial s} = \left(\frac{\partial B^u}{\partial s}\mathbf{e}_u + B^u \mathbf{e}_{us} + \frac{\partial B^v}{\partial s}\mathbf{e}_v + B^v \mathbf{e}_{vs}\right) \cdot \mathbf{e}_u + (B^u \mathbf{e}_u + B^v \mathbf{e}_v) \cdot \mathbf{e}_{us}$$','interpreter','latex')
subplot(3,1,2)
ax2 = gca;

plot(s_plot,term1(s_ind:end-4,u_index,v_index))

title(term1_str,'interpreter','latex')

subplot(3,1,3)
ax3 = gca;

plot(s_plot,(term2(s_ind:end-4,u_index,v_index)))

title(term2_str,'interpreter','latex')


linkaxes([ax1 ax2 ax3],'x')

%% first term of above is the issue
%check it too
%Bu_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4) + dot(BU.*eu + BV.*ev,eus,4);
term1 = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4);
term1_str = '$$[\left(\frac{\partial B^u}{\partial s}\mathbf{e}_u + B^u \mathbf{e}_{us} + \frac{\partial B^v}{\partial s}\mathbf{e}_v + B^v \mathbf{e}_{vs}\right) \cdot \mathbf{e}_u]$$';

t1 = dot( BU_s.*eu,eu,4);
t1_s = '$$\frac{\partial B^u}{\partial s}\mathbf{e}_u \cdot \mathbf{e}_u$$';
t2 = dot(BU.*eus,eu,4);
t2_s = '$$B^u \mathbf{e}_{us}\cdot \mathbf{e}_u$$';
t3 = dot(BV_s.*ev,eu,4);
t3_s = '$$\frac{\partial B^v}{\partial s}\mathbf{e}_v \cdot \mathbf{e}_u$$';
t4 = dot(BV .*evs,eu,4);
t4_s = '$$B^v \mathbf{e}_{vs}\ \cdot \mathbf{e}_u$$';
figure
subplot(5,1,1)
ax1 = gca;
plot(s_plot,term1(s_ind:end-4,u_index,v_index))
title(term1_str,'interpreter','latex','FontSize',fsize)

subplot(5,1,2)
ax2 = gca;

plot(s_plot,(t1(s_ind:end-4,u_index,v_index)))

title(t1_s,'interpreter','latex','FontSize',fsize)
subplot(5,1,5)
ax5 = gca;

plot(s_plot,(t2(s_ind:end-4,u_index,v_index)))

title(t2_s,'interpreter','latex','FontSize',fsize)

subplot(5,1,4)
ax4 = gca;

plot(s_plot,(t3(s_ind:end-4,u_index,v_index)))

title(t3_s,'interpreter','latex','FontSize',fsize)

subplot(5,1,3)
ax3 = gca;

plot(s_plot,(t4(s_ind:end-4,u_index,v_index)))

title(t4_s,'interpreter','latex','FontSize',fsize)


linkaxes([ax1 ax2 ax3 ax4 ax5],'x')

%% check BV_s and eu dot eu
figure
subplot(4,1,1)
ax1 = gca;

plot(s_plot,(t1(s_ind:end-4,u_index,v_index)))

title(t1_s,'interpreter','latex','FontSize',fsize)

subplot(4,1,2)
ax2 = gca;
plot(s_plot,BV(s_ind:end-4,u_index,v_index))
title('$$B^v$$','interpreter','latex','FontSize',fsize)

subplot(4,1,3)
ax3 = gca;
plot(s_plot,BV_s(s_ind:end-4,u_index,v_index))
title('$$\frac{\partial B^v}{\partial s}$$','interpreter','latex','FontSize',fsize)

subplot(4,1,4)
ax4 = gca;
plot(s_plot,g_uu(s_ind:end-4,u_index,v_index))
title('$$g_{uu} = \mathbf{e}_u \cdot \mathbf{e}_u$$','interpreter','latex','FontSize',fsize)

linkaxes([ax1 ax2 ax3 ax4],'x')

