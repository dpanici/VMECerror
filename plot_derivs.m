% close all
% Script for use after running force_error, so suvgrid and real-space
% quantities are already defined

% Check derivative values by plotting derivs and original quantities
% together

% check iota, iotar and chi, chir - Good
% figure()
% plot(data.phi,iota)
% hold on
% plot(data.phi,s_deriv(iota,data))
% legend('iota','iota_s')
% 
% figure()
% plot(data.phi,chi)
% hold on
% plot(data.phi,s_deriv(chi,data))
% hold on
% plot(data.phi,-iota.*s_deriv(Phi,data))
% legend('chi','chi_s','-iota.*Phir')

%% Angular Derivs - All good after realizing sign error in v derivs
% derivplot(R,R_u,u) % good
% derivplot(Z,Z_u,u) % good
% derivplot(R,R_v,v) % good
% derivplot(Z,Z_v,v) % good
% derivplot(L,L_u,u) % good % Matches VMECplot
% derivplot(L,L_v,v) % good % L actually matches VMEC plot too...
%  
% derivplot(R_u,R_uu,u) % good
% derivplot(Z_u,Z_uu,u) % good
% derivplot(R_v,R_vv,v) % good
% derivplot(Z_v,Z_vv,v) % good
% derivplot(L_u,L_uu,u) % good
% derivplot(L_v,L_vv,v) % good
% 
% derivplot(R_u,R_uv,v) % good
% derivplot(Z_u,Z_uv,v) % good
% derivplot(R_v,R_uv,u) % good
% derivplot(Z_v,Z_uv,u) % good
% derivplot(L_u,L_uv,v) % good
% derivplot(L_v,L_uv,u) % good

%% Radial Derivatives
derivplot(R,R_s,s) % good
derivplot(Z,Z_s,s) % good
% 
derivplot(R_s,R_ss,s) % good
derivplot(Z_s,Z_ss,s) % good
% 
% derivplot(R_u,R_su,s) % good
% derivplot(Z_u,Z_su,s) % good
% derivplot(R_v,R_sv,s) % good
% derivplot(Z_v,Z_sv,s) % good
% derivplot(L_u,L_su,s) % good
% derivplot(L_v,L_sv,s) % good


% 
% figure() % dR/ds is multi-valued at s=0... I think dan had said this before
% % will likely just average out the force balance error calculated at s=0,
% % u for each v (unless I think of a better way to handle this)
% u = linspace(0,2*pi,100);
% R_s_s0v0 = zeros([1,length(u)]);
% for i = 1:length(d.xm)
%     R_s_s0v0 = R_s_s0v0 + rsmnc(i,1) .* cos(i.*u);
% end
% plot(u,R_s_s0v0)
% title('dR/ds versus u at s=v=0')
% xlabel('u')
% ylabel('dR/ds')

function foo = derivplot(value,deriv,var_wrt_to)
% plot a value and its deriv against a variable the deriv is wrt to, 
valname = inputname(1);
derivname = inputname(2);
wrtname = inputname(3);
sindex=5;
uindex=5;
vindex=1;
figure()
yyaxis left
if wrtname == 'u'
    plot(var_wrt_to,value(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to,value(:,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(value(sindex,uindex,:),size(var_wrt_to)))
end     
hold on
ylabel(valname)
yyaxis right
if wrtname == 'u'
    plot(var_wrt_to,deriv(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to,deriv(:,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(deriv(sindex,uindex,:),size(var_wrt_to)))
end     
hold on
yline(0,'--')
title(sprintf('%s and %s versus %s',valname,derivname,wrtname))
xlabel(wrtname)
ylabel(derivname)
legend(valname,derivname)
set(gcf, 'Position',  [200, 200, 1000, 800])

% fin dif to approx the deriv, use ismembertol to check that the fin dif is
% similar to the calculated deriv

end