% Calculate Force error


deriv_method_1d = 'finite difference'; %deriv method to calculate dp/ds, phi' and chi'
numerical_covariant_B_derivs = 0; % calculate cov_B derivs analytically (0) or numerically (1)
numerical_contravariant_B_derivs = 0; % calculate contra_B derivs analytically (0) or numerically (1)
use_VMEC_contr_B = 0; % use VMEC contr B components to take derivs of and get cov B derivs and then J
use_VMEC_cov_B = 0;% use VMEC cov B components to take derivs of and get J (supersedes use VMEC contr B)
only_energy=false; % only calculate energy, don't calculate force error

%% constants
mu0 = 4*pi * 1e-7;

dimS = data.ns;
dimU = 60;
if data.nfp > 1
    dimV = 60;
else
     dimV=2;
end
%% Define the (s,u,v) 3D grid on which we are evaluating the force error

Phi = -data.phi./2./pi; % toroidal flux normalized by 2*pi (Hirshman p.3, 
% chi, Phi are actually the normalized fluxes)

data.phi = data.phi/ data.phi(end); % normalize toroidal flux to get s
s = data.phi ; % radial variable
u = linspace(0,2*pi,dimU); % poloidal angle
v = linspace(0,2*pi/data.nfp,dimV); % toroidal angle

suvgrid = ndgrid(s,u,v);


%% Flux, iota,and pressure derivs
iota = data.iotaf;% rotational transform
iotar = repmat((s_deriv_non_refl(iota,data,deriv_method_1d))',1,dimU,dimV); % radial deriv of rotational transform


% deriv of phir is one
Phir = s_deriv_non_refl(Phi,data,deriv_method_1d);
chi = data.chi./2./pi; % poloidal flux
chir = repmat((iota .* Phir)',1,dimU,dimV); % radial deriv of poloidal flux 
Phir = repmat(Phir',1,dimU,dimV); % radial deriv of toroidal flux (constant =1 for Heliotron case)
chirr = iotar .* Phir;

iota = repmat(iota',1,dimU,dimV);

presr = s_deriv_non_refl(data.presf,data,deriv_method_1d);
presr = repmat(presr',1,dimU,dimV);
%% Flux surface locations
R = eval_series(suvgrid,data.rmnc,data,'c');
Z = eval_series(suvgrid,data.zmns,data,'s');
L = eval_series(suvgrid,data.lmns,data,'s');

% get the fourier coefficients for the derivatives of R,Z, and lambda
RZ_derivs
%% get derivs evaluated on the grid
%% angular derivs
R_u = eval_series(suvgrid, rumns, data, 's');
R_v = eval_series(suvgrid, rvmns, data, 's');



R_uu = eval_series(suvgrid, ruumnc, data, 'c');
R_vv = eval_series(suvgrid, rvvmnc, data, 'c');
R_uv = eval_series(suvgrid, ruvmnc, data, 'c');

Z_u = eval_series(suvgrid, zumnc, data, 'c');
Z_v = eval_series(suvgrid, zvmnc, data, 'c');


Z_uu = eval_series(suvgrid, zuumns, data, 's');
Z_vv = eval_series(suvgrid, zvvmns, data, 's');
Z_uv = eval_series(suvgrid, zuvmns, data, 's');

L_u = eval_series(suvgrid, lumnc, data, 'c');
L_v = eval_series(suvgrid, lvmnc, data, 'c');

L_uu = eval_series(suvgrid, luumns, data, 's');
L_vv = eval_series(suvgrid, lvvmns, data, 's');
L_uv = eval_series(suvgrid, luvmns, data, 's');
%% radial derivatives
if ~use_real_space_radial_derivs
R_s = eval_series(suvgrid, rsmnc, data, 'c');
Z_s = eval_series(suvgrid, zsmns, data, 's');
R_sv = eval_series(suvgrid, rsvmns, data, 's');
R_su = eval_series(suvgrid, rsumns, data, 's');
Z_sv = eval_series(suvgrid, zsvmnc, data, 'c');
Z_su = eval_series(suvgrid, zsumnc, data, 'c');

L_su = eval_series(suvgrid, lsumnc, data, 'c');
L_sv = eval_series(suvgrid, lsvmnc, data, 'c');

R_ss = eval_series(suvgrid, rssmnc, data, 'c');
Z_ss = eval_series(suvgrid, zssmns, data, 's');
else
[R_s,R_ss] = real_space_deriv_2(R, s,'c');
[Z_s,Z_ss] = real_space_deriv_2(Z, s, 's');
R_sv = real_space_deriv(R_v, s, 's');
R_su = real_space_deriv(R_u, s, 's');
Z_sv = real_space_deriv(Z_v, s,  'c');
Z_su = real_space_deriv(Z_u, s,  'c');

L_su = real_space_deriv(L_u, s,  'c');
L_sv = real_space_deriv(L_v, s,  'c');
if ~use_2nd_deriv
R_ss = real_space_deriv(R_s, s,  'c');
Z_ss = real_space_deriv(Z_s, s,  's');
end
    
end
% below derivs only used at magnetic axis

% R_sss = eval_series(suvgrid, rsssmnc, data, 'c');
% Z_sss = eval_series(suvgrid, zsssmns, data, 's');
% 
% 
% R_suu = eval_series(suvgrid, rsuumnc, data, 'c');
% Z_suu = eval_series(suvgrid, zsuumns, data, 's');
% R_suv = eval_series(suvgrid, rsuvmnc, data, 'c');
% Z_suv = eval_series(suvgrid, zsuvmns, data, 's');
% R_ssu = eval_series(suvgrid, rssumns, data, 's');
% Z_ssu = eval_series(suvgrid, zssumnc, data, 'c');
% R_ssv = eval_series(suvgrid, rssvmns, data, 's');
% Z_ssv = eval_series(suvgrid, zssvmnc, data, 'c');
% R_svv = eval_series(suvgrid, rsvvmnc, data, 'c');
% Z_svv = eval_series(suvgrid, zsvvmns, data, 's');
% 
% L_ssv = eval_series(suvgrid, lssvmnc, data, 'c');
% L_svv = eval_series(suvgrid, lsvvmns, data, 's');
% 
% R_sssu = eval_series(suvgrid, rsssumns, data, 's');
% Z_sssu = eval_series(suvgrid, zsssumnc, data, 'c');
% R_ssuu = eval_series(suvgrid, rssuumnc, data, 'c');
% Z_ssuu = eval_series(suvgrid, zssuumns, data, 's');
% R_ssuv = eval_series(suvgrid, rssuvmnc, data, 'c');
% Z_ssuv = eval_series(suvgrid, zssuvmns, data, 's');


%% covariant basis vector
es = cat(4,R_s,zeros(size(R_s)),Z_s);
eu = cat(4,R_u, zeros(size(R_s)), Z_u);
ev = cat(4,R_v, R, Z_v);

ess = cat(4,R_ss, zeros(size(R_s)), Z_ss);
esu = cat(4,R_su, zeros(size(R_s)), Z_su);
esv = cat(4,R_sv, zeros(size(R_s)), Z_sv);

eus = esu;
euu = cat(4,R_uu, zeros(size(R_s)), Z_uu);
euv = cat(4,R_uv, zeros(size(R_s)), Z_uv);

evs = cat(4,R_sv, R_s, Z_sv);
evu = cat(4,R_uv, R_u, Z_uv);
evv = cat(4,R_vv, R_v, Z_vv);
% Derivs only used at magnetic axis
% esss = cat(4,R_sss, zeros(size(R_s)), Z_sss);
% euss = cat(4,R_ssu, zeros(size(R_s)), Z_ssu);
% eusu = cat(4,R_suu, zeros(size(R_s)), Z_suu);
% evvs = cat(4,R_svv, R_sv, Z_svv);
% evsv=evvs;
% 
% essu = euss;
% evss = cat(4,R_ssv, R_ss, Z_ssv);
% essv = cat(4,R_ssv, zeros(size(R_s)), Z_ssv);
% evsu = cat(4,R_suv, R_su, Z_suv);
% esuv = cat(4,R_suv, zeros(size(R_s)), Z_suv);
% eusv = esuv;
% 
% eusss = cat(4,R_sssu, zeros(size(R_s)), Z_sssu);
% eussu = cat(4,R_ssuu, zeros(size(R_s)), Z_ssuu);
% eussv = cat(4,R_ssuv, zeros(size(R_s)), Z_ssuv);

%% Jacobian (sqrt(g)) and its derivatives 
g = dot(es,cross(eu,ev,4),4); 
gs = dot(ess,cross(eu,ev,4),4) + dot(es,cross(eus,ev,4),4) + dot(es,cross(eu,evs,4),4);
gu = dot(esu,cross(eu,ev,4),4) + dot(es,cross(euu,ev,4),4) + dot(es,cross(eu,evu,4),4);
gv = dot(esv,cross(eu,ev,4),4) + dot(es,cross(euv,ev,4),4) + dot(es,cross(eu,evv,4),4);

% derivs only used at the axis
% g_ss = dot(ess,cross(eus,ev,4),4) + dot(ess,cross(eus,ev,4),4) + dot(es,cross(euss,ev,4),4) + dot(es,cross(eus,evs,4),4) ...
%      + dot(es,cross(eus,evs,4),4);
% g_su = dot(esu,cross(eus,ev,4),4) + dot(es,cross(eusu,ev,4),4);
% g_sv = dot(esv,cross(eus,ev,4),4) + dot(es,cross(eusv,ev,4),4) + dot(es,cross(eus,evv,4),4);
% 
% g_sss = dot(esss,3*cross(eus,ev,4),4) + dot(ess,3*cross(euss,ev,4) + 6*cross(eus,evs,4),4)...
%      + dot(es,cross(eusss,ev,4) + 3*cross(euss,evs,4) + cross(euss,evss,4) + 2*cross(eus,evss,4),4);
% g_ssu = dot(essv,2*cross(eus,ev,4),4) + dot(ess,2*cross(eusv,ev,4) + 2*cross(eus,evv,4),4)...
%      + dot(esv,cross(euss,ev,4) + 2*cross(eus,evs,4),4)...
%      + dot(es,cross(eussv,ev,4) + cross(euss,evv,4) + 2*cross(eusv,evs,4) + 2*cross(eus,evsu,4),4);
% g_ssv = dot(essv,2*cross(eus,ev,4),4) + dot(ess,2*cross(eusv,ev,4)+2*cross(eus,evv,4),4)...
%     + dot(esv,cross(euss,ev,4)+2*cross(eus,evs,4),4) ...
%     + dot(es,cross(eussv,ev,4)+cross(euss,evv,4) + 2*cross(eusv,evs,4) + 2*cross(eus,evsv,4),4);
 %% contravariant basis vectors
eS = cross(eu,ev,4)./g;
eU = cross(ev,es,4)./g;
eV = cross(es,eu,4)./g;

%% metric tensor components
gSS = dot(eS,eS,4);
% define gss at magnetic axis
temp1= cross(eus,ev,4)./gs;% gotten from limit of gss at s->0
temp2 = dot(temp1,temp1,4);
gSS(1,:,:)=temp2(1,:,:);

gvv = dot(eV,eV,4);
guu = dot(eU,eU,4);
guv = dot(eU,eV,4);

%% contravariant B components
BU = (chir - Phir.*L_v)./g;
% define at magnetic axis
BU(1,:,:) = (chirr(1,:,:) - Phir(1,:,:).*L_sv(1,:,:)) ./ gs(1,:,:);


BV = Phir .* (1 + L_u)./g;
%define at magnetic axis 
BV(1,:,:) = (Phir(1,:,:).* -L_su(1,:,:)) ./ gs(1,:,:); 

if only_energy
    get_energy
    return
end

%% partial derivatives of contravariant B components
BU_s = - gs./g./g .* (chir - Phir.*L_v) + (chirr - Phir .* L_sv)./g;
BU_u = - gu./g./g .* (chir - Phir.*L_v) + ( - Phir .* L_uv)./g;
BU_v = - gv./g./g .* (chir - Phir.*L_v) + ( - Phir .* L_vv)./g;

BV_s = (-gs./g./g .* Phir).*(1 + L_u) + Phir./g .* (L_su);
BV_u = - gu./g./g .*Phir .* (1 + L_u) + Phir./g .* (L_uu);
BV_v = - gv./g./g .*Phir .* (1 + L_u) + Phir./g .* (L_uv);

% Define above at the magnetic axis
% BU_s(1,:,:) = (-g_sss(1,:,:) .*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - g_ss(1,:,:).*(chirr(1,:,:) - Phir(1,:,:).*L_sv(1,:,:)) + gs(1,:,:).*(-Phir(1,:,:).*L_ssv(1,:,:)))...
%  ./ 2 ./ (gs(1,:,:).^2);
% BU_u(1,:,:) = (-g_ssu(1,:,:).*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - 2.*g_su(1,:,:).*(chirr(1,:,:)-Phir(1,:,:).*L_sv(1,:,:)) + gu(1,:,:).*Phir(1,:,:).*L_ssv(1,:,:) )...
%  ./ 2 ./ (gs(1,:,:).^2);
% BU_v(1,:,:) = (-g_ssv(1,:,:).*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - 2.*g_sv(1,:,:).*(chirr(1,:,:)-Phir(1,:,:).*L_sv(1,:,:)) + gv(1,:,:).*Phir(1,:,:).*L_ssv(1,:,:) - g_ss(1,:,:).*Phir(1,:,:).*L_vv(1,:,:) - 2.*gs(1,:,:).*Phir(1,:,:).*L_svv(1,:,:) )...
%  ./ 2 ./ (gs(1,:,:).^2);

% BV_s(1,:,:) = ( -g_sss(1,:,:).*Phir(1,:,:)) ./ 2 ./ (gs(1,:,:).^2);
% BV_u(1,:,:) = (-g_ssu(1,:,:).*Phir(1,:,:)) ./ 2 ./ (gs(1,:,:).^2);
% BV_v(1,:,:) = ( -g_ssv(1,:,:).*Phir(1,:,:) ) ./ 2 ./ (gs(1,:,:).^2);

% save as analytic to compare to numeric later
BU_sa = BU_s;
BU_ua = BU_u;
BU_va = BU_v;
BV_sa = BV_s;
BV_ua = BV_u;
BV_va = BV_v;

% use VMEC BU,BV fourier coeffs, get derivs
if use_VMEC_contr_B
BV_smnc = s_deriv_nyq(data.bsupvmnc,data,deriv_method);
BV_s = eval_series_nyq(suvgrid,BV_smnc,data,'c');
BU_smnc = s_deriv_nyq(data.bsupvmnc,data,deriv_method);
BU_s = eval_series_nyq(suvgrid,BU_smnc,data,'c');
end

%% numerical contravariant B derivatives
if numerical_contravariant_B_derivs
    BU_s = real_space_deriv(BU,s,deriv_method);
    BU_u = real_space_deriv(BU,u,deriv_method);
    BU_v = real_space_deriv(BU,v,deriv_method);
    BV_s = real_space_deriv(BV,s,deriv_method);
    BV_u = real_space_deriv(BV,u,deriv_method);
    BV_v = real_space_deriv(BV,v,deriv_method);
    run('plotting/plot_contr_B_derivs')
end


%% covariant B derivatives
Bu_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4) + dot(BU.*eu + BV.*ev,eus,4);
Bv_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,ev,4) + dot(BU.*eu + BV.*ev,evs,4);

Bs_u = dot( BU_u.*eu + BU.*euu+ BV_u.*ev + BV .*evu,es,4) + dot(BU.*eu + BV.*ev,esu,4);
Bv_u = dot( BU_u.*eu + BU.*euu+ BV_u.*ev + BV .*evu,ev,4) + dot(BU.*eu + BV.*ev,evu,4);

Bs_v = dot( BU_v.*eu + BU.*euv+ BV_v.*ev + BV .*evv,es,4) + dot(BU.*eu + BV.*ev,esv,4);
Bu_v = dot( BU_v.*eu + BU.*euv+ BV_v.*ev + BV .*evv,eu,4) + dot(BU.*eu + BV.*ev,euv,4);

if max(data.presf) < 1e-3
   Bu_u =  dot( BU_u.*eu + BU.*euu+ BV_u.*ev + BV .*evu,eu,4) + dot(BU.*eu + BV.*ev,euu,4);
   Bv_v = dot( BU_v.*eu + BU.*euv+ BV_v.*ev + BV .*evv,ev,4) + dot(BU.*eu + BV.*ev,evv,4);
end


% mark as analytic derivatives, as opposed to numeric like I calc later for
% comparison
Bu_sa = Bu_s;
Bv_sa = Bv_s;

Bs_ua = Bs_u;
Bv_ua = Bv_u;

Bs_va = Bs_v;
Bu_va = Bu_v;

%% covariant B components
g_sv = dot(es,ev,4);

g_us = dot(eu,es,4);
g_uu = dot(eu,eu,4);
g_uv = dot(eu,ev,4);%g_uv = g_vu

g_vs = dot(ev,es,4);
g_vu = dot(ev,eu,4);
g_vv = dot(ev,ev,4);

Bs = BU.*g_us + BV .* g_vs;
Bu = BU.*g_uu + BV .* g_vu;
Bv = BU.*g_uv + BV .* g_vv;

%% numerical covariant B derivatives


if numerical_covariant_B_derivs
    Bs_u = real_space_deriv(Bs,u,deriv_method);
    Bs_v = real_space_deriv(Bs,v,deriv_method);

    Bu_s = real_space_deriv(Bu,s,deriv_method);
    Bu_v = real_space_deriv(Bu,v,deriv_method);

    Bv_s = real_space_deriv(Bv,s,deriv_method);
    Bv_u = real_space_deriv(Bv,u,deriv_method);
    run('plotting/plot_cov_B')
end

if use_VMEC_cov_B
Bv_smnc = s_deriv_nyq(data.bsubvmnc,data,deriv_method);
Bv_s = eval_series_nyq(suvgrid,Bv_smnc,data,'c');
Bv_umns = -data.bsubvmnc .* data.xm_nyq';
Bv_u = eval_series_nyq(suvgrid,Bv_umns,data,'s');

Bu_smnc = s_deriv_nyq(data.bsubumnc,data,deriv_method);
Bu_s = eval_series_nyq(suvgrid,Bu_smnc,data,'c');
Bu_vmns = -data.bsubumnc .* data.xn_nyq';
Bu_v = eval_series_nyq(suvgrid,Bu_vmns,data,'s');

Bs_vmnc = data.bsubsmns .* data.xn_nyq';
Bs_v = eval_series_nyq(suvgrid,Bs_vmnc,data,'c');
Bs_umnc = data.bsubsmns .* data.xm_nyq';
Bs_u = eval_series_nyq(suvgrid,Bs_umnc,data,'c');

end




%% Contravariant current (J) components
JS = (Bv_u - Bu_v) ./ mu0 ./ g;
JU = (Bs_v - Bv_s) ./ mu0 ./ g;
JU(1,:,:) = (Bs_v(1,:,:) - Bv_s(1,:,:)) ./ mu0; %at axis g is zero, cancel it out with the g in the Fs eqn

JV = (Bu_s - Bs_u) ./ mu0 ./ g;
JV(1,:,:) = (Bu_s(1,:,:) - Bs_u(1,:,:)) ./ mu0;


%% Magnitudes of direction vectors
mag_eS = sqrt(gSS);
new_mag_beta = sqrt((BV.^2).*dot(cross(ev,es,4),cross(ev,es,4),4) + (BU.^2).*dot(cross(es,eu,4),cross(es,eu,4),4) - 2.*BV.*BU.*dot(cross(ev,es,4),cross(es,eu,4),4));% just checking that if beta written with g cancelled with the denominators of gss,guu, is the same as the above line (it is)


%% Magnitude of Force error (N/m3)
F_s = g .* (JV.*BU - JU.* BV) + presr;
F_s(1,:,:) = (JV(1,:,:).*BU(1,:,:) - JU(1,:,:).*BV(1,:,:)) + presr(1,:,:);% define at axis by cancelling the g's
F_beta = JS;


F = sqrt((F_s.^2).*gSS + (F_beta.^2).*(new_mag_beta.^2));

% calculate energy
get_energy






