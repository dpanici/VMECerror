% Calculate Force error
clearvars, close all
% file = 'wout_HELIOTRON_32x8x8.nc';
% file = 'wout_HELIOTRON_16x4x4.nc';
file = 'wout_HELIOTRON.nc';
% file = 'wout_DSHAPE.nc';

data = read_vmec(file);

%% constants
mu0 = 4*pi * 1e-7;

dimS = data.ns;
dimU = 50;
dimV = 380;
%% Define the (s,u,v) 3D grid on which we are evaluating the force error
s = linspace(0,1,dimS);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi/data.nfp,dimV);

suvgrid = ndgrid(s,u,v);

volS = 150;
s_vol = linspace(0,1,volS);
volgrid = ndgrid(s_vol,u,v);

%% Flux, iota,and pressure derivs
iota = data.iotaf;% rotational transform
iotar = repmat((s_deriv(iota,data))',1,dimU,dimV); % radial deriv or rotational transform


Phi = data.phi./2./pi; % toroidal flux normalized by 2*pi (Hirshman p.3, 
% chi, Phi are actually the normalized fluxes
Phir = s_deriv(Phi,data);
chi = data.chi./2./pi; % poloidal flux
chir = repmat((-iota .* Phir)',1,dimU,dimV); % radial deriv of poloidal flux 
Phir = repmat(Phir',1,dimU,dimV); % radial deriv of toroidal flux (constant =1 for Heliotron case)
chirr = iotar .* Phir;

iota = repmat(iota',1,dimU,dimV);

presr = s_deriv(data.presf,data);
presr = repmat(presr',1,dimU,dimV);
%% Flux surface locations
R = eval_series(suvgrid,data.rmnc,data,'c');
Z = eval_series(suvgrid,data.zmns,data,'s');
L = eval_series(suvgrid,data.lmns,data,'s');

% get the fourier coefficients for the derivatives of R,Z, and lambda
RZ_derivs
% get derivs evaluated on the grid
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

% below derivs only used at magnetic axis
R_sss = eval_series(suvgrid, rsssmnc, data, 'c');
Z_sss = eval_series(suvgrid, zsssmns, data, 's');


R_suu = eval_series(suvgrid, rsuumnc, data, 'c');
Z_suu = eval_series(suvgrid, zsuumns, data, 's');
R_suv = eval_series(suvgrid, rsuvmnc, data, 'c');
Z_suv = eval_series(suvgrid, zsuvmns, data, 's');
R_ssu = eval_series(suvgrid, rssumns, data, 's');
Z_ssu = eval_series(suvgrid, zssumnc, data, 'c');
R_ssv = eval_series(suvgrid, rssvmns, data, 's');
Z_ssv = eval_series(suvgrid, zssvmnc, data, 'c');
R_svv = eval_series(suvgrid, rsvvmnc, data, 'c');
Z_svv = eval_series(suvgrid, zsvvmns, data, 's');

L_ssv = eval_series(suvgrid, lssvmnc, data, 'c');
L_svv = eval_series(suvgrid, lsvvmns, data, 's');

R_sssu = eval_series(suvgrid, rsssumns, data, 's');
Z_sssu = eval_series(suvgrid, zsssumnc, data, 'c');
R_ssuu = eval_series(suvgrid, rssuumnc, data, 'c');
Z_ssuu = eval_series(suvgrid, zssuumns, data, 's');
R_ssuv = eval_series(suvgrid, rssuvmnc, data, 'c');
Z_ssuv = eval_series(suvgrid, zssuvmns, data, 's');


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
esss = cat(4,R_sss, zeros(size(R_s)), Z_sss);
euss = cat(4,R_ssu, zeros(size(R_s)), Z_ssu);
eusu = cat(4,R_suu, zeros(size(R_s)), Z_suu);
evvs = cat(4,R_svv, R_sv, Z_svv);
evsv=evvs;

essu = euss;
evss = cat(4,R_ssv, R_ss, Z_ssv);
essv = cat(4,R_ssv, zeros(size(R_s)), Z_ssv);
evsu = cat(4,R_suv, R_su, Z_suv);
esuv = cat(4,R_suv, zeros(size(R_s)), Z_suv);
eusv = esuv;

eusss = cat(4,R_sssu, zeros(size(R_s)), Z_sssu);
eussu = cat(4,R_ssuu, zeros(size(R_s)), Z_ssuu);
eussv = cat(4,R_ssuv, zeros(size(R_s)), Z_ssuv);

%% Jacobian (sqrt(g)) and its derivatives 
g = dot(es,cross(eu,ev,4),4); % g is negative... but matches matlabVMEC, this is b/c HELIOTRON input has a left-handed poloidal convention
gs = dot(ess,cross(eu,ev,4),4) + dot(es,cross(eus,ev,4),4) + dot(es,cross(eu,evs,4),4);
gu = dot(esu,cross(eu,ev,4),4) + dot(es,cross(euu,ev,4),4) + dot(es,cross(eu,evu,4),4);
gv = dot(esv,cross(eu,ev,4),4) + dot(es,cross(euv,ev,4),4) + dot(es,cross(eu,evv,4),4);

% derivs only used at the axis
g_ss = dot(ess,cross(eus,ev,4),4) + dot(ess,cross(eus,ev,4),4) + dot(es,cross(euss,ev,4),4) + dot(es,cross(eus,evs,4),4) ...
     + dot(es,cross(eus,evs,4),4);
g_su = dot(esu,cross(eus,ev,4),4) + dot(es,cross(eusu,ev,4),4);
g_sv = dot(esv,cross(eus,ev,4),4) + dot(es,cross(eusv,ev,4),4) + dot(es,cross(eus,evv,4),4);

g_sss = dot(esss,3*cross(eus,ev,4),4) + dot(ess,3*cross(euss,ev,4) + 6*cross(eus,evs,4),4)...
     + dot(es,cross(eusss,ev,4) + 3*cross(euss,evs,4) + cross(euss,evss,4) + 2*cross(eus,evss,4),4);
g_ssu = dot(essv,2*cross(eus,ev,4),4) + dot(ess,2*cross(eusv,ev,4) + 2*cross(eus,evv,4),4)...
     + dot(esv,cross(euss,ev,4) + 2*cross(eus,evs,4),4)...
     + dot(es,cross(eussv,ev,4) + cross(euss,evv,4) + 2*cross(eusv,evs,4) + 2*cross(eus,evsu,4),4);
g_ssv = dot(essv,2*cross(eus,ev,4),4) + dot(ess,2*cross(eusv,ev,4)+2*cross(eus,evv,4),4)...
    + dot(esv,cross(euss,ev,4)+2*cross(eus,evs,4),4) ...
    + dot(es,cross(eussv,ev,4)+cross(euss,evv,4) + 2*cross(eusv,evs,4) + 2*cross(eus,evsv,4),4);
 %% contravariant basis vectors
eS = cross(eu,ev,4)./g;
eU = cross(ev,es,4)./g;
eV = cross(es,eu,4)./g;

%% metric tensor components
gss = dot(eS,eS,4);
% define gss at magnetic axis
temp1= cross(eus,ev,4)./gs;% gotten from limit of gss at s->0
temp2 = dot(temp1,temp1,4);
gss(1,:,:)=temp2(1,:,:);

gvv = dot(eV,eV,4);
guu = dot(eU,eU,4);
guv = dot(eU,eV,4);

%% contravariant B components
BU = (chir - Phir.*L_v)./g;
% define at magnetic axis
BU(1,:,:) = (chirr(1,:,:) - Phir(1,:,:).*L_sv(1,:,:)) ./ gs(1,:,:);

BV = -Phir .* (1 + L_u)./g;
%define at magnetic axis 
BV(1,:,:) = (Phir(1,:,:).* -L_su(1,:,:)) ./ gs(1,:,:); % L_su is NOT zero, we can define it (L_u is maybe zero tho)

%% partial derivatives of contravariant B components
BU_s = - gs./g./g .* (chir - Phir.*L_v) + (chirr - Phir .* L_sv)./g;
BU_u = - gu./g./g .* (chir - Phir.*L_v) + ( - Phir .* L_uv)./g;
BU_v = - gv./g./g .* (chir - Phir.*L_v) + ( - Phir .* L_vv)./g;

BV_s = (-gs./g./g .* Phir).*(1 - L_u) + Phir./g .* (-L_su);
BV_u = - gu./g./g .*Phir .* (1 - L_u) + Phir./g .* (-L_uu);
BV_v = - gv./g./g .*Phir .* (1 - L_u) + Phir./g .* (-L_uv);

% Define above at the magnetic axis
BU_s(1,:,:) = (-g_sss(1,:,:) .*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - g_ss(1,:,:).*(chirr(1,:,:) - Phir(1,:,:).*L_sv(1,:,:)) + gs(1,:,:).*(-Phir(1,:,:).*L_ssv(1,:,:)))...
 ./ 2 ./ (gs(1,:,:).^2);
BU_u(1,:,:) = (-g_ssu(1,:,:).*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - 2.*g_su(1,:,:).*(chirr(1,:,:)-Phir(1,:,:).*L_sv(1,:,:)) + gu(1,:,:).*Phir(1,:,:).*L_ssv(1,:,:) )...
 ./ 2 ./ (gs(1,:,:).^2);
BU_v(1,:,:) = (-g_ssv(1,:,:).*(chir(1,:,:) - Phir(1,:,:).*L_v(1,:,:)) - 2.*g_sv(1,:,:).*(chirr(1,:,:)-Phir(1,:,:).*L_sv(1,:,:)) + gv(1,:,:).*Phir(1,:,:).*L_ssv(1,:,:) - g_ss(1,:,:).*Phir(1,:,:).*L_vv(1,:,:) - 2.*gs(1,:,:).*Phir(1,:,:).*L_svv(1,:,:) )...
 ./ 2 ./ (gs(1,:,:).^2);

BV_s(1,:,:) = ( -g_sss(1,:,:).*Phir(1,:,:)) ./ 2 ./ (gs(1,:,:).^2);
BV_u(1,:,:) = (-g_ssu(1,:,:).*Phir(1,:,:)) ./ 2 ./ (gs(1,:,:).^2);
BV_v(1,:,:) = ( -g_ssv(1,:,:).*Phir(1,:,:) ) ./ 2 ./ (gs(1,:,:).^2);
%% covariant B derivatives
Bu_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,eu,4) + dot(BU.*eu + BV.*ev,eus,4);
Bv_s = dot( BU_s.*eu + BU.*eus+ BV_s.*ev + BV .*evs,ev,4) + dot(BU.*eu + BV.*ev,evs,4);

Bs_u = dot( BU_u.*eu + BU.*euu+ BV_u.*ev + BV .*evu,es,4) + dot(BU.*eu + BV.*ev,esu,4);
Bv_u = dot( BU_u.*eu + BU.*euu+ BV_u.*ev + BV .*evu,ev,4) + dot(BU.*eu + BV.*ev,evu,4);

Bs_v = dot( BU_v.*eu + BU.*euv+ BV_v.*ev + BV .*evv,es,4) + dot(BU.*eu + BV.*ev,esv,4);
Bu_v = dot( BU_v.*eu + BU.*euv+ BV_v.*ev + BV .*evv,eu,4) + dot(BU.*eu + BV.*ev,euv,4);

%% Contravariant current (J) components
JS = (Bv_u - Bu_v) ./ mu0 ./ g;
JU = (Bs_v - Bv_s) ./ mu0 ./ g;
JU(1,:,:) = (Bs_v(1,:,:) - Bv_s(1,:,:)) ./ mu0; %at axis g is zero, cancel it out with the g in the Fs eqn
JV = (Bu_s - Bs_u) ./ mu0 ./ g;
JV(1,:,:) = (Bu_s(1,:,:) - Bs_u(1,:,:)) ./ mu0;


%% Magnitudes of direction vectors
mag_eS = sqrt(gss);
mag_beta = g .* sqrt((BV.^2).*guu + (BU.^2).*gvv - 2.*BV.*BU.*guv);
new_mag_beta = sqrt((BV.^2).*dot(cross(ev,es,4),cross(ev,es,4),4) + (BU.^2).*dot(cross(es,eu,4),cross(es,eu,4),4) - 2.*BV.*BU.*dot(cross(ev,es,4),cross(es,eu,4),4));% just checking that if beta written with g cancelled with the denominators of gss,guu, is the same as the above line (it is)
beta_is_same = min(ismembertol(abs(mag_beta(2:end,:,:)),new_mag_beta(2:end,:,:),1e-2),[],'all'); % can write beta without the magnitude of the Jacobian, tho it is
%zero at the magnetic axis. 
%mag_beta(1,:,:) = BV(1,:,:);
%% Magnitude of Force error (N/m3)
F_s = g .* (JV.*BU - JU.* BV) + presr;
F_s(1,:,:) = (JV(1,:,:).*BU(1,:,:) - JU(1,:,:).*BV(1,:,:)) + presr(1,:,:);% define at axis by cancelling the g's
F_beta = JS;
% beta_comp = F_beta .* mag_beta;
% beta_comp(1,:,:) = 

F = sqrt((F_s.^2).*gss) + (F_beta.^2).*(mag_beta.^2);

%% Plot
plot_force_error
% [s1,u1,v1] = ndgrid(s,u,v);
% [s_vol,u,v] = ndgrid(s_vol,u,v); % don't do this, interp the derivs and R,Z instead
% Rint = griddedInterpolant(s1,u1,v1,R,'linear');
% Zint = griddedInterpolant(s1,u1,v1,Z,'linear');
% Fint = griddedInterpolant(s1,u1,v1,F,'linear');
% 
% R = Rint(s_vol,u,v);
% Z = Zint(s_vol,u,v);
% F = Fint(s_vol,u,v);
% 
% plot_force_error


