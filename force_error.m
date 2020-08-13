% Calculate Force error
file = 'wout_HELIOTRON_16x4x4.nc';
% file = 'wout_HELIOTRON_32x8x8.nc';
data = read_vmec(file);

%% constants
mu0 = 4*pi * 1e-7;

dimS = data.ns;
dimU = 50;
dimV = 10;
%% Define the (s,u,v) 3D grid on which we are evaluating the force error
s = linspace(0,1,dimS);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

suvgrid = ndgrid(s,u,v);


%% Flux, iota,and pressure derivs
iota = data.iotaf;% rotational transform
iotar = repmat((s_deriv(iota,data))',1,dimU,dimV); % radial deriv or rotational transform

Phi = data.phi; % toroidal flux
Phir = s_deriv(Phi,data);
chi = data.chi; % poloidal flux
chir = repmat((iota .* Phir)',1,dimU,dimV); % radial deriv of poloidal flux 
Phir = repmat(Phir',1,dimU,dimV); % radial deriv of toroidal flux (constant =1 for Heliotron case)
chirr = iotar .* Phir;

presr = s_deriv(data.presf,data);
presr = repmat(presr',1,dimU,dimV);
%% Flux surface locations
R = eval_series(suvgrid,data.rmnc,data,'c');
Z = eval_series(suvgrid,data.zmns,data,'s');

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


%% covariant basis vector
es = cat(4,R_s,zeros(size(R_s)),R_v);
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

%% Jacobian (sqrt(g)) and its derivatives 
g = dot(es,cross(eu,ev,4),4);
gs = dot(ess,cross(eu,ev,4),4) + dot(es,cross(eus,ev,4),4) + dot(es,cross(eu,evs,4),4);
gu = dot(esu,cross(eu,ev,4),4) + dot(es,cross(euu,ev,4),4) + dot(es,cross(eu,evu,4),4);
gv = dot(esv,cross(eu,ev,4),4) + dot(es,cross(euv,ev,4),4) + dot(es,cross(eu,evv,4),4);

%% contravariant basis vectors
eS = cross(eu,ev,4)./g;
eU = cross(ev,es,4)./g;
eV = cross(es,eu,4)./g;

%% metric tensor components
gss = dot(eS,eS,4);
gvv = dot(eV,eV,4);
guu = dot(eU,eU,4);
guv = dot(eU,eV,4);

%% contravariant B components
BU = (chir - Phir.*L_v)./g;
BU(1,:,:) = (chirr(1,:,:) - Phir(1,:,:).*L_sv(1,:,:)) ./ gs(1,:,:);

BV = Phir .* (1 - L_u)./g;
BV(1,:,:) = (Phir(1,:,:).* -L_su(1,:,:)) ./ gs(1,:,:);

%% partial derivatives of contravariant B components
BU_s = - gs./g .* (chir - Phir.*L_v) + (chirr - Phir .* L_sv)./g;
BU_u = - gu./g .* (chir - Phir.*L_v) + ( - Phir .* L_uv)./g;
BU_v = - gv./g .* (chir - Phir.*L_v) + ( - Phir .* L_vv)./g;

BV_s = (-gs./g .* Phir).*(1 - L_u) + Phir./g .* (-L_su);
BV_u = - gu./g .*Phir .* (1 - L_u) + Phir./g .* (-L_uu);
BV_v = - gv./g .*Phir .* (1 - L_u) + Phir./g .* (-L_uv);

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
JV = (Bu_s - Bs_u) ./ mu0 ./ g;

%% Magnitudes of direction vectors
mag_eS = sqrt(gss);
mag_beta = g .* sqrt((BV.^2).*guu + (BU.^2).*gvv - 2.*BV.*BU.*guv);

%% Magnitude of Force error (N/m3)
F_s = g .* (JV.*BU - JU.* BV) + presr;
F_beta = JS;
F = sqrt((F_s.^2).*gss + (F_beta.^2).*(mag_beta.^2));

%% Plot
plot_force_error
