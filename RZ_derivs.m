% Create the fourier coeffs for the derivatives of R,Z, and lambda

% d = read_vmec('wout_HELIOTRON_16x4x4.nc');
d = data;
xn = -d.xn; % d.xn returned from read_vmec assumes (mu + nv), but here I am assuming (mu - nv) so must flip the sign
xm = d.xm;


% d.rmnc is fourier coefs for R
% d.zmns is fourier coefs for Z
% indexed first by fourier mode, second by radial coordinate s
% so for example, R(:,1) returns all cos fourier coeffs of the outermost
% flux surface defined by s(1)
% in vmec outputs rmnc, etc, s goes from 1 until just before the magnetic
% axis (i.e. rmnc(:,end) is NOT the magnetic axis fourier coeffs)

% need to first be able to get the derivs wrt u,v at the discrete s
% surfaces we have

% grid we will evaluate errors on will be a 3D array index by (s,u,nfp*v),
% 0 <= s <= 1, 0 <= u <= 2pi, 0 <= v <= 2pi
% if I have grid as nfp*v, I should be able to plot just one field period,
% no need to do multiple field periods. Though I think I do need to keep
% nfp in the derivatives? maybe not, unsure of this

% may need to have a second one where s grid is defined by the discrete
% surfaces given by VMEC and use that one to evaluate the derivatives to
% interpolate onto the larger volumetric grid

%functional form is R = rmnc(s) * cos(m*u - n*v*nfp)
% and               Z = zmns(s) * sin(m*u - n*v*nfp)
% this is the fxn form for the fourier series R and Z for most stellarator
% symmetric shapes(i.e. R is always defined by cos and Z by sin because of
% the symmetries they have), since Z is negative when u = -u and R stays
% the same when u = -u, so R is cos and Z is sin

%denote derivatives with a subscript
% rumns is partial deriv of R wrt u sine coeffs, with same mode ordering in the rows as in xm,xn etc
% matlabVMEC's xn is actually xn*nfp, so do not need to use nfp in doing these derivatives 

%% Angular Derivatives functional forms
% R_u = -rmnc * m * sin(m*u - n*v*nfp)
rumns = -d.rmnc .* xm'; % 2D array of the coeffs of the sine terms that are the deriv of R wrt u, first index gives fourier coeff, second gives s index

assert(isequal(rumns, d.rumns)) % check against VMEC output  which also calculates rumns

% R_v =  rmnc * n*nfp*sin(m*u - n*v*nfp)
rvmns = d.rmnc .* xn';%.* d.nfp;

assert(isequal(rvmns, d.rvmns)) % this assertion was failing until I changed xn = -d.xn

% R_uu = -rmnc * m^2 * cos(m*u - n*v*nfp)
ruumnc = -d.rmnc .* (xm.^2)';

% R_vv = -rmnc * n^2 * nfp^2 * cos(m*u - n*v*nfp)
rvvmnc = -d.rmnc .* (xn.^2)';

% R_uv = rmnc * m * n * nfp * cos(m*u - n*v*nfp)
ruvmnc = d.rmnc .*xm' .* xn';

% Z_u = zmns * m * cos(m*u - n*v*nfp)
zumnc = d.zmns .* xm';

assert(isequal(zumnc,d.zumnc))

% Z_v = -zmns * n * nfp * cos(m*u - n*v*nfp)
zvmnc = -d.zmns .* xn';

assert(isequal(zvmnc,d.zvmnc))

% Z_uu = -zmns * m^2 * sin(m*u - n*v*nfp)
zuumns = -d.zmns .* (xm.^2)';

% Z_vv = -zmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
zvvmns = - d.zmns .* (xn.^2)';

% Z_uv = zmns * m * n * nfp * sin(m*u - n*v*nfp)
zuvmns = d.zmns .* xm' .* xn';

%L_u = lmns * m * cos(m*u-n*v*nfp)
lumnc = d.lmns .* xm';

% L_v = -lmns * n * nfp * cos(m*u - n*v*nfp)
lvmnc = -d.lmns .* xn';

% L_uu = -lmns * m^2 * sin(m*u - n*v*nfp)
luumns = -d.lmns .* (xm.^2)';

% L_vv = -lmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
lvvmns = - d.lmns .* (xn.^2)';

% L_uv = lmns * m * n * nfp * sin(m*u - n*v*nfp)
luvmns = d.lmns .* xm' .* xn';

%% Numerical first radial derivatives
rsmnc = s_deriv(d.rmnc,d);
zsmns = s_deriv(d.zmns,d);
rsvmns = s_deriv(d.rvmns,d);
rsumns = s_deriv(d.rumns,d);
zsvmnc = s_deriv(d.zvmnc,d);
zsumnc = s_deriv(d.zumnc,d);


lsumnc = s_deriv(lumnc,d);
lsvmnc = s_deriv(lvmnc,d);

%% Numerical second radial derivatives
rssmnc = s2_deriv(d.rmnc,d);
zssmns = s2_deriv(d.zmns,d);
