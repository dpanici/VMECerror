% Create the fourier coeffs for the derivatives of R,Z, and lambda

% d = read_vmec('wout_HELIOTRON_16x4x4.nc');
d = data;
xn = d.xn;
xm = d.xm;

% matlabvmec uses (mu+nv) convention, so my fourer coeffs involving v derivs
% will be of the opposite sign to those calculated by matlabVMEC

% d.rmnc is fourier coefs for R
% d.zmns is fourier coefs for Z
% indexed first by fourier mode, second by radial coordinate s
% so for example, R(:,end) returns all cos fourier coeffs of the outermost
% flux surface defined by s=1
% in vmec outputs rmnc, etc, s goes from 1/ns to the 1

%I think raxiscc and zaxiscs are the fourier coeffs of the magnetic axis
%(as the components they have do not have any poloidal modes)



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

% functional form is R = rmnc(s) * cos(m*u - n*v*nfp)
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

assert(isequal(-rvmns, d.rvmns)) % must assert the - of ours is equal bc matlabvmec uses (mu-nvnfp)

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

assert(isequal(-zvmnc,d.zvmnc))

% Z_uu = -zmns * m^2 * sin(m*u - n*v*nfp)
zuumns = -d.zmns .* (xm.^2)';

% Z_vv = -zmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
zvvmns = - d.zmns .* (xn.^2)';

% Z_uv = zmns * m * n * nfp * sin(m*u - n*v*nfp)
zuvmns = d.zmns .* xm' .* xn';

% Convert lambda from the half-grid onto the full grid
% 
% lmns = zeros(size(d.lmns));
% lmns(:,1) = 1.5*d.lmns(:,1) - 0.5*d.lmns(:,2);
% lmns(:,2:d.ns) = 0.5 * (d.lmns(:,1:d.ns-1) + d.lmns(:,2:d.ns));
% lmns(:,end) = 2 * d.lmns(end) - d.lmns(d.ns-2);

%L_u = lmns * m * cos(m*u-n*v*nfp)
lumnc = d.lmns .* xm';

% L_v = -lmns * n * nfp * cos(m*u - n*v*nfp)
lvmnc = -d.lmns .* xn';

% L_uu = -lmns * m^2 * sin(m*u - n*v*nfp)
luumns = -d.lmns .* (xm.^2)';

% L_vv = -lmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
lvvmns = - d.lmns .* (xn.^2)';

% L_uv = lmns * m * n * nfp * sin(m*u - n*v*nfp)
luvmns = d.lmns .* xm' .* xn'; % for some reason sign is flipped or something here

%% Numerical first radial derivatives
rsmnc = s_deriv(d.rmnc,d); %R_s
zsmns = s_deriv(d.zmns,d); %Z_s
rsvmns = s_deriv(rvmns,d); %R_sv
rsumns = s_deriv(d.rumns,d); %R_su
zsvmnc = s_deriv(zvmnc,d); %Z_sv
zsumnc = s_deriv(d.zumnc,d); %Z_su

rsuumnc = s_deriv(ruumnc,d); %R_suu
zsuumns = s_deriv(zuumns,d); %Z_suu
rsuvmnc = s_deriv(ruvmnc,d); %R_suv
zsuvmns = s_deriv(zuvmns,d); %Z_suv

rsvvmnc = s_deriv(rvvmnc,d); %R_svv
zsvvmns = s_deriv(zvvmns,d); %Z_svv



lsumnc = s_deriv(lumnc,d); %L_su
lsvmnc = s_deriv(lvmnc,d); %L_sv
lsvvmns = s_deriv(lvvmns,d); %L_svv
%% Numerical second radial derivatives
rssmnc = s2_deriv(d.rmnc,d); %R_ss
zssmns = s2_deriv(d.zmns,d); %Z_ss

rssumns = s2_deriv(d.rumns,d); %R_ssu
zssumnc = s2_deriv(d.zumnc,d); %Z_ssu
rssvmns = s2_deriv(rvmns,d); %R_ssv
zssvmnc = s2_deriv(zvmnc,d); %Z_ssv

rssuvmnc = s2_deriv(ruvmnc,d); %R_ussv
zssuvmns = s2_deriv(zuvmns,d); %Z_ussv

rssuumnc = s2_deriv(ruumnc,d); %R_ssuu
zssuumns = s2_deriv(zuumns,d); %Z_ssuu

lssvmnc = s2_deriv(lvmnc,d); %L_ssv

%% Numerical third radial derivatives
rsssmnc = s3_deriv(d.rmnc,d); %R_sss
zsssmns = s3_deriv(d.zmns,d); %Z_sss

rsssumns = s3_deriv(d.rumns,d); %R_sssu
zsssumnc = s3_deriv(d.zumnc,d); %Z_sssu
