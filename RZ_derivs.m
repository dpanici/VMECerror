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
% R_u = -rmnc * m * sin(m*u + n*v*nfp)
rumns = -d.rmnc .* xm'; % 2D array of the coeffs of the sine terms that are the deriv of R wrt u, first index gives fourier coeff, second gives s index

assert(isequal(rumns, d.rumns)) % check against VMEC output  which also calculates rumns

% R_v =  -rmnc * n*nfp*sin(m*u + n*v*nfp)
rvmns = -d.rmnc .* xn';%.* d.nfp;

assert(isequal(rvmns, d.rvmns)) 

% R_uu = -rmnc * m^2 * cos(m*u - n*v*nfp)
ruumnc = -d.rmnc .* (xm.^2)';

% R_vv = -rmnc * n^2 * nfp^2 * cos(m*u - n*v*nfp)
rvvmnc = -d.rmnc .* (xn.^2)';

% R_uv = rmnc * m * n * nfp * cos(m*u - n*v*nfp)
ruvmnc = -d.rmnc .*xm' .* xn';

% Z_u = zmns * m * cos(m*u - n*v*nfp)
zumnc = d.zmns .* xm';

assert(isequal(zumnc,d.zumnc))

% Z_v = -zmns * n * nfp * cos(m*u - n*v*nfp)
zvmnc = d.zmns .* xn';

assert(isequal(zvmnc,d.zvmnc))

% Z_uu = -zmns * m^2 * sin(m*u - n*v*nfp)
zuumns = -d.zmns .* (xm.^2)';

% Z_vv = -zmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
zvvmns = - d.zmns .* (xn.^2)';

% Z_uv = zmns * m * n * nfp * sin(m*u - n*v*nfp)
zuvmns = -d.zmns .* xm' .* xn';

% Convert lambda from the half-mesh onto the full mesh
%  just linear interpolation onto the full mesh from the half-mesh
lmns = zeros(size(d.lmns));
% lmns(:,1) = 1.5*d.lmns(:,1) - 0.5*d.lmns(:,2); % should be axis limit for lambda
% lmns(:,2:d.ns) = 0.5 * (d.lmns(:,1:d.ns-1) + d.lmns(:,2:d.ns));
% lmns(:,end) = 3/2 * d.lmns(:,end-1) - 0.5*d.lmns(:,end-2);
% figure
% plot(data.phi,lmns(3,:),'DisplayName','full mesh 1st way')
% hold on

if use_lambda_full_mesh
for i=2:d.ns
    lmns(:,i) = 1/2*(d.lmns(:,i-1)+d.lmns(:,i));
end
lmns(:,end) = 3/2 * d.lmns(:,end-1) - 0.5*d.lmns(:,end-2);
lmns(:,1) = 1.5*d.lmns(:,1) - 0.5*d.lmns(:,2);
lmns(:,1) = 1.5*d.lmns(:,1) - 0.5*d.lmns(:,2); % should be axis limit for lambda

else
lmns = d.lmns;
end
% plot(data.phi+(s(2)-s(1))/2,data.lmns(3,:),'x','DisplayName','half mesh')
% hold on
% plot(data.phi,lmns(3,:),'.','DisplayName','full mesh 2nd way')
% legend

%L_u = lmns * m * cos(m*u+n*v*nfp)
lumnc = lmns .* xm';

% L_v = lmns * n * nfp * cos(m*u + n*v*nfp)
lvmnc = lmns .* xn';

% L_uu = -lmns * m^2 * sin(m*u - n*v*nfp)
luumns = -lmns .* (xm.^2)';

% L_vv = -lmns * n^2 * nfp^2 * sin(m*u - n*v*nfp)
lvvmns = - lmns .* (xn.^2)';

% L_uv = lmns * m * n * nfp * sin(m*u - n*v*nfp)
luvmns = -lmns .* xm' .* xn';

%% Numerical first radial derivatives
if use_piecewise_lsq == false
rsmnc = s_deriv(d.rmnc,d,deriv_method); %R_s
zsmns = s_deriv(d.zmns,d,deriv_method); %Z_s
rsvmns = s_deriv(rvmns,d,deriv_method); %R_sv
rsumns = s_deriv(d.rumns,d,deriv_method); %R_su
zsvmnc = s_deriv(zvmnc,d,deriv_method); %Z_sv
zsumnc = s_deriv(d.zumnc,d,deriv_method); %Z_su

rsuumnc = s_deriv(ruumnc,d,deriv_method); %R_suu
zsuumns = s_deriv(zuumns,d,deriv_method); %Z_suu
rsuvmnc = s_deriv(ruvmnc,d,deriv_method); %R_suv
zsuvmns = s_deriv(zuvmns,d,deriv_method); %Z_suv

rsvvmnc = s_deriv(rvvmnc,d,deriv_method); %R_svv
zsvvmns = s_deriv(zvvmns,d,deriv_method); %Z_svv



lsumnc = s_deriv(lumnc,d,deriv_method); %L_su
lsvmnc = s_deriv(lvmnc,d,deriv_method); %L_sv
lsvvmns = s_deriv(lvvmns,d,deriv_method); %L_svv
end
% rsmnc = rho_deriv(d.rmnc,d,deriv_method); %R_s
% zsmns = rho_deriv(d.zmns,d,deriv_method); %Z_s
% rsvmns = rho_deriv(rvmns,d,deriv_method); %R_sv
% rsumns = rho_deriv(d.rumns,d,deriv_method); %R_su
% zsvmnc = rho_deriv(zvmnc,d,deriv_method); %Z_sv
% zsumnc = rho_deriv(d.zumnc,d,deriv_method); %Z_su
% 
% rsuumnc = rho_deriv(ruumnc,d,deriv_method); %R_suu
% zsuumns = rho_deriv(zuumns,d,deriv_method); %Z_suu
% rsuvmnc = rho_deriv(ruvmnc,d,deriv_method); %R_suv
% zsuvmns = rho_deriv(zuvmns,d,deriv_method); %Z_suv
% 
% rsvvmnc = rho_deriv(rvvmnc,d,deriv_method); %R_svv
% zsvvmns = rho_deriv(zvvmns,d,deriv_method); %Z_svv
% 
% lsumnc = rho_deriv(lumnc,d,deriv_method); %L_su
% lsvmnc = rho_deriv(lvmnc,d,deriv_method); %L_sv
% lsvvmns = rho_deriv(lvvmns,d,deriv_method); %L_svv




%% Numerical second radial derivatives
if use_piecewise_lsq == false
rssmnc = s2_deriv(d.rmnc,d,deriv_method); %R_ss
zssmns = s2_deriv(d.zmns,d,deriv_method); %Z_ss

rssumns = s2_deriv(d.rumns,d,deriv_method); %R_ssu
zssumnc = s2_deriv(d.zumnc,d,deriv_method); %Z_ssu
rssvmns = s2_deriv(rvmns,d,deriv_method); %R_ssv
zssvmnc = s2_deriv(zvmnc,d,deriv_method); %Z_ssv

rssuvmnc = s2_deriv(ruvmnc,d,deriv_method); %R_ussv
zssuvmns = s2_deriv(zuvmns,d,deriv_method); %Z_ussv

rssuumnc = s2_deriv(ruumnc,d,deriv_method); %R_ssuu
zssuumns = s2_deriv(zuumns,d,deriv_method); %Z_ssuu

lssvmnc = s2_deriv(lvmnc,d,deriv_method); %L_ssv
end
% try with taking s deriv of s deriv
% rssmnc = s_deriv(rsmnc,d,deriv_method); %R_ss
% zssmns = s_deriv(zsmns,d,deriv_method); %Z_ss
% 
% rssumns = s_deriv(rsumns,d,deriv_method); %R_ssu
% zssumnc = s_deriv(zsumnc,d,deriv_method); %Z_ssu
% rssvmns = s_deriv(rsvmns,d,deriv_method); %R_ssv
% zssvmnc = s_deriv(zsvmnc,d,deriv_method); %Z_ssv
% 
% rssuvmnc = s_deriv(rsuvmnc,d,deriv_method); %R_ussv
% zssuvmns = s_deriv(zsuvmns,d,deriv_method); %Z_ussv
% 
% rssuumnc = s_deriv(rsuumnc,d,deriv_method); %R_ssuu
% zssuumns = s_deriv(zsuumns,d,deriv_method); %Z_ssuu
% 
% lssvmnc = s_deriv(lsvmnc,d,deriv_method); %L_ssv

%% use piecewise least squares
% save the polynomial coeffs 
global POLY_LSQ_WINDOW_SIZE
POLY_LSQ_WINDOW_SIZE=8; % 16 was good and 6 poly order
global POLY_LSQ_ORDER
POLY_LSQ_ORDER = 5; % polynomial order
if use_piecewise_lsq
    [rsmnc,rssmnc] = least_squares_fit_coeffs(d.rmnc,d,deriv_method); %R_s
    [zsmns,zssmns] = least_squares_fit_coeffs(d.zmns,d,deriv_method); %Z_s
    [rsvmns,rssvmns] = least_squares_fit_coeffs(rvmns,d,deriv_method); %R_sv
    [rsumns,rssumns] = least_squares_fit_coeffs(d.rumns,d,deriv_method); %R_su
    [zsvmnc,zssvmnc] = least_squares_fit_coeffs(zvmnc,d,deriv_method); %Z_sv
    [zsumnc,zssumnc] = least_squares_fit_coeffs(d.zumnc,d,deriv_method); %Z_su

    [rsuumnc,rssuumnc] = least_squares_fit_coeffs(ruumnc,d,deriv_method); %R_suu
    [zsuumns,zssuumns] = least_squares_fit_coeffs(zuumns,d,deriv_method); %Z_suu
    [rsuvmnc,rssuvmnc] = least_squares_fit_coeffs(ruvmnc,d,deriv_method); %R_suv
    [zsuvmns,zssuvmns] = least_squares_fit_coeffs(zuvmns,d,deriv_method); %Z_suv

    [rsvvmnc,rssvvmnc] = least_squares_fit_coeffs(rvvmnc,d,deriv_method); %R_svv
    [zsvvmns,zssvvmns] = least_squares_fit_coeffs(zvvmns,d,deriv_method); %Z_svv

    [lsumnc,lssumnc] = least_squares_fit_coeffs(lumnc,d,deriv_method); %L_su
    [lsvmnc,lssvmnc] = least_squares_fit_coeffs(lvmnc,d,deriv_method); %L_sv
    [lsvvmns,lssvvmns] = least_squares_fit_coeffs(lvvmns,d,deriv_method); %L_svv
end




%% Numerical third radial derivatives
rsssmnc = s3_deriv(d.rmnc,d,deriv_method); %R_sss
zsssmns = s3_deriv(d.zmns,d,deriv_method); %Z_sss

rsssumns = s3_deriv(d.rumns,d,deriv_method); %R_sssu
zsssumnc = s3_deriv(d.zumnc,d,deriv_method); %Z_sssu
