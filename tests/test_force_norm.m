function test_force_norm(source_path, results_path)
% Test the actual norm expression against analytic cylindrical curl and force.
% No wout file, matlabVMEC, numerical differentiation or equilibrium solve.
if nargin < 1
    source_path = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fxns', 'force_error.m');
end
if nargin < 2
    results_path = '';
end
source = fileread(source_path);
first = strfind(source, 'F_beta = JS;');
last = strfind(source, '% calculate energy');
assert(numel(first) == 1 && numel(last) == 1 && last > first);
norm_expression = source(first + length('F_beta = JS;'):last - 1);

% Columns: elongation, shaping, radial angle shear, poloidal orientation.
shapes = [1, 0, 0, 1; 1, 0, 0.8, 1; 1.7, 0.25, 0.6, 1; 1.7, 0.25, 0.6, -1];
labels = {'circle_orthogonal', 'circle_sheared', 'shaped_sheared', 'shaped_reversed_chart'};
% Toroidal flux slope, poloidal flux slope, curvature, pressure slope/curvature.
profiles = [2, 1, 0.2, 0.7, -0.1; -2, 1, 0.2, 0.7, -0.1; ...
            2, -1, -0.2, 0.7, -0.1; 2, 1, 0.2, 0, 0];
fid = -1;
if ~isempty(results_path)
    fid = fopen(results_path, 'w');
    assert(fid >= 0);
    fprintf(fid, 'case,profile,points,max_relative_error,max_absolute_error,max_cartesian_curl_relative_error\n');
end
failed = 0;
for shape = 1:size(shapes, 1)
    for profile = 1:size(profiles, 1)
        [relative_error, absolute_error, points, curl_error] = check_case(norm_expression, shapes(shape,:), profiles(profile,:));
        fprintf('%s profile%d: max_relative=%.12g max_absolute=%.12g points=%d curl=%.12g\n', ...
                labels{shape}, profile, relative_error, absolute_error, points, curl_error);
        if fid >= 0
            fprintf(fid, '%s,%d,%d,%.17g,%.17g,%.17g\n', labels{shape}, profile, points, relative_error, absolute_error, curl_error);
        end
        failed = failed + (relative_error > 3e-12);
    end
end
if fid >= 0
    fclose(fid);
end
assert(failed == 0, '%d analytic force-norm cases failed', failed);
end

function [relative_error, absolute_error, points, curl_error] = check_case(norm_expression, shape, profile)
kappa = shape(1); delta = shape(2); shear = shape(3); orientation = shape(4);
flux_slope = profile(1); chi_slope = profile(2); chi_curvature = profile(3);
pressure_slope = profile(4); pressure_curvature = profile(5);
s_grid = [0.05, 0.2, 0.6, 0.9];
u_grid = (0:23) * (2*pi/24);
ns = length(s_grid); nu = length(u_grid);
es = zeros(ns, nu, 1, 3); eu = es; ev = es; eS = es;
gSS = zeros(ns, nu); new_mag_beta = gSS; F_s = gSS; F_beta = gSS;
BU = gSS; BV = gSS; expected = gSS;
curl_error = 0;
for i = 1:ns
    for j = 1:nu
        s = s_grid(i); angle = orientation*u_grid(j) + shear*s;
        ct = cos(angle); st = sin(angle);
        % Level surfaces s=(x^2+(Z/kappa)^2)/(1+delta*x), x=R-3.
        radius = (s*delta*ct + sqrt(s*s*delta*delta*ct*ct + 4*s))/2;
        denominator = 2*radius - s*delta*ct;
        radius_s = (1 + delta*ct*radius - s*delta*st*radius*shear)/denominator;
        radius_u = -s*delta*st*radius*orientation/denominator;
        R = 3 + radius*ct; Z = -kappa*radius*st;
        rs = [radius_s*ct-radius*st*shear, 0, -kappa*(radius_s*st+radius*ct*shear)];
        ru = [radius_u*ct-radius*st*orientation, 0, -kappa*(radius_u*st+radius*ct*orientation)];
        rv = [0, R, 0];
        jacobian = dot(rs, cross(ru, rv));
        grad_s = cross(ru, rv)/jacobian;
        grad_u = cross(rv, rs)/jacobian;
        grad_v = cross(rs, ru)/jacobian;
        x = R-3; D = 1+delta*x; numerator = x*x+Z*Z/(kappa*kappa);
        s_R = 2*x/D-delta*numerator/(D*D);
        s_Z = 2*Z/(kappa*kappa*D);
        s_RR = 2/D-4*delta*x/(D*D)+2*delta*delta*numerator/(D*D*D);
        s_ZZ = 2/(kappa*kappa*D);
        assert(abs(numerator/D-s) < 2e-14);
        assert(norm(grad_s-[s_R,0,s_Z]) < 2e-13);
        cp = chi_slope+2*chi_curvature*s;
        cpp = 2*chi_curvature;
        pp = pressure_slope+2*pressure_curvature*s;
        % B=grad(v) cross grad(chi)+grad(psi_t) cross grad(u), mu0=1.
        % Bphi follows from differentiating u=orientation*(atan2(-Z/kappa,x)-shear*s).
        bphi = orientation*flux_slope*(2+delta*x)/(kappa*D*D);
        bphi_R = orientation*flux_slope*delta/kappa*(1/(D*D)-2*(2+delta*x)/(D*D*D));
        B = [cp*s_Z/R, bphi, -cp*s_R/R];
        % Independent cylindrical curl, including the rotating-basis Bphi/R term.
        J = [0, (cpp*(s_R*s_R+s_Z*s_Z)+cp*(s_RR+s_ZZ-s_R/R))/R, bphi_R+bphi/R];
        if i == 2 && mod(j-1, 6) == 0
            coarse = cartesian_curl([R,0,Z], shape, profile, 2e-3);
            fine = cartesian_curl([R,0,Z], shape, profile, 1e-3);
            ec = norm(coarse-J); ef = norm(fine-J);
            assert(ef < 1e-5*max(norm(J),1));
            assert(ef < 0.35*ec || ec < 1e-8, 'curl refinement ec=%g ef=%g R=%g Z=%g', ec,ef,R,Z);
            curl_error = max(curl_error, ef/max(norm(J),1));
        end
        grad_p = pp*[s_R,0,s_Z];
        physical_force = grad_p-cross(J,B);
        assert(abs(dot(B, grad_s)) < 3e-13);
        bu = dot(B, grad_u); bv = dot(B, grad_v);
        ju = dot(J, grad_u); jv = dot(J, grad_v); js = dot(J, grad_s);
        beta = cross(B, rs);
        es(i,j,1,:) = rs; eu(i,j,1,:) = ru; ev(i,j,1,:) = rv; eS(i,j,1,:) = grad_s;
        gSS(i,j) = dot(grad_s,grad_s); new_mag_beta(i,j) = norm(beta);
        F_s(i,j) = jacobian*(jv*bu-ju*bv)+pp; F_beta(i,j) = js;
        BU(i,j) = bu; BV(i,j) = bv; expected(i,j) = norm(physical_force);
    end
end

% Exercise the existing first-row reciprocal-gradient limit independently.
% This is a supplied-limit branch control, not a magnetic-axis convergence test.
temp1 = eS;
eS(1,:,:,:) = NaN;
eval(norm_expression);
assert(isequal(size(F), size(expected)) && all(isfinite(F(:))) && isreal(F));
absolute_error = max(abs(F(:)-expected(:)));
relative_error = max(abs(F(:)-expected(:))./max(expected(:),1));
points = numel(expected);
end

function J = cartesian_curl(position, shape, profile, h)
derivatives = zeros(3,3);
for coordinate = 1:3
    offset = zeros(1,3); offset(coordinate) = h;
    derivatives(coordinate,:) = (cartesian_field(position+offset,shape,profile)- ...
                                cartesian_field(position-offset,shape,profile))/(2*h);
end
J = [derivatives(2,3)-derivatives(3,2), derivatives(3,1)-derivatives(1,3), ...
     derivatives(1,2)-derivatives(2,1)];
end

function B = cartesian_field(position, shape, profile)
R = hypot(position(1),position(2)); Z = position(3); x = R-3;
kappa = shape(1); delta = shape(2); orientation = shape(4);
D = 1+delta*x; numerator = x*x+Z*Z/(kappa*kappa); s = numerator/D;
s_R = 2*x/D-delta*numerator/(D*D); s_Z = 2*Z/(kappa*kappa*D);
cp = profile(2)+2*profile(3)*s;
br = cp*s_Z/R; bz = -cp*s_R/R;
bphi = orientation*profile(1)*(2+delta*x)/(kappa*D*D);
cos_phi = position(1)/R; sin_phi = position(2)/R;
B = [br*cos_phi-bphi*sin_phi, br*sin_phi+bphi*cos_phi, bz];
end
