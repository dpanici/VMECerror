clear all
d = read_vmec('wout_HELIOTRON_16x4x4.nc');

% with the grid defined, I can use the fourier coeffs to evaluate R,Z, L
% and its derivs on the grid. 

%then it is as simple as multiplying them together

% But, grid must be 3d (s,u,v), and to get s derivs must interpolate
% between the given flux surfaces from vmec

% Can interpolate by interpolating the fourier coefficients
% R_s can simply be a finite diff, we can define it everywhere with forward
% diff s=1 until we get to s=0, then can do backward difference there (or
% say it is zero?)
%R_s(s=s1, u,v)= (R(s=s2) - R(s=s1)) / (s2-s1)

% the corresponding s values are in the phi array (phi_b = 1, s = phi/phi_b
% = phi)
rsmnc = zeros(size(d.rmnc));

for i=2:d.ns-1
    rsmnc(:,i) = (d.rmnc(:,i+1) - d.rmnc(:,i-1)) / (d.phi(i+1) - d.phi(i-1));
end
rsmnc(:,1) = (d.rmnc(:,2) - d.rmnc(:,1)) / (d.phi(2) - d.phi(1));
rsmnc(:,end) = (d.rmnc(:,end) - d.rmnc(:,end-1)) / (d.phi(end) - d.phi(end-1));
rsmnc = s_deriv(d.rmnc,d);
figure()
R_s_u0v0 = zeros([1,d.ns]);
for i=1:d.ns
    R_s_u0v0(i) = sum(rsmnc(:,i));
end
plot(d.phi,R_s_u0v0)
title('dR/ds versus s at u=v=0')
xlabel('s')
ylabel('dR/ds')

figure() % dR/ds is multi-valued at s=0... I think dan had said this before
% will likely just average out the force balance error calculated at s=0,
% u for each v (unless I think of a better way to handle this)
u = linspace(0,2*pi,100);
R_s_s0v0 = zeros([1,length(u)]);
for i = 1:length(d.xm)
    R_s_s0v0 = R_s_s0v0 + rsmnc(i,1) .* cos(i.*u);
end
plot(u,R_s_s0v0)
title('dR/ds versus u at s=v=0')
xlabel('u')
ylabel('dR/ds')