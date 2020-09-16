function value_s3_deriv=s3_deriv(fourier_coeffs,data)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the read_vmec data structure, returns a 2D array of fourier coeffs that
% give the 3rdd derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
% This method yields a multi-valued derivative at the s=0 magnetic axis
value_s3_deriv = zeros(size(data.rmnc));

h = data.phi(2) - data.phi(1);

for i=3:data.ns-2
    value_s3_deriv(:,i) = (0.5.*fourier_coeffs(:,i+2) - fourier_coeffs(:,i+1) + fourier_coeffs(:,i-1) - 0.5 .* fourier_coeffs(:,i-2)) / (h)^3;
end
for i = 1:2
value_s3_deriv(:,i) = (-fourier_coeffs(:,i) + 3*fourier_coeffs(:,i+1) - 3 * fourier_coeffs(:,i+2) + fourier_coeffs(:,i+3) ) / (h)^3;
end
for i=data.ns-1:data.ns
value_s3_deriv(:,end) = (fourier_coeffs(:,i) - 3*fourier_coeffs(:,i-1) + 3 * fourier_coeffs(:,i-2) - fourier_coeffs(:,i-3) ) / (h)^3;
end
end