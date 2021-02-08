function value_s2_deriv=s2_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the 2nd derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
% This method yields a multi-valued derivative at the s=0 magnetic axis
value_s2_deriv = zeros(size(data.rmnc));
if strcmp(deriv_method,'finite difference')
    for i=2:data.ns-1
        value_s2_deriv(:,i) = (fourier_coeffs(:,i+1) + fourier_coeffs(:,i-1) - 2 * fourier_coeffs(:,i)) / (data.phi(i+1) - data.phi(i))^2;
    end
    value_s2_deriv(:,1) = (fourier_coeffs(:,3) - 2*fourier_coeffs(:,2) + fourier_coeffs(:,1)) / (data.phi(2) - data.phi(1))^2;
    value_s2_deriv(:,end) = (fourier_coeffs(:,end) - 2*fourier_coeffs(:,end-1) + fourier_coeffs(:,end-2))  / (data.phi(end) - data.phi(end-1))^2;

elseif strcmp(deriv_method,'spline')
    spline_fit = spline(data.phi,fourier_coeffs);
    spline_deriv_2 = fnder(spline_fit,2);
    value_s2_deriv = ppval(spline_deriv_2,data.phi);

elseif strcmp(deriv_method,'pchip')
    fit = pchip(data.phi,fourier_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv = ppval(deriv_2,data.phi);
elseif strcmp(deriv_method,'makima')
    fit = makima(data.phi,fourier_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv = ppval(deriv_2,data.phi);

end


end