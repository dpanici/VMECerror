function value_s_deriv=s_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
% This method yields a multi-values derivative at the s=0 magnetic axis
value_s_deriv = zeros(size(fourier_coeffs));


if strcmp(deriv_method,'finite difference')
    for i=2:data.ns-1
        value_s_deriv(:,i) = (fourier_coeffs(:,i+1) - fourier_coeffs(:,i-1)) / (data.phi(i+1) - data.phi(i-1));
    end
    
    value_s_deriv(:,1) = (fourier_coeffs(:,2) - fourier_coeffs(:,1)) / (data.phi(2) - data.phi(1));
    value_s_deriv(:,end) = (fourier_coeffs(:,end) - fourier_coeffs(:,end-1)) / (data.phi(end) - data.phi(end-1));

elseif strcmp(deriv_method,'spline')
    spline_fit = spline(data.phi,fourier_coeffs);
    spline_deriv_1 = fnder(spline_fit,1);
    value_s_deriv = ppval(spline_deriv_1,data.phi);
end

end