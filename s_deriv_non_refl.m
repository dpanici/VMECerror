function value_s_deriv=s_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the derivative wrt s of that value, again on the full mesh
% Does not reflect phi across the axis
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

elseif strcmp(deriv_method,'pchip')
    fit = pchip(data.phi,fourier_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv = ppval(deriv_1,data.phi);
elseif strcmp(deriv_method,'makima')
    fit = makima(data.phi,fourier_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv = ppval(deriv_1,data.phi);

end

end