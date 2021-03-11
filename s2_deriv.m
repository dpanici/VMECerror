function value_s2_deriv=s2_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the 2nd derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
value_s2_deriv = zeros(size(data.rmnc));

refl_phi = -data.phi(data.ns):data.phi(data.ns)/(data.ns-1):data.phi(data.ns);
arr_size = size(fourier_coeffs);
refl_coeffs = zeros([arr_size(1) 2*data.ns-1]);
refl_coeffs(:,data.ns:end) = fourier_coeffs;

for i=1:length(data.xm)
    if mod(data.xm(i),2) == 0
        refl_coeffs(i,1:data.ns-1) = flip(fourier_coeffs(i,2:end));
    else % odd parity, negate the coeff when reflecting about phi = 0
        refl_coeffs(i,1:data.ns-1) = -flip(fourier_coeffs(i,2:end));
    end
end


if strcmp(deriv_method,'finite difference')
    for i=2:data.ns-1
        value_s2_deriv(:,i) = (fourier_coeffs(:,i+1) + fourier_coeffs(:,i-1) - 2 * fourier_coeffs(:,i)) / (data.phi(i+1) - data.phi(i))^2;
    end
    value_s2_deriv(:,1) = (fourier_coeffs(:,3) - 2*fourier_coeffs(:,2) + fourier_coeffs(:,1)) / (data.phi(2) - data.phi(1))^2;
    value_s2_deriv(:,end) = (fourier_coeffs(:,end) - 2*fourier_coeffs(:,end-1) + fourier_coeffs(:,end-2))  / (data.phi(end) - data.phi(end-1))^2;

elseif strcmp(deriv_method,'spline')
    spline_fit = spline(refl_phi,refl_coeffs);
% 	spline_fit = csaps(refl_phi,refl_coeffs,0.999999); % smoothing spline, results in pretty bad looks near axis
    spline_deriv_2 = fnder(spline_fit,2);
    value_s2_deriv_refl = ppval(spline_deriv_2,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'pchip')
    fit = pchip(data.phi,fourier_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv = ppval(deriv_2,data.phi);
elseif strcmp(deriv_method,'makima')
    fit = makima(data.phi,fourier_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv = ppval(deriv_2,data.phi);
% elseif strcmp(deriv_method,'chebfun') % requires chebfun toolbox, gives
% highly oscillatory results
%     for i=1:length(data.xm)
%         fit = chebfun(refl_coeffs(i,:)','equi');
%         deriv_2 = diff(fit,2);
%         value_s_deriv_refl = feval(deriv_2,refl_phi);
%         value_s2_deriv(i,:) = value_s_deriv_refl(data.ns:end);
%     end
end

end