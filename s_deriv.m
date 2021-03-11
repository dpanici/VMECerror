function value_s_deriv=s_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the VMEC data structure, returns a 2D array of fourier coeffs that
% give the derivative wrt s of that value, again on the full mesh 
% reflects fourier coeffs over phi = 0, symmetrically bc they are cos, and 
% then uses central difference/splines to get s derivatives
value_s_deriv = zeros(size(fourier_coeffs));
% reflect phi over axis for smoother derivatives at axis
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
    for i=data.ns:2*data.ns-2
        value_s_deriv(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) / (refl_phi(i+1) - refl_phi(i-1));
    end
    
    value_s_deriv(:,end) = (refl_coeffs(:,end) - refl_coeffs(:,end-1)) / (refl_phi(end) - refl_phi(end-1));

elseif strcmp(deriv_method,'spline')
    spline_fit = spline(refl_phi,refl_coeffs);
    spline_deriv_1 = fnder(spline_fit,1);
    value_s_deriv_refl = ppval(spline_deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'pchip')
    fit = pchip(data.phi,refl_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'makima')
    fit = makima(data.phi,refl_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
% elseif strcmp(deriv_method,'chebfun') %requires the chebfuntoolbox
% installed in local directory, gives highly oscillatory results
%     for i=1:length(data.xm)
%         fit = chebfun(refl_coeffs(i,:)',[refl_phi(1) refl_phi(end)],'equi');
%         deriv_1 = diff(fit);
%         value_s_deriv_refl = feval(deriv_1,refl_phi);
%         value_s_deriv(i,:) = value_s_deriv_refl(data.ns:end);
%     end
end
end