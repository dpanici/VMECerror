function value_s2_deriv=s2_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the 2nd derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
global SPLINE_ORDER_SPAPI
global SMOOTH_FACTOR
global POLYFIT_DEGREE
value_s2_deriv = zeros(size(data.rmnc));

refl_phi = -data.phi(data.ns):data.phi(data.ns)/(data.ns-1):data.phi(data.ns);
arr_size = size(fourier_coeffs);
refl_coeffs = zeros([arr_size(1) 2*data.ns-1]);
refl_coeffs(:,data.ns:end) = fourier_coeffs;
value_s2_deriv_r = zeros(size(refl_coeffs));

for i=1:length(data.xm)
    if mod(data.xm(i),2) == 0
        refl_coeffs(i,1:data.ns-1) = flip(fourier_coeffs(i,2:end));
    else % odd parity, negate the coeff when reflecting about phi = 0
        refl_coeffs(i,1:data.ns-1) = -flip(fourier_coeffs(i,2:end));
    end
end


if strcmp(deriv_method,'finite difference')
    for i=data.ns:2*data.ns-2
        value_s2_deriv_r(:,i) = (refl_coeffs(:,i+1) + refl_coeffs(:,i-1) - 2 * refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i))^2;
    end
    value_s2_deriv_r(:,end) = (fourier_coeffs(:,end) - 2*fourier_coeffs(:,end-1) + fourier_coeffs(:,end-2))  / (data.phi(end) - data.phi(end-1))^2;
    value_s2_deriv = value_s2_deriv_r(:,data.ns:end);
elseif strcmp(deriv_method,'finite difference 4th') % 4th order acc finite difference for 2nd deriv
    for i=data.ns:2*data.ns-3
        value_s2_deriv_r(:,i) = (- refl_coeffs(:,i-2) + 16 * refl_coeffs(:,i-1) - 30 * refl_coeffs(:,i) + 16 * refl_coeffs(:,i+1) - refl_coeffs(:,i+2)) / 12/(refl_phi(i+1) - refl_phi(i))^2;
    end
    % second to last pt with 3rd order accurate
    i = 2*data.ns-2;
    value_s2_deriv_r(:,i) = (11 * refl_coeffs(:,i-3) - 20 * refl_coeffs(:,i-2) + 6 * refl_coeffs(:,i-1) + 4 * refl_coeffs(:,i) - refl_coeffs(:,i+1)) / 12/(refl_phi(i+1) - refl_phi(i))^2;
    % last with first order accurate
%     i = 2*data.ns-1:
    value_s2_deriv_r(:,i) = (fourier_coeffs(:,end) - 2*fourier_coeffs(:,end-1) + fourier_coeffs(:,end-2))  / (data.phi(end) - data.phi(end-1))^2;
    value_s2_deriv = value_s2_deriv_r(:,data.ns:end);
    

elseif strcmp(deriv_method,'spline')
    spline_fit = spline(refl_phi,refl_coeffs);
    spline_deriv_2 = fnder(spline_fit,2);
    value_s2_deriv_refl = ppval(spline_deriv_2,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'smooth_spline')
    [spline_fit,p2] = csaps(refl_phi,refl_coeffs,SMOOTH_FACTOR);
    p2=p2
    spline_deriv_2 = fnder(spline_fit,2);
    value_s2_deriv_refl = ppval(spline_deriv_2,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'spapi') % we can change the order of the spline being used
    k = SPLINE_ORDER_SPAPI; % k-1 is the order of spline used
    fit = spapi(k,refl_phi,refl_coeffs);
    deriv = fnder(fit,2);
    value_s2_deriv_refl = fnval(deriv,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'pchip')
    fit = pchip(data.phi,fourier_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv = ppval(deriv_2,data.phi);
elseif strcmp(deriv_method,'makima')
    fit = makima(refl_phi,refl_coeffs);
    deriv_2 = fnder(fit,2);
    value_s2_deriv_refl = ppval(deriv_2,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'tension_spline') %% here could either fit 4 points at a time naively or use the T-spline surface fit thing
    % kind of learning towards naive 4-point spline fitting though, I know
    % the derivative direction there at least. Maybe use both tho for
    % comparison
%     spline_fit = spline(refl_phi,refl_coeffs);
    spline_fit = spline(refl_phi,refl_coeffs);
    spline_deriv_2 = fnder(spline_fit,2);
    value_s2_deriv_refl = ppval(spline_deriv_2,refl_phi);
    value_s2_deriv = value_s2_deriv_refl(:,data.ns:end);
    % elseif strcmp(deriv_method,'chebfun') % requires chebfun toolbox, gives
% highly oscillatory results
%     for i=1:length(data.xm)
%         fit = chebfun(refl_coeffs(i,:)','equi');
%         deriv_2 = diff(fit,2);
%         value_s_deriv_refl = feval(deriv_2,refl_phi);
%         value_s2_deriv(i,:) = value_s_deriv_refl(data.ns:end);
%     end
elseif strcmp(deriv_method,'poly') % piecewise polyfit
    for i=1:length(data.xm)
    fit = polyfit(refl_phi,refl_coeffs(i,:),POLYFIT_DEGREE);
    deriv_1 = polyder(fit);
    deriv_2 = polyder(deriv_1);
    
    value_s2_deriv_refl = polyval(deriv_2,refl_phi);
    value_s2_deriv(i,:) = value_s2_deriv_refl(data.ns:end);
    end
end

end