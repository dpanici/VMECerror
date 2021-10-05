function value_s2_deriv=s2_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the data structure, returns a 2D array of fourier coeffs that
% give the 2nd derivative wrt s of that value, again on the full mesh
% Uses central difference for central points, and forward/backward diff for
% endpoints
global SPLINE_ORDER_SPAPI
global FACTOR_S
global SMOOTH_FACTOR
value_s2_deriv = zeros(size(data.rmnc));

refl_phi = -data.phi(data.ns):data.phi(data.ns)/(data.ns-1):data.phi(data.ns);
arr_size = size(fourier_coeffs);
refl_coeffs = zeros([arr_size(1) 2*data.ns-1]);
refl_coeffs(:,data.ns:end) = fourier_coeffs;
value_s2_deriv_r = zeros(size(refl_coeffs));
value_s2_deriv_refl_fac = zeros(size(refl_coeffs));

value_s_deriv_refl_fac = zeros(size(refl_coeffs));
value_s_deriv_refl = zeros(size(refl_coeffs));
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
elseif strcmp(deriv_method,'factor difference')
    % factor rho^m then do finite diff
    % but only factor out to do the first few derivatives, not the entire
    % thing?
    % not reflecting
    %% get first factored deriv
    [val,s_ind] = min(abs(refl_phi-FACTOR_S))
    % only do this for until s=0.1? then use non factored ones?
    for i=1:length(data.xm)
        refl_coeffs(i,1:data.ns-1) = refl_coeffs(i,1:data.ns-1) ./ sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i);
        refl_coeffs(i,data.ns+1:s_ind+2) = refl_coeffs(i,data.ns+1:s_ind+2) ./ sqrt(abs(refl_phi(data.ns+1:s_ind+2))).^data.xm(i);
    end
    % centered diff everywhere the first index after rho=0 (do not want to
    % use rho=0 value in any of these diffs)
    % point
    for i=data.ns+2:s_ind
        value_s_deriv_refl_fac(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) ./ (refl_phi(i+1) - refl_phi(i-1));
    end
    i=data.ns; % do rho=0 point centered diff
    value_s_deriv_refl_fac(:,data.ns) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) ./ (refl_phi(i+1) - refl_phi(i-1));
% do point right after rho=0 using forward diff so dont use rho=0 value
    i = data.ns+1;
    value_s_deriv_refl_fac(:,data.ns+1) = (refl_coeffs(:,i+1) - refl_coeffs(:,i)) ./ (refl_phi(i+1) - refl_phi(i));
% do endpoint at rho=rho(s_ind)
%     value_s_deriv_refl_fac(:,s_ind) = (refl_coeffs(:,s_ind) - refl_coeffs(:,s_ind-1)) ./ (refl_phi(s_ind) - refl_phi(s_ind-1));
%% do factored second deriv
    % centered diff everywhere the first index after rho=0 (do not want to
    % use rho=0 value in any of these diffs)
    % point
    %% centered diff O(h^2) acc
%     for i=data.ns+2:s_ind
%         value_s2_deriv_refl_fac(:,i) = (refl_coeffs(:,i+1) + refl_coeffs(:,i-1) - 2 * refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i))^2;
%     end
%     i=data.ns; % do rho=0 point centered diff
%     value_s2_deriv_refl_fac(:,data.ns) = (refl_coeffs(:,i+1) + refl_coeffs(:,i-1) - 2 * refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i))^2;
%     i = data.ns+1;
%     value_s2_deriv_refl_fac(:,data.ns+1) = (refl_coeffs(:,i+1) + refl_coeffs(:,i-1) - 2 * refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i))^2;
%% centered diff O(h^4) acc
    for i=data.ns+2:s_ind
        value_s2_deriv_refl_fac(:,i) = (-refl_coeffs(:,i+2) + 16.*refl_coeffs(:,i+1) - 30 .*refl_coeffs(:,i) + 16 .*refl_coeffs(:,i-1) -refl_coeffs(:,i-2)) / 12 /(refl_phi(i+1) - refl_phi(i))^2;
    end


%% get non factored first deriv
    % multiply back rho^m
    for i=1:length(data.xm)
        refl_coeffs(i,1:data.ns-1) = refl_coeffs(i,1:data.ns-1) .* sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i);
        refl_coeffs(i,data.ns+1:s_ind+2) = refl_coeffs(i,data.ns+1:s_ind+2) .* sqrt(abs(refl_phi(data.ns+1:s_ind+2))).^data.xm(i);
    end
    for i=1:length(data.xm)
        value_s_deriv_refl_fac(i,1:data.ns-1) = value_s_deriv_refl_fac(i,1:data.ns-1) .* sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i) + refl_coeffs(i,1:data.ns-1).*data.xm(i)./abs(refl_phi(1:data.ns-1))./2;
        value_s_deriv_refl_fac(i,data.ns+1:s_ind) = value_s_deriv_refl_fac(i,data.ns+1:s_ind) .* sqrt(abs(refl_phi(data.ns+1:s_ind))).^data.xm(i) + refl_coeffs(i,data.ns+1:s_ind).*data.xm(i)./abs(refl_phi(data.ns+1:s_ind))./2;
    end    
%% finally, get non factored second deriv for s=0 -> s=0.1

%     for i=1:length(data.xm)
%         m = data.xm(i);
%         term1 = data.xm(i)/2 .* abs(refl_phi(1:data.ns-1)).^(m/2-1) .* value_s_deriv_refl_fac(i,1:data.ns-1);
%         term2 = abs(refl_phi(1:data.ns-1)).^(m/2).*value_s2_deriv_refl_fac(i,1:data.ns-1);
%         term3 = m ./2 .*abs(refl_phi(1:data.ns-1)).^(-1) .* value_s_deriv_refl(i,1:data.ns-1);
%         term4 =  - m./2 .*abs(refl_phi(1:data.ns-1)).^(-2) .*  refl_coeffs(i,1:data.ns-1);
%         value_s2_deriv_r(i,1:data.ns-1) = term1    +     term2     +    term3  + term4;  
%         term1 = data.xm(i)./2 .* abs(refl_phi(data.ns+1:end)).^(m/2-1) .* value_s_deriv_refl_fac(i,data.ns+1:end);
%         term2 = abs(refl_phi(data.ns+1:end)).^(m/2).*value_s2_deriv_refl_fac(i,data.ns+1:end);
%         term3 = m ./2 .*abs(refl_phi(data.ns+1:end)).^(-1) .* value_s_deriv_refl(i,data.ns+1:end);
%         term4 = - m ./2 .*abs(refl_phi(data.ns+1:end)).^(-2) .*  refl_coeffs(i,data.ns+1:end);
%         value_s2_deriv_r(i,data.ns+1:end) =  term1  +     term2     +     term3     + term4;  
%     end    
% try other way
    for i=1:length(data.xm)
        m = data.xm(i);
        term1 = (m/2).*(m/2-1) .* abs(refl_phi(1:data.ns-1)).^(-2) .* refl_coeffs(i,1:data.ns-1);
        term2 = (m/2).*abs(refl_phi(1:data.ns-1)).^(m/2-1).*value_s_deriv_refl_fac(i,1:data.ns-1);
        term3 = m ./2 .*abs(refl_phi(1:data.ns-1)).^(m/2-1) .* value_s_deriv_refl(i,1:data.ns-1);
        term4 =  m./2 .*abs(refl_phi(1:data.ns-1)).^(m/2) .*  value_s2_deriv_refl_fac(i,1:data.ns-1);
        value_s2_deriv_r(i,1:data.ns-1) = term1    +     term2     +    term3  + term4;  
        term1 = (m/2).*(m/2-1) .* abs(refl_phi(data.ns+1:s_ind)).^(-2) .* refl_coeffs(i,data.ns+1:s_ind);
        term2 = (m/2).*abs(refl_phi(data.ns+1:s_ind)).^(m/2-1).*value_s_deriv_refl_fac(i,data.ns+1:s_ind);
        term3 = m ./2 .*abs(refl_phi(data.ns+1:s_ind)).^(m/2-1) .* value_s_deriv_refl(i,data.ns+1:s_ind);
        term4 =  m./2 .*abs(refl_phi(data.ns+1:s_ind)).^(m/2) .*  value_s2_deriv_refl_fac(i,data.ns+1:s_ind);
        value_s2_deriv_r(i,data.ns+1:s_ind) =  term1  +     term2     +     term3     + term4;  
    end    
    %% overleaf deriv
    for i=1:length(data.xm)
        m = data.xm(i);
        ss = abs(refl_phi);
        term1 = (m/2).*(m/2-1) .* ss(1:data.ns-1).^(-2) .* refl_coeffs(i,1:data.ns-1);
        term2 = (m).*ss(1:data.ns-1).^(m/2-1) .*value_s_deriv_refl_fac(i,1:data.ns-1);
%         term3 = m ./2 .*abs(refl_phi(1:data.ns-1)).^(m/2-1) .* value_s_deriv_refl(i,1:data.ns-1);
        term4 =  abs(refl_phi(1:data.ns-1)).^(m/2) .*  value_s2_deriv_refl_fac(i,1:data.ns-1);
        value_s2_deriv_r(i,1:data.ns-1) = term1    +     term2     + term4;  
        term1 = (m/2).*(m/2-1) .* ss(data.ns+1:s_ind).^(-2) .* refl_coeffs(i,data.ns+1:s_ind);
        term2 = (m).*ss(data.ns+1:s_ind).^(m/2-1) .*value_s_deriv_refl_fac(i,data.ns+1:s_ind);
%         term3 = m ./2 .*abs(refl_phi(1:data.ns-1)).^(m/2-1) .* value_s_deriv_refl(i,1:data.ns-1);
        term4 =  abs(refl_phi(data.ns+1:s_ind)).^(m/2) .*  value_s2_deriv_refl_fac(i,data.ns+1:s_ind);

        value_s2_deriv_r(i,data.ns+1:s_ind) =  term1  +     term2     +     term3     + term4;  
    end    
 %% get second deriv with centered for s=0.1 -> s=1
    for i=s_ind+2:2*data.ns-2
        value_s2_deriv_r(:,i) = (refl_coeffs(:,i+1) + refl_coeffs(:,i-1) - 2 * refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i))^2;
    end
    value_s2_deriv_r(:,end) = (fourier_coeffs(:,end) - 2*fourier_coeffs(:,end-1) + fourier_coeffs(:,end-2))  / (data.phi(end) - data.phi(end-1))^2;
 
    value_s2_deriv = value_s2_deriv_r(:,data.ns:end);
%%
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
end

end