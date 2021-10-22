function value_s_deriv=s_deriv(fourier_coeffs,data,deriv_method)
% accepts the 2D fourier coeffs of a value that is evaluated on the full
% mesh and the VMEC data structure, returns a 2D array of fourier coeffs that
% give the derivative wrt s of that value, again on the full mesh 
% reflects fourier coeffs over phi = 0, symmetrically bc they are cos, and 
% then uses central difference/splines to get s derivatives
global SPLINE_ORDER_SPAPI
global SPLINE_TENSION
global FACTOR_S
global SMOOTH_FACTOR
global POLYFIT_DEGREE
value_s_deriv = zeros(size(fourier_coeffs));
% reflect phi over axis for smoother derivatives at axis
refl_phi = -data.phi(data.ns):data.phi(data.ns)/(data.ns-1):data.phi(data.ns);

arr_size = size(fourier_coeffs);
refl_coeffs = zeros([arr_size(1) 2*data.ns-1]);
refl_coeffs(:,data.ns:end) = fourier_coeffs;
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
        value_s_deriv_refl(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) / (refl_phi(i+1) - refl_phi(i-1));
    end
    
    value_s_deriv_refl(:,end) = (refl_coeffs(:,end) - refl_coeffs(:,end-1)) / (refl_phi(end) - refl_phi(end-1));
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'finite difference 1st')% only use forward or backward diff
    for i=data.ns:2*data.ns-2
        value_s_deriv_refl(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i)) / (refl_phi(i+1) - refl_phi(i));
    end
    
    value_s_deriv_refl(:,end) = (refl_coeffs(:,end) - refl_coeffs(:,end-1)) / (refl_phi(end) - refl_phi(end-1));
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
elseif strcmp(deriv_method,'finite difference 4th') % 4th order accurate central diff using 5-point formula
    for i=data.ns:2*data.ns-3
        value_s_deriv_refl(:,i) = (refl_coeffs(:,i-2) + 8*refl_coeffs(:,i+1) - 8*refl_coeffs(:,i-1) - refl_coeffs(:,i+2)) / 12 / (refl_phi(i+1) - refl_phi(i));
    end
    % get deriv of second to last point using central diff (2nd order
    % accurate)
    i=2*data.ns-2;
    value_s_deriv_refl(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) / (refl_phi(i+1) - refl_phi(i-1));
    %get deriv of last point using backwards diff (first order acc)
    value_s_deriv_refl(:,end) = (refl_coeffs(:,end) - refl_coeffs(:,end-1)) / (refl_phi(end) - refl_phi(end-1));
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);

%%
elseif strcmp(deriv_method,'factor difference')
    % factor rho^m then do finite diff
    % but only factor out to do the first few derivatives, not the entire
    % thing?
    % not reflecting

    [val,s_ind] = min(abs(refl_phi-FACTOR_S));

    % only do this for until s=0.1? then use non factored ones?
    for i=1:length(data.xm)
        refl_coeffs(i,1:data.ns-1) = refl_coeffs(i,1:data.ns-1) ./ sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i);
        refl_coeffs(i,data.ns+1:s_ind) = refl_coeffs(i,data.ns+1:s_ind) ./ sqrt(abs(refl_phi(data.ns+1:s_ind))).^data.xm(i);
    end
    % centered diff everywhere the first index after rho=0 (do not want to
    % use rho=0 value in any of these diffs)
    % point
    for i=data.ns+2:s_ind
        value_s_deriv_refl(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) ./ (refl_phi(i+1) - refl_phi(i-1));
    end
    i=data.ns; % do rho=0 point centered diff
    value_s_deriv_refl(:,data.ns) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) ./ (refl_phi(i+1) - refl_phi(i-1));
% do point right after rho=0 using forward diff so dont use rho=0 value
    i = data.ns+1;
    value_s_deriv_refl(:,data.ns+1) = (refl_coeffs(:,i+1) - refl_coeffs(:,i)) ./ (refl_phi(i+1) - refl_phi(i));
% do endpoint at rho=rho(s_ind)
    value_s_deriv_refl(:,s_ind) = (refl_coeffs(:,s_ind) - refl_coeffs(:,s_ind-1)) ./ (refl_phi(s_ind) - refl_phi(s_ind-1));
    % multiply back rho^m
    for i=1:length(data.xm)
        refl_coeffs(i,1:data.ns-1) = refl_coeffs(i,1:data.ns-1) .* sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i);
        refl_coeffs(i,data.ns+1:s_ind) = refl_coeffs(i,data.ns+1:s_ind) .* sqrt(abs(refl_phi(data.ns+1:s_ind))).^data.xm(i);
    end
    for i=1:length(data.xm)
        value_s_deriv_refl(i,1:data.ns-1) = value_s_deriv_refl(i,1:data.ns-1) .* sqrt(abs(refl_phi(1:data.ns-1))).^data.xm(i) + refl_coeffs(i,1:data.ns-1).*data.xm(i)./abs(refl_phi(1:data.ns-1))./2;
        value_s_deriv_refl(i,data.ns+1:s_ind) = value_s_deriv_refl(i,data.ns+1:s_ind) .* sqrt(abs(refl_phi(data.ns+1:s_ind))).^data.xm(i) + refl_coeffs(i,data.ns+1:s_ind).*data.xm(i)./abs(refl_phi(data.ns+1:s_ind))./2;
    end    

    % now we have values for s=0 -> s=0.1, lets get the rest of the values
    % like normal from finite diff
    for i=s_ind+1:2*data.ns-2
        value_s_deriv_refl(:,i) = (refl_coeffs(:,i+1) - refl_coeffs(:,i-1)) / (refl_phi(i+1) - refl_phi(i-1));
    end
    
    value_s_deriv_refl(:,end) = (refl_coeffs(:,end) - refl_coeffs(:,end-1)) / (refl_phi(end) - refl_phi(end-1));
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
%%
elseif strcmp(deriv_method,'spline')
    fit = spline(refl_phi,refl_coeffs);
    spline_deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(spline_deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
%%
elseif strcmp(deriv_method,'factor_spline')
    % factor out rho^m from data (except for rho=0) before spline fitting
    for i=1:length(data.xm)
    refl_coeffs(:,1:data.ns-1) = refl_coeffs(:,1:data.ns-1) / refl_phi(1:data.ns-1).^data.xm(i);
    refl_coeffs(:,data.ns+1:end) = refl_coeffs(:,data.ns+1:end) / refl_phi(data.ns+1:end).^data.xm(i);
    end
    % above are the R tilde ones
    fit = spline(refl_phi,refl_coeffs);

    spline_deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(spline_deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
%%    
elseif strcmp(deriv_method,'smooth_spline')
%     spline_fit = spline(refl_phi,refl_coeffs);
    [fit,p1] = csaps(refl_phi,refl_coeffs,SMOOTH_FACTOR); % smoothing spline, results in pretty bad looks near axis
    p1=p1
    spline_deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(spline_deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
%%
elseif strcmp(deriv_method,'tension_spline') %% here could either fit 4 points at a time naively or use the T-spline surface fit thing
    % kind of learning towards naive 4-point spline fitting though, I know
    % the derivative direction there at least. Maybe use both tho for
    % comparison
%     spline_fit = spline(refl_phi,refl_coeffs);
    first_deriv= true;
    second_deriv = false;
    t = SPLINE_TENSION;
    for i=1:length(data.xm)
        [~,deriv_out,~] = spline1d(refl_phi,refl_phi,refl_coeffs(i,1:end),[],[],t,first_deriv,second_deriv); % smoothing spline, results in pretty bad looks near axis
        value_s_deriv(i,1:end) = deriv_out(data.ns:end);
    end

elseif strcmp(deriv_method,'spapi') % we can change the order of the spline being used
    k = SPLINE_ORDER_SPAPI; % k-1 is the order of spline used
    fit = spapi(k,refl_phi,refl_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv_refl = fnval(deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);

    
elseif strcmp(deriv_method,'pchip')
    fit = pchip(refl_phi,refl_coeffs);
    deriv_1 = fnder(fit,1);
    value_s_deriv_refl = ppval(deriv_1,refl_phi);
    value_s_deriv = value_s_deriv_refl(:,data.ns:end);
    
   
elseif strcmp(deriv_method,'makima')
    fit = makima(refl_phi,refl_coeffs);
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
elseif strcmp(deriv_method,'poly') % global polyfit, but maybe a local piecewise would be better, going 2^n at a time
    for i=1:length(data.xm)
    fit = polyfit(refl_phi,refl_coeffs(i,:),POLYFIT_DEGREE);
    deriv_1 = polyder(fit);
    value_s_deriv_refl = polyval(deriv_1,refl_phi);
    value_s_deriv(i,:) = value_s_deriv_refl(data.ns:end);
    end
end
% for i=1:5
% f=figure;
% scatter(refl_phi,refl_coeffs(i,:))
% hold on
% plot(refl_phi,ppval(fit,refl_phi),'DisplayName',sprintf('%s Fit',deriv_method))
% xlabel('Phi')
% ylabel('Fourier Coefficient')
% title(sprintf('Fourier coefficient m = %s n = %s Fit',data.xm(i),data.xn(i)))
% uiwait(f)
% end
end