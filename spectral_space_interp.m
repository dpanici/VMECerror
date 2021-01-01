function interp_coeffs = spectral_space_interp(coeffs,inputArg2)
%SPECTRAL_SPACE_INTERP accepts the fourier coeffs for a value (R,Z,R_u,
%etc) at the number of flux surfaces s from VMEC, and interpolates them for
%more flux surfaces
size_val = size(value);
space_dim = size(size_val);
deriv = zeros(size(value));
dim = length(var_wrt_to);
var = var_wrt_to;


if space_dim(2) ~= 3 % means dimV is one, and case is axisymmetric (s,u only, and d/dv = 0)
    if inputname(2) == 'v'
        return
    
    elseif inputname(2) =='s'
        dvar = var(2) - var(1);
        deriv(1,:) = (value(2,:) - value(1,:)) / (dvar);
        deriv(2:dim-1,:) = (value(3:dim,:) - value(1:dim-2,:)) /  (2*dvar);
        deriv(end,:) = (value(dim,:) - value(dim-1,:)) / (2*dvar);
        return
    
    elseif inputname(2) =='u'
        dvar = var(2) - var(1);
        deriv(:,1) = (value(:,2) - value(:,1)) / (dvar);
        deriv(:,2:dim-1) = (value(:,3:dim) - value(:,1:dim-2)) /  (2*dvar);
        deriv(:,dim) = (value(:,dim) - value(:,dim-1)) / (2*dvar);
        return
    end
end



if strcmp(deriv_method,'finite difference')
    deriv(1,:) = (value(2,:) - value(1,:)) / (dvar);
    deriv(2:dim-1,:) = (value(3:dim,:) - value(1:dim-2,:)) /  (2*dvar);
    deriv(end,:) = (value(dim,:) - value(dim-1,:)) / (2*dvar);
end

return

end


