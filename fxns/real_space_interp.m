function interp_val = real_space_interp(value,orig_grid,fine_grid,deriv_method)
%REAL_SPACE_INTERP Take a value in real space (defined on s,u,v grid) and
%return it interpolated onto a finer (in s) grid

size_val = size(value);
space_dim = size(size_val);
deriv = zeros(size(value));
dim = length(fine_grid);
var = fine_grid;


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

