function [y_interp] = Eval_Spline_One_Sided_Basis(x,x_data,coeffs,deriv_order)
%EVAL_SPLINE_ONE_SIDED_BASIS Takes coeffs of a one sided spline basis and
%evaluates the interpolated values (ie sums over coeffs * basis)
% assumes only simple knots were used (multiplicity of 1)
%   x is the points to evaluate the one sided basis at
%   x_data is the origial data points used (needed to get the knots)
%   coeffs are the coeffs found by spline_fit_one_sided.m
%   deriv_order is either 0,1,2,3

x_knots = x_data(2:end-1);
k = length(x_knots);
m = length(coeffs)-k; % m+k - k = m 
y_interp = zeros(size(x));
for p=1:m+k
    y_interp = y_interp + coeffs(p).*spl_1sided_basis_fxn(x,x_knots,p,m,deriv_order);
end


end

