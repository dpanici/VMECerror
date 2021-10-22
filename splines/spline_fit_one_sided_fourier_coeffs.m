function [s_deriv,s2_deriv] = spline_fit_one_sided_fourier_coeffs(coeffs,data,~)
%SPLINE_FIT fits given x,y data using interior points as knots
% with a cubic spline one sided basis using not a knot or natural end conditions 
% currently only handles cubic splines
%   x_points is the x data points, the interior points of which will be
%       knots, this must be at least 4 data points(I think)
%   y_points are the y data points
%   end condition is either 'natural' or 'not-a-knot'
%       'natural' makes 2nd derivs of end points equal to zero
%       'not-a-knot' makes 3rd derivs of end points equal to zero
global MY_SPLINE_END_CONDITION

s_deriv = zeros(size(coeffs));
s2_deriv = zeros(size(coeffs));

for mode_ind=1:length(data.xm)
    y = coeffs(mode_ind,:);
    
    basis_coeffs= spline_fit_one_sided(data.phi, y, MY_SPLINE_END_CONDITION); 

    s_deriv(mode_ind,:)=Eval_Spline_One_Sided_Basis(data.phi,data.phi,basis_coeffs,1);

    s2_deriv(mode_ind,:)=Eval_Spline_One_Sided_Basis(data.phi,data.phi,basis_coeffs,2);
    
end



end

