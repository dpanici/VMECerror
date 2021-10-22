function [coeffs] = spline_fit_one_sided(x_points,y_points,end_condition)
%SPLINE_FIT fits given x,y data using interior points as knots
% with a cubic spline one sided basis using not a knot or natural end conditions 
% currently only handles cubic splines
%   x_points is the x data points, the interior points of which will be
%       knots, this must be at least 4 data points(I think)
%   y_points are the y data points
%   end condition is either 'natural' or 'not-a-knot'
%       'natural' makes 2nd derivs of end points equal to zero
%       'not-a-knot' makes 3rd derivs of end points equal to zero
x_knots = x_points(2:end-1);

N = length(x_points);
m = 4;
k = length(x_knots); % k knots, the first point here tho is an endpoint
A = zeros(m+k);

if strcmp(end_condition,'natural')
    deriv_order=2;
elseif strcmp(end_condition,'not-a-knot')
    deriv_order=3;
end

for p=1:m+k
    A(1:N,p) = spl_1sided_basis_fxn(x_points,x_knots,p,m,0);
    A(N+1:N+2,p) = spl_1sided_basis_fxn([x_points(1); x_points(end)],x_knots,p,m,deriv_order); % % do end conditions
end

coeffs = A\[y_points,0,0]'; 

end

