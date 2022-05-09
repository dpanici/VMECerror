function [y] = spl_1sided_basis_pluses(x,x_knots,i,j,deriv_order)
%SPL_1SIDED_BASIS returns value of the one sided basis function
% for the one sided spline basis of order m (order being degree + 1)
%  for the ith knot with jth exponent (jth really meaning m-j but j=1 for
%  simple knots)
% x is vector (k+2 points) at which one sided basis is defined
% x_knots is the k+1 knots (the first k+1 points in interval, excludes last
% points basically
% i is the index of knot being evaluated with (basis fxn is zero for x <
% x_knots(i)
% j is the degree of the polynomial (will be always = m-1 for simple knots)
% deriv order is which derivative to return, 0 for value, 1 for first, etc

bool_x_greater_than_x_i = x>=x_knots(i); % multiply final y by this to enforce being defined
% only for x > x_i
if deriv_order == 0
y = (x-x_knots(i)).^j;

elseif deriv_order==1
y = (j).*(x-x_knots(i)).^(j-1);

elseif deriv_order==2
y = (j)*(j-1).*(x-x_knots(i)).^(j-2);

elseif deriv_order==3
y = (j)*(j-1)*(j-2).*(x-x_knots(i)).^(j-3);
end
y = y .* bool_x_greater_than_x_i;

end

