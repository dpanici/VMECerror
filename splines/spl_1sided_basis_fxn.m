function [y] = spl_1sided_basis_fxn(x,x_knots,p,spl_order,deriv_order)
%SPL_1SIDED_BASIS returns value of the pth one sided basis function
% where p = 1... m+k
% p=1...m are 1,x-a,...(x-a).^(m-1)
% p = m+1...m+k are (x-a)_+ .^(m-1)
% for the one sided spline basis fxns of order m (order being degree + 1)
%  for the ith knot with jth exponent (jth really meaning m-j but j=1 for
%  simple knots)
% x is vector (k+2 points) at which one sided basis is defined
% x_knots is the k+1 knots (the first k+1 points in interval, excludes last
% points basically
% i is the index of knot being evaluated with (basis fxn is zero for x <
% x_knots(i)
% j is the degree of the polynomial (will be always = m-1 for simple knots)
% deriv order is which derivative to return, 0 for value, 1 for first, etc

if p > spl_order
y=spl_1sided_basis_pluses(x,x_knots,p-spl_order,spl_order-1,deriv_order);
else
    if deriv_order == 0
        y = (x-x(1)).^(p-1);
    
    elseif deriv_order ==1
        is_linear_at_least = (p - 1) > 0;
        y = (p-1).*(x-x(1)).^(p-2);
        y = y .*is_linear_at_least; % ensure only linear term or above is 
        y(isnan(y))=0;
        
    elseif deriv_order ==2
        is_quadratic_at_least = (p - 2) > 0;
        y = (p-1)*(p-2).*(x-x(1)).^(p-3);
        y = y .*is_quadratic_at_least; % ensure only quadratic term or above is 
        y(isnan(y))=0;
        
    elseif deriv_order ==3
        is_cubic_at_least = (p - 3) > 0;
        y = (p-1)*(p-2)*(p-3).*(x-x(1)).^(p-4);
        y = y .*is_cubic_at_least; % ensure only cubic term or above is 
        y(isnan(y))=0;
    end
end

