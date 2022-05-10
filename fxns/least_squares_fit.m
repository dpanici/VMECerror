function c = least_squares_fit(x,y,k)
%LEAST_SQUARES_FIT Summary of this function goes here
%   takes in points x that you want to fit with poly of order k
%   to fit points y
%   returns c coefficients of length k
   %in descending power (so how polyval takes them)
N=length(x);
% choose to fit with a k order polynomial, 
% k coefficients to find
% now, form a vandermonde matrix for our fitting function
% we want Ac = y, where
% y is [N x 1] is our data
% c is [k x 1] are our coefficients (in ascending order let's say)
% A is [N x k] is our vandermonde matrix
% The ith row is the polynomial evaluated at the ith x value
% the jth column is the jth part of our basis polynomial
% i.e. A_11 is the first part (constant 1) evaluated at the first x (0) = 1
% polynomial is :
%f(x) = c(1) * 1 + c(2) * x + c(3)* x^2 + c(4) * x^3 + ... + c(k) * x^(k-1)


A = zeros([N,k]);
for j=1:k
        A(:,j) = x'.^(j-1);
end

c = A\y';

c = flip(c);

end

