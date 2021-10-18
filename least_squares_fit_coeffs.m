function [s_deriv,s2_deriv] = least_squares_fit_coeffs(coeffs,data,~)
%LEAST_SQUARES_FIT_COEFFS Summary of this function goes here
%   takes in fourier coeff matrix from VMEC coeffs
% data from vmec output
% window_size, the number of data points at a time you want to fit with poly of order k
%   returns the radial deriv of the coeffs, same size as input coeffs
% if return_2nd_deriv is true, it also returns the 2nd radial deriv of the
% coeffs
global POLY_LSQ_WINDOW_SIZE
global POLY_LSQ_ORDER
window_size=POLY_LSQ_WINDOW_SIZE;
k = POLY_LSQ_ORDER;
s_deriv = zeros(size(coeffs));
s2_deriv = zeros(size(coeffs));

m = data.ns / window_size;
N = data.ns;

x = linspace(0,1,N);
for mode_ind=1:length(data.xm)
    y = coeffs(mode_ind,:);
    cs = zeros([m,k]);
    for i=1:m
    
    inds = (1+(i-1)*N/m:i*N/m);
    cs(i,:)= least_squares_fit(x(inds), y(inds), k)'; 
    s_deriv_c = polyder(cs(i,:));
    s_deriv(mode_ind,inds)=polyval(s_deriv_c,x(inds));

    s2_deriv_c = polyder(s_deriv_c);
    s2_deriv(mode_ind,inds)=polyval(s2_deriv_c,x(inds));
    
    end


end
end
