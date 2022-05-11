function evaluated_value = eval_series(suvgrid, coeffs, data, s_or_c)
% accepts a 3d array of (s,u,v) coordinates suvgrid, the fourier coeffs of
% the value to be evaluated, the read_vmec output data struct, and a string
% 's' or 'c' for if the coeffs are a sin or cos series
% the second dim of coeffs and the first dim of suvgrid should be the same
% (i.e. there should be one set of fourier coeffs for each s value in the
% suvgrid)
% vectorized for speed
% credit to matlabVMEC for implementation inspiration
xn = data.xn;
evaluated_value = zeros(size(suvgrid));
dimS = size(suvgrid,1);
dimU = size(suvgrid,2);
dimV = size(suvgrid,3);

u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi/data.nfp,dimV);
% [u,v] = ndgrid(u,v);

is_sin = s_or_c=='s';
is_cos = s_or_c=='c';

% mode angle arrays
mu = data.xm'*u;
nv = data.xn'*v;
% evaluated trig fxn arrays
cosmu = cos(mu);
sinmu = sin(mu);
cosnv = cos(nv);
sinnv = sin(nv);


for is=1:dimS
    coeffs_rep = repmat(coeffs(:,is),[1 dimU]); % replicate the fourier coeffs at this flux surface across the values of u
    sin_term = (coeffs_rep.*sinmu)'*cosnv + (coeffs_rep.*cosmu)'*sinnv;
    cos_term = (coeffs_rep.*cosmu)'*cosnv - (coeffs_rep.*sinmu)'*sinnv;
    evaluated_value(is,:,:)= is_sin.*sin_term + is_cos.*cos_term;

end

end