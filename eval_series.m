function evaluated_value = eval_series(suvgrid, coeffs, data, s_or_c)
% accepts a 3d array of (s,u,v) coordinates suvgrid, the fourier coeffs of
% the value to be evaluated, the read_vmec output data struct, and a string
% 's' or 'c' for if the coeffs are a sin or cos series
% the second dim of coeffs and the first dim of suvgrid should be the same
% (i.e. there should be one set of fourier coeffs for each s value in the
% suvgrid)

evaluated_value = zeros(size(suvgrid));
dimS = size(suvgrid,1);
dimU = size(suvgrid,2);
dimV = size(suvgrid,3);

u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

is_sin = s_or_c=='s';
is_cos = s_or_c=='c';

for is = 1:dimS% is will be the second dim of the fourier series coeffs called
    for iu = 1:dimU
        for iv = 1:dimV
            for i=1:length(data.xm)
                sin_term = sin(data.xm(i)*u(iu) - data.nfp*data.xn(i).*v(iv));
                cos_term = cos(data.xm(i)*u(iu) - data.nfp*data.xn(i).*v(iv));
                evaluated_value(is,iu,iv) = evaluated_value(is,iu,iv) ...
                    + coeffs(i,is) * (is_sin*sin_term + is_cos*cos_term);
            end
        end
    end
end

end