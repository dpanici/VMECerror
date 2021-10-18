% find and plot rational surfaces in s
% close all % for now
s_dens = linspace(0,1,50000);
iota_dense = polyval(flip(data.ai),s_dens); % needs to be iotabar, not just iota
% figure
% plot(s,polyval(flip(data.ai),s))
% hold on
% plot(s,data.iotaf,'--','DisplayName','iota')
% hold on
rationals = [];
rat_s = [];
ind = 1;
for m=1:max(data.xm_nyq)
    m
    for n=1:max(data.xn_nyq/data.nfp)
        i_rational = n/m;
        [val,s_rat_ind] = min(abs(iota_dense-i_rational));
%         xline(s_dens(s_rat_ind))
%         hold on
        if  and((s_dens(s_rat_ind) > 0.1), (s_dens(s_rat_ind) < 0.99))
        rationals(ind)=i_rational;
        rat_s(ind) = s_dens(s_rat_ind);
        ind = ind + 1;
        end
    end
end