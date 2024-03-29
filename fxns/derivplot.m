function foo = derivplot(value,deriv,var_wrt_to)
% plot a value and its deriv against a variable the deriv is wrt to, 
% e.g. derivplot(R,R_s,'s') to plot R, R_s vs s
% used for debugging
valname = inputname(1);
derivname = inputname(2);
wrtname = inputname(3);
sindex=2;
uindex=6;
vindex=1;
figure()
yyaxis left
if wrtname == 'u'
    plot(var_wrt_to,value(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to(sindex:end),value(sindex:end,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(value(sindex,uindex,:),size(var_wrt_to)))
end     
hold on

ylabel(valname)
yyaxis right
if wrtname == 'u'
    plot(var_wrt_to,deriv(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to(sindex:end),deriv(sindex:end,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(deriv(sindex,uindex,:),size(var_wrt_to)))
end     
hold on
yline(0,'--')
title(sprintf('%s and %s versus %s',valname,derivname,wrtname))
xlabel(wrtname)
ylabel(derivname)
legend(valname,derivname)


end