function foo = derivplot(value,deriv,var_wrt_to)
% plot a value and its deriv against a variable the deriv is wrt to, 
valname = inputname(1);
derivname = inputname(2);
wrtname = inputname(3);
sindex=10;
uindex=6;
vindex=1;
figure()
yyaxis left
if wrtname == 'u'
    plot(var_wrt_to,value(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to,value(:,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(value(sindex,uindex,:),size(var_wrt_to)))
end     
hold on

ylabel(valname)
yyaxis right
if wrtname == 'u'
    plot(var_wrt_to,deriv(sindex,:,vindex))
elseif wrtname == 's'
    plot(var_wrt_to,deriv(:,uindex,vindex))
elseif wrtname == 'v'
    plot(var_wrt_to,reshape(deriv(sindex,uindex,:),size(var_wrt_to)))
end     
hold on
yline(0,'--')
title(sprintf('%s and %s versus %s',valname,derivname,wrtname))
xlabel(wrtname)
ylabel(derivname)
legend(valname,derivname)
% set(gcf, 'Position',  [200, 200, 900, 700])

% fin dif to approx the deriv, use ismembertol to check that the fin dif is
% similar to the calculated deriv

end