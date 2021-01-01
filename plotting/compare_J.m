if exist('JU_vmec','var') == false
JU_vmec = eval_series_nyq(suvgrid,data.currumnc,data,'c');
end
if exist('JV_vmec','var') == false
JV_vmec = eval_series_nyq(suvgrid,data.currvmnc,data,'c');
end

figure()

plot(data.phi(s_index:end),JU(s_index:end,u_index,v_nfp_index))
title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u')

hold on
plot(data.phi(s_index:end),JU_vmec(s_index:end,u_index,v_nfp_index))
title(sprintf('J^u VMEC vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u VMEC')
legend('My Calc','VMEC')

figure()
plot(data.phi(s_index:end),JU(s_index:end,u_index,v_nfp_index))
title(sprintf('J^u vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^u')
ylim([0,1e6])
hold on
plot(data.phi(s_index:end),JU_vmec(s_index:end,u_index,v_nfp_index))
legend('mine','VMEC')


figure()
plot(data.phi(s_index:end),JV(s_index:end,u_index,v_nfp_index))
title(sprintf('J^v vs s at u=%f, nfp*phi=%f',u(u_index),v(v_nfp_index)))
xlabel('s')
ylabel('J^v')
hold on
plot(data.phi(s_index:end),JV_vmec(s_index:end,u_index,v_nfp_index))
legend('mine','VMEC')


