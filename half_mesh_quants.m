% get R,Z derivs on half mesh for use in finding g near axis
ns = data.ns;
rumns_h = f2h(rumns,ns);
zumnc_h = f2h(zumnc,ns);
rsmnc_h = f2h(rsmnc,ns);
zsmns_h = f2h(zsmns,ns);
rvmns_h = f2h(rvmns,ns);
zvmnc_h = f2h(zvmnc,ns);
rmnc_h = f2h(data.rmnc,ns);
gmnc_h = rmnc_h.*rsmnc_h.*zumnc_h + rmnc_h.*rumns_h.*zsmns_h;

% just need to make suvgrid that is a half grid on s, to pass to
% eval_series

s_h = linspace(0.5*s(2),1,dimS-1);
u = linspace(0,2*pi,dimU);
v = linspace(0,2*pi,dimV);

suvgrid_h = ndgrid(s_h,u,v);

R_h = eval_series(suvgrid_h,rmnc_h,data,'c');
R_u_h = eval_series(suvgrid_h,rumns_h,data,'s');
Z_u_h = eval_series(suvgrid_h,zumnc_h,data,'c');
R_v_h = eval_series(suvgrid_h,rvmns_h,data,'s');
Z_v_h = eval_series(suvgrid_h,zvmnc_h,data,'c');
R_s_h = eval_series(suvgrid_h,rsmnc_h,data,'c');
Z_s_h = eval_series(suvgrid_h,zsmns_h,data,'s');

es_h = cat(4,R_s_h,zeros(size(R_s_h)),Z_s_h);
eu_h = cat(4,R_u_h, zeros(size(R_s_h)), Z_u_h);
ev_h = cat(4,R_v_h, R_h, Z_v_h);

% g_h = R_h .* (R_s_h .* Z_u_h + R_u_h .* Z_s_h);
g_h = dot(es_h,cross(eu_h,ev_h,4),4);
g_for_h = eval_series(suvgrid_h,gmnc_h,data,'c');