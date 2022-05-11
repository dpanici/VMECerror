% Run this script after running force_error.m to get the MHD energy
%% Energy

BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
magB_vmec = eval_series_nyq(suvgrid,data.bmnc,data,'c');
magB_vmec_sq = magB_vmec.^2;
abs_g_vmec = abs(eval_series_nyq(suvgrid,data.gmnc,data,'c'));

% W = W_magnetic + W_pressure
% W_magnetic = |B|^2 /2/mu0 , W_pressure = p / (gamma-1)
gamma=data.gamma;
presf = repmat(data.presf',1,dimU,dimV);
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4) + (BU.*BV).*dot(eu,ev,4) + (BU.*BV).*dot(ev,eu,4));
magB_sq = magB.^2;

abs_g = abs(g);


if dimV >1
    W_B = trapz(v,trapz(u,trapz(s,magB_sq.*abs_g))) / 2/mu0;
    VMEC_W_B = trapz(v,trapz(u,trapz(s,magB_vmec_sq.*abs_g_vmec))) / 2/mu0;
else
    W_B = trapz(u,trapz(s,magB_sq.*abs_g)) / 2/mu0 *2*pi;
    VMEC_W_B = trapz(u,trapz(s,magB_vmec_sq.*abs_g_vmec)) / 2/mu0 *2*pi;
end
 % bc we only have v going around 1 field period
W_B = W_B * data.nfp
VMEC_W_B = VMEC_W_B * data.nfp


if dimV > 1
    W_p = data.nfp*trapz(v,trapz(u,trapz(s,abs(g).*presf./(gamma-1))));
else
    W_p = trapz(u,trapz(s,abs(g).*presf./(gamma-1)));
end
W = W_B + W_p;
