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
% magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4) + (BU.*BV).*dot(eu,ev,4) + (BU.*BV).*dot(ev,eu,4));
magB_sq = magB.^2;

abs_g = abs(g);


% W_B = trapz(v,trapz(u,trapz(s,magB.*g))) / 2/mu0
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


% %Use LHS
% VMEC_WB_LHS=0;
% if dimV > 1
%     NV = dimV
% else
%     NV = 2;
% end
% for is=1:dimS-1
%     ds = s(is+1)-s(is);
%     for iu = 1:dimU-1
%         du = u(iu+1) - u(iu);
%         for iv = 1:NV-1
%             if dimV > 1
%                 dv = v(iv+1)-v(iv);
%             else
%                 dv = 2*pi;
%             end
%             VMEC_WB_LHS = VMEC_WB_LHS + magB_vmec_sq(is,iu,iv) * abs(abs_g_vmec(is,iu,iv))*ds*du*dv;
%         end
%     end
% end
% VMEC_W_B_LHS = VMEC_WB_LHS/2/mu0
% 
% %Use RHS
% if dimV > 1
%     NV = dimV;
% else
%     NV = 2;
% end
% VMEC_WB_RHS=0;
% for is=2:dimS
%     ds = s(is)-s(is-1);
%     for iu = 2:dimU
%         du = u(iu) - u(iu-1);
%         for iv = 2:NV
%             if dimV > 1
%                 dv = v(iv)-v(iv-1);
%             else
%                 dv = 2*pi;
%                 iv = 1;
%             end
%             VMEC_WB_RHS = VMEC_WB_RHS + magB_vmec_sq(is,iu,iv) * abs(abs_g_vmec(is,iu,iv))*ds*du*dv;
%         end
%     end
% end
% VMEC_W_B_RHS = VMEC_WB_RHS/2/mu0
if dimV > 1
    W_p = trapz(v,trapz(u,trapz(s,g.*presf./(gamma-1))));
else
    W_p = trapz(u,trapz(s,g.*presf./(gamma-1)));
end
W = W_B + W_p;
