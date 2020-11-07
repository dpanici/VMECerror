% Run this script after running force_error.m to get the MHD energy
%% Energy

BU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
BV_vmec = eval_series_nyq(suvgrid,data.bsupvmnc,data,'c');
% JU_vmec = eval_series_nyq(suvgrid,data.bsupumnc,data,'c');
magB_vmec = (BU_vmec.^2).*dot(eu,eu,4) + (BV_vmec.^2).*dot(ev,ev,4);
g_vmec = eval_series_nyq(suvgrid,data.gmnc,data,'c');

% W = W_magnetic + W_pressure
%W_magnetic = |B|^2 /2/mu0 , W_pressure = p / (gamma-1)
gamma=data.gamma;
presf = repmat(data.presf',1,dimU,dimV);
magB = (BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4);

magB_axis = mean(magB(4,:,:),'all');
normed_magB = magB./magB_axis;

[s1,u1,v1] = meshgrid(s,u,v);

% W_B = trapz(v,trapz(u,trapz(s,magB.*g))) / 2/mu0 
W_B = trapz(u,trapz(s,magB.*g)) / 2/mu0 *2*pi;
VMEC_W_B = trapz(u,trapz(s,magB_vmec.*g)) / 2/mu0 *2*pi;
% W_B_norm = trapz(u,trapz(s,normed_magB.*g)) / 2/mu0 *2*pi;
% W_B = mean(W_B_norm,'all') 

%Use LHS
VMEC_WB_LHS=0
for is=1:dimS-1
    for iu = 1:dimU-1
        for iv = 1:dimV
            ds = s(is+1)-s(is);
            du = u(iu+1) - u(iu);
            %dv = 2*pi;
            VMEC_WB_LHS = VMEC_WB_LHS + magB_vmec(is,iu,iv) * g_vmec(is,iu,iv)*ds*du;
        end
    end
end
VMEC_W_B_LHS = VMEC_WB_LHS/2/mu0 * 2 * pi

%Use RHS
VMEC_WB_RHS=0
for is=2:dimS
    for iu = 2:dimU
        for iv = 1:dimV
            ds = s(is)-s(is-1);
            du = u(iu) - u(iu-1);
            %dv = 2*pi;
            VMEC_WB_RHS = VMEC_WB_RHS + magB_vmec(is,iu,iv) * g_vmec(is,iu,iv)*ds*du;
        end
    end
end
VMEC_W_B_RHS = VMEC_WB_RHS/2/mu0* 2 * pi

W_p = trapz(u,trapz(s,g.*presf./(gamma-1)));

W = W_B + W_p;
