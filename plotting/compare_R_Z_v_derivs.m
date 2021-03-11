% compare R_v, Z_v to actual R, Z

%% numerical derivs for R_v, Z_v
if exist('R_v_n','var') == false   
R_v_n = real_space_deriv(R,v,deriv_method);
end
if exist('Z_v_n','var') == false   
Z_v_n = real_space_deriv(Z,v,deriv_method);
end



%% plot Raxis vs v
figure()

yyaxis left
plot(reshape(v,[190 1]),reshape(R(s_index,u_index,:),[190 1]),'b')
xlabel('v')
ylabel('R of axis (m)')

hold on
yyaxis right
plot(reshape(v,[190 1]),reshape(R_v(s_index,u_index,:),[190 1]))
hold on 
plot(reshape(v,[190 1]),reshape(R_v_n(s_index,u_index,:),[190 1]),'g--')
legend('R','R_v', 'numerical R_V')


%% plot Zaxis vs v
figure()
for s_index = 25:30
yyaxis left
plot(reshape(v,[190 1]),reshape(Z(s_index,u_index,:),[190 1]))
xlabel('v')
ylabel('Z of axis (m)')

hold on
yyaxis right
plot(reshape(v,[190 1]),reshape(Z_v(s_index,u_index,:),[190 1]))
hold on
plot(reshape(v,[190 1]),reshape(Z_v_n(s_index,u_index,:),[190 1]),'g--')
legend('Z','Z_v', 'numerical Z_v')
end