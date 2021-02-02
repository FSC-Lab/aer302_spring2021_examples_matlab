% cj1tt.m - thrust  at sea level
% 
% created on: 24-Sep-00
% updated on:
%

% run data file first
cj1data;

% altitude
h = 6705.6;
[T_h,p_h,rho_h] = stdatm(h);

% from thrust available to calculate speed
c_s = [1 -(F_s/S)/C_d_0 K*WS^2/C_d_0];
q_s = roots(c_s);
V_s = sqrt(2* q_s / rho_s);
Vstall_s = sqrt(2*W / (rho_s * S * C_L_max));

F_h = rho_h / rho_s * F_s;
c_h = [1 -(F_h/S)/C_d_0 K*WS^2/C_d_0];
q_h = roots(c_h);
V_h = sqrt(2* q_h / rho_h);
Vstall_h = sqrt(2*W / (rho_h * S * C_L_max)); % assume C_L_max independent of the altitude



% thrut required and available
for i = 1:400
    % at sea level
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_d_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i); %TR
    Fs(i) = F_s; %TA
    % at altitude (varies V = \sqrt{sigma} Ve)
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_d_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i); %TR
    Fh(i) = (rho_h / rho_s) * Fs(i); %TA
end


% minimum Thrust
[TRmin_s,iMin_s] = min(Ds);
VTRmin_s = 10 + iMin_s *1;
[TRmin_h,iMin_h] = min(Dh);
VTRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h);

figure(1)
%plot(Vs,Ds,'-',Vs,Fs,'-',Vh,Dh,'--',Vh,Fh,'--',V,D,':',V,F,':')
plot(Vs,Ds,'-',Vs,Fs,'-',Vh,Dh,'--',Vh,Fh,'--',VTRmin_s,TRmin_s,'o', VTRmin_h,TRmin_h,'*')
title('Thrust Required and Available');
xlabel(' velocity (m/s)');
ylabel(' thrust (N)');
legend('T_R: sea-level','T_A: sea-level','T_R: altitude','T_A: altitude')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
% plotdlg;

% CL/CD progile, illustration @ sea level
LDs = C_Ls ./ C_Ds;
LDmax_s = max(LDs);
CLmax_s = sqrt(C_d_0 / K);
% [LDmax_s,iMax_s] = max(LDs);
% VLDmax_s = 10+iMax_s*1;

figure(2)
plot(C_Ls, LDs,CLmax_s,LDmax_s,'o')
axis([0 5 0 20])
grid
title('Lift to Drag Ratio')
xlabel(' C_L')
ylabel(' C_L / C_D')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
% plotdlg;

figure(3)
plot(C_L, C_D)
axis([0 5 0 1.5])
grid
title('Drag Polar')
xlabel(' C_L')
ylabel(' C_D')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);

  