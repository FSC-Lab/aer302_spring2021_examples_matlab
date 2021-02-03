% cj1tt.m - thrust  at sea level
% 
% created on: 24-Sep-00
% updated on:
%

%% basic data for CP1, P45 lecture 08
cj1data;

% some data at altitude h=6.7km
h = 6705.6;
[T_h,p_h,rho_h] = stdatm(h);

%% compute speed for given thrust available (P47 lecture 08)
% at sea level
c_s = [1 -(F_s/S)/C_D_0 K*WS^2/C_D_0];
q_s = roots(c_s);
V_s = sqrt(2* q_s / rho_s);
Vstall_s = sqrt(2*W / (rho_s * S * C_L_max));

% at altitude h=6.7km
F_h = rho_h / rho_s * F_s;              % thrust available at h is different from sea  level
c_h = [1 -(F_h/S)/C_D_0 K*WS^2/C_D_0];
q_h = roots(c_h);
V_h = sqrt(2* q_h / rho_h);
Vstall_h = sqrt(2*W / (rho_h * S * C_L_max));  % assume C_L_max independent of the altitude

%% thrut required and thrust available plots
for i = 1:400
    % at sea level
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i);              % thrust required
    Fs(i) = F_s;                                % thrust available 
    
    % at altitude h=6.7km
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i);                         % thrust required
    Fh(i) = (rho_h / rho_s) * Fs(i);             % thrust available 
end

% minimum Thrust
[TRmin_s,iMin_s] = min(Ds);
VTRmin_s = 10 + iMin_s *1;
[TRmin_h,iMin_h] = min(Dh);
VTRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h);

% plot
figure(1)
plot(Vs,Ds,'-',Vs,Fs,'-',Vh,Dh,'--',Vh,Fh,'--',VTRmin_s,TRmin_s,'o', VTRmin_h,TRmin_h,'*');
grid on;
title('Thrust Required and Available');
xlabel(' velocity (m/s)');
ylabel(' thrust (N)');
legend('T_R: sea-level','T_A: sea-level','T_R: altitude','T_A: altitude')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);

%% CL/CD progile, illustration @ sea level
LDs = C_Ls ./ C_Ds;
LDmax_s = max(LDs);
CLmax_s = sqrt(C_D_0 / K);

% lift to drag ratio vs. drag
figure(2)
plot(C_Ls, LDs,CLmax_s,LDmax_s,'o')
axis([0 5 0 20])
grid;
title('Lift to Drag Ratio')
xlabel(' C_L')
ylabel(' C_L / C_D')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);

% drag pola, lift coefficient vs. drag coefficient
figure(3)
plot(C_L, C_D)
axis([0 5 0 1.5])
grid on;
title('Drag Polar')
xlabel(' C_L')
ylabel(' C_D')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);

  