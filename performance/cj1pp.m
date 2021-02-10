% cj1pp.m - power required and available
% 
% created on: 25-Sep-00
% updated on: 27-Sep-00
%

%% basic data for CP1, P15 lecture 08
cj1data;

%% power required at sea-level
V0 = 152.4;                                     % velocity for minimum power required
C_L = W / (.5 * rho_s * V0^2 * S);
C_D = C_D_0 + K*C_L^2;
D0 = (.5 * rho_s * V0^2 * S)*C_D;

% because V0 is the velocity for minimum power required, RR0=PR1
PR0 = D0*V0;
PR1 = W * sqrt(2*WS/rho_s) / (C_L^(3/2)/C_D);

%% power required, power available plots
% at sea level
for i = 1:400
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i);
    Fs(i) = F_s;
    PR_s(i) = Ds(i) * Vs(i);
    PA_s(i) = Fs(i) * Vs(i);
end

% altitude h=6.7km
h = 6705.6;
[T_h,p_h,rho_h] = stdatm(h);

for i = 1:400
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i);
    Fh(i) = (rho_h / rho_s) * Fs(i);
    PR_h(i) = Dh(i) * Vh(i);
    PA_h(i) = PA_s(i) * sqrt(rho_h / rho_s);
end

% minimum power
[PRmin_s,iMin_s] = min(PR_s);
VPRmin_s = Vs(iMin_s);
[PRmin_h,iMin_h] = min(PR_h);
VPRmin_h = Vh(iMin_h);

% minimum thrust
[TRmin_s,iMin_s] = min(Ds);
VTRmin_s = Vs(iMin_s);

% power required and power available at different altitude
figure(1)
plot(Vs,PR_s,'-',Vs,PA_s,'-',Vh,PR_h,'--',Vh,PA_h,'--',VPRmin_s,PRmin_s,'o',VPRmin_h,PRmin_h,'*')
grid on;
title('Power Required and Available');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s)');
legend('PR at sea level','PA at sea level',['PR at altitude ' num2str(h) ' m'],['PA at altitude ' num2str(h) ' m']);
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

% power required and thrust available at sea-level
figure(2) 
plot(Vs,PR_s/100,'-',Vs,Ds,'--',VPRmin_s,PRmin_s/100,'o',VTRmin_s,TRmin_s,'*');
grid on;
title('Power Required and Thrust Required at sea leval');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s) or thrust (N)');
legend('PR at sea level (scaled by 100)','TR at sea level');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
