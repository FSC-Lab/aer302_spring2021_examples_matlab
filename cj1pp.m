% cj1pp.m - power required and available
% 
% created on: 25-Sep-00
% updated on: 27-Sep-00
%

% run data file first
cj1data;

% sea-level [example ch4]
V0 = 152.4;
C_L = W / (.5 * rho_s * V0^2 * S);
C_D = C_D_0 + K*C_L^2;
D0 = (.5 * rho_s * V0^2 * S)*C_D;
PR0 = D0*V0;
PR1 = W * sqrt(2*WS/rho_s) / (C_L^(3/2)/C_D); % for verification == PR0


% altitude
h = 6705.6;
[T_h,p_h,rho_h] = stdatm(h);

% power required
for i = 1:400
    % at sea level
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i);
    Fs(i) = F_s;
    PR_s(i) = Ds(i) * Vs(i);
    PA_s(i) = Fs(i) * Vs(i);
    % at altitude (varies V = \sqrt{sigma} Ve)
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i);
    Fh(i) = (rho_h / rho_s) * Fs(i);
    PR_h(i) = Dh(i) * Vh(i);
%     PR_h(i) = PR_s(i) * sqrt(rho_s / rho_h);
    PA_h(i) = PA_s(i) * sqrt(rho_h / rho_s);
%     % at altitude (varies V = Ve)
%     V(i) = Vs(i);
%     qS = 0.5 * rho_h * V(i)^2 * S;
%     C_L = W / qS;
%     C_D = C_d_0 + K * C_L^2;
%     D(i) = qS * C_D;
%     F(i) = Fh(i);
%     PR(i) = D(i) * V(i);
%     PA(i) = F(i) * V(i);
end

% minimum power
[PRmin_s,iMin_s] = min(PR_s);
VPRmin_s = 10+iMin_s*1;
[PRmin_h,iMin_h] = min(PR_h);
VPRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h);

% minimum thrust (repeat from cj1tt)
[TRmin_s,iMin_s] = min(Ds);
VTRmin_s = 10 + iMin_s *1;



figure(1) %PR-PA Plot with different altitudes
plot(Vs,PR_s,'-',Vs,PA_s,'-',Vh,PR_h,'--',Vh,PA_h,'--',VPRmin_s,PRmin_s,'o',VPRmin_h,PRmin_h,'*')
% plot(V,Pr_s,'-',V,Pa_s,'-',Vh,Pr_h,'--',Vh,Pa_h,'--')
%plot(V,Pr_s,V,Pa_s,Vh,Pr_h,Vh,Pa_h)
title('Power Required and Available');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s)');
legend('PR at sea level','PA at sea level',['PR at altitude ' num2str(h) ' m'],['PA at altitude ' num2str(h) ' m']);
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;

figure(2) % PR and TR plot together at sea-level
plot(Vs,PR_s/100,'-',Vs,Ds,'--',VPRmin_s,PRmin_s/100,'o',VTRmin_s,TRmin_s,'*')
% plot(V,Pr_s,'-',V,Pa_s,'-',Vh,Pr_h,'--',Vh,Pa_h,'--')
%plot(V,Pr_s,V,Pa_s,Vh,Pr_h,Vh,Pa_h)
title('Power Required and Thrust Required at sea leval');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s) or thrust (N)');
legend('PR at sea level (scaled by 100)','TR at sea level');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;
