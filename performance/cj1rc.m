% cj1rc.m - rate of climbing
% 
% created on: 03-Oct-00
% updated on:
%
% Revision History


% run data file first
cj1data;

% altitude
h = 6705.6;
[T_h,p_h,rho_h] = stdatm(h);
F_h = (rho_h / rho_s) * F_s; %TA
    

% RC at sea level - assume gamma (flight path angle) << 1
for i = 1:400
% at sea level
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i); %TR
    Fs(i) = F_s; %TA
    RC_s(i) = (Fs(i) * Vs(i) - Ds(i) * Vs(i)) / W;
    gamma_s(i) = asin(RC_s(i) / Vs(i)) * 180/pi;
    % at altitude 
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i); %TR
    
    Fh(i) = F_h; %TA
    RC_h(i) = (Fh(i) * Vh(i) - Dh(i) * Vh(i)) / W;
    gamma_h(i) = asin(RC_h(i) / Vh(i)) * 180/pi;
    
    % compared to power curve
    PR_s(i) = Ds(i) * Vs(i);
    PA_s(i) = F_s * Vs(i);

    % hodograph (approximation)
    Vv_s(i) = RC_s(i);
    % Vh_s(i) = sqrt(V(i)^2 - Vv_s(i)^2);
    Vh_s(i) = Vs(i) * cos(gamma_s(i)*pi/180);
    Vv_h(i) = RC_h(i);
    Vh_h(i) = Vh(i) * cos(gamma_h(i)*pi/180);
 end

% compare maximum values
[RCmax_s,iRCmax_s] = max(RC_s);
Vrcmax_s = Vs(iRCmax_s);
[gammamax_s,igammamax_s] = max(real(gamma_s));
Vgammamax_s = Vs(igammamax_s);

% analytical calculation of maximum values
TWs = F_s / W;
Vrcmax_s_calc = sqrt( (WS*TWs)/(3*rho_s*C_D_0) *(1+sqrt(1+12*K*C_D_0 / TWs^2)) );
% alternative calculation use V* as base
Vstar = sqrt(2*WS/(rho_s *sqrt(C_D_0/K)));
z = sqrt(TWs^2 / (12*K*C_D_0) +1);
Vrcmax_s_calc2 = Vstar * sqrt(1/(2*sqrt(3)))*(sqrt(z+1) + sqrt(z-1));
Vgammamax_s_calc = sqrt(2*WS/(rho_s *sqrt(C_D_0/K)));


%Eqn = [-rho_s^2 * S^2 * C_d_0, rho_s * S * F_s, 0, 2*K*W^2];
%Vcal = roots(Eqn);

figure(1)
subplot(211)
plot(Vs,PA_s,'-',Vs,PR_s,'--');
title(' Power available and required at Sea-Level');
xlabel(' velocity (m/s)');
ylabel('Power (N.m/s)');
legend('PA','PR');
subplot(212)
plot(Vs,RC_s,'-',Vrcmax_s,RCmax_s,'o')
axis([0 400 -100 50]);
title('Rate of Climb at Sea-Level');
xlabel(' velocity (m/s)');
ylabel('R/C (m/s)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;

figure(2)
subplot(211)
plot(Vs,PA_s,'-',Vs,PR_s,'--');
title(' Power available and required at Sea-Level');
xlabel(' velocity (m/s)');
ylabel('Power (N.m/s)');
legend('PA','PR');
subplot(212)
plot(Vs,gamma_s,'-',Vgammamax_s,gammamax_s,'o')
axis([0 400 -100 50]);
title('Climbing Angle at Sea-Level');
xlabel(' velocity (m/s)');
ylabel('\gamma (degree)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;

figure(3)
plot(Vh_s,Vv_s,'-',Vh_h,Vv_h,'--');
title(' Hodograph for [CJ1]');
legend(['T_A =',num2str(F_s)], ['T_A =',num2str(F_h)])
axis([0 400 0 50]); 
xlabel(' Vh (m/s)');
ylabel('Vv (m/s)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
plotdlg;  