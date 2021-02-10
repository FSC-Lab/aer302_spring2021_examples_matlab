% cp1pp.m - power required and available, including climbing performance
% 
% created on: 10.Feb.2020
% updated on: 
%

%% basic data for CP1, P16 lecture 08
cp1data;

%% power required at sea-level P14 lecture 08
V0 = 60.96;                                     % velocity for minimum power required
C_L = W / (.5 * rho_s * V0^2 * S);
C_D = C_D_0 + K*C_L^2;
D0 = (.5 * rho_s * V0^2 * S)*C_D;

% because V0 is the velocity for minimum power required, RR0=PR1
PR0 = D0*V0;                                                    
PR1 = W * sqrt(2*WS/rho_s) / (C_L^(3/2)/C_D); 

%% power required and power available plots
% at sea level
for i = 1:200
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i);
    PR_s(i) = Ds(i) * Vs(i);
    PA_s(i) = PAs;
    RC_s(i) = (PA_s(i) - PR_s(i))/W;                        % climb rate
    gamma_s(i) = asin(RC_s(i) / Vs(i)) * 180/pi;      % flight path angle
end

% at h=6km
h = 6000;
[T_h,p_h,rho_h] = stdatm(h);
PAh = PAs * sqrt(rho_h / rho_s);
for i = 1:200
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);                    % Vh(i) = Vs(i);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i);
    PR_h(i) = Dh(i) * Vh(i);
    PA_h(i) = PAh;
    RC_h(i) = (PA_h(i) - PR_h(i))/W;
    gamma_h(i) = asin(RC_h(i) / Vs(i)) * 180/pi;
end

% minimum power from graphics
[PRmin_s,iMin_s] = min(PR_s);
VPRmin_s = 10+iMin_s*1;
[PRmin_h,iMin_h] = min(PR_h);
VPRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h);

% minimum power from calculation (analytical)
VPRmin_s_calc = sqrt( 2*WS / rho_s / sqrt(3*C_D_0/K) );
VPRmin_h_calc = sqrt( 2*WS / rho_h / sqrt(3*C_D_0/K) );

% plot
figure(1)
plot(Vs,PR_s,'-',Vs,PA_s,'-',Vh,PR_h,'--',Vh,PA_h,'--',VPRmin_s,PRmin_s,'o',VPRmin_h,PRmin_h,'*')
grid on;
axis([0 250 0 2e06])
title('Power Required and Available for CP1');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s)');
legend('PR at sea level','PA at sea level',['PR at altitude ' num2str(h) ' m'],['PA at altitude ' num2str(h) ' m']);
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

%% hodograph (approximation)
Vv_s = RC_s;
Vh_s = Vs .* cos(gamma_s*pi/180);
Vv_h = RC_h;
Vh_h = Vh .* cos(gamma_h*pi/180);

figure(3)
plot(Vh_s,Vv_s,'-',Vh_h,Vv_h,'--');
grid on;
axis([0 100 0 10]); 
title(' Hodograph for [CP1]');
legend('sea-level', '6000m')
xlabel(' Vh (m/s)');
ylabel('Vv (m/s)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

%% other performance computations
% maximum rate of climb
VRC_max_s_calc = VPRmin_s_calc;          % sea level 
VRC_max_h_calc = VPRmin_h_calc;         % h=6km
C_D_rcmax = 4*C_D_0;                            % C_D_0 = 1/3 K CL^2;
RCmax_s = (PAs - PRmin_s)/W;                                    % sea level
RCmax_h = (PAs * sqrt(rho_h / rho_s) - PRmin_h)/W;    % h=6km

% maximum angle of climb: 1) solve for the velocity for maximum gamma
poly_s = [1/2 *rho_s^2* C_D_0/WS 0 0 1/2*rho_s*PAs/W -2*K*WS];
poly_h = [1/2 *rho_h^2* C_D_0/WS 0 0 1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W -2*K*WS];
Vgamma_s = roots(poly_s);
Vgamma_h = roots(poly_h);
Vgammamax_s_est = 2*K*WS / (1/2*rho_s*PAs/W);
Vgammamax_h_est = 2*K*WS / (1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W);

% maximum angle of climb: 2) compute gamma from the velocity (sea level)
q_gamma_s = 1/2*rho_s*Vgammamax_s_est^2;
CL_gamma_s = W / (q_gamma_s *S);
CD_gamma_s = C_D_0 + K*CL_gamma_s^2;
D_gamma_s = q_gamma_s*S*CD_gamma_s;
RC_gamma_s = (PAs - D_gamma_s*Vgammamax_s_est)/W;
gammamax_s = asin(RC_gamma_s / Vgammamax_s_est) * 180/pi;

% maximum angle of climb: 2) compute gamma from the velocity (h=6km)
q_gamma_h = 1/2*rho_s*Vgammamax_h_est^2;
CL_gamma_h = W / (q_gamma_h *S);
CD_gamma_h = C_D_0 + K*CL_gamma_h^2;
D_gamma_h = q_gamma_h*S*CD_gamma_h;
RC_gamma_h = (PAs * sqrt(rho_h / rho_s) - D_gamma_h*Vgammamax_h_est)/W;
gammamax_h = asin(RC_gamma_h / Vgammamax_h_est) * 180/pi;

%% power available to sustain given speeds
P = PAs * sqrt(rho_h / rho_s);
v = 0.1:1:100;
p1 = [.5*rho_h*S*C_D_0 0 0 -P K*S*WS^2/(.5*rho_h)];
f1 = polyval(p1,v);

figure(2)
plot(v,f1); grid on;
title(['speed profile under P_A =', num2str(P)])
xlabel('velocity (m/s)');
ylabel('P_A (N m/s)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

  