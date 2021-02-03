% cp1pp.m - power required and available, including climbing performance
% 
% created on: 10.Feb.2020
% updated on: 
%

% run data file first
cp1data;

% sea-level [example ch4]
V0 = 60.96;
C_L = W / (.5 * rho_s * V0^2 * S);
C_D = C_D_0 + K*C_L^2;
D0 = (.5 * rho_s * V0^2 * S)*C_D;
PR0 = D0*V0;
PR1 = W * sqrt(2*WS/rho_s) / (C_L^(3/2)/C_D); % for verification == PR0

% altitude
h = 6000;
[T_h,p_h,rho_h] = stdatm(h);
PAh = PAs * sqrt(rho_h / rho_s);

% power required
for i = 1:200
    % at sea level
    Vs(i) = 10 + 1 * i;
    qsS(i) = 0.5 * rho_s * Vs(i)^2 * S;
    C_Ls(i) = W / qsS(i);
    C_Ds(i) = C_D_0 + K * C_Ls(i)^2;
    Ds(i) = qsS(i) * C_Ds(i);
    PR_s(i) = Ds(i) * Vs(i);
    PA_s(i) = PAs;
    RC_s(i) = (PA_s(i) - PR_s(i))/W;
    gamma_s(i) = asin(RC_s(i) / Vs(i)) * 180/pi;
    % at altitude (varies V = \sqrt{sigma} Ve)
    Vh(i) = Vs(i) * sqrt(rho_s / rho_h);
    qS(i) = 0.5 * rho_h * Vh(i)^2 * S;
    C_L(i) = W / qS(i);
    C_D(i) = C_D_0 + K * C_L(i)^2;
    Dh(i) = qS(i) * C_D(i);
    PR_h(i) = Dh(i) * Vh(i);
%     PR_h(i) = PR_s(i) * sqrt(rho_s / rho_h);
    PA_h(i) = PAh;
    RC_h(i) = (PA_h(i) - PR_h(i))/W;
    gamma_h(i) = asin(RC_h(i) / Vs(i)) * 180/pi;
    
    % hodograph (approximation)
    Vv_s(i) = RC_s(i);
    % Vh_s(i) = sqrt(V(i)^2 - Vv_s(i)^2);
    Vh_s(i) = Vs(i) * cos(gamma_s(i)*pi/180);
    Vv_h(i) = RC_h(i);
    Vh_h(i) = Vh(i) * cos(gamma_h(i)*pi/180);
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

% minimum power from graphics
[PRmin_s,iMin_s] = min(PR_s);
VPRmin_s = 10+iMin_s*1;
[PRmin_h,iMin_h] = min(PR_h);
VPRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h);

% minimum power from calculation (analytical)
VPRmin_s_calc = sqrt( 2*WS / rho_s / sqrt(3*C_D_0/K) );
VPRmin_h_calc = sqrt( 2*WS / rho_h / sqrt(3*C_D_0/K) );

% maximum rate of climb
VRC_max_s_calc = VPRmin_s_calc;
VRC_max_h_calc = VPRmin_h_calc;
C_D_rcmax = 4*C_D_0; % C_D_0 = 1/3 K CL^2;
RCmax_s = (PAs - PRmin_s)/W;
RCmax_h = (PAs * sqrt(rho_h / rho_s) - PRmin_h)/W;


% maximum angle of climb
poly_s = [1/2 *rho_s^2* C_D_0/WS 0 0 1/2*rho_s*PAs/W -2*K*WS];
poly_h = [1/2 *rho_h^2* C_D_0/WS 0 0 1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W -2*K*WS];
Vgamma_s = roots(poly_s);
Vgamma_h = roots(poly_h);
Vgammamax_s_est = 2*K*WS / (1/2*rho_s*PAs/W);
Vgammamax_h_est = 2*K*WS / (1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W);

q_gamma_s = 1/2*rho_s*Vgammamax_s_est^2;
CL_gamma_s = W / (q_gamma_s *S);
CD_gamma_s = C_D_0 + K*CL_gamma_s^2;
D_gamma_s = q_gamma_s*S*CD_gamma_s;
RC_gamma_s = (PAs - D_gamma_s*Vgammamax_s_est)/W;
gammamax_s = asin(RC_gamma_s / Vgammamax_s_est) * 180/pi;;

q_gamma_h = 1/2*rho_s*Vgammamax_h_est^2;
CL_gamma_h = W / (q_gamma_h *S);
CD_gamma_h = C_D_0 + K*CL_gamma_h^2;
D_gamma_h = q_gamma_h*S*CD_gamma_h;
RC_gamma_h = (PAs * sqrt(rho_h / rho_s) - D_gamma_h*Vgammamax_h_est)/W;
gammamax_h = asin(RC_gamma_h / Vgammamax_h_est) * 180/pi;;




figure(1)
plot(Vs,PR_s,'-',Vs,PA_s,'-',Vh,PR_h,'--',Vh,PA_h,'--',VPRmin_s,PRmin_s,'o',VPRmin_h,PRmin_h,'*')
axis([0 250 0 2e06])
% plot(V,Pr_s,'-',V,Pa_s,'-',Vh,Pr_h,'--',Vh,Pa_h,'--')
%plot(V,Pr_s,V,Pa_s,Vh,Pr_h,Vh,Pa_h)
title('Power Required and Available for CP1');
xlabel(' velocity (m/s)');
ylabel(' power (N m/s)');
legend('PR at sea level','PA at sea level',['PR at altitude ' num2str(h) ' m'],['PA at altitude ' num2str(h) ' m']);
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;

% get a sense of corresponding speed with respect to P_A
P = PAs * sqrt(rho_h / rho_s);
v = 0.1:1:100;
p1 = [.5*rho_h*S*C_D_0 0 0 -P K*S*WS^2/(.5*rho_h)];
% p2 = [.5*rho_h*S*C_D_0 0 -P 0 K*S*WS^2/(.5*rho_h)];
f1 = polyval(p1,v);
% f2 = polyval(p2,v);

figure(2)
plot(v,f1)
title(['speed profile under P_A =', num2str(P)])
% plot(v,f1,v,f2)
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

figure(3)
plot(Vh_s,Vv_s,'-',Vh_h,Vv_h,'--');
title(' Hodograph for [CP1]');
legend('sea-level', '6000m')
axis([0 100 0 10]); 
xlabel(' Vh (m/s)');
ylabel('Vv (m/s)');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);
% plotdlg;  

  