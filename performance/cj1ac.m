% cj1ac.m - absolute ceiling and time to climb
% 
% created on: 25-Sep-00
% updated on: 27-Sep-00
%

% run data file first
cj1data;

for i = 1:500
    h(i) = 50 * (i-1); % altitude 0 - 25000 m
    [T,p,rho] = stdatm(h(i));
    sigma = rho / rho_s; % density ratio 
    % calculate thrust load
    TW = F_s * sigma / W;
    % V for RC max
    Vrcmax = sqrt(WS*TW/(3*rho*C_D_0)*(1+sqrt(1+12*K*C_D_0/TW^2)));
    % max rate of climb
    q_rcmax = 1/2 * rho * Vrcmax^2;
    RCmax_h(i) = Vrcmax * (TW - q_rcmax * C_D_0/WS - K * WS / q_rcmax);
end
    
 
figure(1) 
plot(RCmax_h,h,'-');
grid()
axis([0 45 0 20000])
title('ceiling');
xlabel(' RC max (m/s)');
ylabel(' altitude (m)');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
% plotdlg;  


% plot(h,InvMaxRC,'-')
% title('time to climb');
% ylabel(' 1 / RC max (m/s)');
% xlabel(' altitude (m)');
% axis([0 6100 0 0.1]);
% plotdlg;
%   