% cp1ac.m - absolute and service ceilings 
% 
% created on: 10.Feb.2020
% updated on: 
%

% run data file first
cp1data;

for i = 1:500
    h(i) = 50 * (i-1); % altitude 0 - 25000 m
    [T,p,rho] = stdatm(h(i));
    sigma = sqrt(rho / rho_s); % density ratio 
    % calculate power
    PA_h = PA_s * sigma;
    % min PR
    PRmin = W * sqrt(K*C_D_0) * 4/sqrt(3) * sqrt(2*WS/sqrt(3*C_D_0/K)/rho);
    % max RC
    RCmax_h(i) = PA_h / W - PRmin/W;
end
    
figure(1) 
plot(RCmax_h,h,'-')
axis([0 10 0 15000])
title('ceiling');
xlabel(' RC max (m/s)');
ylabel(' altitude (m)');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
plotdlg;  

% plot(h,InvMaxRC,'-')
% title('time to climb');
% ylabel(' 1 / RC max (m/s)');
% xlabel(' altitude (m)');
% axis([0 6100 0 0.1]);
% plotdlg;
  