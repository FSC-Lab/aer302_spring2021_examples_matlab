% cj1er.m - endurance and range
% 
% created on: 24-Sep-00
% updated on:
%

% run data file first
cj1data;

% altitude
h = 6705.6; % (m)
[T,p,rho] = stdatm(h);

% lift-to-drag ratio
for i = 1:400
    V(i) = 10 + 1 * i;
    % at sea level
    C_L_s = W / (0.5 * rho_s * V(i)^2 * S);
    C_D_s = C_D_0 + K * C_L_s^2;
    LD_s(i) = C_L_s / C_D_s;
    LsD_s(i) = sqrt(C_L_s)  / C_D_s;
    % at altitude
    C_L = W / (0.5 * rho * V(i)^2 * S);
    C_D = C_D_0 + K * C_L^2;
    LD(i) = C_L / C_D;
    LsD(i) = sqrt(C_L)  / C_D;
end

% maximum ratio value
LD_sMax = max(LD_s);
LsD_sMax = max(LsD_s);
LDMax = max(LD);
LsDMax = max(LsD);

% calculation of endurance and range
W0 = W;
W1 = W - Wf; % empty tank
Es = 1/c_t * LD_sMax * log(W0/W1) / 3600; % hour
E  = 1/c_t * LDMax * log(W0/W1) / 3600; % hour
Rs = 2 * sqrt(2 / (rho*S)) * 1 / c_t * LsD_sMax * (sqrt(W0) - sqrt(W1));
R = 2 * sqrt(2 / (rho*S)) * 1 / c_t * LsDMax * (sqrt(W0) - sqrt(W1));

% presentation of results
plot(V,LD_s,'-',V,LsD_s,'-',V,LD,'--',V,LsD,'--');
legend('C_L/C_D sea level', '\sqrt(C_L)/C_D sea level', 'C_L/C_D h', '\sqrt(C_L)/C_D h')
grid;
title('C_L to C_D ratio');
xlabel(' velocity (m/s)');
ylabel(' ');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
% plotprn;
%   
% disp([' the maximum C_L/C_D calcualted from sea level is:       ', num2str(LD_sMax)]);
% disp([' the endurance is:                                       ', num2str(Es), '  (hour)']);
% disp([' the maximum C_L^1/2/C_D calcualted from sea level is:   ' num2str(LsD_sMax)]);
% disp([' the range is:                                           ' num2str(Rs) ' (m)']);
% disp([' the maximum C_L/C_D calcualted from altitude  is:       ', num2str(LDMax)]);
% disp([' the endurance is:                                       ', num2str(E), '  (hour)']);
% disp([' the maximum C_L^1/2/C_D calcualted from altitude  is:   ' num2str(LsDMax)]);
% disp([' the range is:                                           ' num2str(R) ' (m)']);
