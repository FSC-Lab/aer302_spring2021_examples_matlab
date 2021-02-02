%	Paper Airplane Flight Path
%	Copyright 2005 by Robert Stengel
%	August 23, 2005
%   modified by Hugh Liu
%   updated on: 09 December 2008
%   revised on: 02 Oct 2018

    close all
    clear all
	global CL CD S m g rho	
	S		=	0.017;			% Reference Area, m^2
	AR		=	0.86;			% Wing Aspect Ratio
	e		=	0.9;			% Oswald Efficiency Factor;
	m		=	0.003;			% Mass, kg
	g		=	9.8;			% Gravitational acceleration, m/s^2
	rho		=	1.225;			% Air density at Sea Level, kg/m^3	
	CLa		=	3.141592 * AR/(1 + sqrt(1 + (AR / 2)^2));
                            	% Lift-Coefficient Slope, per rad
	CDo		=	0.02;			% Zero-Lift Drag Coefficient
	kappa	=	1 / (3.141592 * e * AR);% Induced Drag Factor
    
%   optimal Condition at (CL/CD)max    
	CL_star 	=	sqrt(CDo / kappa);        % CL for Maximum Lift/Drag Ratio
	CD_star     =	CDo + kappa * CL_star^2;  % Corresponding CD
	LD_star     =	CL_star / CD_star;        % Maximum Lift/Drag Ratio
	Gamma_star	=	-atan(1 / LD_star);               % Corresponding Flight Path Angle, rad
    Gamma_starDeg = Gamma_star * 180/pi;
	V_star      =	sqrt(2 * m * g * cos(Gamma_star)/(rho * S * (CL_star)));
    %Ve      =	sqrt(2 * m * g /(rho * S * (CLe * cos(Gammae) - CDe * sin(Gammae))));
                                            % Corresponding Velocity, m/s
	Alpha_star  =	CL_star / CLa;                  % Corresponding Angle of Attack, rad
    Alpha_starDeg = Alpha_star * 180/pi;
    
%   equilibrium condition  under a fixed alpha  
	Alphae   =   8 * pi/180;       % assume constant angle of attack, rad
	%Alpha   =   Alphae + 10 * pi/180;       % angle of attack, rad
    CLe      =   CLa * Alphae;
    CDe      =   CDo + kappa * CLe^2;
    CL = CLe;
    CD = CDe;
    Gammae	=	-atan(1 / (CLe/CDe));  
    Ve      =	sqrt(2 * m * g * cos(Gammae)/(rho * S * (CLe)));

% initial conditions    
    Vo      =   Ve;         % initial velocity, m/s
	H		=	2;			% Initial Height, m
	R		=	0;			% Initial Range, m
    Gammao  =   Gammae;     % Initial flight path angle, rad
	to		=	0;			% Initial Time, sec
	tf		=	6;			% Final Time, sec
	tspan	=	[to tf];

%	a) Equilibrium Glide at Maximum Lift/Drag Ratio
	xo		=	[Ve;Gammae;H;R];
	[ta,xa]	=	ode23('EqMotion',tspan,xo);
	
%	b) Oscillating Glide due to Zero Initial Flight Path Angle
	xo		=	[Ve;0;H;R];
	[tb,xb]	=	ode23('EqMotion',tspan,xo);

%	c) Effect of Increased Initial Velocity
	xo		=	[1.5*Ve;0;H;R];
	[tc,xc]	=	ode23('EqMotion',tspan,xo);

%	d) Effect of Further Increase in Initial Velocity
	xo		=	[3*Ve;0;H;R];
	[td,xd]	=	ode23('EqMotion',tspan,xo);
    
    %   Linearized system around equilibrium
    A = [-rho*Ve*CD*S/m -g*cos(Gammae);
       rho*CL*S/m g*sin(Gammae)/Ve];
    B = zeros(2,1);
    C = eye(2);
    D = zeros(2,1);
    Sys1 = ss(A,B,C,D);
    X0a = [Ve,Gammae]';
    X0b = [Ve,0]';
    X0c = [1.5*Ve,0]';
    X0d = [3*Ve,0]';
        
    tl = to:0.01:tf;
    ya = initial(Sys1,X0a,tl);
    yb = initial(Sys1,X0b,tl);
    yc = initial(Sys1,X0c,tl);
    yd = initial(Sys1,X0d,tl);
    
    Ahat = [2*g*sin(Gammae)/Ve        -g*cos(Gammae);
            2*g*cos(Gammae)/(Ve^2)     g*sin(Gammae)/Ve];
    Bhat = [2*g*cos(Gammae) * kappa * CLa
            g*cos(Gammae) / (Ve*Alphae)];
    Chat = eye(2);
    Dhat = zeros(2,1);
    Sys2 = ss(Ahat,Bhat,Chat,Dhat);
    deltaAlpha = 10*pi/180;
    y = step(Sys2,tl) .* deltaAlpha;
    Vl = y(:,1) + Ve;
    Gammal = y(:,2) + Gammae;

	
   
        
        
   
    figure(1)
	plot(xa(:,4),xa(:,3),xb(:,4),xb(:,3),xc(:,4),xc(:,3),xd(:,4),xd(:,3))
    legend('V_e,\gamma_e', 'V_e, \gamma_0=0', '1.5V_e, \gamma_0=0','3V_e, \gamma_0=0')
	xlabel('Distance, m'), ylabel('Height, m'), grid,
    title('Longitudinal Trajectory')
    htype = findobj(gcf,'type','line');
    set(htype,'linewidth',2);
%     plotdlg

	figure(2)
	subplot(2,2,1)
	plot(ta,xa(:,1),tb,xb(:,1),tc,xc(:,1),td,xd(:,1))
    legend('V_e,\gamma_e', 'V_e, \gamma_0=0', '1.5V_e, \gamma_0=0','3V_e, \gamma_0=0')
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
	subplot(2,2,2)
	plot(ta,xa(:,2),tb,xb(:,2),tc,xc(:,2),td,xd(:,2))
	xlabel('Time, s'), ylabel('Flight Path Angle, rad'), grid
	subplot(2,2,3)
	plot(ta,xa(:,3),tb,xb(:,3),tc,xc(:,3),td,xd(:,3))
	xlabel('Time, s'), ylabel('Altitude, m'), grid
	subplot(2,2,4)
	plot(ta,xa(:,4),tb,xb(:,4),tc,xc(:,4),td,xd(:,4))
	xlabel('Time, s'), ylabel('Range, m'), grid
    htype = findobj(gcf,'type','line');
    set(htype,'linewidth',2);

    figure(3)
	subplot(2,1,1)
	plot(ta,xa(:,1),tl,Vl)
	%plot(tl,Vl)
    title('V_e \gamma_e')
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
    legend('no control','linearized control, \Delta\alpha = 10 deg')
	subplot(2,1,2)
	plot(ta,xa(:,2),tl,Gammal)
	%plot(tl,Gammal)
	xlabel('Time, s'), ylabel('Flight Path Angle, rad'), grid
    htype = findobj(gcf,'type','line');
    set(htype,'linewidth',2);

    figure(4)
    title('nonlinear vs linear, initial response')
	subplot(2,2,1)
	plot(ta,xa(:,1),tl,ya(:,1))
    legend('numerical','linearized')
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
    title('initial: V_e \gamma_e')
	subplot(2,2,2)
	plot(tb,xb(:,1),tl,yb(:,1))
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
    title('initial: V_e \gamma_0=0')
    subplot(2,2,3)
	plot(tc,xc(:,1),tl,yc(:,1))
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
    title('initial: 1.5V_e \gamma_0=0')
    subplot(2,2,4)
	plot(td,xd(:,1),tl,yd(:,1))
	xlabel('Time, s'), ylabel('Velocity, m/s'), grid
    title('initial: 3V_e \gamma_0=0')
    htype = findobj(gcf,'type','line');
    set(htype,'linewidth',2);

    figure(5)
    title('nonlinear vs linear, initial response')
	subplot(2,2,1)
	plot(ta,xa(:,2),tl,ya(:,2))
    legend('numerical','linearized')
	xlabel('Time, s'), ylabel('\gamma'), grid
    title('initial: V_e \gamma_e')
	subplot(2,2,2)
	plot(tb,xb(:,2),tl,yb(:,2))
	xlabel('Time, s'), ylabel('\gamma'), grid
    title('initial: V_e \gamma_0=0')
    subplot(2,2,3)
	plot(tc,xc(:,2),tl,yc(:,2))
	xlabel('Time, s'), ylabel('\gamma'), grid
    title('initial: 1.5V_e \gamma_0=0')
    subplot(2,2,4)
	plot(td,xd(:,2),tl,yd(:,2))
	xlabel('Time, s'), ylabel('\gamma'), grid
    title('initial: 3V_e \gamma_0=0')
    htype = findobj(gcf,'type','line');
    set(htype,'linewidth',2);
