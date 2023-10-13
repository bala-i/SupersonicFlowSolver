% Program calculates the properties of a supersonic flow field over a cone given ...
% freestream mach number and oblique shock angle
% By Bala Iganus
% 4/13/20
% Whole world is shutdown due to covid-19. literally

clear all
clc

% Input values
M_1 = 2;                    % upstream mach number
gamma = 1.4;                  % gamma for diatomic gas air
theta_s = 36 * pi/180;        % shock wave angle


% Initial conditions
v = (  (2/(gamma -1)) * (1/M_1)^2 +1) ^ (-1/2); % velocity

t=[];
t(1) = 29 *pi/180;    % theta
h = -0.001 *pi/180; % increment for runge kutta

y1= [];               % V_r
y2= [];               % V theta
y1(1) = v *cos(theta_s);
y2(1) = -v*sin(theta_s) * ( ( (gamma-1)*M_1^2 *sin(theta_s)^2 +2) ...
    / ( (gamma +1) * M_1^2 * sin(theta_s)^2));

% Differential equations
f1 = @(t, y1, y2) y2;
f2 = @(t, y1, y2) ( y2^2 * y1 - ((gamma-1)/2)* (1-y1^2 - y2^2) * ...
    ((2*y1) + (y2*cot(t)))) / ( ((gamma-1)/2) * (1- y1^2 -y2^2) - y2^2);
count = 1;

vector = [29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11] .* (pi/180); % indexing purposes
vtheta = zeros(1,length(vector));
vr = zeros(1,length(vector));
i=1;
% Runge Kutta
while y2< 0
    k1 = h*f1(t(i), y1(i), y2(i));
    m1 = h*f2(t(i), y1(i), y2(i));
    k2 = h*f1(t(i)+0.5*h, y1(i)+ 0.5*k1, y2(i)+ 0.5*m1);
    m2 = h*f2(t(i)+0.5*h, y1(i)+ 0.5*k1, y2(i)+0.5*m1);
    k3 = h*f1(t(i)+0.5*h, y1(i)+ 0.5*k2, y2(i)+0.5*m2);
    m3 = h*f2(t(i)+0.5*h, y1(i)+ 0.5*k2, y2(i)+0.5*m2);
    k4 = h*f1(t(i)+h, y1(i)+k3, y2(i)+m3);
    m4 = h*f2(t(i)+h, y1(i)+k3, y2(i)+m3);

    t(i+1) = t(i)+h;
    y1(i+1) = y1(i) + (k1+2*k2+2*k3+k4)/6;
    y2(i+1) = y2(i) + (m1+2*m2+2*m3+m4)/6;
   
    i=i+1;
    count= count+1;
    
end

theta_c = t(end)* 180/pi;                      % cone angle
v_c = sqrt(y1(end)^2 + y2(end)^2);             % cone surface velocity

m_n1 = M_1 * sin(theta_s);                     % M_n1

m_n2 = sqrt((m_n1^2 + (2/(gamma - 1) )) / ( (2*gamma / (gamma-1)) ...
    * m_n1^2 -1));                             % M_n2

m_c = 1/ sqrt( ( (gamma - 1)*((1/v_c^2) -1) ) / 2);    % mach number at cone surface


pressure_sw = 1+ ( 2*gamma*(M_1^2 *sin(theta_s)^2 -1))/(gamma+1); %Pressure ratio p2/p1
density_sw = ( (gamma+1) * M_1^2 * sin(theta_s)^2) / ( (gamma-1)*M_1^2 *sin(theta_s)^2 +2); % rho2/rho1
temperature_sw = pressure_sw * (1/density_sw);                % T2/T1
theta = atan(2* cot(theta_s) * ( (M_1^2 * sin(theta_s)^2 -1) ...
    / (M_1^2 *(gamma + cos(2*theta_s)) +2)));              % theta from theta beta M eqtn
M_2 = m_n2/sin(theta_s - theta);                          % M_2


% index wanted values for vr, vtheta, and thetaconical flow
l=1;
p=1;
thetaconicalflow= zeros(1,19);   

for op=1:18
    thetaconicalflow(p) = t(l);
    vr(p) = y1(l);
    vtheta(p) = y2(l);
    l=l+999;          % best rn is at 999 at h=-.001 interval of one degree
    p=p+1;
end

thetaconicalflow(end)= t(end);    % matches vector to first decimal point, excellent approx
vr(end) = y1(end);
vtheta(end) = y2(end);




%%% before shock properties
theta1 = (60:-10:30)*pi/180;                     % theta 
M1 = 2.2 * ones(1,length(theta1));               % mach number
v1 = (  (2/(gamma -1)) .* (1./M1).^2 +1) .^ (-1/2); % velocity
vtheta1= -v1.*sin(theta1);                        
vr1= v1.*cos(theta1);
%%% end before shock properties


% combine full flow vectors
vr=[vr1,vr];
vtheta=[vtheta1,vtheta];
fullthetarange = [theta1, thetaconicalflow]*180/pi;

% flow inside cone
vconical = sqrt(vr.^2 + vtheta.^2);   % velocity along rays
mconical = 1./ sqrt( ( (gamma - 1).*((1./vconical.^2) -1) ) ./ 2); % mach # along rays

%%% more full flow vector(s)
machno= [M1,mconical];                    % complete flow mach numbers up and downstream
%%% end more full flow vector(s)

% Pressure
ptot2 = 1./ ((1+ ((gamma-1)./2).* mconical.^2).^(gamma/(gamma-1))); % p/po2
p2ptot2 = 1./ ((1+ ((gamma-1)./2).* M_2.^2).^(gamma/(gamma-1)));    % P2/Po2

pp1 = ( ptot2./ p2ptot2) .* pressure_sw;                           % p/p_1


% Temperature
Ttot2 = 1./ (1+ ( (gamma-1)/2 .*mconical.^2));      % T/To2
T2Ttot2 = 1./ (1+ ( (gamma-1)/2 *M_2^2));           % T2/To2

TT1 = ( Ttot2./T2Ttot2) .* temperature_sw;          % T/T_1

% Density
rhotot2 = 1./ (1+ ((gamma-1)/2).* mconical.^2).^(1/(gamma-1)); %Rho/Rho_o2
rho2rhotot2 = 1./ (1+ ((gamma-1)/2).* M_2.^2).^(1/(gamma-1));  %Rho_2 / Rho_o2

rhorho1 = ( rhotot2 ./ rho2rhotot2) .* density_sw;             %Rho/Rhp_1




%%%%% Index flow properties between 12 deg and theta_c
% all flow properties are between 12 deg and theta_c for table
% step size of 0.1 deg

x=17001;
d=1;
thetalil = zeros(1,10);      %theta
vrlil = zeros(1,10);         %V_r
vthetalil = zeros(1,10);     %V_theta

for yg=1:9
    thetalil(d)=t(x);
    vrlil(d) = y1(x);
    vthetalil(d)= y2(x);
    x=x+97;
    d=d+1;
end

thetalil(end) = t(end);
vrlil(end) = y1(end);
vthetalil(end) = y2(end);

vconicallil = sqrt(vrlil.^2 + vthetalil.^2);   % velocity along rays
mconicallil = 1./ sqrt( ( (gamma - 1).*((1./vconicallil.^2) -1) ) ./ 2); % mach # along rays
ptot2lil = 1./ ((1+ ((gamma-1)./2).* mconicallil.^2).^(gamma/(gamma-1))); % p/po2
p2ptot2lil = 1./ ((1+ ((gamma-1)./2).* M_2.^2).^(gamma/(gamma-1)));    % P2/Po2

pp1lil = ( ptot2lil./ p2ptot2lil) .* pressure_sw;                           % p/p_1
Ttot2lil = 1./ (1+ ( (gamma-1)/2 .*mconicallil.^2));      % T/To2
T2Ttot2lil = 1./ (1+ ( (gamma-1)/2 *M_2^2));           % T2/To2

TT1lil = ( Ttot2lil./T2Ttot2lil) .* temperature_sw;          % T/T_1
rhotot2lil = 1./ (1+ ((gamma-1)/2).* mconicallil.^2).^(1/(gamma-1)); %Rho/Rho_o2
rho2rhotot2lil = 1./ (1+ ((gamma-1)/2).* M_2.^2).^(1/(gamma-1));  %Rho_2 / Rho_o2

rhorho1lil = ( rhotot2lil ./ rho2rhotot2lil) .* density_sw;             %Rho/Rhp_1

%%%%%%% end flow properties between 12 deg and theta_c



% plots

figure(1)
plot ((fullthetarange), pp1,(fullthetarange),TT1,...
(fullthetarange), rhorho1 );
set(gca, 'XDir','reverse')
grid on
ylim([0.8 1.6])
xlabel 'theta'
legend ('p/p_1', 'rho/rho_1', 'T/T_1')
title ' Flow profiles p/p_1, T/T_1, and rho/rho_1'



figure(2) 
plot((fullthetarange),mconical);
set(gca, 'XDir','reverse')
grid on
xlabel 'theta', ylabel 'M', title 'Mach number vs theta'
ylim([1.9 2.3])

 
figure(3)
plot((fullthetarange), vtheta);
set(gca, 'XDir','reverse')
grid on
xlabel 'theta', ylabel 'V_t_h_e_t_a', title 'V_t_h_e_t_a/V_m_a_x vs theta'
ylim([-0.7 0])

