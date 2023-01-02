                        %% Created by Mo7aMeD Adel %%
                     %% Computitional Fluid Dynamics %%
                            %% 1 / 5 / 2016 %%
clc
clear all
close all

% Notes:
% 1) That this code Solves Joukowski airfoil ONLY using Potintial Method.

%% Airfoil Parameters
t_c = 0.06;     % Max. Thickness to Chord ratio
C_c = 0.05;     % Max. Camber to Chord ratio
Chord = 1;      % Chord length

V_inf = 100;    % Free streem velocity in m/sec.
Alpha = 10*pi./180;      % AOA in degrees
%% Circle Parameters
b = Chord/4;
e = t_c/1.3;
B = 2*C_c;
a = b*(1+e)/cos(B);
xo = -b*e;
yo = a*B;
%% Airfoil Coordinates
[x_upper,y_upper,x_lower,y_lower,theta,r] = JFoil(t_c,C_c,Chord);

%% Solution
[X,Y] = pol2cart(theta,r);
x_dash = X-xo;
y_dash = Y-yo;
[theta_dash,r_dash] = cart2pol(x_dash,y_dash);
vr_dash = V_inf.*(1-a^2./r_dash.^2).*cos(theta_dash-Alpha);
vt_dash = -V_inf.*(sin(theta_dash-Alpha).*(1+a^2./r_dash.^2)+2.*(a./r_dash).*sin(Alpha+B));
AA = vr_dash.*cos(theta_dash)-vt_dash.*sin(theta_dash);
BB = -(vr_dash.*sin(theta_dash)+vt_dash.*cos(theta_dash));
CC = 1-b^2./r.^2.*cos(2.*theta);
DD = b^2./r.^2.*sin(2.*theta);
V = sqrt((AA.^2+BB.^2)./(CC.^2+DD.^2));
Cp = 1-(V./V_inf).^2;
u = (AA.*CC+BB.*DD)./(CC.^2+DD.^2);
v = (BB.*CC-AA.*DD)./(CC.^2+DD.^2);
V_u = V(100:length(x_upper));
u_u = u(100:length(x_upper));
v_u = v(100:length(x_upper));
Cp_u = 1-(V_u./V_inf).^2;
V_l = V(length(x_upper)+1:end-100);
u_l = u(length(x_upper)+1:end-100);
v_l = v(length(x_upper)+1:end-100);
Cp_l = 1-(V_l./V_inf).^2;

%% Plots
% Airfoil
figure
hold on 
grid on
plot(x_upper,y_upper,'r')
plot(x_lower,y_lower,'g')
axis equal
title(['Joukowski Airfoil with t_m_a_x/Chord = ',   num2str(t_c),',     Camber_m_a_x/Chord = ',   num2str(C_c)])
xlabel('x')
ylabel('y')
% Circle
figure
hold on
grid on
% polar(theta,r,'r')
polar(theta_dash,r_dash,'b')
axis equal

% Velocity Distribution on airfoil
figure
hold on
grid on
plot(x_upper,y_upper,x_lower,y_lower,'LineWidth',2,'color','b')
plot(x_upper(100:end),V_u./V_inf,'LineWidth',1.5,'color','r')
plot(x_lower(1:end-100),V_l./V_inf,'LineWidth',1.5,'color','g')
title('Velocity Distribution on the Joukowski Airfoil Surface')
xlabel('x')
ylabel('V/V_\infty')
legend('Airfoil surface','Airfoil surface','Upper Surface','Lower Surface')
% Pressure Distribution on airfoil
figure
hold on
grid on
plot(x_upper,y_upper,x_lower,y_lower,'LineWidth',2,'color','b')
plot(x_upper(100:end),Cp_u,'LineWidth',1.5,'color','r')
plot(x_lower(1:end-100),Cp_l,'LineWidth',1.5,'color','g')
title('Pressure Distribution on the Joukowski Airfoil Surface')
xlabel('x')
ylabel('C_P')
legend('Airfoil surface','Airfoil surface','Upper Surface','Lower Surface')
