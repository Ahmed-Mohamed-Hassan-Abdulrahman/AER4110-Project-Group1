clear all 
close all 
clc

%% Initialization (Inputs)

% The Free Stream Velocity

Vinf=100;

% Grid parameters 

i_max = 50;              % eta_1 grid size 
j_max = 50;                % eta_2 grid size

% Enter the Joukowski Airfoil Parameters

c=1;            % Chord
C_max_c=0.04;     % Maximum Camber/Chord Percentage
t_max_c=0.05;     % Maximum Thickness/Chord Percentage
AoA=4*pi/180;          % Angle of Attack of flow Percentage

%  Joukowski Transformation Parameters.

c = 1;
b=c/4;                  
r=2*b;                  % Radius of the circle
e=t_max_c/1.3;          
beta=2*C_max_c;
a=b*(1+e)/cos(beta);

% Circle Shift in Joukowski Z Plane
x0 = -b*e;
y0 = a*beta;

%% Z_Dash Plane 

% angle gradient with step Delta theta. and size with the eta_1 grid size

D_theta = 2*pi()/(i_max -1);

% theta vector in z_dahs Plane 

Theta_dash_vec = linspace(0,2*pi(),50);
r_dash= a;

% x-y coords in z_dahs plane 

x_dash = r_dash*cos(Theta_dash_vec);
y_dash = r_dash*sin(Theta_dash_vec);

%% Z Plane 

x = x_dash + x0;
y = y_dash + y0;
r_ = sqrt(x(1)^2+ y(1)^2);

theta = atan2(y,x);

%% Z1 Plane , the Airfoil Plane
x1 = x.*(1+(b^2)./(x.^2+y.^2));
y1 = y.*(1-(b^2)./(x.^2+y.^2));
r1 = sqrt(x1(1)^2+ y1(1)^2),

theta1 = atan2(y1,x1);

%% airforil Coordinates 
% x1 = 2*b*cos(Theta_dash_vec);
% y1 = 2*b*e*(1-cos(Theta_dash_vec)).*sin(Theta_dash_vec)+2*b*beta*sin(Theta_dash_vec).^2;

%% Joukowski Analytical solution
V_rDash = Vinf*(1-(a/r_dash)^2).*cos(Theta_dash_vec-AoA);
V_thetaDash = -Vinf *(sin(Theta_dash_vec-AoA).*(1-(a/r_dash)^2)+2*(a/r_dash)*sin(AoA+beta));

% Velocity Magnitude over the Airfoil in Z1 Plane
V1 = 0.01*sqrt(abs(V_thetaDash)./(1-2*(b/r_)^2.*cos(2*theta)+(b/r_)^4))
% pressure Coefficient
C_p = 1-(V1/Vinf)
%% Plot commands 

figure()
      plot(x_dash, y_dash)
      grid on 
   axis([-0.5 0.5 -0.5 0.5])
   axis('equal')
   xlabel('$x$', 'interpreter', 'latex')
figure()
      plot(x, y)
      grid on 
   axis([-0.5 0.5 -0.5 0.5])
   axis('equal')
   xlabel('$x$', 'interpreter', 'latex')
figure()
      plot(x1, y1)
      grid on
   axis('equal')
   xlabel('$x$', 'interpreter', 'latex')
hold on
      plot(x1, V1)
      grid on
   xlabel('$x$', 'interpreter', 'latex')