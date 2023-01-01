clear all 
close all 
clc

%% Initialization (Inputs)

% The Free Stream Velocity

Vinf=100;

% Enter the Joukowski Airfoil Parameters

c=1;              % Chord
C_max_c=0.04;     % Maximum Camber/Chord Percentage
t_max_c=0.05;     % Maximum Thickness/Chord Percentage
AoA=4*pi/180;     % Angle of Attack of flow Percentage

% Grid parameters 

i_max = 100;       % eta_1 grid size 
j_max = 100;       % eta_2 grid size

%  Joukowski Transformation Parameters.

%% circle Parameters 

b=c/4;                  
e=t_max_c/1.3;          
beta=2*C_max_c;
a=b*(1+e)/cos(beta);
                  
% Circle Shift in Joukowski Z Plane
x0 = -b*e;
y0 = a*beta;

%% Z' Plane 


D_theta = 2*pi()/(i_max-1);                % angle step size in the domain  Delta theta = 2*pi/(imax-1) 
Theta_dash_vec = 0:D_theta:2*pi();       % theta vector in z_dahs Plane 
r_dash = a*ones(1,100);                                % Radius of the circle

% x'-y' coords in z_dahs plane 

x_dash = r_dash.*cos(Theta_dash_vec);
y_dash = r_dash.*sin(Theta_dash_vec);

%% Z Plane 

% x-y coords in z_dahs plane 

x = x_dash + x0;
y = y_dash + y0;

r_ = sqrt(x.^2+ y.^2);                  % the radious of the circle in Z plane 

theta_vec = atan2(y,x);

%% Z1 Plane , the Airfoil Plane

% x1-y1 coords in z_dahs plane 

x1 = x.*(1+(b^2)./(x.^2+y.^2));
y1 = y.*(1-(b^2)./(x.^2+y.^2));

r1 = sqrt(x1.^2+ y1.^2);                % the radious of the circle in Z1 plane

theta1_vec = atan2(y1,x1);

%% Joukowski Analytical solution
V_rDash = Vinf*(1-(a./r_dash).^2).*cos(Theta_dash_vec-AoA);
V_thetaDash = -Vinf *(sin(Theta_dash_vec-AoA).*(1+(a./r_dash).^2)+2*(a./r_dash)*sin(AoA+beta));

% Velocity Magnitude over the Airfoil in Z1 Plane
V1 = sqrt(V_thetaDash.^2./(1-2*(b./r_).^2.*cos(2*theta_vec)+(b./r_).^4));
% pressure Coefficient
C_p = 1-(V1/Vinf).^2;

%% airforil Coordinates with Formula  
 X = 2*b*cos(Theta_dash_vec);
 Y = 2*b*e*(1-cos(Theta_dash_vec)).*sin(Theta_dash_vec)+2*b*beta*sin(Theta_dash_vec).^2;

%% Plot commands 

% plot the Z Circle and the Circle shift

figure('Name', 'Joukowski Circles' )
        plot(x_dash, y_dash)
        grid on 
        axis([-0.5 0.5 -0.5 0.5])
        axis('equal')
        xlabel('$x$', 'interpreter', 'latex')
        ylabel('$y$', 'interpreter', 'latex')
hold on
        plot(x, y)
        grid on 
        axis([-0.5 0.5 -0.5 0.5])
        axis('equal')
        xlabel('$x''$', 'interpreter', 'latex')
        ylabel('$y''$', 'interpreter', 'latex')


figure('Name', 'Joukowski Airfoil comparision' )

      plot(x1, y1)
      grid on
hold on
      plot(X, Y, '*')
      grid on
      axis('equal')


figure('Name', 'Velocity over the Airfoil' )
      plot(x1, V1, '-')
      grid on
figure('Name', 'Velocity over the Airfoil' )
      plot(x1, C_p, '-')
      grid on


