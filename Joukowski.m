%% Flow Over Joukowski Airfoil. 
% This Fucion Calculates the Velocity and pressure Distribution over an
% Airforil by using the analysical joukowski transformaiotn.
% Funcion Argumentd : (joukowski(Vinf, AoA, c, C_max_c, t_max_c) )
%
%   V_inf: the free streem Velocity.
%     AoA: angle of attack
%       c: cord line length
% C_max_c: max camber (% of cord)
% t_max_c: max thickness (% of cord)
%   i_max: number of points  
%     
% Published by: Mohamed Tarek Mohamed Amien 
% publishing year: 1st January /2023

function [V1, C_p] = Joukowski(Vinf, AoA, c, C_max_c, t_max_c,  i_max ) 
%% Initialization (Inputs)

% The Free Stream Velocity

Vinf=100;

% Enter the Joukowski Airfoil Parameters

% c=1;              % Chord
% C_max_c=0.04;     % Maximum Camber/Chord Percentage
% t_max_c=0.05;     % Maximum Thickness/Chord Percentage
   
% Grid parameters 

%  i_max = 10;       % eta_1 grid size 
% j_max = 100;       % eta_2 grid size

%  Joukowski Transformation Parameters.

%% joukowski circle Parameters 

b=c/4;                  
e=t_max_c/1.3;          
beta=2*C_max_c;
a=b*(1+e)/cos(beta);
                  
% Circle Shift in Joukowski Z Plane
x0 = -b*e;
y0 = a*beta;

%% Z' Plane 

D_theta = 2*pi()/(i_max-1);              % angle step size in the domain  Delta theta = 2*pi/(imax-1) 
Theta_dash_vec = 0:D_theta:2*pi();       % theta vector in z_dahs Plane 
r_dash = a*ones(1,i_max);                  % Radius of the circle

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

V1_x(1:1:100,:) = V1.*cos(theta1_vec).*ones(100,1);
V1_y(1:1:100,:) = V1.*sin(theta1_vec).*ones(100,1);
% pressure Coefficient
C_p = 1-(V1/Vinf).^2;

%% airforil Coordinates with Formula  
 X = 2*b*cos(Theta_dash_vec);
 Y = 2*b*e*(1-cos(Theta_dash_vec)).*sin(Theta_dash_vec)+2*b*beta*sin(Theta_dash_vec).^2;

 

%% Plot commands 

% plot the Z Circle and the Circle shift

% figure('Name', 'Joukowski Circles' )
%         plot(x_dash, y_dash)
%         grid on 
%         axis([-0.5 0.5 -0.5 0.5])
%         axis('equal')
%         xlabel('$x$', 'interpreter', 'latex')
%         ylabel('$y$', 'interpreter', 'latex')
% hold on
%         plot(x, y)
%         grid on 
%         axis([-0.5 0.5 -0.5 0.5])
%         axis('equal')
%         xlabel('$x''$', 'interpreter', 'latex')
%         ylabel('$y''$', 'interpreter', 'latex')


figure('Name', 'Joukowski Airfoil check' )

      plot(x1, y1)
      grid on
hold on
      plot(X, Y, '.')
      grid on
      axis('equal')


figure('Name', 'Velocity and Cp distribution over the Airfoil' )
tiledlayout(2,1);
nexttile 
      hold on 
      plot(x1, V1, '-' ,'LineWidth',1.5,'color','red')
      plot(x1, 700*y1,'LineWidth',0.5,'color','black')
      fill(x1, 700*y1,'cyan')
      grid on
      grid on
      xlabel('$x_1$', 'interpreter', 'latex')
      ylabel('$V_1$', 'interpreter', 'latex')
      legend('Velocity distribution','Airforil')
      title('Velocity distribution', 'FontName','lm roman 9')
nexttile
      hold on
      fill(x1, 20*y1,'cyan')
      plot(x1, 20*y1,'LineWidth',0.5,'color','black')
      plot(x1, C_p, '-','LineWidth',1.5,'color','blue')
      grid on
      xlabel('$x_1$', 'interpreter', 'latex')
      ylabel('$C_p$', 'interpreter', 'latex')
      legend('Cp distribution','Airforil')
      title('Pressure Coefficient distribution', 'FontName','lm roman 9')


 end
