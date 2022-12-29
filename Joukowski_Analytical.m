clear all 
close all 
clc


% Grid parameters 

i_max = 50;              % eta_1 grid size 
j_max                   % eta_2 grid size


%  Joukowski Transformation Parameters.

c = 1;
b=c/4;                  
r=2*b;                  % Radius of the circle
e=t_max_c/1.3;          
beta=2*C_max_c;
a=b*(1+e)/cos(beta);
x0 = -b*e;
y0 = a*beta;

% angle gradient with step Delta theta. and size with the eta_1 grid size
D_theta = 2*pi()/(i_mas -1);
Theta_vec = 0:D_theta:2*pi();

% Airfoil coordinates 

dtheta = 2*pi()/(i_max-1); 
theta_vec = 0:dtheta:2*pi();

% New edits
%% Joukowski Analytical solution

% Circle Shift in Joukowski Z Plane

x0 = -b*e;              
y0 = a*beta;


theta_dash_vec = 0:Delta_theta:2*pi;
r_dash = a;
x_dath 
    
                                                                        