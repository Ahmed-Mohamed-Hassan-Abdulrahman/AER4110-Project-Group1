%-------------------------------------------------------------------------
%                     JOUKOWSKI TRANSFORMATION
%-------------------------------------------------------------------------
% Created by:                   Dario Isola, Politecnico di Milano, Italy
% Date:                                                 29, october, 2005
%-------------------------------------------------------------------------
clear all
close all
clc
disp('-------------------------------------------------------------------')
disp('  Joukowski Transformation Input Manager                           ')
disp('-------------------------------------------------------------------')
v_inf = input('  Asymptotic Speed Modulus [m/s]: ');
v = v_inf/v_inf;
theta = input('  Asymptotic Speed Angle [deg]: ');
theta = theta*pi/180;
disp('-------------------------------------------------------------------')
s_x = input('  Circle Origin, X_0 [m]: ');
s_y = input('  Circle Origin, Y_0 [m]: ');
s = s_x + i*s_y;
r = input('  Radius [m]: ');
disp('-------------------------------------------------------------------')
disp('  If Solution visualization is uncorrect try modify Tolerance TOLL')
disp('-------------------------------------------------------------------')
% FLUID PARAMETER
rho = 1.225;
% TRANSFORMATION PARAMETER
lambda = r-s;
% CIRCULATION
beta = (theta);
k = 2*r*v*sin(beta);
Gamma = k/(2*pi); %CIRCULATION
%COMPLEX ASYMPTOTIC SPEED 
w = v * exp(i*theta);
%TOLLERANCE
toll = +5e-2;
% GENERATING MESH
x = meshgrid(-5:.1:5);
y = x';
% COMPLEX PLANE
z = x + i*y;
% Inside-circle points are Excluded!
for a = 1:length(x)
    for b = 1:length(y)
        if abs(z(a,b)-s) <=  r - toll
            z(a,b) = NaN;
        end
    end
end
% AERODYNAMIC POTENTIAL
f = w*(z) + (v*exp(-i*theta)*r^2)./(z-s) + i*k*log(z);
% JOUKOWSKI TRANSFORMATION, 
J = z+lambda^2./z;
%GRAPHIC - Circle and Joukowski Airfoil
angle = 0:.1:2*pi;
z_circle = r*(cos(angle)+i*sin(angle)) + s;
z_airfoil = z_circle+lambda^2./z_circle;
% KUTTA JOUKOWSKI THEOREM
L = v_inf*rho*Gamma;
L_str = num2str(L);
%PLOTTING SOLUTION
figure(1)
hold on
contour(real(z),imag(z),imag(f),[-5:.2:5])
fill(real(z_circle),imag(z_circle),'y')
axis equal
axis([-5 5 -5 5])
title(strcat('Flow Around a Circle.   Lift:  ',L_str,'  [N/m]'));
figure(2)
hold on
contour(real(J),imag(J),imag(f),[-5:.2:5])
fill(real(z_airfoil),imag(z_airfoil),'y')
axis equal
axis([-5 5 -5 5])
title(strcat('Flow Around the Corresponding Airfoil.   Lift:  ',L_str,'  [N/m]'));
