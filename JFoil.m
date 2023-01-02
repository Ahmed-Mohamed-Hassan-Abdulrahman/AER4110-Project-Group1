function [x_upper,y_upper,x_lower,y_lower,theta,r] = JFoil(t_c,C_c,Chord)
%% Defenetions 
%
% Outputs:-
% x_upper : x Coordinates of upper surface
% y_upper : y Coordinates of upper surface
% x_lower : x Coordinates of lower surface
% y_lower : y Coordinates of lower surface
% r       : Radius of the circle
% 
% Inputs:-
% t_c     : Max. Thickness to Chord ratio
% C_c     : Max. Camber to Chord ratio
% Chord   : Chord Length
% 

    disp(['Joukowski Airfoil:  Thickness/Chord   =  ',  num2str(t_c)])
    disp(['                    Camber/Chord      =  ',  num2str(C_c)])
    disp(['                    Chord             =  ',  num2str(Chord)])
%% Circle Parameters
b = Chord/4;
e = t_c/1.3;
B = 2*C_c;
a = b*(1+e)/cos(B);
xo = -b*e;
yo = a*B;
%% Airfoil Coordinates
theta = 0:0.0001:2*pi;
for i = 1:length(theta)
    x1(i) = 2*b*cos(theta(i));
    y1(i) = 2*b*e*(1-cos(theta(i)))*sin(theta(i))+2*b*B*sin(theta(i))^2;
    r(i)  = b*(1+e*(1-cos(theta(i)))+B*sin(theta(i)));
end
x = x1+Chord/2;
y = y1;
for m = 1:length(x)-1
    if y(m)>=0 && y(m+1)<=0
        le = m;
    end
end
x_upper = x(1:le);
x_upper(1) = Chord;
x_upper(end) = 0;
y_upper = y(1:le);
y_upper(1) = 0;
y_upper(end) = 0;
x_lower = x(le+1:end);
x_lower(1) = 0;
x_lower(end) = Chord;
y_lower = y(le+1:end);
y_lower(1) = 0;
y_lower(end) = 0;
end