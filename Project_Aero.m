% O grid of airfoil using transformations
% Project made by Ahmed Mohamed Hassan Abdulrahman

clc
clearvars
close all

%% Initialization (Inputs)

% Choose the maximum mesh size:

i_max=49; j_max=50;

% Enter the Joukowski Airfoil Parameters

c=1;            % Chord
C_max_c=0.04;     % Maximum Camber/Chord Percentage
t_max_c=0.05;     % Maximum Thickness/Chord Percentage
AoA=4;          % Angle of Attack of flow Percentage

% Drawing Parameters

airfoil_points=500;     % Number of points on the airfoil; More points = Better accuracy
R=5*c;                  % Far field radius assumption

% The Transformation Parameters

eta1_max=1;     eta1_min=0;
eta2_max=1;     eta2_min=0;

%% Generating Airfoil Coordinates

% Joukowski airfoil parameters calculations

b=c/4;                  
r=2*b;                  % Radius of the circle
e=t_max_c/1.3;          
beta=2*C_max_c;
a=b*(1+e)/cos(beta);

% x coordinates of the points on the upper and lower surfaces of the airfoil:
x_airfoil_coord=linspace(-r,r,airfoil_points);

% y coordinates equation on the upper surface of the airfoil:
y_upper= @(x) 2*b*e*(1-x/2/b).*(sqrt(1-(x/2/b).^2))+2*b*beta*(1-(x/2/b).^2);
% y coordinates equation on the lower surface of the airfoil:
y_lower= @(x) 2*b*e*(1-x/2/b).*(-sqrt(1-(x/2/b).^2))+2*b*beta*(1-(x/2/b).^2);

airfoil_plot=plot(x_airfoil_coord,y_upper(x_airfoil_coord),x_airfoil_coord,y_lower(x_airfoil_coord));
hold on

%% Generating Circle around airfoil

 % We must divide the circumference of the circle by the number of
 % solving points (in i)
Delta_theta=2*pi/(i_max-1);
 % to plot a circle we need points in x and y, this can be
 % achieved by using x=r*cos(theta) and y=r*sin(theta)

for i=1:i_max
    x_circle_plot(i)=r*cos(Delta_theta*(i-1));  % x coordinates of circle
    y_circle_plot(i)=r*sin(Delta_theta*(i-1));  % y coordinates of circle
end

 % Plotting the plain circle
plot(x_circle_plot,y_circle_plot)

hold on

%% Plotting i_max lines from the center to the circle of the airfoil

for i=1:i_max
    plot([0 x_circle_plot(i)],[0 y_circle_plot(i)])
end
hold on

%% Projecting the intersection points
for i=1:i_max
    if i<=i_max/2
        airfoil_proj(i)=y_upper(r*cos(Delta_theta*(i-1)));  % Projected coordinate
        plot([x_circle_plot(i) x_circle_plot(i)],[y_circle_plot(i) airfoil_proj(i)])
    elseif i>=i_max/2
        airfoil_proj(i)=y_lower(r*cos(Delta_theta*(i-1)));  % Projected coordinate
        plot([x_circle_plot(i) x_circle_plot(i)],[y_circle_plot(i) airfoil_proj(i)])
    end
end

%% Plotting Far field and its lines

for i=1:i_max
    x_circleR_plot(i)=R*cos(Delta_theta*(i-1));  % x coordinates of circle
    y_circleR_plot(i)=R*sin(Delta_theta*(i-1));  % y coordinates of circle

    plot([0 x_circleR_plot(i)],[0 y_circleR_plot(i)])
end

 % Plotting the far field circle
plot(x_circleR_plot,y_circleR_plot)

%% Discretizing the Domain

% Discretizing the domain in delta x and delta y
Delta_x=(x_circleR_plot-x_circle_plot)/(j_max-1);   % The discretized change in the x direction
Delta_y=(y_circleR_plot-airfoil_proj)/(j_max-1);    % The discretized change in the y direction

%% Plotting the whole discretized domain

figure(2)
plot(x_airfoil_coord,y_upper(x_airfoil_coord),x_airfoil_coord,y_lower(x_airfoil_coord))
hold on
plot(x_circleR_plot,y_circleR_plot)
hold on
for i=1:i_max
    plot([x_circle_plot(i) x_circleR_plot(i)],[airfoil_proj(i) y_circleR_plot(i)])
end
hold on

for j=1:j_max-2
    for i=1:i_max
        x_coords(j,i)=x_circle_plot(i)+Delta_x(i)*j;
        y_coords(j,i)=airfoil_proj(i)+Delta_y(i)*j;
    end
    plot(x_coords(j,:),y_coords(j,:))
end

%% Transforming into the computational domain

Delta_eta1=(eta1_max-eta1_min)/(i_max-1);
Delta_eta2=(eta2_max-eta2_min)/(j_max-1);











