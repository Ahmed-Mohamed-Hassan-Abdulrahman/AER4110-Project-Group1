% O grid of airfoil using transformations
% Project made by Ahmed Mohamed Hassan Abdulrahman

clc
clearvars
close all

%% Initialization (Inputs)

% The Free Stream Velocity

Vinf=100;

% Choose the maximum mesh size:

i_max=50; j_max=50;

% Enter the Joukowski Airfoil Parameters

c=1;            % Chord
C_max_c=0.04;     % Maximum Camber/Chord Percentage
t_max_c=0.05;     % Maximum Thickness/Chord Percentage
AoA=4*pi/180;          % Angle of Attack of flow Percentage

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

x_circle_plot=zeros(1,i_max);
y_circle_plot=zeros(1,i_max);

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
axis equal

%% Projecting the intersection points

airfoil_proj=zeros(1,i_max);

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

x_circleR_plot=zeros(1,i_max);
y_circleR_plot=zeros(1,i_max);

for i=1:i_max
    x_circleR_plot(i)=R*cos(Delta_theta*(i-1));  % x coordinates of farfield circle
    y_circleR_plot(i)=R*sin(Delta_theta*(i-1));  % y coordinates of farfield circle

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
axis equal

plot(x_circleR_plot,y_circleR_plot)
for i=1:i_max
    plot([x_circle_plot(i) x_circleR_plot(i)],[airfoil_proj(i) y_circleR_plot(i)])
end

% initialiazing variables used in plotting
x_coords=zeros(j_max,i_max);
y_coords=zeros(j_max,i_max);

% Setting the coordinates of the airfoil and the farfield into the
% variables

x_coords(1,:)=x_circle_plot;
y_coords(1,:)=airfoil_proj;
x_coords(j_max,:)=x_circleR_plot;
y_coords(j_max,:)=y_circleR_plot;

% plotting the lines

for j=2:j_max-1
    for i=1:i_max
        x_coords(j,i)=x_circle_plot(i)+Delta_x(i)*(j-1);
        y_coords(j,i)=airfoil_proj(i)+Delta_y(i)*(j-1);
    end
    plot(x_coords(j,:),y_coords(j,:))
end

%% Transforming into the computational domain

Delta_eta1=(eta1_max-eta1_min)/(i_max-1);
Delta_eta2=(eta2_max-eta2_min)/(j_max-1);

% initializing the eta domain

eta1_coords=zeros(j_max,i_max);
eta2_coords=zeros(j_max,i_max);

% Getting the eta1 and eta2 coordinates

for j=1:j_max
    for i=1:i_max
        eta1_coords(j,i)=Delta_eta1*(i-1);
        eta2_coords(j,i)=Delta_eta2*(j-1);
    end
end

% Calculating computational domain derivatives

% Initialization
x_eta1=zeros(j_max,i_max); y_eta1=zeros(j_max,i_max);
x_eta2=zeros(j_max,i_max); y_eta2=zeros(j_max,i_max);

x_eta1(:,1)=(-1*x_coords(:,i_max-1)+x_coords(:,2))./(2*Delta_eta1);
x_eta2(1,:)=(-3*x_coords(1,:)+4*x_coords(2,:)-x_coords(3,:))./(2*Delta_eta2);
y_eta1(:,1)=(-1*y_coords(:,i_max-1)+y_coords(:,2))./(2*Delta_eta1);
y_eta2(1,:)=(-3*y_coords(1,:)+4*y_coords(2,:)-y_coords(3,:))./(2*Delta_eta2);

x_eta1(:,i_max)=(-1*x_coords(:,i_max-1)+x_coords(:,2))./(2*Delta_eta1);
x_eta2(j_max,:)=(3*x_coords(j_max,:)-4*x_coords(j_max-1,:)+x_coords(j_max-2,:))./(2*Delta_eta2);
y_eta1(:,i_max)=(-1*y_coords(:,i_max-1)+y_coords(:,2))./(2*Delta_eta1);
y_eta2(j_max,:)=(3*y_coords(j_max,:)-4*y_coords(j_max-1,:)+y_coords(j_max-2,:))./(2*Delta_eta2);

for i=2:i_max-1
    x_eta1(:,i)=(-1*x_coords(:,i-1)+x_coords(:,i+1))/(2*Delta_eta1);
    y_eta1(:,i)=(-1*y_coords(:,i-1)+y_coords(:,i+1))/(2*Delta_eta1);
end

for j=2:j_max-1
    x_eta2(j,:)=(-1*x_coords(j-1,:)+x_coords(j+1,:))/(2*Delta_eta2);
    y_eta2(j,:)=(-1*y_coords(j-1,:)+y_coords(j+1,:))/(2*Delta_eta2);
end

J=x_eta1.*y_eta2-x_eta2.*y_eta1;
c11=(x_eta2.^2+y_eta2.^2)./J;
c12=-1*(x_eta1.*x_eta2+y_eta1.*y_eta2)./J;
c22=(x_eta1.^2+y_eta1.^1)./J;


%% Calculating Boundary Conditions

% at j=jmax, psi=0
% Calculting velocity components at the farfield

uinf=Vinf*cos(AoA);
vinf=Vinf*cos(AoA);

% Calculating delta psi at the farfield

Delta_psi=uinf*Delta_y-vinf*Delta_x;












