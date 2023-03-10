% O grid of airfoil using transformations
% Project made by Ahmed Mohamed Hassan Abdulrahman

clc
clearvars
close all

%% Initialization (Inputs)

% The Free Stream Velocity

Vinf=100;

% Choose the maximum mesh size:

i_max=100; j_max=100;

% Enter the Joukowski Airfoil Parameters

c=1;                    % Chord
C_max_c=0.04;           % Maximum Camber/Chord Percentage
t_max_c=0.05;            % Maximum Thickness/Chord Percentage
AoA=4*pi/180;           % Angle of Attack of flow Percentage

% Drawing Parameters

airfoil_points=500;     % Number of points on the airfoil; More points = Better accuracy
R=5*c;                  % Far field radius assumption
Contour_Detail=100;     % How fine the contour plots would be

% The Transformation Parameters

eta1_max=1;     eta1_min=0;
eta2_max=1;     eta2_min=0;

% Solution Limits

RMS_limit=1e-5;

%% Generating Airfoil Coordinates

% Joukowski airfoil parameters calculations

b=c/4;                  
r=2*b;                  % Radius of the circle
e=t_max_c/1.3;          % eccentricity of the circle
beta=2*C_max_c;
a=b*(1+e)/cos(beta);

% x coordinates of the points on the upper and lower surfaces of the airfoil:
x_airfoil_coord=linspace(-r,r,airfoil_points);

% y coordinates equation on the upper surface of the airfoil:
y_upper= @(x) 2*b*e*(1-x/2/b).*(sqrt(1-(x/2/b).^2))+2*b*beta*(1-(x/2/b).^2);
% y coordinates equation on the lower surface of the airfoil:
y_lower= @(x) 2*b*e*(1-x/2/b).*(-sqrt(1-(x/2/b).^2))+2*b*beta*(1-(x/2/b).^2);

% Plotting the Airfoil on its own

figure('Name','The Airfoil on its own')
plot(x_airfoil_coord,y_upper(x_airfoil_coord),x_airfoil_coord,y_lower(x_airfoil_coord))
title('The Airfoil on its own')

% plotting the Airfoil along with the discretizing circles
figure('Name','Projections over the Airfoil')
hold on
plot(x_airfoil_coord,y_upper(x_airfoil_coord),x_airfoil_coord,y_lower(x_airfoil_coord))
title('Projections over the Airfoil')

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
hold on

%% Discretizing the Domain

% Discretizing the domain in delta x and delta y
Delta_x=(x_circleR_plot-x_circle_plot)/(j_max-1);   % The discretized change in the x direction
Delta_y=(y_circleR_plot-airfoil_proj)/(j_max-1);    % The discretized change in the y direction

%% Plotting the whole discretized domain

figure('Name','The Whole Discretized Domain')
plot(x_airfoil_coord,y_upper(x_airfoil_coord),x_airfoil_coord,y_lower(x_airfoil_coord))
title('The Whole Discretized Domain')
hold on

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
c22=(x_eta1.^2+y_eta1.^2)./J;

%% Calculating psi at the zero condition

% Calculating velocity components

uinf=Vinf*cos(AoA);
vinf=Vinf*sin(AoA);

% Initialization

psi=zeros(j_max,i_max);
Delta_x_farfield=zeros(1,i_max-1);
Delta_y_farfield=zeros(1,i_max-1);
Delta_psi=zeros(1,i_max-1);

% Calculation of the zero iteration and boundary condition iteration

% assuming that psi is zero on the airfoil
for i=2:i_max
    Delta_x_farfield(i-1)=x_circleR_plot(i)-x_circleR_plot(i-1);
    Delta_y_farfield(i-1)=y_circleR_plot(i)-y_circleR_plot(i-1);
    Delta_psi(i-1)=uinf*Delta_y_farfield(i-1)-vinf*Delta_x_farfield(i-1);
    psi(j_max,i)=psi(j_max,i-1)+Delta_psi(i-1);
    psi(:,i)=linspace(psi(1,i),psi(j_max,i),j_max);
end

psi(:,1)=psi(:,i_max);
if C_max_c==0
    psi(1,:)=psi(2,1);
elseif C_max_c>0
    psi(1,:)=psi(2,i_max-ceil(i_max*C_max_c/2));
elseif C_max_c<0
    psi(1,:)=psi(2,1+floor(i_max*C_max_c/2));
end

%% Calculating the half coefficients used in iterations

%iph= i+0.5     %jph= j+0.5
%imh= i-0.5     %jmh= j-0.5
c11_ih=zeros(j_max,i_max);
c22_jh=zeros(j_max,i_max);

for i=1:i_max-1
    c11_ih(:,i+1)=(c11(:,i)+c11(:,i+1))/2; % from i=1.5 to i=i_max-0.5
end
c11_ih(:,1)=c11_ih(:,i_max); % value at i=0.5 is the same at i=i_max-0.5

% since our iterations will start from j=2
for j=1:j_max-1
    c22_jh(j,:)=(c22(j,:)+c22(j+1,:))/2; % from j=1.5 to j=j_max-0.5
end

%% Coefficients used in iterations

for j=2:j_max-1
    for i=1:i_max-1
        S_i_j(j-1,i)=(c11_ih(j,i+1)+c11_ih(j,i))/Delta_eta1^2 + ...
            (c22_jh(j,i)+c22_jh(j-1,i))/Delta_eta2^2;

        S_ip1_j(j-1,i)=c11_ih(j,i+1)/Delta_eta1^2;
        S_im1_j(j-1,i)=c11_ih(j,i)/Delta_eta1^2;

        S_i_jp1(j-1,i)=c22_jh(j,i)/Delta_eta2^2;
        S_i_jm1(j-1,i)=c22_jh(j-1,i)/Delta_eta2^2;

        S_ip1_jp1(j-1,i)=(c12(j,i+1)+c12(j+1,i))/(4*Delta_eta1*Delta_eta2);
        S_ip1_jm1(j-1,i)=(-c12(j,i+1)-c12(j-1,i))/(4*Delta_eta1*Delta_eta2);
        if i==1
            S_im1_jp1(j-1,i)=(-c12(j,i_max-1)-c12(j+1,i))/(4*Delta_eta1*Delta_eta2);
            S_im1_jm1(j-1,i)=(c12(j,i_max-1)+c12(j-1,i))/(4*Delta_eta1*Delta_eta2);
        else
        S_im1_jp1(j-1,i)=(-c12(j,i-1)-c12(j+1,i))/(4*Delta_eta1*Delta_eta2);
        S_im1_jm1(j-1,i)=(c12(j,i-1)+c12(j-1,i))/(4*Delta_eta1*Delta_eta2);
        end     
    end
end

%% Iterations

iteration_No=0;
psi_intitial=psi;
psi_iteration=psi_intitial;

RMS=RMS_limit*1e5;

while RMS>RMS_limit
    for j=2:j_max-1
        for i=1:i_max-1
            if i==1
                psi(j,i)=(psi_iteration(j,i+1)*S_ip1_j(j-1,i)+psi_iteration(j,i_max-1)*S_im1_j(j-1,i)+...
                psi_iteration(j+1,i+1)*S_ip1_jp1(j-1,i)+psi_iteration(j-1,i+1)*S_ip1_jm1(j-1,i)+...
                psi_iteration(j+1,i_max-1)*S_im1_jp1(j-1,i)+psi_iteration(j-1,i_max-1)*S_im1_jm1(j-1,i)+...
                psi_iteration(j+1,i)*S_i_jp1(j-1,i)+psi_iteration(j-1,i)*S_i_jm1(j-1,i))/S_i_j(j-1,i);
            else
                psi(j,i)=(psi_iteration(j,i+1)*S_ip1_j(j-1,i)+psi_iteration(j,i-1)*S_im1_j(j-1,i)+...
                psi_iteration(j+1,i+1)*S_ip1_jp1(j-1,i)+psi_iteration(j-1,i+1)*S_ip1_jm1(j-1,i)+...
                psi_iteration(j+1,i-1)*S_im1_jp1(j-1,i)+psi_iteration(j-1,i-1)*S_im1_jm1(j-1,i)+...
                psi_iteration(j+1,i)*S_i_jp1(j-1,i)+psi_iteration(j-1,i)*S_i_jm1(j-1,i))./S_i_j(j-1,i);
            end
        end
    end
    psi(:,i_max)=psi(:,1);

    if C_max_c==0
        psi(1,:)=psi(2,1);
    elseif C_max_c>0
        psi(1,:)=psi(2,i_max-1);
    elseif C_max_c<0
        psi(1,:)=psi(2,1+floor(i_max*C_max_c/2));
    end

    RMS=sqrt(sum(sum((psi-psi_iteration).^2))/((i_max-1)*(j_max-1)));

    psi_iteration=psi;

    iteration_No=iteration_No+1;
end

%% Velocity Calculation

% since at i=1 and i=i_max is the same point
% we calculate psi using central difference

psi_eta1(:,1)=(psi(:,2)-psi(:,i_max-1))/(2*Delta_eta1);
psi_eta1(:,i_max)=psi_eta1(:,1);

% calculating at j=1 using forward difference
psi_eta2(1,:)=(-3*psi(1,:)+4*psi(2,:)-1*psi(3,:))/(2*Delta_eta2);

% calculating at j=j_max using backward difference
psi_eta2(j_max,:)=(3*psi(j_max,:)-4*psi(j_max-1,:)+psi(j_max-2,:))/(2*Delta_eta2);

for j=2:j_max-1
    for i=1:i_max
        psi_eta2(j,i)=(psi(j+1,i)-psi(j-1,i))/(2*Delta_eta2);
    end
end
for j=2:j_max
    for i=2:i_max-1
        psi_eta1(j,i)=(psi(j,i+1)-psi(j,i-1))/(2*Delta_eta1);
    end
end

eta1_x=y_eta2./J;
eta1_y=-x_eta2./J;
eta2_x=-y_eta1./J;
eta2_y=x_eta1./J;

u=psi_eta1.*eta1_y+psi_eta2.*eta2_y;
v=-psi_eta1.*eta1_x-psi_eta2.*eta2_x;

%% Cp Calculation

V=sqrt(u.^2+v.^2);
[j_ind,i_ind]=find(V>=3*mean(mean(V)));
V(j_ind,i_ind)=(V(j_ind,i_ind-1)+V(j_ind,i_ind+1))/2;
V(j_ind:j_ind,i_ind-1:i_ind+1)=linspace(V(j_ind,i_ind-1),V(j_ind,i_ind+1),length(i_ind));

Cp=1-(V/Vinf).^2;

%% Analytical solution using Joukowski Airfoil

[V_analytical, Cp_analytical, x_coords_analytical]=Joukowski(Vinf,AoA,c,C_max_c,t_max_c,i_max/2);


%% Results Graphs and plots

% Numerical pressure and velocity distribution over the Airfoil

figure('Name', 'Numerical Pressure and Velocity Distribution over the Airfoil')
     tiledlayout(2,1);
        nexttile 
        hold on
            plot(x_circle_plot,V(1,:),'-' ,'LineWidth',1.5,'color','red')      
            plot(x_circle_plot,800*airfoil_proj,'LineWidth',0.5,'color','black')
            fill(x_circle_plot,800*airfoil_proj,'cyan')             
            grid on
            xlabel('$x$', 'interpreter', 'latex')
            ylabel('$V$', 'interpreter', 'latex')
            title('Velocity distribution', 'FontName','lm roman 9')
            xlim([-0.6 0.6])
        nexttile 
        hold on
            plot(x_circle_plot,Cp(1,:),'-' ,'LineWidth',1.5,'color','blue')
            plot(x_circle_plot,15*airfoil_proj+min(Cp_analytical),'LineWidth',0.5,'color','black')
            fill(x_circle_plot,15*airfoil_proj+min(Cp_analytical),'cyan') 
            grid on
            xlabel('$x$', 'interpreter', 'latex')
            ylabel('$C_p$', 'interpreter', 'latex')
            title('pressure coefficient distribution', 'FontName','lm roman 9')
            xlim([-0.6 0.6])

% Comparison Between the Joukowski Analytical and Numerical Solution

figure('Name', 'Comparison Between the Joukowski Analytical and Numerical Solution')
    tiledlayout(2,1);
        nexttile 
            hold on
            % velocity results
            plot(x_circle_plot,V(1,:),'-' ,'LineWidth',1.5,'color','red')
            plot(x_coords_analytical,V_analytical,'-' ,'LineWidth',1.5,'color','blue')
            % the Airfoil 
            plot(x_circle_plot,800*airfoil_proj,'LineWidth',0.5,'color','black')
            fill(x_circle_plot,800*airfoil_proj,'cyan') 
            % figure settings
            grid on
            xlabel('$x$', 'interpreter', 'latex')
            ylabel('$V$', 'interpreter', 'latex')
            legend('Numerical','Analytical','Airfoil')
            title('Velocity distribution', 'FontName','lm roman 9')
        nexttile 
            hold on
            % velocity results
            plot(x_circle_plot,Cp(1,:),'-' ,'LineWidth',1.5,'color','red')
            plot(x_coords_analytical,Cp_analytical,'-' ,'LineWidth',1.5,'color','blue')
            % the Airfoil 
            plot(x_circle_plot,20*airfoil_proj+min(Cp_analytical),'LineWidth',0.5,'color','black')
            fill(x_circle_plot,20*airfoil_proj+min(Cp_analytical),'cyan') 
            % figure settings
            grid on
            xlabel('$x$', 'interpreter', 'latex')
            ylabel('$C_p$', 'interpreter', 'latex')
            legend('Numerical','Analytical','Airfoil')
            title('Pressure coefficient distribution', 'FontName','lm roman 9')

% Stream Lines Over the Airfoil
%Zoomed out
figure('Name', 'Streamlines Over the Airfoil')
        hold on
            plot(x_circle_plot,airfoil_proj)
            fill(x_circle_plot,airfoil_proj,'cyan')
            contour(x_coords,y_coords,psi,linspace(min(min(psi)),max(max(psi)),Contour_Detail));
            axis equal
            axis off
            title('Stream Lines Over the Airfoil', 'FontName','lm roman 12')
%Zoomed in
figure('Name', 'Streamlines Over the Airfoil')
        hold on
            plot(x_circle_plot,airfoil_proj)
            fill(x_circle_plot,airfoil_proj,'cyan')
            contour(x_coords,y_coords,psi,linspace(min(min(psi)),max(max(psi)),Contour_Detail));
            axis([-1.5 1.5 -1 1 ])
            axis off
            title('Stream Lines Over the Airfoil', 'FontName','lm roman 12')

% Velocity and pressure Distribution Contour Around the Airfoil
% zoomed out
figure('Name', 'Velocity contours')
        hold on
            plot(x_circle_plot,airfoil_proj,'LineWidth',1.5)
            title('Velocity Contour Around the Airfoil')
            contourf(x_coords,y_coords,V,Contour_Detail,'edgecolor','none')
            axis off
            xlabel('x','FontSize',16)
            ylabel('y','FontSize',16)
            c = colorbar;
            c.Label.String = 'Velocity';
            c.Label.FontSize = 16;
            axis equal

figure('Name', 'Pressure Coefficient contours')
            hold on
            plot(x_circle_plot,airfoil_proj,'LineWidth',1.5)
            title('Pressure Contour Around the Airfoil')
            contourf(x_coords,y_coords,Cp,Contour_Detail,'edgecolor','none')
            axis off
            xlabel('x','FontSize',16)
            ylabel('y','FontSize',16)
            c = colorbar;
            c.Label.String = 'Coefficient of Pressure';
            c.Label.FontSize = 16;
            axis equal

% zoomed in
figure('Name', 'Velocity and Pressure Coefficient contours')
    tiledlayout(1,2)
        nexttile
        hold on
            plot(x_circle_plot,airfoil_proj,'LineWidth',1.5)
            title('Velocity Contour Around the Airfoil')
            contourf(x_coords,y_coords,V,Contour_Detail,'edgecolor','none')
            axis([-1.5 1.5 -1.5 1.5 ])
            axis off
            xlabel('x','FontSize',16)
            ylabel('y','FontSize',16)
            c = colorbar;
            c.Label.String = 'Velocity';
            c.Label.FontSize = 16;
        nexttile
        hold on
            plot(x_circle_plot,airfoil_proj,'LineWidth',1.5)
            title('Pressure Contour Around the Airfoil')
            contourf(x_coords,y_coords,Cp,Contour_Detail,'edgecolor','none')
            axis([-1.5 1.5 -1.5 1.5 ])
            axis off
            xlabel('x','FontSize',16)
            ylabel('y','FontSize',16)
            c = colorbar;
            c.Label.String = 'Coefficient of Pressure';
            c.Label.FontSize = 16;
