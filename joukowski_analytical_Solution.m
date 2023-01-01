clc
clear vars
close all
%% Givens:
V_inf=100 ; %%%Free stream velocity in m/s
chord=1 ;   %%%chord Length in meter
maxthickness=10/100; %%maximum thickness to chord ratio
maxcamber=5/100;    %%maximum camber to chord ratio
alpha= 8;           %% angle of attack on degrees
%% Circle Parameters:
b=chord/4 ;    %% intrsection of the shifted circle with the positve x-axis
e=maxthickness/1.3;
beta=(2*maxcamber)/chord; 
beta_degree=beta*(180/pi);
a=(b*(1+e))/(cos(beta)); %%radius of the shifted circle in Z_1 plane
x_0=-b*e;   %% x corrdinate of the shifted circle
y_0=a*beta; %% y  corrdinate of the shifted circle
%% jouwkoski Airfoil Calculations:
theta=linspace(1,360-1,2*180); %% the angle theta for any point on the airfoil in Z_1 plane and on circle in Z dash plane
x_1=2*b*cosd(theta);           
y_1=2*b*e*(1-cosd(theta)).*sind(theta)+2*b*beta*sind(theta).^2;
r=b*(1+e*(1-cosd(theta))+beta*sind(theta));
%% jouwkoski airfoil plot:
figure
plot(x_1,y_1,'*');
title('The General Joukowski airfoil ')
xlabel('airfoil X axis')
ylabel('airfoil Y axis')
axis([-0.5 0.5 -0.2 0.2])
%% stream lines plot:
r_dash1=linspace(a,3*a,360);
theta_dash1=linspace(0,2*pi,360);
alpha_free_stream=pi/180*5; %% freestream angle in Rad
Gamma=4*pi*V_inf*a*sin(alpha_free_stream+beta); %% vortex Strength Caculation
[Theta_dash1,R_dash]= meshgrid(theta_dash1,r_dash1);
Z_dash_streamline=R_dash.*exp(1i*Theta_dash1)+x_0+1i*y_0;
Z1=Z_dash_streamline+b^2./Z_dash_streamline;
W = V_inf .*exp(-1i.*alpha_free_stream(length(alpha_free_stream))).*Z_dash_streamline+V_inf*a^2*exp( 1i *alpha_free_stream(length(alpha_free_stream)))./Z_dash_streamline+1i*log(Z_dash_streamline)*Gamma/2/pi;
psi=imag(W); %% Lines of stream line
psi1=psi+b^2./psi; %% cnvorting the stream lines from Z_dash plane to Z_1_plane
figure
hold on
axis equal
plot(Z1(1,:),'b');
fill(real(Z1(1,:)),imag(Z1(1,:)),'r')
box on;
contour(real(Z1),imag(Z1),psi1,'levelstep',2,'color','b');
xlim([-2/3*chord,2/3*chord]);
ylim([-a,a]);
title('Stream lines on the joukowski Airfoil')
%% calculations in the z-plane:
x=r.*cosd(theta);
y=r.*sind(theta);
%% calculations in z_dash plane:
x_dash=x-x_0;  %% calculation of x points on the z_dash plane where the axis on the center of the circle
y_dash=y-y_0;  %% calculation of y points on the z_dash plane where the axis on the center of the circle
z_dash=x_dash+1i*y_dash; %% making the points on the real- imaginary plane
r_dash=abs(z_dash);
theta_dash=angle(z_dash);
%% calculation of velocities in Z_dash plane:
v_r_dash= V_inf*(1-(a^2./r_dash.^2)).*cos(theta_dash-alpha*(pi/180));
v_theta_dash=-V_inf*((1+(a^2./r_dash.^2)).* sin(theta_dash-alpha*(pi/180))+2*(a./r_dash)*sin((alpha*pi/180)+beta));
%% calculation of velocities in z_1 plane:
V_1=(v_theta_dash.^2./(1-2*(b^2./r.^2).*cosd(2*theta)+(b^4./r.^4))).^0.5;
%% plot the velocity over the airfoil
figure
plot(x_1,V_1)
title('velocit along the joukowskki airfoil')
xlabel('x_1')
ylabel('V_1')
%% plot the flow velocity over the airfoil:
figure 
plot(x_1,V_1./V_inf)
title('V_1/V_infinity along the joukowskki airfoil')
xlabel('x_1')
ylabel('V_1/V_infinity')
%% C_p calculation and plot:
C_p=1-(V_1./V_inf).^2;
figure 
plot(x_1,C_p)
title('C_p over the airfoil')
xlabel('x_1')
ylabel('C_p')
axis([-0.5 0.5 -6 1])
%% C_l calculations:
angle_attack= -5:10 ;
C_l=2*pi*(1+e).*sind(angle_attack+beta_degree);
%% c_l plot with angle of attacks:
figure
plot(angle_attack,C_l)
title('C_l Vs alpha(-5 to 10)')
xlabel('alpha')
ylabel('C_L')
%% Cm calculations:
angel_attack_rad=(pi/180)*angle_attack;
for ii=1:length(angle_attack)
    
    v_theta_dash=-V_inf*((1+(a^2./r_dash.^2)).* sin(theta_dash-angel_attack_rad(ii))+2*(a./r_dash).*sin((angel_attack_rad(ii))+beta));
    V1=(v_theta_dash.^2./(1-2*(b^2./r.^2).*cos(2*theta)+(b^4./r.^4))).^0.5;
    Cp = 1-(V1/V_inf).^2;
    C_m_LE(ii)=1/chord*(trapz(-Cp.*(x_1+chord/2),x_1+chord/2) + trapz(-Cp.*y_1,y_1) );
end
%% cm plot:
figure
plot(angle_attack,C_m_LE)
title('C_m Vs alpha(-5 to 10)')
xlabel('alpha')
ylabel('C_m')
