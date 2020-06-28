clear all;
close all;

dt = 0.1;                                    
T_eph = 49500 + 18 + 3*3600;
X = .311761962891E+07;                     
Y = .158781806641E+08;
Z = .196852387695E+08;
Vx = -.179497814178E+04;
Vy = -.202221393585E+04;
Vz = .192774677277E+04;
Ax = -.186264514923E-05;
Ay = -.186264514923E-05;
Az = -.186264514923E-05;

coordpotr = [55.756657, 37.703288 190];           %55°45'24.0"N 37°42'11.8"E
n = fix(12*3600/dt);                          
n_eph = fix((T_eph-(12+3)*3600)/dt);     
t12h = (12+3)*60*60 + dt.*(1:1:n);          
coordPZ90 = zeros(n,3);                        
coordECI = zeros(n,6);
w_e = 7.292115E-05;
t_G0 = 9*3600+18*60+10.5009;           %9:18:10.5009; 
t_G = t_G0 + w_e*(T_eph - 3*3600);
Yn = zeros(6,1);

Yn(1) = X*cos(t_G) - Y*sin(t_G); 
Yn(2) = X*sin(t_G) + Y*cos(t_G);
Yn(3) = Z;
Yn(4) = Vx*cos(t_G) - Vy*sin(t_G) - w_e*Yn(2);
Yn(5) = Vx*sin(t_G) + Vy*cos(t_G) + w_e*Yn(1);
Yn(6) = Vz;

coordECI(n_eph,:) = Yn;
% [t, Yn1] = ode45('difury', T_eph:-dt:t_12h(1), Yn);
% coordECI(1:n_eph,:) = Yn1(end:-1:1,1:3);
% 
% [t, Yn1] = ode45('difury', T_eph:dt:t_12h(n), Yn);
% coordECI(n_eph:end,:) = Yn1(1:end,1:3);

for i = n_eph+1:1:length(t12h)
    K1 = difury(t12h(i),coordECI(i-1,:));
    K2 = difury(t12h(i)+0.5*dt, coordECI(i-1,:)+0.5*dt*K1');
    K3 = difury(t12h(i)+0.5*dt, coordECI(i-1,:)+0.5*dt*K2');
    K4 = difury(t12h(i)+dt, coordECI(i-1,:)+dt*K3');
    coordECI(i,:) = coordECI(i-1,:) + dt*(K1' + 2*K2'+2*K3'+K4')/6;
end
for i = n_eph-1:-1:1
    K1 = difury(t12h(i),coordECI(i+1,:));
    K2 = difury(t12h(i)-0.5*dt, coordECI(i+1,:)-0.5*dt*K1');
    K3 = difury(t12h(i)-0.5*dt, coordECI(i+1,:)-0.5*dt*K2');
    K4 = difury(t12h(i)-dt, coordECI(i+1,:)-dt*K3');
    coordECI(i,:) = coordECI(i+1,:) - dt*(K1'+2*K2'+2*K3'+K4')/6;
end

t_G = t_G0 + w_e*(t12h' - 3*3600);
coordPZ90(:,1) = coordECI(:,1).*cos(t_G) + coordECI(:,2).*sin(t_G);
coordPZ90(:,2) = -coordECI(:,1).*sin(t_G) + coordECI(:,2).*cos(t_G);
coordPZ90(:,3) = coordECI(:,3);

X = zeros(n,1);
Y = zeros(n,1);
Z = zeros(n,1);
r = zeros(n,1);
teta = zeros(n,1);
phi = zeros(n,1);
parfor i = 1:length(coordPZ90(:,1))
    [X(i), Y(i), Z(i)] = ecef2enu(coordPZ90(i,1), coordPZ90(i,2), coordPZ90(i,3), coordpotr(1), coordpotr(2), coordpotr(3), wgs84Ellipsoid);
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(Z(i)/r(i));
        phi(i) = -atan2(Y(i),X(i))+pi/2;
    else
        teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end

[X_sf, Y_sf, Z_sf] = sphere(25);
cutoff_grads = pi/180.*(1:359)';
cutoff_angle = 85.*ones(359,1);
R_z = 6371;

% ECI & ПЗ-90
plot3(coordPZ90(:,1)/1000, coordPZ90(:,2)/1000, coordPZ90(:,3)/1000,'r');
hold on
plot3(coordECI(:,1)/1000, coordECI(:,2)/1000, coordECI(:,3)/1000,'b');
grid on
title('Траектория движения спутника ГЛОНАСС №22 в СК ПЗ-90 и ECI');
xlabel('OX, км');
ylabel('OY, км');
zlabel('OZ, км');
axis('square');
axis('equal');
surf(X_sf*R_z, Y_sf*R_z, Z_sf*R_z);
hold off

% SkyView
figure;
axes = polaraxes;
hold on
polarplot(axes,phi,teta*180/pi,'r')
polarplot(axes,cutoff_grads,cutoff_angle,'b')
hold off
axes.ThetaDir = 'clockwise';
axes.ThetaZeroLocation = 'top';
title('SkyView спутника ГЛОНАСС №22')

