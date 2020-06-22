% эфемериды, полученные на первом этапе
% 22 20  2 10 13 45  0.0 -.134212896228E-04 -.181898940355E-11  .495000000000E+05
%      .311761962891E+04 -.179497814178E+01 -.186264514923E-08  .000000000000E+00
%      .158781806641E+05 -.202221393585E+01 -.186264514923E-08 -.300000000000E+01
%      .196852387695E+05  .192774677277E+01 -.186264514923E-08  .000000000000E+00

clear all;
close all;

tic;
% Константы
dt = 1.0E-01;               % шаг изменения времени, с
w_e = 7.292115E-05;      % угловая скорость вращения Земли
mu = 3.986004418E+14;   % константа гравитационного поля Земли
R_z = 6371;              % радиус Земли, км
% Эфемериды
T_eph = 49500 + 18 + 3*3600;  % время задания эфемерид
% Координаты
X = .311761962891E+07;                     
Y = .158781806641E+08;
Z = .196852387695E+08;
% Составляющие скорости
Vx = -.179497814178E+04;
Vy = -.202221393585E+04;
Vz = .192774677277E+04;
% Составляющие ускорения от прочих небесных тел (Луна и Солнце) 
Ax = -.186264514923E-05;
Ay = -.186264514923E-05;
Az = -.186264514923E-05;

coordpotr = [55.756657, 37.703288 190];           % координаты потребителя 55°45'24.0"N 37°42'11.8"E
n = fix(12*3600/dt);                          % Количество отсчетов за 12 часовой интервал расчета
n_eph = fix((T_eph-(12+3)*3600)/dt);     % номер отсчета с эфемеридными данными
t_12h = (12+3)*60*60 + dt.*(1:1:n);          % Вектор отсчетов времени
coordPZ90 = zeros(n,3);                        % Вектора координат с нулевым заполнением
coordECI = zeros(n,3);
t_G0 = 9*3600+18*60+10.5009;           %9:18:10.5009; истинное звездное время в гринвичскую полночь даты задания tэ
t_G = t_G0 + w_e*(T_eph - 3*3600);

% Пересчет координат из ПЗ-90 в ECI
Xa = X*cos(t_G) - Y*sin(t_G); 
Ya = X*sin(t_G) + Y*cos(t_G);
Za = Z;
Vxa = Vx*cos(t_G) - Vy*sin(t_G) - w_e*Ya;
Vya = Vx*sin(t_G) + Vy*cos(t_G) + w_e*Xa;
Vza = Vz;
Yn = [Xa Ya Za Vxa Vya Vza];

% Интегрирование методом Рунге-Кутты
% Расчет координат от начала диапазона до момента времени эфемерид
[t, Yn1] = ode45('difury', T_eph:-dt:t_12h(1), Yn);
coordECI(1:n_eph,:) = Yn1(end:-1:1,1:3);
% Расчет координат от момента времени эфемерид до конца исследуемого
% диапазона
[t, Yn1] = ode45('difury', T_eph:dt:t_12h(n), Yn);
coordECI(n_eph:end,:) = Yn1(1:end,1:3);


% Пересчет полученных координат из ECI в ПЗ-90
for i = 1:n
    t_G = t_G0 + w_e*(t_12h(i)- 3*3600);
    cos_t_G = cos(t_G);
    sin_t_G = sin(t_G);
    coordPZ90(i,1) = coordECI(i,1)*cos_t_G + coordECI(i,2)*sin_t_G;
    coordPZ90(i,2) = -coordECI(i,1)*sin_t_G + coordECI(i,2)*cos_t_G;
    coordPZ90(i,3) = coordECI(i,3);
end
% Пересчет координат из ПЗ-90 в WGS84 (для получения SkyView)
ppb = 1e-9;
mas = 1e-3/206264.8; % [рад]
M_WGS84 = [-3*ppb -353*mas -4*mas;
            353*mas -3*ppb 19*mas;
            4*mas -19*mas -3*ppb];

coordWGS84 = coordPZ90.'; % Переход к вектору-столбцу
for i = 1:length(coordWGS84(1,:))
    coordWGS84(:,i) =  coordWGS84(:,i) + M_WGS84 * coordWGS84(:,i) + [0.07; -0; -0.77];
end
coordWGS84 = coordWGS84.'; % Переход к вектору-строке
% Пересчет координат из WGS84 в SkyView
X = zeros(n,1);
Y = zeros(n,1);
Z = zeros(n,1);
r = zeros(n,1);
teta = zeros(n,1);
phi = zeros(n,1);
for i = 1:length(coordWGS84(:,1))
    [X(i), Y(i), Z(i)] = ecef2enu(coordWGS84(i,1),coordWGS84(i,2),coordWGS84(i,3),coordpotr(1),coordpotr(2),coordpotr(3),wgs84Ellipsoid);
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(Z(i)/r(i));
        if X(i) > 0
            phi(i) = -atan(Y(i)/X(i))+pi/2;
        elseif (X(i)<0)&&(Y(i)>0)
            phi(i) = -atan(Y(i)/X(i))+3*pi/2;
        elseif (X(i)<0)&&(Y(i)<0)
            phi(i) = -atan(Y(i)/X(i))-pi/2;
        end
    else
        teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end
% Построение графиков

% Расчет сферы, изображения земли
[X_sf, Y_sf, Z_sf] = sphere(25);
% Окружность, для изображения угла отсечки
cutoff_grads = pi/180.*(1:359)';
cutoff_angle = 85.*ones(359,1);

% График координат в ПЗ-90
plot3(coordPZ90(:,1)/1000, coordPZ90(:,2)/1000, coordPZ90(:,3)/1000);
hold on
grid on
title('Положение спутника в СК ПЗ-90');
xlabel('OX, км');
ylabel('OY, км');
zlabel('OZ, км');
axis('square');
axis('equal');
surf(X_sf*R_z, Y_sf*R_z, Z_sf*R_z);
hold off

% График координат в ECI
figure;
plot3(coordECI(:,1)/1000, coordECI(:,2)/1000, coordECI(:,3)/1000);
hold on
grid on
title('Положение спутника в СК ECI');
xlabel('OX, км');
ylabel('OY, км');
zlabel('OZ, км');
axis('square');
axis('equal');
surf(X_sf*R_z, Y_sf*R_z, Z_sf*R_z);
hold off

% SkyView спутника (угол-место) относительно корпуса Е
figure;
axes = polaraxes;
hold on
polarplot(axes,phi,teta*180/pi,'r')
polarplot(axes,cutoff_grads,cutoff_angle,'b')
hold off
axes.ThetaDir = 'clockwise';
axes.ThetaZeroLocation = 'top';
title('SkyView ГЛОНАСС №22')

toc;


