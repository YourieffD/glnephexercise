clc
clear
% Точный алгоритм пересчёта координат и составляющих вектора 
% скорости центра масс НКА на заданный момент времени шкалы МВД
% Расчет времени формата ГЛОНАСС
N4=floor((2020-1996)/4)+1;
Nt=365*(2020-1996-4*(N4-1))+31+25+1;
tb=13*60*60 + 45*60 + 18 + 10800;
Toe = (12+3)*60*60;
Tof = (24+3)*60*60;
Ts = 1;
ti = Toe:Ts:Tof;

% Эфемериды
X = 3117619.63; 
Y = 15878180.66;
Z = 19685238.77;
VX = -1794.97814;
VY = -2022.21934;
VZ = 1927.74677;
AX = -0.0000019;
AY = -0.0000019;
AZ = -0.0000019;

% Константы
ae = 6378136; 
omega_z = 7.2921151467e-5;

% Юлианская дата
JD0 = 1461*(N4 - 1) + Nt + 2450008.5 - (Nt -3)/25; 
% Среднее время по Гринвичу
T_delta = (JD0-2451545)/36525;

% Угол поворота Земли
ERA = 2*pi*(0.7790572732640 + 1.00273781191135448*(JD0-2451545));

GMST = ERA+0.0000000703270726+0.0223603658710194*T_delta+...
    +0.0000067465784654*T_delta^2-0.0000000000021332*T_delta^3+...
    - 0.0000000001452308*T_delta^4-0.0000000000001784*T_delta^5;

% Перевод в ГСК
S = GMST + omega_z*(tb - 3*60*60);
Xate = X*cos(S)-Y*sin(S);
Yate = X*sin(S)+Y*cos(S);
Zate = Z;

Vxate = VX*cos(S)-VY*sin(S)-omega_z*Yate;
Vyate = VX*sin(S)+VY*cos(S)+omega_z*Xate;
Vzate = VZ;

Axte = AX*cos(S)-AY*sin(S);
Ayte = AX*sin(S)+AY*cos(S);
Azte = AZ;


% Меод Рунге-Кутты
res0 = [Xate Yate Zate Vxate Vyate Vzate];
[t, Y] = ode45('diffs', tb:-Ts:ti(1), res0);
res1 = Y(end:-1:2,:);
t1 = t(end:-1:2,:);
[t, Y] = ode45('diffs', tb:Ts:ti(end), res0);
res1 = [res1; Y];
t1 = [t1;t];

% Учет ускорений
tau1 = t1 - tb;
AXTE = AX*(tau1.^2)/2;
AYTE = AY*(tau1.^2)/2;
AZTE = AZ*(tau1.^2)/2;

delta_VX = AX*tau1;
delta_VY = AY*tau1;
delta_VZ = AZ*tau1;

delta_A = [AXTE AYTE AZTE delta_VX delta_VY delta_VZ];


res1 = res1 + delta_A;
% Пересчет координат в ПЗ-90
S = GMST + omega_z*(t1 - 3*60*60);
pz90(:,1) = res1(:,1).*cos(S) + res1(:,2).*sin(S);
pz90(:,2) = -res1(:,1).*sin(S) + res1(:,2).*cos(S);
pz90(:,3) = res1(:,3);

% Корпус Е
N = 55.756687916667*pi/180;
E = 37.703077194444*pi/180;
H = 500;
cord_E = [N E H];


% Скайплот
for i = 1:length(pz90(:,1))
[X(i) Y(i) Z(i)] = ecef2enu(pz90(i,1),pz90(i,2),pz90(i,3),N,E,H,wgs84Ellipsoid,'radians');
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(Z(i)/r(i));

        if X(i) > 0
            phi(i) = -atan(Y(i)/X(i))+pi/2;
        elseif (X(i)<0)&&(Y(i)>0)
            phi(i) = -atan(Y(i)/X(i))+3*pi/2;
        elseif (X(i)<0)&&(res2(i)<0)
            phi(i) = -atan(res2(i)/X(i))-pi/2;
        end
    else teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end

figure(1)
[Xx,Yy,Zz]=sphere(50);
Rz=6371000;
surf(Rz*Xx,Rz*Yy,Rz*Zz)
hold on
grid on
plot3(res1(:,1), res1(:,2), res1(:,3), 'b')
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
surf(Rz*Xx,Rz*Yy,Rz*Zz)
grid on
plot3(pz90(:,1),pz90(:,2),pz90(:,3),'r')
title({'Траектория движения спутника ГЛОНАСС №22'})
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
hold off


figure (2)
pax = polaraxes;
polarplot(pax,phi,teta*180/pi,'r')
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
title('SkyView спутника ГЛОНАСС №22')
th = hours(t1./3600-3);

figure(3);
grid on
hold on
plot(th,(-teta*180/pi+90),'DurationTickFormat','hh:mm:ss')
title('Угол места')
xlabel('Время в МДВ')
ylabel('Угол места спутника ГЛОНАСС №22, град')





