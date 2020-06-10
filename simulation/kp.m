clc
clear

N4 = floor((2020-1996)/4)+1;
Nt=365*(2020-1996-4*(N4-1))+31+25+1;
% Эфемериды спутника ГЛОНАСС №22
X0 = 3117619.63;
Y0 = 15878180.66;
Z0 = 19685238.77;

VX = -1794.97814;
VY = -2022.21934;
VZ = 1927.74677;

AX = -0.0000019;
AY = -0.0000019;
AZ = -0.0000019;

tau = 13421.3;
gamma = -0.0018;

% Дата и время 2020/02/25 13:45:18
Toe = 13*60*60+45*60+18;
Ti = 12*60*60;

% Время по Гринвичу
Omega = 7.2921151467e-5;
TG0 = 16*60*60+45*60+18;
TGr = TG0+Omega*(Toe-3*3600);

% Геоцетрическая система ускорений
Xate = X0*cos(TGr) - Y0*sin(TGr);
Yate = X0*sin(TGr) + Y0*cos(TGr);
Zate = Z0;

% Пересчет координат
VXate = VX*cos(TGr)-VY*sin(TGr)-Omega*Yate;
VYate = VX*sin(TGr)+VY*cos(TGr)+Omega*Xate;
VZate = VZ;

x = AX*cos(TGr)-AY*sin(TGr)-Omega*Yate;
y = VX*sin(TGr)+VY*cos(TGr)+Omega*Xate;
z = AZ;

% const
J02 = 1082625.75 * 10^-9;
ae = 6378136;
GE = 398600441.8*10^6;

r = sqrt(Xate^2 + Yate^2 + Zate^2);
GE1 = GE/r^2;

xs = Xate/r;
ys = Yate/r;
zs = Zate/r;
p = ae/r;

% Система ур-й движения спутника
dx0 = VXate;
dy0 = VYate;
dz0 = VYate;

dVXate = -GE1*xs - 1.5*J02*GE1*xs*(p^2)*(1-5*zs^2)+ AX;
dVYate = -GE1*ys - 1.5*J02*GE1*ys*(p^2)*(1-5*zs^2)+ AY;
dVZate = -GE1*zs - 1.5*J02*GE1*zs*(p^2)*(3-5*zs^2)+ AZ;

