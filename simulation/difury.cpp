#include<difury.h>
#include<cmath>
#include <iostream>

void difury(double t, double  *y, double  *dy){
    double mu = 3.986004418E+14;
    double Rz = 6378136;
    double w_e = 7.292115E-05;
    double C20 = -1082.62575E-06;
    double t_G0 = (9*3600+18*60+10.5009+3*3600);
    double t_G = t_G0 + w_e*(t - 3*3600);
    double Ax = -.186264514923E-05;
    double Ay = -.186264514923E-05;
    double Az = -.186264514923E-05;
    double JX = Ax*cos(t_G)-Ay*sin(t_G);
    double JY = Ax*sin(t_G)+Ay*cos(t_G);
    double JZ = Az;
    double r = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
    mu = mu/(r*r*r);
    Rz = Rz/r;
    dy[0] = y[3];
    dy[1] = y[4];
    dy[2] = y[5];
    dy[3] = -mu*y[0] + 3.0/2*C20*mu*y[0]*Rz*Rz*(1-5*pow(y[2]/r,2)) + JX;
    dy[4] = -mu*y[1] + 3.0/2*C20*mu*y[1]*Rz*Rz*(1-5*pow(y[2]/r,2)) + JY;
    dy[5] = -mu*y[2] + 3.0/2*C20*mu*y[2]*Rz*Rz*(3-5*pow(y[2]/r,2)) + JZ;
}
