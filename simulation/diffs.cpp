#include <cmath>
#include <include\libglnsvpos\diffs.h>

// 22 20  2 10 13 45  0.0 -.134212896228E-04 -.181898940355E-11  .495000000000E+05
//      .311761962891E+04 -.179497814178E+01 -.186264514923E-08  .000000000000E+00
//      .158781806641E+05 -.202221393585E+01 -.186264514923E-08 -.300000000000E+01
//      .196852387695E+05  .192774677277E+01 -.186264514923E-08  .000000000000E+00

void diffs(double t, double *koord, double *dif)
{
    double mu = 3.986004418E+14; // конствнтва гравитационного поля Земли
    double R_z = 6378136; // экваториальный радиус Земли
    double omega_e = 7.292115E-05; // угловая скорость вращения Земли
    double C_20 = -1082.62575E-06;
    double t_g0 = (9*3600 + 18*60 + 10.5009);
    double t_g = t_g0 + omega_e*(t - 3*3600);

    // Ускорения, принимаем постоянными на всем интервале расчета
    double Ax = -.186264514923E-05;
    double Ay = -.186264514923E-05;
    double Az = -.186264514923E-05;
    double Jsum_x = Ax*cos(t_g) - Ay*sin(t_g);
    double Jsum_y = Ax*sin(t_g) + Ay*cos(t_g);
    double Jsum_z = Az;

    dif[0] = koord[3];
    dif[1] = koord[4];
    dif[2] = koord[5];
    double r = sqrt(koord[0]*koord[0] + koord[1]*koord[1] + koord[2]*koord[2]);
    double mu_ = mu/(r*r);
    double x_ = koord[0]/r;
    double y_ = koord[1]/r;
    double z_ = koord[2]/r;
    double ro = R_z/r;
    dif[3] = -mu_*x_ + 3/2*C_20*mu_*x_*ro*ro*(1-5*z_*z_) + Jsum_x;
    dif[4] = -mu_*y_ + 3/2*C_20*mu_*y_*ro*ro*(1-5*z_*z_) + Jsum_y;
    dif[5] = -mu_*z_ + 3/2*C_20*mu_*z_*ro*ro*(3-5*z_*z_) + Jsum_z;
}
