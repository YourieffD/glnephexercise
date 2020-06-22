#include <include\libglnsvpos\glnsvpos.h>
#include <include\libglnsvpos\rungekutta.h>

#include <iostream>
#include <cmath>
#include <ostream>

using namespace std;

// 22 20  2 10 13 45  0.0 -.134212896228E-04 -.181898940355E-11  .495000000000E+05
//      .311761962891E+04 -.179497814178E+01 -.186264514923E-08  .000000000000E+00
//      .158781806641E+05 -.202221393585E+01 -.186264514923E-08 -.300000000000E+01
//      .196852387695E+05  .192774677277E+01 -.186264514923E-08  .000000000000E+00

void koordinate_GLONASS(double **koord_n)
{
    double t_0 = 49500 + 18 + 3*3600;
    double t_nachalo = (12 + 3)*3600;
    double t_konez = t_nachalo + 12*3600;
    double del_t = 1E-01; //шаг по времени
    double omega_e = 7.292115E-05; // угловая скорость вращения Земли
    int num = (int) abs((t_konez - t_nachalo)/del_t);
    double *t = new double [num];
    double *koord_prom = new double [6];
    for (int i = 0; i < num; i++)
    {
        t[i] = t_nachalo + i*del_t;
    }
    int num_eph = (int) (t_0 - t_nachalo)/del_t;

    // начальные условия
    koord_prom[0] = .311761962891E+07;
    koord_prom[1] = .158781806641E+08;
    koord_prom[2] = .196852387695E+08;
    koord_prom[3] = -.179497814178E+04;
    koord_prom[4] = -.202221393585E+04;
    koord_prom[5] = .192774677277E+04;

    double t_G0 = (9*3600+18*60+10.5009); //9:18:10.5009; Истинное звездное время на гринвичевскую полночь текущей даты
    double t_G = t_G0 + omega_e*(t[num_eph-1]- 3*3600);
    double cos_t = cos(t_G);
    double sin_t = sin(t_G);
    // пересчет координат из ПЗ-90 в ECI
    koord_n[num_eph-1][0] = koord_prom[0]*cos_t - koord_prom[1]*sin_t;
    koord_n[num_eph-1][1] = koord_prom[0]*sin_t + koord_prom[1]*cos_t;
    koord_n[num_eph-1][2] = koord_prom[2];
    koord_n[num_eph-1][3] = koord_prom[3]*cos_t - koord_prom[4]*sin_t - omega_e*koord_n[num_eph-1][1];
    koord_n[num_eph-1][4] = koord_prom[3]*sin_t + koord_prom[4]*cos_t + omega_e*koord_n[num_eph-1][0];
    koord_n[num_eph-1][5] = koord_prom[5];
    delete[] koord_prom;
    koord_prom = nullptr;

    for (int i = num_eph-1; i > 0; i--)
    {
        RungeKutta(t[i], t[i-1], koord_n[i], koord_n[i-1]);
    }
    for (int i = num_eph-1; i < num-1 ; i++)
    {
        RungeKutta(t[i], t[i+1], koord_n[i], koord_n[i+1]);
    }
    delete[] t;
    t = nullptr;
    delete[] koord_prom;
    koord_prom = nullptr;
}
