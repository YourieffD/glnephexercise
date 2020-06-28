#include<difury.h>
#include<Runge_kutta.h>
#include<cmath>
#include <iostream>

void RungeKutta(double t_start, double dt, double **&Yy, int size, int size_s){
    double *K1 = new double[6];
    double *K2 = new double[6];
    double *K3 = new double[6];
    double *K4 = new double[6];
    double *y_pr1 = new double[6];
    double *y_pr2 = new double[6];
    double *y_pr3 = new double[6];

    for (int i = 0; i < size-1; i++)
    {
        double t_tek = t_start + i*dt;
        difury(t_tek, Yy[i], K1);
        for (int j = 0; j < size_s; j++){
            y_pr1[j] = Yy[i][j] + 0.5*dt*K1[j];
        }

        difury(t_tek + 0.5*dt, y_pr1, K2);
        for (int j = 0; j < size_s; j++){
            y_pr2[j] = Yy[i][j] + 0.5*dt*K2[j];
        }

        difury(t_tek + 0.5*dt, y_pr2, K3);
        for (int j = 0; j < size_s; j++){
            y_pr3[j] = Yy[i][j] + dt*K3[j];
        }

        difury(t_tek + dt, y_pr3, K4);
        for (int j = 0; j < size_s; j++){
            Yy[i+1][j] = Yy[i][j] + dt/6 * (K1[j]+2*K2[j]+2*K3[j]+K4[j]);
        }
    }
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] y_pr1;
    delete[] y_pr2;
    delete[] y_pr3;
}
