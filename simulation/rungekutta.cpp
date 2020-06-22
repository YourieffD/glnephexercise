#include <include\libglnsvpos\rungekutta.h>
#include <cmath>
#include <include\libglnsvpos\diffs.h>
#include <iostream>
#include <fstream>

using namespace std;

void RungeKutta(double t0, double tn, double *koord0, double *koord){
    double del_t = tn-t0;
    double *proizv = new double [6];
    double *K_1 = new double [6];
    double *K_2 = new double [6];
    double *K_3 = new double [6];
    double *K_4 = new double [6];
    double *vremen = new double [6];
// Rschet K1
    diffs(t0, koord0, proizv);
    for (int j = 0; j<=5; j++){
        K_1[j] = del_t*proizv[j];
        vremen[j] = koord0[j]+K_1[j]/2.0;
    }
// Rschet K2
    diffs(t0+del_t/2.0, vremen, proizv);
    for (int j = 0; j<=5; j++){
        K_2[j] = del_t*proizv[j];
        vremen[j] = koord0[j]+K_2[j]/2.0;
    }
// Rschet K3
    diffs(t0+del_t/2.0, vremen, proizv);
    for (int j = 0; j<=5; j++){
        K_3[j] = del_t*proizv[j];
        vremen[j] = koord0[j]+K_3[j];
    }
// Rschet K4
    diffs(tn, vremen, proizv);
    for (int j = 0; j<=5; j++){
        K_4[j] = del_t*proizv[j];
        koord[j] = koord0[j] + (K_1[j]+2*K_2[j]+2*K_3[j]+K_4[j])/6.0;
    }
    delete[] K_1;
    K_1 = nullptr;
    delete[] K_2;
    K_2 = nullptr;
    delete[] K_3;
    K_3 = nullptr;
    delete[] K_4;
    K_4 = nullptr;
    delete[] proizv;
    proizv = nullptr;
    delete[] vremen;
    vremen = nullptr;
}
