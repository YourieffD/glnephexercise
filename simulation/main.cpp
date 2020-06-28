#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <Runge_kutta.h>

using namespace std;

int main()
{
    time_t start, end;
    double dt = 0.1, T_eph = 49500 + 18 + 3*3600;
    int n = (int) 12*3600/dt,  n_eph = (int) (T_eph-15*3600)/dt, i_delta = 0;
    double t_G0 = 9*3600+18*60+10.5009,  w_e = 7.292115E-05,  t_G = t_G0 + w_e*(T_eph - 3*3600);
    double cosTg = cos(t_G),  sinTg = sin(t_G);
    double X = .311761962891E+07,  Y = .158781806641E+08,  Z = .196852387695E+08;
    double Vx = -.179497814178E+04,  Vy = -.202221393585E+04,  Vz = .192774677277E+04;
    double *Y0 = new double [6] {X*cosTg - Y*sinTg, X*sinTg + Y*cosTg, Z, Vx*cosTg - Vy*sinTg - w_e*(X*sinTg + Y*cosTg), Vx*sinTg + Vy*cosTg + w_e*(X*cosTg - Y*sinTg),Vz};
    double **Yy = new double *[n];
    double *t_12h = new double [n];
    for (int i = 0; i < n; i++){
        Yy[i] = new double [6];
        t_12h[i] = 15*3600+i*dt;
    }
    double **Yy1 = new double *[n_eph];
    for (int i = 0; i < n_eph; i++){
        Yy1[i] = new double [6];
    }
    double **Yy2 = new double *[n-n_eph+1];
    for (int i = 0; i < n-n_eph+1; i++){
        Yy2[i] = new double [6];
    }
    for (int j = 0; j < 6; j++){
        Yy1[0][j] = Y0[j];
        Yy2[0][j] = Y0[j];
    }
    time(&start);
    RungeKutta(T_eph, -dt, Yy1, n_eph,6);
    RungeKutta(T_eph, dt, Yy2, n-n_eph+1,6);
    for (int i =0; i < n; i++)
    {
        if (i <= n_eph-1){
            for (int j = 0; j < 6; j++)
            {
                Yy[i][j] = Yy1[n_eph-i-1][j];
            }
        } else {
            for (int j = 0; j < 6; j++)
            {
                Yy[i][j] = Yy2[i-n_eph+1][j];
            }
        }
    }
    time(&end);
    for (int i = 0; i < n_eph; i++){
        delete [] Yy1[i];
    }
    for (int i = 0; i < n-n_eph+1; i++){
        delete [] Yy2[i];
    }
    delete [] Yy1;
    delete [] Yy2;

    double *Yy_matlab = new double[3];
    double max_delta = 0;
    ifstream input("D:\\MATLAB.txt");
    if (!input){
        cout << "Cant open file, check location of file" << endl;
    } else {
        cout << "OK, Matlab file opened" << endl;
    }
    ofstream output;
    output.open("D:\\CPP.txt");
    for (int i = 0; i < n; i++){
        input >> Yy_matlab[0] >> Yy_matlab[1] >> Yy_matlab[2];
        string YY_str1 = to_string(Yy[i][0]), YY_str2 = to_string(Yy[i][1]), YY_str3 = to_string(Yy[i][2]);
        output << YY_str1 << "\t" << YY_str2 << "\t" << YY_str3 << endl;
        for (int j = 0; j < 3; j++){
            if (abs(Yy[i][j] - Yy_matlab[j]) > max_delta){
                max_delta = abs(Yy[i][j] - Yy_matlab[j]);
                i_delta = i;
            }
        }
    }
    time(&end);
    input.close();
    output.close();
    delete [] t_12h;
    delete [] Y0;
    for (int i = 0; i < n; i++){
        delete [] Yy[i];
    }
    delete [] Yy;
    delete [] Yy_matlab;
    double time_RK = difftime(end, start);
    string time_RK1 = to_string(time_RK*1000000/n), max_delta1 = to_string(max_delta), idelta = to_string(i_delta);
    cout << "time of work of one cicle~: " << time_RK1 << " mcs" << endl;
    cout << "max delta of coords: " << max_delta1 << " m" << endl;
    cout << "number of max delta: " << idelta << endl;
    return 0;
}
