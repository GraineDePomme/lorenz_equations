#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 10
#define R 28
#define B 8.0 / 3.0

#define T0 0.0
#define TF 50.0
#define STEPS_COUNT 5000
#define H (TF - T0) / STEPS_COUNT

double dx_dt(double x, double y, double z) {
    return A * (y - x);
}

double dy_dt(double x, double y, double z) {
    return R * x - y - x * z;
}

double dz_dt(double x, double y, double z) {
    return x * y - B * z;
}

void next_step_rk4(double* x, double* y, double* z) {
    double k1_x = H * dx_dt(*x, *y, *z);
    double k1_y = H * dy_dt(*x, *y, *z);
    double k1_z = H * dz_dt(*x, *y, *z);

    double k2_x = H * dx_dt(*x + k1_x / 2, *y + k1_y / 2, *z + k1_z / 2);
    double k2_y = H * dy_dt(*x + k1_x / 2, *y + k1_y / 2, *z + k1_z / 2);
    double k2_z = H * dz_dt(*x + k1_x / 2, *y + k1_y / 2, *z + k1_z / 2);

    double k3_x = H * dx_dt(*x + k2_x / 2, *y + k2_y / 2, *z + k2_z / 2);
    double k3_y = H * dy_dt(*x + k2_x / 2, *y + k2_y / 2, *z + k2_z / 2);
    double k3_z = H * dz_dt(*x + k2_x / 2, *y + k2_y / 2, *z + k2_z / 2);

    double k4_x = H * dx_dt(*x + k3_x, *y + k3_y, *z + k3_z);
    double k4_y = H * dy_dt(*x + k3_x, *y + k3_y, *z + k3_z);
    double k4_z = H * dz_dt(*x + k3_x, *y + k3_y, *z + k3_z);

    *x += (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6;
    *y += (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6;
    *z += (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6;
}

int main() {
    double t = T0;
    double x = 0;
    double y = 1;
    double z = 0;

    printf("t,x,y,z\n");

    for (t = T0 ; t <= TF+H ; t += H) {
        printf("%f,%f,%f,%f\n", t, x, y, z);
        next_step_rk4(&x, &y, &z);
    }

    return 0;
}