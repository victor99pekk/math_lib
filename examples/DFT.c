#include <stdio.h>
#include "../DFT.h"

static double f_true(double t, double T) {
    return cos(2.0*M_PI*3*t/T) + 0.3*sin(2.0*M_PI*5*t/T);
}

int main(void) {
    int N = 10;         // number of sample points
    double dt = 0.001;     // arbitrary sampling frequency; T = N*dt
    double T = N * dt;
    double T_function = T;// / 2; // Period of the true function

    double *y = (double*)malloc((size_t)N * sizeof(double));
    for (int n = 0; n < N; ++n) {
        double t = n * dt;
        y[n] = f_true(t, T_function);
    }

    FourierFunction F = construct_dft(y, N, dt);

    // Evaluate at dense grid
    double gamma = 1; // gamma has to be at least 1 to avoid aliasing
    for (int i = 0; i < 10; ++i) {
        double t_eval = gamma * i * dt;
        double y_hat = fourier_approximation(F, t_eval);
        printf("t=%.6f  f_true=%.9f fhat=%.9f\n", t_eval, f_true(t_eval, T_function), y_hat);
    }

    fourier_free(&F);
    free(y);
    return 0;
}