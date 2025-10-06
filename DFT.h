

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


typedef struct {
    int N;
    double dt; //sampling freq
    double T; //period
    double complex *X; //dft coefficients
    int is_real;
} FourierFunction;

void dft_complex(const double* y, int N, double complex* X){
    for (int k = 0; k < N; ++k){
        double complex sum = 0.0 + 0.0*I;
        for (int n = 0; n < N; ++n){
            double theta = -2.0 * M_PI * (double)k * (double)n / (double)N;
            sum += y[n] * cexp(I*theta);//(cos(theta) + I*sin(theta));
        }
        X[k] = sum;
    }
}

FourierFunction construct_dft(const double *y, int N, double dt) {
    assert(N > 0);
    FourierFunction F;
    F.N = N;
    F.dt = dt;
    F.T = (double)N * dt;
    F.X = (double complex*)malloc((size_t)N * sizeof(double complex));
    assert(F.X != NULL);
    dft_complex(y, N, F.X);
    F.is_real = 1;
    return F;
}

double fourier_approximation(FourierFunction fourierfunction, double const t){
    double complex f_hat = 0.0 + 0.0*I;
    for (int i = 0; i < fourierfunction.N; ++i){
        double complex theta = 2 * M_PI * (double)i * (t / fourierfunction.T);;
        f_hat += fourierfunction.X[i] * cexp(I*theta);//(cos(theta) + I*sin(theta));
    }

    return creal(f_hat) / (double)fourierfunction.N;
}

void frequency_filter(FourierFunction* fourierfunction, int low, int high){
    for (int i = 0; i < fourierfunction->N; ++i){
        if (i < low || high < i){
            fourierfunction->X[i] = 0.0 + 0.0*I;
        }
    }
}

void fourier_free(FourierFunction *F) {
    if (F && F->X) {
        free(F->X);
        F->X = NULL;
    }
}