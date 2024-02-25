
#include <stddef.h>
#include <assert.h>  // for assert
#include <stdio.h>
#include <time.h>
#include <math.h>

/**
 * @brief Compute the factorial of a number
 * 
*/
float newton_raphson(float (*f)(float), float (*f_prime)(float), float x0, float tol, int max_iter) {
    float x = x0;

    for (int i = 0; i < max_iter; i++) {
        float fx = f(x);
        float fpx = f_prime(x);
        float x_new = x - fx / fpx;
        if (fabs(x_new - x) < tol) {
            return x_new;
        }
        x = x_new;
    }
    return x;
} 


/**
 * @brief Compute the factorial of a number
 * 
*/
float bisection(float (*f)(float), float a, float b, float tol, int max_iter) {
    float fa = f(a);
    float fb = f(b);
    assert(fa * fb < 0);
    for (int i = 0; i < max_iter; i++) {
        float c = (a + b) / 2;
        float fc = f(c);
        if (fc == 0 || (b - a) / 2 < tol) {
            return c;
        }
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return (a + b) / 2;
}