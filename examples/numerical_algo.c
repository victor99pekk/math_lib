#include <stdio.h>
#include "../numerical_algorithms.h"

// Test function: f(x) = x^2
float test_func(float x) {
    return x * x;
}

int main(void) {
    float interval[2] = {0.0, 1.0};
    int nbr_points = 11; // Simpson's rule requires odd number of points
    float exact = 1.0f / 3.0f;

    float simpson_result = simpson(interval, nbr_points, test_func);
    printf("Simpson's rule: %.8f (error: %.8f)\n", simpson_result, fabsf(simpson_result - exact));
    // Simpson's rule: 0.31650004 (error: 0.01683331)

    float epsilon = 1e-6f;
    int max_depth = 20;
    float adaptive_result = adaptive_Simpson(test_func, interval[0], interval[1], epsilon, max_depth);
    printf("Adaptive Simpson's rule: %.8f (error: %.8f)\n", adaptive_result, fabsf(adaptive_result - exact));
    // Adaptive Simpson's rule: 0.33333331 (error: 0.00000003)
    return 0;
}