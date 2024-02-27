
#include <stddef.h>
#include <assert.h>  // for assert
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>  // for malloc and free

/**
 * @brief A tuple of two floats.
 * @param x constant 1.
 * @param y constant 1.
 */
typedef struct Point
{
    float x;
    float y;
} Point;


/**
 * @brief generates h number of equidistant x values in the interval [a, b].
 * @param interval array with the start and end of an interval the interval [a, b].
 * @param h the step size.
 * @return an array of x values.
*/
float* generateX_array(float interval[], float h){
    int nbr_points = (interval[1] - interval[0]) / h;
    float *x = (float*)malloc(nbr_points * sizeof(float));
    for(int index = 0; index < nbr_points; index++){
        x[index] = interval[0] + index * h;
    }
    return x;
}


/**
 * @brief The Simpson rule is a method for approximating the definite integral of a function.
 * @param interval array with the start and end of an interval the interval [a, b].
 * @param nbr_points the number of equidistant points to use in the approximation.
 * @param f the function to integrate.
 * @return the integral of the function.
 */
float simpson(float interval[], int nbr_points, float (*f)(float)) {
    float h = (interval[1] - interval[0]) / (nbr_points-1);
    float integral = 0;

    float* x = generateX_array(interval, h);
    //  Ym = a + mh

    for(int index = 1; index < nbr_points; index++) {
        float x0 = interval[0] + (index-1) * h;
        float xn = x0 + h;
        float middle = (xn + x0) / 2;
        integral += (h/6) * (f(x[index + 1]) + 4 * f(middle) + f(x[index]));
    }

    /* for(int index = 1; index < nbr_points; index++) {
        // float middle = x[index-1] + (h/2);
        float middle = (x[index-1] + x[index]) / 2;
        integral += (h/6) * (f(x[index + 1]) + 4 * f(middle) + f(x[index]));
    } */
    free(x);  // don't forget to free the memory when you're done with it

    return integral;
}
