#include "matrix.h"
#include <assert.h>
#include <stdio.h>

#define ASSERT(x, msg) if(!(x)){printf("%s\n", msg); assert(x);}
//#define printFourier(f) printf("Fourier: nbr points: %d, degree: %d\n", f.nbr_points, f.degree); MAT_PRINT(f.matrix); MAT_PRINT(f.y_values);
#define printFourier(f) printf("\nFourier: nbr points: %d, degree: %d\n", f.nbr_points, f.degree); fourier_print(f, "Matrix", "Y values", 4);


typedef struct Fourier
{
    Mat matrix;
    Mat y_values;
    int degree;
    int nbr_points;
} Fourier;

Fourier fourier_alloc(int n, int degree)
{
    Fourier f;
    f.matrix = mat_alloc(n, degree);
    f.y_values = mat_alloc(n, 1);
    f.nbr_points = n;
    f.degree = degree;
    return f;
}

void fourier_print(Fourier fourier, const char *name, const char *y, int padding)
{
    int extraPadd = fourier.degree * 7.5;
    printf("%*s%s: ", (int)padding, "", name);
    printf("        %*s%s: \n", (int)extraPadd, "", y);        
    for(size_t i = 0; i < fourier.matrix.rows; i++){
        printf("%*s    ", (int)padding, "");
        for(size_t j = 0; j < fourier.matrix.cols; j++){
            printf("%f ", MAT_AT(fourier.matrix, i, j));
        }
        printf("        %f\n", MAT_AT(fourier.y_values, i, 0));
    }
    printf("\n");
}


float P(float x)
{
    float modulo = fmod(x, M_PI);

    if (modulo == 0) return 0;
    else return modulo/2;
}

float random_float() {
    // Seed the random number generator only once
    static int seeded = 0;
    if (!seeded) {
        srand((unsigned int)time(NULL));
        seeded = 1;
    }

    // Generate a random integer between 0 and 100
    int randInt = rand() % 101;

    // Convert to a floating-point number between 0 and 1 with two decimal places
    float randomValue = (float)randInt / 100.0;

    return randomValue;
}

void produce_randX(float start, float end, int n, float x[], float y[])
{
    //produce random data points from function P
    float interval = (end - start);

    for(int i = 0; i < n; i++){
        x[i] = start + interval*random_float();
        y[i] = P(x[i]);
    }
}

void x_values(float start, float end, int n, float x[], float y[])
{
    for(int i = 1; i <= n; i++){
        x[i-1] = (-M_PI + (i-1)*(2*M_PI/(n)));
        y[i-1] = P(x[i-1]);
    }
}

void print_dataPoints(float x[], float y[], int n)
{
    for(int i = 0; i < n; i++){
        printf("x: %f, y: %f\n", x[i], y[i]);
    }
}

void fill_fourierMatrix(Mat x, Fourier f, int nbr_points, int degree)
{
    Mat fourierMatrix = f.matrix;
    for(int row = 0; row < nbr_points; row++){
        float x_value = MAT_AT(x, row, 0);
        for(int col = 0; col < degree; col++){
            MAT_AT(fourierMatrix, row, col) = sin(col*x_value);   
        }
    }
}

Fourier create_fourier(int n, int degree, float x[], float y[])
{
    Fourier f = fourier_alloc(n, degree);
    Mat fourierX = mat_alloc(n, 1);
    Mat fourierY = mat_alloc(n, 1);
    fill_with_array(f.y_values, y);
    fill_with_array(fourierX, x);
    fill_fourierMatrix(fourierX, f, n, degree);
    return f;
}

int plot(const char *filename, const char *outputFilename, float x[], float y[], const int n) {
    FILE *dataFile = fopen(filename, "w");

    if (dataFile == 0){
        perror("Error opening file");
        return 1;
    }

    // Generate data for the function (e.g., sine function)
    for(int i = 0; i < n; i++){
        fprintf(dataFile, "%f %f\n", x[i], y[i]);
    }

    fclose(dataFile);

    // Use the plot command to visualize the data
    char plotCommand[100]; // Allocate enough space for the command
    int formatted = snprintf(plotCommand, sizeof(plotCommand), "gnuplot -e \"set term png; set output '%s'; plot '%s' with points\"", outputFilename, filename);

    if(formatted < 0 || formatted >= 100){
        // Handle the error, e.g., print an error message or return an error code
        return -1;
    }

    int result = system(plotCommand);

    return result;
}

