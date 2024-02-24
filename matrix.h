#ifndef matrix_H_
#define matrix_H_

#include <stddef.h>
#include <assert.h>  // for assert
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifndef matrix_MALLOC
#include <stdlib.h>  // for malloc
#define matrix_MALLOC malloc
#endif //NN_MALLOC

#ifndef matrix_ASSERT
#include <assert.h>
#define matrix_ASSERT assert
#endif // NN_ASSERT
#include <stdbool.h>

typedef struct{
    size_t rows;
    size_t cols;
    size_t stride;
    float *es;
} Mat;

#define MAT_AT(m, i , j) (m).es[(i)*(m).stride + (j)]
#define MAT_PRINT(m) mat_print(m, #m, 0);


float rand_float(void)
{
    return (float)rand() / (float) RAND_MAX;
}

Mat mat_alloc(size_t rows, size_t cols)
{
    Mat m;
    m.rows = rows;
    m.cols = cols;
    m.stride = cols;
    m.es = (float*)matrix_MALLOC(sizeof(*m.es)*rows*cols);
    matrix_ASSERT(m.es != NULL);
    return m;
}

float sigmoidf(float x)
{
    return 1.f / (1.f + expf(-x));
}

void mat_sig(Mat m)
{
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(m, i, j) = sigmoidf(MAT_AT(m, i, j));
        }
    }
}

void mat_dot(Mat dst, Mat a, Mat b)
{
    matrix_ASSERT(a.cols == b.rows);
    matrix_ASSERT(dst.rows == a.rows);
    matrix_ASSERT(dst.cols == b.cols);
    size_t n = a.cols;
    
    for(size_t i = 0; i < dst.rows; i++){
        for(size_t j = 0; j < dst.cols; j++){
            MAT_AT(dst, i, j) = 0;
            for(size_t k = 0; k < n; k++){
                MAT_AT(dst, i, j) += MAT_AT(a, i, k) * MAT_AT(b, k, j);
            }
        }
    }
}

Mat mat_row(Mat m, size_t row)
{
    return (Mat){
        .rows = 1,
        .cols = m.cols,
        .stride = m.stride,
        .es = &MAT_AT(m, row, 0)
    };
}

Mat mat_col(Mat m, size_t col)
{
    Mat result = mat_alloc(m.rows, 1);
    for(size_t i = 0; i < m.rows; i++){
        MAT_AT(result, i, 0) = MAT_AT(m, i, col);
    }
    return result;
}

void mat_copy(Mat dst, Mat src)
{
    matrix_ASSERT(dst.rows == src.rows);
    matrix_ASSERT(dst.cols == src.cols);
    for(size_t i = 0; i < dst.rows; i++){
        for(size_t j = 0; j < dst.cols; j++){
            MAT_AT(dst, i, j) = MAT_AT(src, i, j);
        }
    }
}

void mat_sum(Mat dst, Mat a)
{
    matrix_ASSERT(dst.rows == a.rows);
    matrix_ASSERT(dst.cols == a.cols);
    for(size_t i = 0; i < dst.rows; i++){
        for(size_t j = 0; j < dst.cols; j++){
            MAT_AT(dst, i, j) += MAT_AT(a, i, j);
        }
    }
}

void mat_print(Mat m, const char *name, int padding)
{
    printf("%*s%s: \n", (int)padding, "", name);
    for(size_t i = 0; i < m.rows; i++){
        printf("%*s    ", (int)padding, "");
        for(size_t j = 0; j < m.cols; j++){
            printf("%f ", MAT_AT(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void mat_fill(Mat m, float x)
{
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(m, i, j) = x;
        }
    }
}

void mat_rand(Mat m, float low, float high)
{
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(m, i, j) = rand_float()*(high-low)+low;
        }
    }
}

void mat_symmetric_rand(Mat m, float low, float high)
{
    for(size_t i = 0; i < m.cols; i++){
        float number = rand_float()*(high-low)+low;
        for(size_t j = 0; j < m.rows; j++){
            MAT_AT(m, j, i) = number;
        }
    }
}

void diagonalDominant(Mat m)
{
    for(size_t i = 0; i < m.rows; i++){
        float sum = 0;
        for(size_t j = 0; j < m.cols; j++){
            if(i != j){
                sum += fabs(MAT_AT(m, i, j));
            }
        }
        MAT_AT(m, i, i) = sum + 1;
    }
}

void fill_matrix(Mat m, float *data)
{
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(m, i, j) = data[i*m.cols + j];
        }
    }
}

void fill_with_array(Mat m, float array[])
{
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(m, i, j) = array[i*m.cols + j];
        }
    }

}

void fill_with_array_col(Mat m, float array[], size_t col)
{
    for(size_t i = 0; i < m.rows; i++){
        MAT_AT(m, i, col) = array[i];
    }
}

void copy_col(Mat dst, Mat src, int destCol, int srcCol)
{
    for(size_t i = 0; i < src.rows; i++){
        MAT_AT(dst, i, destCol) = MAT_AT(src, i, srcCol);
    }
}




void switchRow(Mat m, int row1, int row2) {
    for (int i = 0; i < m.cols; ++i) {
        float temp = MAT_AT(m, row1, i);
        MAT_AT(m, row1, i) = MAT_AT(m, row2, i);
        MAT_AT(m, row2, i) = temp;
    }
}

Mat fill_identity(Mat m)
{
    m = mat_alloc(m.rows, m.cols);
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            if(i == j){
                MAT_AT(m, i, j) = 1;
            } else {
                MAT_AT(m, i, j) = 0;
            }
        }
    }
    return m;
}

Mat subtract_cols(Mat a, Mat b, int colA, int colB)
{
    Mat result = mat_alloc(a.rows, 1);
    for(size_t i = 0; i < a.rows; i++){
        MAT_AT(result, i, 0) = MAT_AT(a, i, colA) - MAT_AT(b, i, colB);
    }
    return result;
}

void prepare(Mat m, Mat x, Mat identity, int index, int ndim) {
    matrix_ASSERT(m.rows == x.rows);

    if(index >= m.cols || index >= m.rows){
        return;
    }

    for (int row = index; row <= ndim; ++row) {
        bool largest = false;

        while (!largest) {
            largest = true;

            for (int i = index + 1; i <= ndim; ++i) {
                if (MAT_AT(m, row, index) < MAT_AT(m, i, index)) {
                    largest = false;
                    switchRow(m, row, i);
                    switchRow(identity, row, i);
                    switchRow(x, row, i);
                }
            }
        }
    }
}

void eliminate(Mat m, Mat x, Mat identity)
{
    for(size_t i = 0; i < m.cols; i++){
        //prepare(m, x, identity, i, m.rows - 1);
        for(size_t j = i + 1; j < m.rows; j++){
            if(MAT_AT(m, j, i) == 0 || MAT_AT(m, i, i) == 0){
                continue;
            } 
            float ratio = MAT_AT(m, j, i) / MAT_AT(m, i, i);
            for(size_t k = i; k < m.cols; k++){
                MAT_AT(m, j, k) -= ratio * MAT_AT(m, i, k);
            }
            MAT_AT(x, j, 0) -= ratio * MAT_AT(x, i, 0);
        }
    }
}

Mat back_sub(Mat m, Mat x)
{

    matrix_ASSERT(m.rows == x.rows);
    Mat result = mat_alloc(m.rows, 1);
    mat_copy(result, x);
    for(int i = m.rows - 1; i >= 0; i--){
        for(int j = m.cols; j > (i); j--){
            MAT_AT(result, i, 0) -= MAT_AT(m, i, j) * MAT_AT(result, j, 0);
        }
        if(MAT_AT(m, i, i) != 0){
            MAT_AT(result, i, 0) /= MAT_AT(m, i, i);
        }else{
            MAT_AT(result, i, 0) = 0;
        }
    }
    return result;
}

Mat mat_transpose(Mat m)
{
    Mat result = mat_alloc(m.cols, m.rows);
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.cols; j++){
            MAT_AT(result, j, i) = MAT_AT(m, i, j);
        }
    }
    return result;
}

Mat mat_mul(Mat a, Mat b)
{
    matrix_ASSERT(a.cols == b.rows);
    Mat result = mat_alloc(a.rows, b.cols);
    mat_dot(result, a, b);
    return result;
}

Mat solve(Mat m, Mat x)
{
    Mat id = fill_identity(m);
    eliminate(m, x, id);
    return back_sub(m, x);
}

Mat solve_leastSquares(Mat m, Mat x)
{
    Mat transposed = mat_transpose(m);
    Mat mtm = mat_alloc(transposed.rows, m.cols);
    mat_dot(mtm, transposed, m);
    Mat right = mat_alloc(transposed.rows, x.cols);
    mat_dot(right, transposed, x);
    return solve(mtm, right);
}

void print_result(Mat m, Mat x)
{
    Mat result = solve(m, x);
    for(size_t i = 0; i < m.rows; i++){
        printf("x%zu = %f\n", i, MAT_AT(result, i, 0));
    }
}

void print_resultMatrix(Mat x)
{
    for(size_t i = 0; i < x.rows; i++){
        printf("x%zu = %f\n", i, MAT_AT(x, i, 0));
    }
}

Mat compute_gramSchmidt(Mat m)
{

    Mat R = mat_alloc(m.cols, m.cols);
    mat_fill(R, 0);
    Mat Q = mat_alloc(m.rows, m.rows);
    mat_copy(Q, m);

    for(size_t i = 0; i < m.cols; i++){
        Mat vi = mat_col(m, i);
        float k;
        
        for(size_t j = 0; j < i; j++){
            Mat qi = mat_col(Q, j);
            Mat transpose = mat_transpose(qi);
            k = MAT_AT(mat_mul(transpose, vi), 0, 0);
            MAT_AT(R, j, i) = k;
            //MAT_AT(Q, i, 0) = (vi, qi, i, 0);
            for(size_t h = 0; h < m.rows; h++){
                MAT_AT(Q, h, i) -= MAT_AT(qi, h, 0) * k;
            }
        }

        float norm = 0;
        for(size_t j = 0; j < Q.rows; j++){
            norm += MAT_AT(Q, j, i) * MAT_AT(Q, j, i);
        }
        norm = sqrt(norm);
        MAT_AT(R, i, i) = norm;
        MAT_PRINT(R);
        if(norm == 0) continue;

        for(size_t j = 0; j < Q.rows; j++){
            MAT_AT(Q, j, i) = MAT_AT(Q, j, i) / MAT_AT(R, i, i);
        }
    }
    MAT_PRINT(Q);
    MAT_PRINT(R);
    return mat_mul(Q, R);
}

#endif //NN_H_

