#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"

/*--------------------------------------------------------------------------------------------------------------------*/
void validate_input(int condition) {
    if (condition == FALSE) {
        printf("Invalid Input!");
        exit(1);
    }
}

void handle_error(int error) {
    if (error == TRUE) {
        printf("An Error Has Occurred");
        exit(1);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
void *create_array(int n_items, int item_size) {
    void *arr;

    arr = malloc(n_items * item_size);
    handle_error(arr == NULL);

    return arr;
}

double **create_matrix(int rows, int cols) {
    double **mat;
    int row;

    mat = (double **) create_array(rows, sizeof(double *));

    for (row = 0; row < rows; row++)
        mat[row] = (double *) create_array(cols, sizeof(double));

    return mat;
}

double **create_square_matrix(int order) {
    return create_matrix(order, order);
}

double **create_identity_matrix(int order) {
    double **mat;
    int row, col;

    mat = create_square_matrix(order);

    for (row = 0; row < order; row++) {
        for (col = 0; col < order; col++)
            mat[row][col] = (row == col);
    }

    return mat;
}

void free_matrix(int rows, double **mat) {
    int row;

    for (row = 0; row < rows; row++)
        free(mat[row]);

    free(mat);
}

void print_matrix(int rows, int cols, double **mat, int jacobi) {
    int row, col;
    double val;
    char suffix;

    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            /* jacobi eigenvalues are supposed to be positive */
            val = mat[row][col];
            if (jacobi == TRUE && row == 0 && fabs(val) < 0.0001)
                val = 0;

            suffix = (col < cols - 1) ? ',' : '\n';

            printf("%.4f%c", val, suffix);
        }
    }
}

void copy_matrix(int rows, int cols, double **dest, double **src) {
    int row;

    for (row = 0; row < rows; row++)
        memcpy(dest[row], src[row], cols * sizeof(double));
}

void copy_square_matrix(int order, double **dest, double **src) {
    copy_matrix(order, order, dest, src);
}

void transpose_square_matrix(int order, double **mat) {
    int row, col;
    double tmp;

    for (row = 1; row < order; row++) {
        for (col = 0; col < row; col++) {
            tmp = mat[row][col];
            mat[row][col] = mat[col][row];
            mat[col][row] = tmp;
        }
    }
}