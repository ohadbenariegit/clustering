#include "global.h"

void handle_error(int);

void validate_input(int);

void *create_array(int, int);

double **create_matrix(int, int);

double **create_square_matrix(int);

double **create_identity_matrix(int order);

void free_matrix(int, double **);

void print_matrix(int, int, double **, int);

void copy_square_matrix(int, double **, double **);

void transpose_square_matrix(int, double **);