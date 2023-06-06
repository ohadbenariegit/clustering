#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "common/helper_methods.h"

/*--------------------------------------------------------------------------------------------------------------------*/
#define MAX_ITER 100
#define EPS 1.0 * pow(10, -5)

/*--------------------------------------------------------------------------------------------------------------------*/
double **datapoints;
int n, m, k;

/*--------------------------------------------------------------------------------------------------------------------*/
/* read data points from input file */
void setup(char *filename) {
    FILE *file;

    char ch, *line;
    int row, col, first_line = TRUE, row_len = 1, max_row_len = 0;

    file = fopen(filename, "r");
    handle_error(file == NULL);

    n = 0, m = 1;

    /* determine number of data points, vector length, and max row length */
    while ((ch = (char) fgetc(file)) != EOF) {
        row_len++;

        if (first_line && ch == ',')
            m++;
        else if (ch == '\n') {
            if (row_len > max_row_len)
                max_row_len = row_len;

            first_line = FALSE;
            row_len = 1;
            n++;
        }
    }

    if (ferror(file) != 0)
        handle_error(TRUE);

    rewind(file);

    datapoints = create_matrix(n, m);
    line = (char *) create_array(max_row_len, sizeof(char));

    /* fill datapoints matrix */
    for (row = 0; row < n; row++) {
        handle_error(fgets(line, max_row_len, file) == NULL);
        line[strlen(line) - 1] = 0;

        for (col = 0; col < m; col++) {
            datapoints[row][col] = atof(strtok((col == 0) ? line : NULL, ","));
        }
    }

    handle_error(fclose(file) == EOF);

    free(line);
}

/*--------------------------------------------------------------------------------------------------------------------*/
typedef struct rotation_matrix {
    double c, s;
    int i, j;
} rotation_matrix;

typedef struct wam_and_ddg {
    double **wam;
    double **ddg;
} wam_and_ddg;

double get_weight(int i, int j) {
    double weight = 0;
    int l;

    for (l = 0; l < m; l++)
        weight += pow(datapoints[i][l] - datapoints[j][l], 2);

    weight = exp(-(sqrt(weight) / 2));

    return weight;
}

double **get_wam(void) {
    double **wam;

    int i, j;

    wam = create_square_matrix(n);

    for (i = 0; i < n; i++) {
        wam[i][i] = 0;

        for (j = i + 1; j < n; j++)
            wam[i][j] = wam[j][i] = get_weight(i, j);
    }

    return wam;
}

double get_degree(int i, double **wam) {
    double degree = 0;
    int j;

    for (j = 0; j < n; j++)
        degree += wam[i][j];

    return degree;
}

wam_and_ddg get_wam_and_ddg(void) {
    wam_and_ddg wamAndDdg;

    int i, j;

    wamAndDdg.wam = get_wam();
    wamAndDdg.ddg = create_square_matrix(n);

    for (i = 0; i < n; i++) {
        wamAndDdg.ddg[i][i] = get_degree(i, wamAndDdg.wam);

        for (j = i + 1; j < n; j++)
            wamAndDdg.ddg[i][j] = wamAndDdg.ddg[j][i] = 0;
    }

    return wamAndDdg;
}

double **get_ddg(void) {
    wam_and_ddg wamAndDdg = get_wam_and_ddg();

    free_matrix(n, wamAndDdg.wam);

    return wamAndDdg.ddg;
}

double **get_lnorm(void) {
    double **lnorm;
    wam_and_ddg wamAndDdg;
    double **wam, **ddg;

    int i, j;

    wamAndDdg = get_wam_and_ddg();
    wam = wamAndDdg.wam;
    ddg = wamAndDdg.ddg;
    lnorm = create_square_matrix(n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            lnorm[i][j] = (i == j) - (wam[i][j] / sqrt(ddg[i][i] * ddg[j][j]));
    }

    free_matrix(n, wam);
    free_matrix(n, ddg);

    return lnorm;
}

/* set rotation matrix by calculating c, s, i, and j */
int set_rotation_mat(double **A, rotation_matrix *P) {
    double abs_val, max = 0, theta, t;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            abs_val = fabs(A[i][j]);
            if (abs_val > max) {
                max = abs_val;
                P->i = i;
                P->j = j;
            }
        }
    }

    /* A is diagonal already */
    if (max == 0)
        return -1;

    theta = (A[P->j][P->j] - A[P->i][P->i]) / (2 * A[P->i][P->j]);
    t = ((theta >= 0) - (theta < 0)) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    P->c = 1 / sqrt(pow(t, 2) + 1);
    P->s = t * P->c;

    return 0;
}

/* set A_tag according to rotation matrix and A matrix */
void set_A_tag(double **A, rotation_matrix *P, double **A_tag) {
    int r;

    for (r = 0; r < n; r++) {
        if (r != P->i && r != P->j) {
            A_tag[P->i][r] = A_tag[r][P->i] = (P->c * A[r][P->i]) - (P->s * A[r][P->j]);
            A_tag[P->j][r] = A_tag[r][P->j] = (P->c * A[r][P->j]) + (P->s * A[r][P->i]);
        }
    }

    A_tag[P->i][P->i] =
            (pow(P->c, 2) * A[P->i][P->i]) + (pow(P->s, 2) * A[P->j][P->j]) - (2 * P->s * P->c * A[P->i][P->j]);
    A_tag[P->j][P->j] =
            (pow(P->s, 2) * A[P->i][P->i]) + (pow(P->c, 2) * A[P->j][P->j]) + (2 * P->s * P->c * A[P->i][P->j]);
    A_tag[P->i][P->j] = A_tag[P->j][P->i] = 0;
}

double get_off_squared(double **mat) {
    double off = 0;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++)
            off += pow(mat[i][j], 2);
    }

    off *= 2;

    return off;
}

/* update V by multiplying it by P */
void update_V(double **V, rotation_matrix *P) {
    double temp;
    int l;

    for (l = 0; l < n; l++) {
        temp = V[l][P->i];
        V[l][P->i] = P->c * temp - P->s * V[l][P->j];
        V[l][P->j] = P->s * temp + P->c * V[l][P->j];
    }
}

double **get_jacobi(double **mat) {
    double **jacobi;
    double **A, **A_tag, **V;
    rotation_matrix P;

    int count = 0, i;
    double off_squared_diff;

    A = create_square_matrix(n);
    copy_square_matrix(n, A, mat);
    A_tag = create_square_matrix(n);
    copy_square_matrix(n, A_tag, mat);
    V = create_identity_matrix(n);
    jacobi = create_matrix(n + 1, n);

    do {
        if (set_rotation_mat(A, &P) == -1)
            break;

        set_A_tag(A, &P, A_tag);

        update_V(V, &P);

        off_squared_diff = get_off_squared(A) - get_off_squared(A_tag);

        copy_square_matrix(n, A, A_tag);

        count++;
    } while (count < MAX_ITER && off_squared_diff > EPS);

    for (i = 0; i < n; i++)
        jacobi[0][i] = A[i][i];

    copy_square_matrix(n, jacobi + 1, V);

    free_matrix(n, A);
    free_matrix(n, A_tag);
    free_matrix(n, V);

    return jacobi;
}

/*--------------------------------------------------------------------------------------------------------------------*/
typedef struct eigen {
    double eigenvalue;
    double *eigenvector;
} eigen;

/* comparing eigen objects */
int compare_eigens(const void *e1, const void *e2) {
    double diff = (((eigen *) e1)->eigenvalue - ((eigen *) e2)->eigenvalue);

    return (diff > 0) - (diff < 0);
}

/* return sorted list of eigen objects based on jacobi */
eigen *get_sorted_eigens(double **jacobi) {
    eigen *eigens;

    int i;

    eigens = (eigen *) create_array(n, sizeof(eigen));
    transpose_square_matrix(n, jacobi + 1);

    for (i = 0; i < n; i++) {
        eigens[i].eigenvalue = jacobi[0][i];
        eigens[i].eigenvector = jacobi[1 + i];
    }

    qsort(eigens, n, sizeof(eigen), compare_eigens);

    return eigens;
}

/* set k using heuristic */
void set_k(eigen *eigens) {
    double val, max = -1;
    int i;

    eigen *eigens_shift = eigens - 1;

    for (i = 1; i <= n / 2; i++) {
        val = fabs(eigens_shift[i].eigenvalue - eigens_shift[i + 1].eigenvalue);

        if (val > max) {
            max = val;
            k = i;
        }
    }
}

/* return the data points to be used for spk */
double **get_spk_datapoints() {
    double **spk_dps;
    double **lnorm;
    double **jacobi;
    eigen *eigens;

    double row_l2_norm;
    int i, j;

    lnorm = get_lnorm();
    jacobi = get_jacobi(lnorm);
    eigens = get_sorted_eigens(jacobi);

    if (k == 0) {
        set_k(eigens);
        handle_error(k <= 1);
    } else {
        validate_input(k > 1 && k < n);
    }

    spk_dps = create_matrix(n, k);

    for (i = 0; i < n; i++) {
        row_l2_norm = 0;

        for (j = 0; j < k; j++) {
            spk_dps[i][j] = eigens[j].eigenvector[i];
            row_l2_norm += pow(spk_dps[i][j], 2);
        }

        if (row_l2_norm != 0) {
            row_l2_norm = sqrt(row_l2_norm);

            for (j = 0; j < k; j++)
                spk_dps[i][j] /= row_l2_norm;
        }
    }

    free_matrix(n, lnorm);
    free_matrix(n + 1, jacobi);
    free(eigens);

    return spk_dps;
}

/*--------------------------------------------------------------------------------------------------------------------*/
matrix output_goal(int input_k, char *goal, char *filename, int with_spk) {
    matrix output;

    k = input_k;

    setup(filename);

    if (strcmp(goal, "wam") == 0) {
        output.mat = get_wam();
    } else if (strcmp(goal, "ddg") == 0) {
        output.mat = get_ddg();
    } else if (strcmp(goal, "lnorm") == 0) {
        output.mat = get_lnorm();
    } else if (strcmp(goal, "jacobi") == 0) {
        output.mat = get_jacobi(datapoints);
    } else if (with_spk == TRUE && strcmp(goal, "spk") == 0) {
        output.mat = get_spk_datapoints();
    } else {
        validate_input(FALSE);
    }

    output.rows = ((strcmp(goal, "jacobi") == 0) ? n + 1 : n);
    output.cols = ((strcmp(goal, "spk") == 0) ? k : n);

    free_matrix(n, datapoints);

    return output;
}

/*--------------------------------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]) {
    matrix output;

    validate_input(argc == 3);

    output = output_goal(FALSE, argv[1], argv[2], FALSE);

    print_matrix(output.rows, output.cols, output.mat, strcmp(argv[1], "jacobi") == 0);

    free_matrix(output.rows, output.mat);

    return 0;
}