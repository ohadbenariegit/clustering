#include <Python.h>
#include "common/helper_methods.h"
#include "kmeans/kmeans.h"
#include "spkmeans.h"

double **python_matrix_to_c_matrix(int rows, int cols, PyObject *py_mat) {
    double **c_mat;

    int row, col;
    PyObject *py_list, *item;

    c_mat = create_matrix(rows, cols);

    for (row = 0; row < rows; row++) {
        py_list = PyList_GetItem(py_mat, row);

        for (col = 0; col < cols; col++) {
            item = PyList_GetItem(py_list, col);
            c_mat[row][col] = PyFloat_AsDouble(item);
        }
    }

    return c_mat;
}

double **python_square_matrix_to_c_square_matrix(int order, PyObject *py_mat) {
    return python_matrix_to_c_matrix(order, order, py_mat);
}

PyObject *c_matrix_to_python_matrix(int rows, int cols, double **c_mat) {
    PyObject *py_mat;

    int row, col;
    PyObject *py_list;

    py_mat = PyList_New(rows);

    for (row = 0; row < rows; row++) {
        py_list = PyList_New(cols);

        for (col = 0; col < cols; col++) {
            PyList_SetItem(py_list, col, Py_BuildValue("f", c_mat[row][col]));
        }

        PyList_SetItem(py_mat, row, py_list);
    }

    return py_mat;
}

PyObject *c_square_matrix_to_python_square_matrix(int order, double **c_mat) {
    return c_matrix_to_python_matrix(order, order, c_mat);
}

static PyObject *fit_goal(PyObject *self, PyObject *args) {
    int k;
    char *goal, *filename;

    matrix output;

    if (!PyArg_ParseTuple(args, "iss", &k, &goal, &filename))
        return NULL;

    output = output_goal(k, goal, filename, TRUE);

    return c_matrix_to_python_matrix(output.rows, output.cols, output.mat);
}

static PyObject *fit_kmeans(PyObject *self, PyObject *args) {
    int n, k;
    PyObject *datapoints_py, *centroids_py;
    double **datapoints, **centroids;

    if (!PyArg_ParseTuple(args, "iiOO", &n, &k, &datapoints_py, &centroids_py))
        return NULL;

    datapoints = python_matrix_to_c_matrix(n, k, datapoints_py);
    centroids = python_square_matrix_to_c_square_matrix(k, centroids_py);

    set_centroids(n, k, datapoints, centroids);

    return c_square_matrix_to_python_square_matrix(k, centroids);
}

static PyMethodDef spkmeansmoduleMethods[] = {
        {"fit_goal",   (PyCFunction) fit_goal,   METH_VARARGS, NULL},
        {"fit_kmeans", (PyCFunction) fit_kmeans, METH_VARARGS, NULL},
        {NULL, NULL, 0,                                        NULL}
};

static struct PyModuleDef spkmeansmodule = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        spkmeansmoduleMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void) {
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m)
        return NULL;
    return m;
}