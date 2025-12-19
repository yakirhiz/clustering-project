#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sc_header.h"

/* This module holds the main part of the Spectural Clustering Algorithm, as well as the Python Api*/

static PyObject *specClust_CAPI(PyObject *self, PyObject *args);
static double *createMatrix(PyObject *list, int len, int dim);
PyObject* matrix_to_pylist(double *matrix, int rows, int cols);

// input: 2d python list of size (dim)x(len), returns: matrix of the same size, implemented as 1d array
static double *createMatrix(PyObject *list, int len, int dim){
    double *matrix;
    int i, j;
    PyObject *vector;
    PyObject *num;

    matrix = (double *) calloc(len*dim, sizeof(double));
    if (matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Unable to allocate memory for matrix");
        return NULL;
    }

    for (i=0; i<len; i++) {
        vector = PyList_GetItem(list, i); /* Return value: Borrowed reference */
        if (!PyList_Check(vector)){
            PyErr_SetString(PyExc_TypeError, "Each item in the list must be a Python list");
            free(matrix);
            return NULL;
        }
        for (j=0; j<dim; j++){
            num = PyList_GetItem(vector, j);
            if (!PyFloat_Check(num)){
                PyErr_SetString(PyExc_TypeError, "Each item in the list must be a float");
                free(matrix);
                return NULL;
            }
            matrix[i*dim + j] = PyFloat_AsDouble(num); /* Convert a Python float object to double */
            if (matrix[i*dim + j] == -1.0 && PyErr_Occurred()){
                PyErr_SetString(PyExc_TypeError, "Error converting Python float to C double");
                free(matrix);
                return NULL;
            }
        }
    }

    return matrix;
}

// This is the C API function that gets called from python, and calls the rest of the C code.
static PyObject *specClust_CAPI(PyObject *self, PyObject *args){
    int i, d, n, k;
    PyObject *N_obs_floats;
    double *points;

	int k_from_user;
    /* Initial values for optional vars. */
    /*const int k_from_user = 0; // if user doesn't provide a k, we find it using eigengap heuristic*/
    /* The mapping of the variables, ending with a NULL sentinel */
    /*static char *kwlist[] = {"k_from_user", NULL};*/

    /* This parses the Python arguments into a int (ii)  variable named z and int (O) variable named n*/
    if(!PyArg_ParseTuple(args, "Oiii", &N_obs_floats, &n, &d, &k_from_user)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    /* check if N_obs_floats is a list*/
    if (!PyList_Check(N_obs_floats)){
        PyErr_SetString(PyExc_TypeError, "N_obs_floats must be a Python list");
        return NULL;
    }

    points = createMatrix(N_obs_floats, n, d);
    if (points == NULL) {
        return NULL;
    }

    double **Tk_tup = specClust(points, n, d, k_from_user); // if k_from_user == 0 --> we use eigengap heuristic 
    k = (int) *Tk_tup[1]; // the dimensions of Tk_tup[0] are nxk

    // create a PyList that is a python list and return it to Python.
    // the list will be a 2d array len==n*k
    PyObject *pyList = matrix_to_pylist(Tk_tup[0], n, k);
    if (pyList == NULL){ // checking if pyList_new worked
        free(points);
        free(Tk_tup[0]);
        free(Tk_tup); 
        return NULL; // returning NULL after a failed list init
    }

    free(points);
    free(Tk_tup[0]);
    free(Tk_tup);
    //Tk_tup[1] is just a double, k

    return pyList;
}

PyObject* matrix_to_pylist(double *matrix, int rows, int cols) {
    PyObject *outer = PyList_New(rows);
    if (outer == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Unable to create outer list");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        PyObject *row = PyList_New(cols);
        if (row == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Unable to create inner list");
            Py_DECREF(outer);
            return NULL;
        }

        for (int j = 0; j < cols; j++) {
            PyObject *num = PyFloat_FromDouble(matrix[i * cols + j]);
            if (num == NULL) {
                PyErr_SetString(PyExc_MemoryError, "Unable to create float");
                Py_DECREF(row);
                Py_DECREF(outer);
                return NULL;
            }

            PyList_SetItem(row, j, num);  // steals reference
        }

        PyList_SetItem(outer, i, row);  // steals reference
    }

    return outer;
}

/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[] = {
        {"specClust_c",                        /* the Python method name that will be used */
            (PyCFunction) specClust_CAPI,  /* the C-function that implements the Python function and returns static PyObject*  */
            METH_VARARGS,                  /* flags indicating parametersaccepted for this function */
            PyDoc_STR("Runs Spectral Clustering Algorithm until step 5 \ninput is : (n_obs, n, d)\nNote that n_obs is a 2d (python) lists")}, /*  The docstring for the function */
        {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef mySpecClust = {
        PyModuleDef_HEAD_INIT,
        "mySpecClust", /* name of module */
        "a module that implements Spectral Clustering", /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mySpecClust(void)
{
    PyObject *m;
    m = PyModule_Create(&mySpecClust);
    if (!m) {
        return NULL;
    }
    return m;
}