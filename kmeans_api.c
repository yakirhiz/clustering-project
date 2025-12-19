#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "km_header.h"

static PyObject* kmeanspp_capi(PyObject *self, PyObject *args);
static double *createMatrix(PyObject *list, int len, int dim);
PyObject* matrix_to_pylist(double *matrix, int rows, int cols);

static PyObject* kmeanspp_capi(PyObject *self, PyObject *args){
    int K, N, d, MAX_ITER, i;
    PyObject *N_obs_floats, *cents_floats;
    double *observations, *cents;
    /* This parses the Python arguments into a int (iiii)  variable named z and int (O) variable named n*/
    if(!PyArg_ParseTuple(args, "iiiiOO", &K, &N, &d, &MAX_ITER, &N_obs_floats, &cents_floats)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    /* check if N_obs_floats is a list*/
    if (!PyList_Check(N_obs_floats) || !(PyList_Check(cents_floats))){
        PyErr_SetString(PyExc_TypeError, "N_obs_floats and cents_floats must be Python lists");
        return NULL;
    }

    observations = createMatrix(N_obs_floats, N, d);
    if (observations == NULL) {
        return NULL;
    }

    cents = createMatrix(cents_floats, K, d);
    if (cents == NULL) {
        free(observations);
        return NULL;
    }

    // uses N_obs and cents and then frees them
    double *newCents = kmeans(K, N, d, MAX_ITER, observations, cents);

    // create a PyList that is a python list and return it to Python.
    // the list will be a 2d-array of size d*k
    PyObject *pyList = matrix_to_pylist(newCents, K, d);
    if (pyList == NULL){ // checking if pyList_new worked
        free(newCents); 
        free(observations);
        return NULL; // returning NULL after a failed list init
    }

    free(newCents);
    free(observations);

    return pyList;
}


// input: 2d python list of size (dim)x(len), returns: matrix of the same size, implemented as 1d array
static double *createMatrix(PyObject *list, int len, int dim){
    double *matrix;
    int i, j;
    PyObject *vector;
    PyObject *num;

    matrix = (double *) malloc(sizeof(double) * len*dim);
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
        {"kmeanspp_c",                   /* the Python method name that will be used */
                (PyCFunction) kmeanspp_capi,  /* the C-function that implements the Python function and returns static PyObject*  */
                     METH_VARARGS,                 /* flags indicating parametersaccepted for this function */
                        PyDoc_STR("Calculates K clusters\ninput is : (k,n,d,maxiter,n_obs,cents)\nNote that n_obs & cents are 2d (python) lists")}, /*  The docstring for the function */
        {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", /* name of module */
        "a module that implements kmeans", /* module documentation, may be NULL */
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
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}
