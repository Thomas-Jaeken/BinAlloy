#define PY_SSIZE_T_CLEAN
#include </Users/thomasjaeken/miniconda3/include/python3.9/Python.h>
#include "functions.h"
#include <time.h>

static PyObject *
OrderParam(PyObject *self, PyObject *args)
{
    // Processing user input
    double T_step, T_start; // Kelvin
    int N1D, metrosteps, n_temps;
    if (!PyArg_ParseTuple(args, "iiidd",
                          &N1D, &metrosteps, &n_temps, &T_step, &T_start))
    {
        printf("please provide a dimension for the lattice, a number of steps for the metropolis algorithm to take,\
        a number of temperature steps, and the start and end of the temperature range you are interested in.");
    }
    int N_total = N1D * N1D * N1D;
    // allocating memory

    double *T_range = malloc(sizeof(double) * n_temps);
    for (int i = 0; i < n_temps; i++)
    {
        T_range[i] = i * T_step + T_start;
    }
    double *P_save = malloc(sizeof(double) * n_temps);
    double *E_save = malloc(sizeof(double) * n_temps);
    double *R_save = malloc(sizeof(double) * n_temps);
    double *var_E_save = malloc(sizeof(double) * n_temps);
    double *var_P_save = malloc(sizeof(double) * n_temps);
    double *var_R_save = malloc(sizeof(double) * n_temps);
    double *ineff_E = malloc(sizeof(double) * n_temps);
    double *ineff_P = malloc(sizeof(double) * n_temps);
    double *ineff_R = malloc(sizeof(double) * n_temps);
    int *Cu = malloc(sizeof(int) * N_total);
    int *Zn = malloc(sizeof(int) * N_total);
    for (int i = 0; i < N_total; i++)
    {
        Cu[i] = 1;
        Zn[i] = 0;
    }
    int A = 0, B = 0;

    // calculation
    srand(time(NULL)); // random seed for the metropolis algorithm

    // randomising the initial state to speed up convergence
    for (int i = 0; i < N_total; i++)
    {
        trail_change(Cu, Zn, N1D, &A, &B);
    }

    // Run through the temperature list and perform MC
    for (int i = 0; i < n_temps; i++)
    {
        // allocate temporary memory to save the generated data points
        double *P_temp = malloc(sizeof(double) * metrosteps);
        double *E_temp = malloc(sizeof(double) * metrosteps);
        double *R_temp = malloc(sizeof(double) * metrosteps);
        double *time_temp = malloc(sizeof(double) * metrosteps);
        for (int j = 0; j < metrosteps; j++)
        {
            time_temp[j] = j;
        };
        int N_equil = 5e5 * 300 / T_range[i]; // Initial guess

        metropolis(N_equil, Cu, Zn, N1D, N_total, T_range[i], metrosteps, P_temp, E_temp, R_temp);

        for (int j = 0; j < metrosteps; j++)
        {
            P_save[i] += P_temp[j] / metrosteps;
            E_save[i] += E_temp[j] / metrosteps;
            R_save[i] += R_temp[j] / metrosteps;
        }
        var_E_save[i] = variance_given_avg(E_temp, metrosteps, E_save[i]);
        var_P_save[i] = variance_given_avg(P_temp, metrosteps, P_save[i]);
        var_R_save[i] = variance_given_avg(R_temp, metrosteps, R_save[i]);
        ineff_E[i] = inefficiency(E_temp, metrosteps);
        ineff_P[i] = inefficiency(P_temp, metrosteps);
        ineff_R[i] = inefficiency(R_temp, metrosteps);

        // clear allocated memory
        free(P_temp);
        P_temp = NULL;
        free(E_temp);
        E_temp = NULL;
        free(R_temp);
        E_temp = NULL;
    }

    // Conversion of result to python types
    PyObject *T_export = PyList_New(n_temps);
    PyObject *P_export = PyList_New(n_temps);
    PyObject *E_export = PyList_New(n_temps);
    PyObject *R_export = PyList_New(n_temps);
    PyObject *P_Var_export = PyList_New(n_temps);
    PyObject *E_Var_export = PyList_New(n_temps);
    PyObject *R_Var_export = PyList_New(n_temps);
    PyObject *P_Ineff_export = PyList_New(n_temps);
    PyObject *E_Ineff_export = PyList_New(n_temps);
    PyObject *R_Ineff_export = PyList_New(n_temps);

    for (int i = 0; i < n_temps; i++)
    {
        PyList_SetItem(T_export, i, PyLong_FromDouble(T_range[i]));
        PyList_SetItem(P_export, i, PyLong_FromDouble(P_save[i]));
        PyList_SetItem(E_export, i, PyLong_FromDouble(E_save[i]));
        PyList_SetItem(R_export, i, PyLong_FromDouble(R_save[i]));
        PyList_SetItem(P_Var_export, i, PyLong_FromDouble(var_P_save[i]));
        PyList_SetItem(E_Var_export, i, PyLong_FromDouble(var_E_save[i]));
        PyList_SetItem(R_Var_export, i, PyLong_FromDouble(var_R_save[i]));
        PyList_SetItem(P_Ineff_export, i, PyLong_FromDouble(ineff_P[i]));
        PyList_SetItem(E_Ineff_export, i, PyLong_FromDouble(ineff_E[i]));
        PyList_SetItem(R_Ineff_export, i, PyLong_FromDouble(ineff_R[i]));
    }

    // clear all C allocated memory
    free(Cu);
    Cu = NULL;
    free(Zn);
    Zn = NULL;
    free(P_save);
    P_save = NULL;
    free(E_save);
    E_save = NULL;
    free(R_save);
    R_save = NULL;
    free(var_E_save);
    var_E_save = NULL;
    free(var_P_save);
    var_P_save = NULL;
    free(var_R_save);
    var_R_save = NULL;
    free(ineff_E);
    ineff_E = NULL;
    free(ineff_P);
    ineff_P = NULL;
    free(ineff_R);
    ineff_R = NULL;
    free(T_range);
    T_range = NULL;

    // return the dictionary holding the results
    return Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O}",
                         "temperature range", T_export,
                         "order parameterr", P_export,
                         "energy", E_export,
                         "short-range order", R_export,
                         "order variance", P_Var_export,
                         "energy variance", E_Var_export,
                         "short order variancne", R_Var_export,
                         "order inefficiency", P_Ineff_export,
                         "energy inefficiency", E_Ineff_export,
                         "short order inefficiency", R_Ineff_export);
}

static PyMethodDef BinAlloyMethods[] = {
    // lists the methods and points to their implementation.
    {"run_monte_carlo", OrderParam, METH_VARARGS,
     "calculate the temperature dependecy of the order parameter."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef binalloymodule = {
    // defines the module name and points to its methods.
    PyModuleDef_HEAD_INIT,
    "BinAlloy", /* name of module */
    NULL,       /* module documentation, may be NULL */
    -1,         /* size of per-interpreter state of the module,
                   or -1 if the module keeps state in global variables. */
    BinAlloyMethods};

PyMODINIT_FUNC
PyInit_BinAlloy(void)
{
    // the entry point for the compiler, this points to the module struct.
    return PyModule_Create(&binalloymodule);
}

int main(int argc, char *argv[])
{
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL)
    {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }

    /* Add a built-in module, before Py_Initialize */
    if (PyImport_AppendInittab("BinAlloy", PyInit_BinAlloy) == -1)
    {
        fprintf(stderr, "Error: could not extend in-built modules table\n");
        exit(1);
    }

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(program);

    /* Initialize the Python interpreter.  Required.
       If this step fails, it will be a fatal error. */
    Py_Initialize();

    /* Optionally import the module; alternatively,
       import can be deferred until the embedded script
       imports it. */
    PyObject *pmodule = PyImport_ImportModule("BinAlloy");
    if (!pmodule)
    {
        PyErr_Print();
        fprintf(stderr, "Error: could not import module 'BinAlloy'\n");
    }

    PyMem_RawFree(program);
    return 0;
}