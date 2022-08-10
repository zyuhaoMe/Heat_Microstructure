//
// Created by zyuhao on 2022-01-18.
// This file contains commonly-used functions about dynamic array allocation
// and deallocation

// common function used both in CA and heat transfer calculation

#include "CA_heat.h"

double *Allocate_1D_Double(int size, double value) {
    double *v;
    int i;

    v = (double *) malloc(size * sizeof (double));

    for (i = 0; i < size; i++) {
        v[i] = value;
    }
    return (v);
}

void Free_1D_Double(double *a, int size) {
    free((void *) a);
}

double **Allocate_2D_Double(int size1, int size2, double value) {
    double** v;
    int i, j;
    v = (double**)malloc(size1 * sizeof(double*));

    for (i = 0; i < size1; i++) {
        v[i] = (double*)malloc(size2 * sizeof(double));
    }
    // initialize
    for (j = 0; j < size2; j++) {
        for (i = 0; i < size1; i++) {
            v[i][j] = value;
        }
    }
    return (v);
}

void Free_2D_Double(double** a, int size1, int size2) {
    int i;
    for (i = 0; i < size1; i++) {
        free((void*) a[i]);
    }
    free(a);
}

double ***Allocate_3D_Double(int size1, int size2, int size3, double value) {
    double ***v;
    int i, j, k;
    v = (double ***) malloc(size1 * sizeof (double **));

    for (i = 0; i < size1; i++) {
        v[i] = (double **) malloc(size2 * sizeof (double *));
        for (j = 0; j < size2; j++) {
            v[i][j] = (double *) malloc(size3 * sizeof (double));
        }
    }
    // initialize
    for (k = 0; k < size3; k++) {
        for (j = 0; j < size2; j++) {
            for (i = 0; i < size1; i++) {
                v[i][j][k] = value;
            }
        }
    }
    return (v);
}

void Free_3D_Double(double*** a, int size1, int size2, int size3) {
    int i, j;
    for (i = 0; i < size1; i++) {
        for (j = 0; j < size2; j++) {
            free((void*) a[i][j]);
        }
    }
    free(a);
}

double ****Allocate_4D_Double(int size1, int size2, int size3, int size4, double value) {
    double**** v;
    int i, j, k, l;
    v = (double****)malloc(size1 * sizeof(double***));

    for (i = 0; i < size1; i++) {
        v[i] = (double***)malloc(size2 * sizeof(double**));
        for (j = 0; j < size2; j++) {
            v[i][j] = (double**)malloc(size3 * sizeof(double*));
            for (k = 0; k < size3; k++) {
                v[i][j][k] = (double*)malloc(size4 * sizeof(double));
            }
        }
    }
    // initialize
    for (l = 0; l < size4; l++) {
        for (k = 0; k < size3; k++) {
            for (j = 0; j < size2; j++) {
                for (i = 0; i < size1; i++) {
                    v[i][j][k][l] = value;
                }
            }
        }
    }
    return (v);
}

void Free_4D_Double(double**** a, int size1, int size2, int size3, int size4) {
    int i, j, k;
    for (i = 0; i < size1; i++) {
        for (j = 0; j < size2; j++) {
            for (k = 0; k < size3; k++) {
                free((void*) a[i][j][k]);
            }
        }
    }
    free(a);
}