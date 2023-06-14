#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
// #define K 0.00008617343 // eV/K
// #define EAA -0.436      // eV
// #define EBB -0.113      // Ev
// #define EAB -0.294      // Ev

double bond(int A, int B);
double hetero(int A, int B);
double eval_energy(int *A, int *B, int N1D);
double eval_delta_energy(int *A, int *B, int *trialA, int *trialB, int N1D, int atomA, int atomB);
double eval_order(int *A, int *B, int N_total);
double short_range_order(int *A, int *B, int N1D);
void write_to_file(char *fname, double *x,
                   double *y, int n_points);
void trail_change(int *A, int *B, int N1D, int *atomA, int *atomB);
void restore_change(int *A, int *B, int N1D, int *atomA, int *atomB);
void metropolis(int N_equil, int *A, int *B, int N1D, int N_total, double T, int steps, double *P, double *E, double *R);
double variance_given_avg(double *energy, int N, double avg);
double corr_abs(double *xarray, int n_points, int k);
double av_sq(double *xarray, int n_points);
void blockavs(double *array, double *blockavis, int n_points, int b);
double variance(double *array, int n_points);
int inefficiency(double *data, int N);