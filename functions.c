#include "functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#define K 0.00008617343 // eV/K
#define EAA -0.436      // eV
#define EBB -0.113      // Ev
#define EAB -0.294      // Ev

double bond(int A, int B)
{
    if (A == B && B == 1)
    {
        return EAA;
    }
    else if (A == B)
    {
        return EBB;
    }
    else
    {
        return EAB;
    }
}
double hetero(int A, int B)
{
    if (A == B)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

double eval_energy(int *A, int *B, int N1D)
{
    double energy = 0;
    for (int i = 0; i < N1D; i++)
    {
        for (int j = 0; j < N1D; j++)
        {
            for (int k = 0; k < N1D; k++)
            {
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k]);
                int k_next = k + 1;
                if (k == N1D - 1)
                {
                    k_next = 0;
                }
                int i_next = i + 1;
                if (i == N1D - 1)
                {
                    i_next = 0;
                }
                int j_next = j + 1;
                if (j == N1D - 1)
                {
                    j_next = 0;
                }

                energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k_next]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k_next]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k_next]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k]);
                energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k_next]);
            }
        }
    }
    return energy; // we don't normalize the energy!!
}

double eval_delta_energy(int *A, int *B, int *trialA, int *trialB, int N1D, int atomA, int atomB)
{
    double energy = 0;
    int i = atomA / (N1D * N1D);
    int j = (atomA % (N1D * N1D)) / N1D;
    int k = atomA % N1D;
    // energy from atom A
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + j * N1D + k]);
    int k_next = k + 1;
    if (k == N1D - 1)
    {
        k_next = 0;
    }
    int i_next = i + 1;
    if (i == N1D - 1)
    {
        i_next = 0;
    }
    int j_next = j + 1;
    if (j == N1D - 1)
    {
        j_next = 0;
    }
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k_next]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + j * N1D + k_next]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + (j_next)*N1D + k]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k_next]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + (j_next)*N1D + k_next]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[(i_next)*N1D * N1D + j * N1D + k]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k_next]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[(i_next)*N1D * N1D + j * N1D + k_next]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[(i_next)*N1D * N1D + (j_next)*N1D + k]);
    energy += bond(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k_next]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[(i_next)*N1D * N1D + (j_next)*N1D + k_next]);
    // energy from atom B
    i = atomB / (N1D * N1D);
    j = (atomB % (N1D * N1D)) / N1D;
    k = atomB % N1D;

    energy += bond(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + j * N1D + k]);
    k_next = k - 1;
    if (k == 0)
    {
        k_next = N1D - 1;
    }
    i_next = i - 1;
    if (i == 0)
    {
        i_next = N1D - 1;
    }
    j_next = j - 1;
    if (j == 0)
    {
        j_next = N1D - 1;
    }
    energy += bond(A[i * N1D * N1D + j * N1D + k_next], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + j * N1D + k_next], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[i * N1D * N1D + (j_next)*N1D + k], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + (j_next)*N1D + k], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[i * N1D * N1D + (j_next)*N1D + k_next], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[i * N1D * N1D + (j_next)*N1D + k_next], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[(i_next)*N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[(i_next)*N1D * N1D + j * N1D + k], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[(i_next)*N1D * N1D + j * N1D + k_next], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[(i_next)*N1D * N1D + j * N1D + k_next], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[(i_next)*N1D * N1D + (j_next)*N1D + k], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[(i_next)*N1D * N1D + (j_next)*N1D + k], trialB[i * N1D * N1D + j * N1D + k]);
    energy += bond(A[(i_next)*N1D * N1D + (j_next)*N1D + k_next], B[i * N1D * N1D + j * N1D + k]) - bond(trialA[(i_next)*N1D * N1D + (j_next)*N1D + k_next], trialB[i * N1D * N1D + j * N1D + k]);

    return -energy; // we don't normalize the energy!!/(2*N1D*N1D*N1D);// per atom!!
}

double eval_order(int *A, int *B, int N_total)
{
    // we represent A atoms with 1 and B with 0
    int A_on_a = 0;
    for (int i = 0; i < N_total; i++)
    {
        A_on_a += A[i];
        if (A[i] < 0)
        {
            printf("ERROR: there is a negative value in the lattice.\n");
        }
    }
    return 2. / N_total * A_on_a - 1;
}

double short_range_order(int *A, int *B, int N1D)
{
    double order = 0;
    for (int i = 0; i < N1D; i++)
    {
        for (int j = 0; j < N1D; j++)
        {
            for (int k = 0; k < N1D; k++)
            {
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k]);
                int k_next = k + 1;
                if (k == N1D - 1)
                {
                    k_next = 0;
                }
                int i_next = i + 1;
                if (i == N1D - 1)
                {
                    i_next = 0;
                }
                int j_next = j + 1;
                if (j == N1D - 1)
                {
                    j_next = 0;
                }

                order += hetero(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + j * N1D + k_next]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[i * N1D * N1D + (j_next)*N1D + k_next]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + j * N1D + k_next]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k]);
                order += hetero(A[i * N1D * N1D + j * N1D + k], B[(i_next)*N1D * N1D + (j_next)*N1D + k_next]);
            }
        }
    }
    int N_tot = N1D * N1D * N1D;
    return 1. / (4 * N_tot) * (order - 4 * N_tot);
}

void write_to_file(char *fname, double *x,
                   double *y, int n_points)
{
    FILE *fp = fopen(fname, "w");
    for (int i = 0; i < n_points; i++)
    {
        fprintf(fp, "%lf, %lf\n", x[i], y[i]);
    }
    fclose(fp);
}

void trail_change(int *A, int *B, int N1D, int *atomA, int *atomB)
{
    int N_total = N1D * N1D * N1D;
    *atomA = ((double)rand() / (double)RAND_MAX) * N_total; // decide an atom in A
    *atomB = ((double)rand() / (double)RAND_MAX) * N_total; // decide a neigbour in b
    // now invert them both, because the total amount of atoms of each kind should be kept constant.
    int tempA = A[*atomA];
    int tempB = B[*atomB];
    while (tempA == tempB)
    {
        *atomA = ((double)rand() / (double)RAND_MAX) * N_total; // decide an atom in A
        *atomB = ((double)rand() / (double)RAND_MAX) * N_total; // decide a neigbour in b
        // now invert them both, because the total amount of atoms of each kind should be kept constant.
        tempA = A[*atomA];
        tempB = B[*atomB];
    }
    A[*atomA] = tempB;
    B[*atomB] = tempA;
    // printf("switched %d and %d",atomA,atomB);
}

void restore_change(int *A, int *B, int N1D, int *atomA, int *atomB)
{
    int changes = 1 * 1 * 1;
    for (int i = 0; i < changes; i++)
    {
        // now invert them both, because the total amount of atoms of each kind should be kept constant.
        int tempA = A[*atomA];
        int tempB = B[*atomB];
        A[*atomA] = tempB;
        B[*atomB] = tempA;
        // printf("switched %d and %d",atomA,atomB);
    }
}

void metropolis(int N_equil, int *A, int *B, int N1D, int N_total, double T, int steps, double *P, double *E, double *R)
{
    int acceptance_rate = 0;
    int *trialA = malloc(sizeof(int) * N_total);
    int *trialB = malloc(sizeof(int) * N_total);
    for (int i = 0; i < N_total; i++)
    {
        trialA[i] = A[i];
        trialB[i] = B[i];
    }

    E[0] = eval_energy(A, B, N1D);
    P[0] = eval_order(A, B, N_total);
    R[0] = short_range_order(A, B, N1D);

    for (int i = 0; i < steps + N_equil; i++)
    {
        int atomA = 0;
        int atomB = 0;
        double deltaE = 0;
        if (i <= N_equil)
        { // equilibration phase
            trail_change(trialA, trialB, N1D, &atomA, &atomB);
            deltaE = eval_delta_energy(A, B, trialA, trialB, N1D, atomA, atomB);

            if (deltaE < 0 || exp(-(deltaE) / K / T) >= ((double)rand() / (double)RAND_MAX))
            {
                restore_change(A, B, N1D, &atomA, &atomB);
                E[0] += deltaE;
                P[0] = eval_order(A, B, N_total);
                R[0] = short_range_order(A, B, N1D);
            }
            else
            { // reject
                restore_change(trialA, trialB, N1D, &atomA, &atomB);
            }
        }
        else
        { // start saving
            trail_change(trialA, trialB, N1D, &atomA, &atomB);
            deltaE = eval_delta_energy(A, B, trialA, trialB, N1D, atomA, atomB);

            if (deltaE < 0 || exp(-(deltaE) / K / T) >= ((double)rand() / (double)RAND_MAX))
            {
                restore_change(A, B, N1D, &atomA, &atomB);
                acceptance_rate++;
                E[i - N_equil] += deltaE;
                P[i - N_equil] = eval_order(A, B, N_total);
                R[i - N_equil] = short_range_order(A, B, N1D);
            }
            else
            { // reject
                restore_change(trialA, trialB, N1D, &atomA, &atomB);
                E[i - N_equil] = E[i - N_equil - 1];
                P[i - N_equil] = P[i - N_equil - 1];
                R[i - N_equil] = R[i - N_equil - 1];
            }
        }
    }
    free(trialA);
    trialA = NULL;
    free(trialB);
    trialB = NULL;
}

double variance_given_avg(double *energy, int N, double avg)
{
    double result = 0; // Var = (<E^2> - <E>^2)
    for (int i = 0; i < N; i++)
    {
        result += energy[i] * energy[i];
    }
    return result / N - avg * avg;
}

double corr_abs(double *xarray, int n_points, int k)
{
    double average = 0;
    for (int i = 0; i < n_points - k; i++)
    {
        average += xarray[i + k] * xarray[i];
    }
    return average / (n_points - k);
}
double av_sq(double *xarray, int n_points)
{
    double average = 0;
    for (int i = 0; i < n_points; i++)
    {
        average += xarray[i];
    }
    return pow(average / (n_points), 2);
}
void blockavs(double *array, double *blockavis, int n_points, int b)
{
    for (int i = 0; i < floor(n_points / b); i++)
    {
        blockavis[i] = 0;
        for (int j = 0; j < b; j++)
        {
            blockavis[i] += array[b * i + j];
        }
        blockavis[i] /= b;
    }
}
double variance(double *array, int n_points)
{
    double average = 0;
    for (int i = 0; i < n_points; i++)
    {
        average += array[i];
    }
    average /= n_points;
    double var = 0;
    for (int i = 0; i < n_points; i++)
    {
        var += pow((array[i] - average), 2);
    }
    return var / n_points;
}

int inefficiency(double *data, int N)
{
    int k_values = N;
    // centralise the values
    double avergi = 0;
    for (int i = 0; i < N; i++)
    {
        avergi += data[i];
    }
    avergi /= N;
    for (int i = 0; i < N; i++)
    {
        data[i] -= avergi;
    }

    double *auto_corr = malloc(sizeof(double) * k_values);
    double fi_sin = corr_abs(data, N, 0);
    double fi_sout = av_sq(data, N);
    double *k_array = malloc(sizeof(double) * k_values);
    for (int i = 0; i < k_values; i++)
    {
        auto_corr[i] = (corr_abs(data, N, i) - fi_sout) / (fi_sin - fi_sout);
        k_array[i] = (double)i;
        if (auto_corr[i] <= exp(-2.))
        {
            return i;
        }
    }

    printf("error: the statistical inefficiency was not found");
    return k_values; // note for cases where it was not found, the variance is likely zero, so it doesnt matter that much.
                     // if it's just very low but not zero, a high
}
