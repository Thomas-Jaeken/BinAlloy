#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define K 0.00008617343 //eV/K
#define EAA -0.436 //eV
#define EBB -0.113 //Ev
#define EAB -0.294 //Ev


void entropy(double *P, double *S, int n_points, int N){
    for(int i=0; i< n_points;i++){
    
        S[i] = K *(2*log(2) - (1+P[i])*log(1+P[i]) - (1-P[i])*log(1-P[i]));
    }
}
void write_to_file(char *fname, double *x,
		   double *y, int n_points)
{
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < n_points; i++){
        fprintf(fp, "%lf, %lf\n",x[i],y[i]);
    }
    fclose(fp);
}

void potential_energy(double *P, double *U, int n_points, int N){
    for(int i=0;i<n_points;i++){
        U[i] = (2*(EAA + EBB + 2*EAB) - 2*P[i]*P[i]*(EAA+EBB-2*EAB));
    }
}

int minim(double *S,double *U,double T,int n_points,double *P_points){
    int result =0;
    double min = 0;
    for(int i=0;i<n_points;i++){

        if(U[i]-T*S[i]<min){
            min=U[i]-T*S[i];
            result = i;          
            }
    }
    return result;
}


int main(){
    int N = 50; // # of atoms
    int n_temps=2000;
    int P_points=100000;
    double P_step = 2./P_points;
    double T_step = 1.;
    double T_start = 0.;//Kelvin
    double *T_range=malloc(sizeof(double)*n_temps);
    double *F_list=malloc(sizeof(double)*n_temps);
    double *C_list=malloc(sizeof(double)*(n_temps));
    double *order=malloc(sizeof(double)*n_temps);
    double *P_list = malloc(sizeof(double)*P_points);
    for(int i=0;i<P_points;i++){
        P_list[i] = i*P_step-1;}
    double *S_list = malloc(sizeof(double)*P_points);
    entropy(P_list,S_list,P_points,N);
    double *U_list = malloc(sizeof(double)*P_points);
    potential_energy(P_list,U_list,P_points,N);

    for(int i=0;i<n_temps;i++){
        T_range[i] = i*T_step + T_start;
        int d=minim(S_list,U_list,T_range[i],P_points,P_list);
        
        order[i] = P_list[d];
        F_list[i] = U_list[d];
        if(i<1 ){;}
        else{C_list[i-1] = (F_list[i] - F_list[i-1])/2/T_step;}
        }

    write_to_file("mean_field_approx.csv",T_range,order,n_temps);
    write_to_file("energy_temp.csv",T_range,F_list,n_temps);
    write_to_file("capacity_temp.csv",T_range,C_list,n_temps-1);
    write_to_file("U.csv",P_list,U_list,P_points);
    write_to_file("S.csv",P_list,S_list,P_points);
}