#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>


//#include "system_equation.h"
#include "datatypes.h"

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <cvodes/cvodes_impl.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#define IJth(A,length,i,j) A[i*length+j] 
#define NV_DATA_S(v) ( NV_CONTENT_S(v)->data )

void WrongArgs(char *name);
void ProcessArgs(int argc, char *argv[], FILE **input, FILE **input2, FILE **output, FILE **parfile, FILE **settings, FILE **initialvaluesfile);
void ProcessArgsML(int argc, char *argv[], FILE **input, FILE **input2, FILE **output, FILE **parfile, FILE **settings, FILE **initialvaluesfile);
int LineCount(FILE *f);
int ExtractColumn(FILE *f, int number_of_columns, int specified_column, int datalength, realtype *data);
realtype  LinearInterpolation(realtype *X, realtype *Y, realtype T, int length);
int Integrate(void *cvode_mem, UserData data, realtype y0[], realtype **output, realtype Rtol, realtype Atol, realtype T0, realtype T1);
int Sensitivity(void *cvode_mem, UserData data, realtype y0[], realtype **output, realtype **outputsol, realtype Rtol, realtype Atol, realtype T0, realtype T1);
int spline(realtype *x, realtype *y, int datalength, realtype yp1, realtype ypn, realtype *y2);
void splint(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *y);
int splintdot(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *ydot);
int splintdotequalspace(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *ydot);
void splintequal(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *y);
//int PrintSensToFile(FILE *outputfile, realtype *current_solution, UserData data);
int printoutput(FILE *outputfile, realtype *time, int datalength, realtype *output);
int AfterSolve(realtype *current_solution, UserData data);
int PrintToFile(FILE *outputfile, realtype *current_solution, UserData data);
//void splintdot(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *ydot);
//int likelihood(realtype *array1, realtype *array2, int data_length, realtype *result, realtype sigma, int stepno10s);
int likelihood(realtype *array1, realtype *array2, realtype error_value, realtype *result, realtype sigma);
realtype error(realtype *array1, realtype *array2, int data_length, int stepno10s);
int generate_proposed_pars(realtype *current_par, realtype *proposed_par, realtype *sigma, realtype *minimum_value, realtype *maximum_value, int no_of_pars, gsl_rng *r);
int next_mc_link(realtype *proposed, realtype *current_parameter, realtype *proposed_error, realtype *current_error, int no_of_pars, gsl_rng *r, realtype sigma);
//int PrintParametersToFile(FILE *outputfile, realtype *current_pars, int number_of_parameters, realtype *likelihood);
int PrintParametersToFile(FILE *outputfile, realtype *current_pars, int number_of_parameters, realtype error_value);
int calculate_hr_from_potential(realtype *result, realtype *times, realtype *potential, realtype datalength);
int calculate_hr_from_times(realtype *result, realtype *times, realtype datalength);
realtype  fast_interp(realtype *X, realtype *Y, realtype T, int length);
//int Sensitivity(void *cvode_mem, UserData data, realtype y0[], realtype *output, realtype *outputsol, realtype Rtol, realtype Atol, realtype T0, realtype T1);
void PrintOutputS(N_Vector *uS, UserData data);
realtype **matrix(long nrl, long nrh, long ncl, long nch);
void printerror(char error_text[]);
int PrintSensToFile(FILE *outputfile, realtype **current_solution, UserData data);
void free_matrix(realtype **m, long nrl, long nrh, long ncl, long nch);
int CalcRelSens(realtype **senssum, realtype **solution, realtype **sens, UserData data);
int AddSens(realtype **senssum, realtype **sens, UserData data);
int PotToHr(realtype **pot, realtype **solution, realtype **HR, UserData data);
int Scale(realtype **sens, realtype scale, UserData data);
