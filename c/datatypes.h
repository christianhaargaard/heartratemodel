/* In this file we define different datatypes to be used */

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */

#ifndef USERDATA_GUARD
#define USERDATA_GUARD
typedef struct {
/* Type : UserData */
  int no_equations;
	int no_parameters;      /* Number of problem parameters */
  realtype *p;           /* problem parameters */
	int inputsize;        /* Contains the number of lines in data*/
	realtype *t;            /* Time series for data */
  realtype *bp;           /* blood pressure data */
  int resp_length;        /* Number of data points for respiration */
  realtype *resp_t;       /* Interpolation should be done in python, and fed to the model with steady time steps */
  realtype *resp_freq;
	} *UserData;
#endif

#define Ith(v,i)    NV_Ith_S(v,i)       /* i-th vector component i=1..NEQ */
#define IJth(A,length,i,j) A[i*length+j]
#define DoubleArrayAccess(a, NEQ, i, j) ((a)[(i) * (NEQ) + (j)])
