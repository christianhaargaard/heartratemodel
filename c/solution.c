#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

/* Sundials libraries */
#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#define Ith(v,i) NV_Ith_S(v,i)       /* i-th vector component i=0..NEQ-1 */
#define DoubleArrayAccess(a, NEQ, i, j) ((a)[(i) * (NEQ) + (j)])

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
  realtype *tiltangle;         /* blood pressure data */
  realtype *resp;
  realtype *bp_spline_coef;
  realtype *tiltangle_spline_coef;
  realtype *resp_spline_coef;
  realtype *flags;
	} *UserData;
#endif

/* Since I will be calling these functions before I define them, input and output of them are prescribed here */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int check_flag(void *flagvalue, char *funcname, int opt);
static int set_initial_values(realtype t, N_Vector y, void *f_data);
static realtype wall_strain(realtype p, void *f_data);
/*static realtype fast_interp(realtype *X, realtype *Y, realtype T, int length);*/
static realtype mean_start_pressure(void *f_data, realtype time);
static realtype firing_rate(N_Vector y, realtype ew, void *f_data);
static realtype parasymp_baro(realtype firing, void *f_data);
static realtype symp_baro(realtype firing, void *f_data);
static realtype resp_weight(realtype angle,void *f_data);
static realtype baro_weight(realtype angle,void *f_data);
static realtype heart_rate(N_Vector y, void *f_data);
/*static int spline(realtype *x, realtype *y, int datalength, realtype yp1, realtype ypn, realtype *y2);*/
static void splint(realtype *xa, realtype *ya, realtype *y2a, int datalength, realtype x, realtype *y);
static realtype Tresp(realtype fresp, void *f_data);
static realtype min(realtype a, realtype b);

/* Function to set initial conditions for states, based on the assumption that they are in steady state at the beginning*/
static int set_initial_values(realtype t, N_Vector y, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype p, AM0, a1, a2, b1, b2, td, tA, qp, kiN, tN, qs, mu, KA, KN;

  /*p0 = data->p[0]; //Needs to call function to calculate p0 from input pressure data. Maybe make it dependent on if it is given as a parameter.*/

  // Load rest of parameters
  AM0=data->p[2];
  a1=data->p[3];
  a2=data->p[4];
  b1=data->p[5];
  b2=data->p[6];
  td=data->p[24];
  tA=data->p[25];
  qp=data->p[26];
  kiN=data->p[27]; 
  tN=data->p[28];
  qs=data->p[29];
  mu=data->p[31];
  KA=data->p[32];
  KN=data->p[35];
  
  /* Interpolate pressure data */
  splint(data->t, data->bp,data->bp_spline_coef,data->inputsize,t, &p);

  Ith(y,0) = p;
  /*printf("Time for initial point= %2.2e\n",t);*/

  // temporary parameters, that cannot be varied ever
  realtype kappa;
  kappa = (a1*b2+a2*b1)/(a1*b2+a2*b1+b1*b2);
  
  // Set s1 and s2 if they are 0
  if(data->p[7]==0){ //s1
    data->p[7] = RCONST(1.0)/(RCONST(1.0)+kappa)*sqrt(AM0)/(sqrt(AM0)-RCONST(1.0));
    }
  if(data->p[8]==0){ //s1
    data->p[8] = kappa/(RCONST(1.0)+kappa);
    }

  // Calculate inintial wall strain
  realtype ew = wall_strain(p,data);

  // Set initial conditions for voigt-model
  Ith(y,1) = kappa*ew; // Should be inintial value of voigt2
  Ith(y,2) = a2/b2*(1.0-kappa)*ew; // Should be inintial value of voigt1

  // Calculate firing
  realtype f, Tpbr, Tsbr;

  f = firing_rate(y, ew, data);

  // Set fs and fp equal to initial firing
  if(data->p[12]==0){
    data->p[12] = f;
    }
  if(data->p[16]==0){
    data->p[16] = f;
    }

  Tpbr = parasymp_baro(f,data);
  Tsbr = symp_baro(f,data);


  // Constants for delayed sympathetic signal 
  int rho;
  realtype scale_constant, delay_alpha, factorial;
  rho = 5;
  factorial = 4*3*2;

  realtype angle, fresp, Tp, Ts;

  if(td <= 0){
    Ith(y,3) = 0;
    Ith(y,4) = 0;
    Ith(y,5) = 0;
    Ith(y,6) = 0;
    Ith(y,7) = 0;
    Ts = Tsbr;
    }
  else{
    delay_alpha = ((realtype) rho)/td;
    scale_constant = pow(delay_alpha,rho)/factorial;

    Ith(y,3) = Tsbr/delay_alpha;
    Ith(y,4) = ((realtype) rho - 4)*Ith(y,3)/delay_alpha;
    Ith(y,5) = ((realtype) rho - 3)*Ith(y,4)/delay_alpha;
    Ith(y,6) = ((realtype) rho - 2)*Ith(y,5)/delay_alpha;
    Ith(y,7) = ((realtype) rho - 1)*Ith(y,6)/delay_alpha;
    Ts = scale_constant*Ith(y,7);
  }

  /*angle = fast_interp(data->t, data->tiltangle,t, data->inputsize);*/
  /*Tresp = fast_interp(data->t, data->resp,t, data->inputsize);*/
  splint(data->t, data->tiltangle,data->tiltangle_spline_coef,data->inputsize,t, &angle);
  splint(data->t, data->resp,data->resp_spline_coef,data->inputsize,t, &fresp);


  Tp = baro_weight(angle,data)*Tpbr + resp_weight(angle,data)*Tresp(fresp,data);

  realtype CA, CN, CAS, CNS, CAF;

  CN = qs*Ts*tN;
  CA = (qp - kiN*CN)*Tp*tA;
  Ith(y,8) = CA;
  Ith(y,9) = CN;

  /* If KA and KN are set to zero, set equal to initial value of CA and CN */
  if(data->p[32]==0){
    data->p[32]=CA;
    KA = data->p[32];
    }
  if(data->p[35]==0){
    data->p[35]=CN;
    KN = data->p[35];
    }


  CAS = (1-mu)*pow(CA,2)/(pow(CA,2) + pow(KA,2));
  CNS = pow(CN,2)/(pow(CN,2) + pow(KN,2));
  Ith(y,10) = CAS;
  Ith(y,11) = CNS;
  CAF = mu*pow(CA,2)/(pow(CA,2) + pow(KA,2));
  Ith(y,12) = 0;

  //Set H0 from initial values
  if(data->p[36]==0){
    realtype hm, hM, hrmean;
    hm=data->p[37]; 
    hM=data->p[38];
    hrmean=data->p[39];

    realtype A,B,C, roota, rootb, h0;
    A = 1;
    B = hM*CNS - hm*(CAS+CAF) - hrmean;  //%max heartbeats per minute * conc. of noradren. at
    C = -hM*hM*CNS*CAF;

    roota = (-B+sqrt(pow(B,2)-4*A*C))/(2*A);
    rootb = (-B-sqrt(pow(B,2)-4*A*C))/(2*A);

    if(-min(-roota,-rootb)<0){
      printf("Error. Calculated h0 <0. Setting h0=100.");
      h0 = RCONST(100.0);
      }
    else{
      if(min(roota,rootb)>0){
        h0 = min(roota,rootb);
        }
      else{
        h0 = -min(-roota,-rootb);
        }
      }
      data->p[36] = h0;

    }

  return(0);
}

/* f routine. Compute f(t,y).  */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data) {
  UserData data;
  data = (UserData) f_data;

  /*// Load rest of parameters*/
  /*realtype p0, k, AM0, a1, a2, b1, b2, s1, s2, Tpm, TpM, xi, fP,  Tsm, TsM, eta, fS, angle_threshold, am, aM, ak, bm, bM, bk, td, tA, qp, kiN, tN, qs, tAF, mu, KA, tAS, tNS, KN, h0, hM, hm;*/
  /*k=data->p[1]               ; AM0=data->p[2]              ; a1=data->p[3]   ; a2=data->p[4]   ; b1=data->p[5]   ;*/
  /*b2=data->p[6]              ; s1=data->p[7]               ; s2=data->p[8]   ; Tpm=data->p[9]  ; TpM=data->p[10] ;*/
  /*xi=data->p[11]             ; fP=data->p[12]              ; Tsm=data->p[13] ; TsM=data->p[14] ; eta=data->p[15] ;*/
  /*fS=data->p[16]             ; angle_threshold=data->p[17] ; am=data->p[18]  ; aM=data->p[19]  ; ak=data->p[20]  ;*/
  /*bm=data->p[21]             ; bM=data->p[22]              ; bk=data->p[23]  ; td=data->p[24]  ; tA=data->p[25]  ;*/
  /*qp=data->p[26]             ; kiN=data->p[27]             ; tN=data->p[28]  ; qs=data->p[29]  ; tAF=data->p[30] ;*/
  /*mu=data->p[31]             ; KA=data->p[32]              ; tAS=data->p[33] ; tNS=data->p[34] ; KN=data->p[35]  ;*/
  /*h0=data->p[36]             ; hm=data->p[37]              ; hM=data->p[38]  ;*/
  
  // Load rest of parameters - here with unused parameters removed
  realtype a1, a2, b1, b2, td, tA, qp, kiN, tN, qs, mu, KA, tAS, tNS, KN;
  a1=data->p[3];
  a2=data->p[4];
  b1=data->p[5];
  b2=data->p[6];
  td=data->p[24]; 
  tA=data->p[25];
  qp=data->p[26];
  kiN=data->p[27];
  tN=data->p[28];
  qs=data->p[29];
  mu=data->p[31];
  KA=data->p[32];
  tAS=data->p[33];
  tNS=data->p[34];
  KN=data->p[35];

  realtype e1, e2, f, Tpbr, Tsbr, fresp, angle, Tp, Ts, p, pbar;

  /* Interpolate pressure data */
  splint(data->t, data->bp,data->bp_spline_coef,data->inputsize,t, &p);

  pbar = Ith(y,0);
  Ith(ydot,0) = p-pbar;

  /* Calculate wall_strain from pressure */
  realtype ew;
  if(data->flags[0]>0){
    ew = wall_strain(pbar,data);
    }
  else{
    ew = wall_strain(p,data);
    }


  e1 = Ith(y,1);
  e2 = Ith(y,2);

  Ith(ydot,1) = -(a1 + a2 + b1)*e1 + (b1 - b2)*e2 + (a1 + a2)*ew;
  Ith(ydot,2) = -a2*e1 - b2*e2 + a2*ew;

  f = firing_rate(y, ew, data);

  Tpbr = parasymp_baro(f,data);
  Tsbr = symp_baro(f,data);

  // Constants for delayed sympathetic signal 
  int rho;
  realtype scale_constant, delay_alpha, factorial;
  rho = 5;
  factorial = 4*3*2;


  // If the delay is 0, don't use the delay implementation.
  if(td <=0){
    Ith(ydot,3) = 0;
    Ith(ydot,4) = 0;
    Ith(ydot,5) = 0;
    Ith(ydot,6) = 0;
    Ith(ydot,7) = 0;
    Ts = Tsbr;
    }
  else{
    delay_alpha = ((realtype) rho)/td;
    scale_constant = pow(delay_alpha,rho)/factorial;
    /*printf("k = %2.2e\n",scale_constant);*/

    Ith(ydot,3) = Tsbr - delay_alpha*Ith(y,3);
    Ith(ydot,4) = ((realtype) rho - 4.0)*Ith(y,3) - delay_alpha*Ith(y,4);
    Ith(ydot,5) = ((realtype) rho - 3.0)*Ith(y,4) - delay_alpha*Ith(y,5);
    Ith(ydot,6) = ((realtype) rho - 2.0)*Ith(y,5) - delay_alpha*Ith(y,6);
    Ith(ydot,7) = ((realtype) rho - 1.0)*Ith(y,6) - delay_alpha*Ith(y,7);
    Ts = scale_constant*Ith(y,7);
    }

  /* Interpolate tiltangle and respiratory data */
  splint(data->t, data->tiltangle,data->tiltangle_spline_coef,data->inputsize,t, &angle);
  splint(data->t, data->resp,data->resp_spline_coef,data->inputsize,t, &fresp);

  if(angle < 0){
    angle = 0;
    }

  // Function call to calculate Tresp

  Tp = baro_weight(angle,data)*Tpbr + resp_weight(angle,data)*Tresp(fresp,data);


  realtype CA, CN, CAS, CNS, CAF, CAT;
  CA = Ith(y,8);
  CN = Ith(y,9);

  Ith(ydot,8) = -CA/tA + (qp - kiN*CN)*Tp;
  Ith(ydot,9) = -CN/tN + qs*Ts;

  CAS = Ith(y,10);
  CNS = Ith(y,11);
  CAF = mu*pow(CA,2)/(pow(CA,2)+pow(KA,2));
  CAT = CAF + CAS;

  Ith(ydot,10) = ((1-mu)*pow(CA,2)/(pow(CA,2)+pow(KA,2)) - CAS)/tAS;
  Ith(ydot,11) = (pow(CN,2)/(pow(CN,2)+pow(KN,2)) - CNS)/tNS;

  realtype h0, hm, hM;

  h0=data->p[36]; hm=data->p[37]; hM=data->p[38];
  Ith(ydot,12) = h0 + hM*CNS - hm*CAT - 1/h0*hM*hm*CNS*CAF;

  return(0);
}

/* Sundials function check if necessary memory is being allocated correctly */
int check_flag(void *flagvalue, char *funcname, int opt){
        int *errflag;

        /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
        if (opt == 0 && flagvalue == NULL) {
                fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
                return(1); 
                }

        /* Check if flag < 0 */
        else if (opt == 1) {
                errflag = (int *) flagvalue;
                if (*errflag < 0) {
                        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
                        return(1); 
                        }
                }

        /* Check if function returned NULL pointer - no memory allocated */
        else if (opt == 2 && flagvalue == NULL) {
                fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
                return(1); 
                }

        return(0);
        }

/* Function to make a fast and dirty interpolation under the assumption that the data is meassured with constant frequency */
realtype fast_interp(realtype *X, realtype *Y, realtype T, int length){
        
        /* We need some variables */
        int j;
        realtype data_step, dt, dp, p;

        /* Find the step size - Assume similar distribution of data throughought the data */
        data_step = X[1]-X[0];

        /* First we check if the time is too large to be interpolated to */
        if(X[length-1]<T){
                printf("# WARNING: The value of T=%f is outside the range of the data (function interp) \n",T);
                return Y[length-1]+(Y[length-1]-Y[length-2])/data_step*(T-X[length-1]);
                }

        /* If T equals the last element of X, return last element of Y */
        if(T==X[length-1]){
                return Y[length-1];
                }

        /* If T equals the first element of X, return first element of Y */
        if(T==X[0]){
                return Y[0];
                }



        /* Find the largest index where the X value is smaller than T */
        j = floorf(T/data_step);

        /* Calculate the distance from the floor value */
        dt = T-j*data_step;

        dp = (Y[j+1]-Y[j])/(X[j+1]-X[j]);

        /* calculate the interpolated value of p */
        p = Y[j] + dt*dp;
        return p;
}

int spline(realtype *x, realtype *y, int datalength, realtype yp1, realtype ypn, realtype *y2){
  /*y2 should have the same length as the datalength - I think... */
  /* Implemented from the recipe in nummerical recipes */
  int i,k;
  
  realtype p,qn,sig,un,*u;

  u = malloc(sizeof(realtype)*datalength);
  
  //for(i=0;i<datalength;i++){
  //  printf("%5.5e\t%5.5e\n", x[i], y[i]);
  //}

  /* The lower boundary condition is set either to be "natural" or else to have specified first derivative */
  if (yp1 > RCONST(0.99e30)) {
    y2[0]=RCONST(0.0);
    u[0]=RCONST(0.0);
  }
  else {
    y2[0] = -RCONST(0.5);
    u[0] = (RCONST(3.0)/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  /* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors */
  for (i=1; i<=datalength-2; i++){
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+RCONST(2.0);
    y2[i] = (sig-RCONST(1.0))/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (RCONST(6.0)*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  /* The upper boundary condition is set either to be "natural or else to have a specified first derivative. */
  if (ypn > RCONST(0.99e30)) {
    qn = un = RCONST(0.0);
  }
  else {
    qn = RCONST(0.5);
    un = (RCONST(3.0)/(x[datalength-1]-x[datalength-2]))*(ypn-(y[datalength-1]-y[datalength-2])/(x[datalength-1]-x[datalength-2]));
  }
  y2[datalength-1] = (un - qn*u[datalength-2])/(qn*y2[datalength-2]+RCONST(1.0));

  /* This is the backsubstitution loop of the tridiagonal algorithm */
  for (k=datalength-2; k>=0; k--){
    y2[k] = y2[k]*y2[k+1]+u[k];
  }
  free(u);
  return(0);
}

void splint(realtype *xa, realtype *ya, realtype *y2a, int datalength, realtype x, realtype *y){
  /* Implemented from numerical recipes */
  
  int klo, khi, k;
  realtype h, b, a;

  /* We will find the right place in the table by means of bisection.  This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous values of klo and khi and test if they remain appropriate on the next call. */
  klo = 0;
  khi = datalength-1;
  
  while(khi-klo>1){
    k=(khi+klo) >> 1;
    if (xa[k] > x){
      khi = k;
    }
    else {
      klo = k;
    }
  }
 
  //printf("klo = %i and khi = %i\n", klo, khi); 
  /* klo and khi now bracket the input value of x */
  h=xa[khi]-xa[klo];

  /* The xa's must be distinct */
  if(h == RCONST(0.0)){
    printf("Bad XA input to routine splint");
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;

  /* Cubic spline polynomial is now evaluated */
  //y[0] = a*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
  *y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
 // printf("\n\n%5.5e\n\n",a*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0));
}

/* Function to calculate mean value of input pressure for some initial period */
realtype mean_start_pressure(void *f_data, realtype time){
  UserData data;
  data = (UserData) f_data;
  realtype sum=0;
  int i,n;
  realtype dt = data->t[1]-data->t[0];
  n = (int) time/dt;
  for(i=0;i<n;i++){
  sum = sum + data->bp[i];
  }
  return sum/( (realtype) n);
}

/* Function to calculate wall strain from pressure */
realtype wall_strain(realtype p, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype p0, k, AM0;
  p0 = data->p[0]; //Needs to call function to calculate p0 from input pressure data. Maybe make it dependent on if it is given as a parameter.
  k = data->p[1];
  AM0 = data->p[2];

  return(1.0-pow((pow(p0,k) +pow(p,k))/(pow(p0,k) + AM0*pow(p,k)),(.5)));
  }

/* Function to calculate firing rate from baroreceptor strain*/
realtype firing_rate(N_Vector y, realtype ew, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype s1 = data->p[7], s2 = data->p[8];
  return s1*(ew-Ith(y,1))+s2;
}

/* Function to calculate Parasympathetic tone from baroreceptors */
realtype parasymp_baro(realtype firing, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype Tpm = data->p[9], xi = data->p[11], fP = data->p[12]; 
  return Tpm + (1.0 -Tpm)*pow(firing,xi)/(pow(firing,xi) + pow(fP,xi));
}

/* Function to calculate Sympathetic tone from baroreceptors */
realtype symp_baro(realtype firing, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype Tsm = data->p[13], eta = data->p[15], fS = data->p[16];
  return 1.0 - (1.0 - Tsm)*pow(firing,eta)/(pow(firing,eta) + pow(fS,eta));
}

/* Function to calculate weight of parasympathetic tone from baroreceptors */
realtype baro_weight(realtype angle, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype angle_threshold, am, aM, ak;
  angle_threshold=data->p[17] ; am=data->p[18]  ; aM=data->p[19]  ; ak=data->p[20]  ;
  return am + (aM-am)*pow(sin(angle),ak)/(pow(sin(angle),ak) + pow(sin(angle_threshold),ak));
}

/* Function to calculate weight of parasympathetic tone from respiration */
realtype resp_weight(realtype angle, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype angle_threshold, bm, bM, bk;
  angle_threshold=data->p[17] ; bm=data->p[21] ; bM=data->p[22]; bk=data->p[23];
  return bM - (bM-bm)*pow(sin(angle),bk)/(pow(sin(angle),bk) + pow(sin(angle_threshold),bk));
}

/* Function to calculate heart rate */
realtype heart_rate(N_Vector y, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype h0, hm, hM, mu, KA, CA, CAF, CAS, CAT, CNT;

  h0=data->p[36]; hm=data->p[37]; hM=data->p[38];
  mu=data->p[31]; KA=data->p[32];
  CA = Ith(y,8);
  CAS = Ith(y,10);
  CAF = mu*pow(CA,2)/(pow(CA,2)+pow(KA,2));
  CAT = CAF + CAS;
  CNT = Ith(y,11);
  return h0 + hM*CNT - hm*CAT - 1/h0*hM*hm*CNT*CAF;
}

/* Function to calculate T_resp */
realtype Tresp(realtype fresp, void *f_data){
  UserData data;
  data = (UserData) f_data;
  realtype fresp_thres = data->p[40];
  realtype fresp_k = data->p[41];
  return RCONST(1.0) - pow(fresp,fresp_k)/(pow(fresp,fresp_k)+pow(fresp_thres,fresp_k));
}

/* The model function */
int model(int input_size, double * input_time , double * input_bp, double * input_tiltangle, double * input_respiration, int output_size, double * output_time, double * output_hr, int parameter_size, double * parameter, double * flags){

  void *cvode_mem;
  UserData data;
  realtype t, tout;
  N_Vector y;
  int flag;

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;

  int NEQ = 13;
  realtype p;

  // Load data into data structure
  data = (UserData) malloc(sizeof *data);
  data->no_parameters = parameter_size;
  data->p = parameter;
  data->inputsize= input_size;
  data->t = input_time;
  data->bp = input_bp;
  /*data->bp_spline_coef = input_bp_spline_coef;*/
  /*data->tiltangle_spline_coef = input_tiltangle_spline_coef;*/
  /*data->resp_spline_coef = input_respiration_spline_coef;*/
  data->bp_spline_coef = malloc(input_size*sizeof(realtype));
  data->tiltangle_spline_coef = malloc(input_size*sizeof(realtype));
  data->resp_spline_coef = malloc(input_size*sizeof(realtype));
  data->flags = malloc(1*sizeof(realtype));
  data->tiltangle = input_tiltangle;
  data->resp = input_respiration;

  if(data->p[0]==0){
    data->p[0] = mean_start_pressure(data, 10.0);
    }


  spline(data->t,data->bp,data->inputsize,RCONST(1e30),RCONST(1e30),data->bp_spline_coef);
  spline(data->t,data->tiltangle,data->inputsize,RCONST(1e30),RCONST(1e30),data->tiltangle_spline_coef);
  spline(data->t,data->resp,data->inputsize,RCONST(1e30),RCONST(1e30),data->resp_spline_coef);

  /*data->p[0] = input_bp[0];*/
  realtype Rtol=1e-6;
  realtype Atol=1e-6;

  // Initialize state variables
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  int i,j;

  for(i=0;i<1;i++){
    data->flags[i] = flags[i];
    }

  // Set initial values.
  // The function calculates values that would be at steady state with the initial pressure input.
  set_initial_values(0, y, data);

  /*[> Create CVODES object <]*/
  cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /*[> Allocate space for CVODES <]*/
  flag = CVodeInit(cvode_mem, f, RCONST(0.0), y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /*[> Use private function to compute error weights <]*/
  flag = CVodeSStolerances(cvode_mem, Rtol, Atol);
  if (check_flag(&flag, "CVodeSetEwtFn", 1)) return(1);

  /*[> Attach user data <]*/
  flag = CVodeSetUserData(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /*[> Attach linear solver <]*/
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "cvdense", 1)) return(1);


  /* We save the initial values to the solution array, and the initial timepoint as well. The +2 is there because we want to store time as well as the state solution, and also make room for one more time series calculated post ode. */
  output_hr[0] = Ith(y,0);
  output_hr[0+output_size] = Ith(y,1);
  output_time[0] = input_time[0];
  for(j = 0; j<NEQ;j++){
    output_hr[j*output_size] = Ith(y,j);
  }
  /*output_hr[NEQ*output_size] = heart_rate(y, data);*/

  /* Loop over output points, call CVode, save results, test for error */
  for (i=1; i< output_size; i++){
          /* Call the integrator, and check the flag */
          tout = output_time[i]; //data->t[i];
          flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
          //if (check_flag(&flag, "CVode", 1)) break;
          if (check_flag(&flag, "CVode", 1)) return(1);

          /* We save the time, and the results for each equation to the output array */
          output_time[i]=t;
          for(j=0; j<NEQ;j++){
            output_hr[i+j*output_size]=Ith(y,j);
            }
          /*output_hr[i+NEQ*output_size]=heart_rate(y,data);*/

          splint(data->t, data->bp,data->bp_spline_coef,data->inputsize,t, &p);
          output_hr[i+(NEQ)*output_size]=wall_strain(p,data);
  }


N_VDestroy_Serial(y);
CVodeFree(&cvode_mem);
free(data->bp_spline_coef);
free(data->tiltangle_spline_coef);
free(data->resp_spline_coef);
free(data);
return(0);
}

/* Function to calculate sensitivities */
int sens(int input_size, double * input_time , double * input_bp, double * input_tiltangle, double * input_respiration,  int output_size, double * output_time, double * output_hr, int parameter_size, double * parameter, double * output_sensitivities, double * flags){

  void *cvode_mem;
  UserData data;
  realtype t, tout;
  N_Vector y;
  int flag;

  realtype *pbar;
  N_Vector *yS;
  int is, *plist, i, j, k;
  booleantype sensi, err_con;
  int sensi_meth;

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
  yS        = NULL;
  plist     = NULL;
  pbar      = NULL;
  sensi     = TRUE;

  /* Some settings */
  err_con = TRUE;
  sensi_meth = CV_SIMULTANEOUS;

  int NEQ = 13;

  // Load data into data structure
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2)) return(1);

  data->no_parameters = parameter_size;
  data->p = parameter;
  data->inputsize= input_size;
  data->t = input_time;
  data->bp = input_bp;
  data->bp_spline_coef = malloc(input_size*sizeof(realtype));
  data->tiltangle_spline_coef = malloc(input_size*sizeof(realtype));
  data->resp_spline_coef = malloc(input_size*sizeof(realtype));
  data->tiltangle = input_tiltangle;
  data->resp = input_respiration;
  data->flags = malloc(1*sizeof(realtype));

  data->p[0] = mean_start_pressure(data, 10.0);

  spline(data->t,data->bp,data->inputsize,RCONST(1e30),RCONST(1e30),data->bp_spline_coef);
  spline(data->t,data->tiltangle,data->inputsize,RCONST(1e30),RCONST(1e30),data->tiltangle_spline_coef);
  spline(data->t,data->resp,data->inputsize,RCONST(1e30),RCONST(1e30),data->resp_spline_coef);

  /*data->p[0] = input_bp[0];*/
  realtype Rtol=1e-8;
  realtype Atol=1e-8;

  realtype *sens_tol;
  sens_tol = malloc(sizeof(realtype)*data->no_parameters);
  for(i=0;i<data->no_parameters;i++){
    sens_tol[i]=Atol;
  }

  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  // Save flags. First one is if input should be smoothed (1) or no (0)
  for(i=0;i<1;i++){
    data->flags[i] = flags[i];
    }

  // Set initial values.
  // The function calculates values that would be at steady state with the initial pressure input.
  set_initial_values(0, y, data);

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Allocate space for CVODES */
  flag = CVodeInit(cvode_mem, f, RCONST(0.0), y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Use private function to compute error weights */
  flag = CVodeSStolerances(cvode_mem, Rtol, Atol);
  if (check_flag(&flag, "CVodeSetEwtFn", 1)) return(1);

  /* Attach user data */
  flag = CVodeSetUserData(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  realtype *sdata;

  if(sensi){


  plist = (int *) malloc(data->no_parameters * sizeof(int));
  if(check_flag((void *)plist, "malloc", 2)) return(1);
  for(is=0; is<data->no_parameters; is++) plist[is] = is;

  pbar  = (realtype *) malloc(data->no_parameters * sizeof(realtype));
  if(check_flag((void *)pbar, "malloc", 2)) return(1);
  for(is=0; is<data->no_parameters; is++) pbar[is] = data->p[plist[is]];

  /*for(j=0; j<data->no_parameters;j++){*/
    /*printf("%2.2f\n",data->p[j]);*/
  /*}*/

  yS = N_VCloneVectorArray_Serial(data->no_parameters, y);
  if(check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
  for(is=0;is<data->no_parameters;is++)
    N_VConst(RCONST(0.0), yS[is]);

  flag = CVodeSensInit1(cvode_mem, data->no_parameters, sensi_meth, NULL, yS);
  if(check_flag(&flag, "CVodeSensInit1", 1)) return(1);

  flag = CVodeSensEEtolerances(cvode_mem);
  if(check_flag(&flag, "CVodeSensEEtolerances", 1)) return(1);

  flag = CVodeSetSensErrCon(cvode_mem, err_con);
  if(check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);

  flag = CVodeSensSStolerances(cvode_mem, Rtol, sens_tol);
  if (check_flag(&flag, "CVodeSensSStolerances", 1)) return(1);

  flag = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, RCONST(0.0));
  if(check_flag(&flag, "CVodeSetSensDQMethod", 1)) return(1);

  flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
  if(check_flag(&flag, "CVodeSetSensParams", 1)) return(1);
  }



printf("Ready to go.\n");


  /*We save the initial values to the solution array, and the initial timepoint as well. The +2 is there because we want to store time as well as the state solution, and also make room for one more time series calculated post ode. */
  output_hr[0] = Ith(y,0);
  output_hr[0+output_size] = Ith(y,1);
  output_time[0] = input_time[0];
  for(j = 0; j<NEQ;j++){
    output_hr[j*output_size] = Ith(y,j);
  }
  /*output_hr[NEQ*output_size] = heart_rate(y, data);*/

   /*Loop over output points, call CVode, save results, test for error */
  for (i=1; i< output_size; i++){
           /*Call the integrator, and check the flag */
          tout = output_time[i]; //data->t[i];
          flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
          //if (check_flag(&flag, "CVode", 1)) break;
          if (check_flag(&flag, "CVode", 1)) return(1);

          output_time[i]=t;
          for(j=0; j<NEQ;j++){
            output_hr[i+j*output_size]=Ith(y,j);
            }
            

          if(sensi){
            flag = CVodeGetSens(cvode_mem, &t, yS);
            if(check_flag(&flag, "CVodeGetSens", 1)) break;

            /* Save sensitivities */
            for(j=0;j<data->no_parameters;j++){
              sdata = NV_DATA_S(yS[j]);
              for(k=0;k<NEQ;k++){
                output_sensitivities[i+j*output_size+k*parameter_size*output_size] = sdata[k];
                }
              }
          }

          /*flag = CVodeGetSens(cvode_mem, &t, yS);*/
          /*if (check_flag(&flag, "CVodeGetSens", 1)) return(1);*/

          /* Here I need to take the sensitivities out */
  }

printf("Finished.\n");

/*N_VDestroy_Parallel(y);*/
N_VDestroy_Serial(y);
if(sensi){
  N_VDestroyVectorArray_Serial(yS, parameter_size);
  free(plist);
  free(pbar);
}
CVodeFree(&cvode_mem);
free(data);
return(0);
}

int spline_evaluation(realtype *x, realtype *y, realtype *coef, int in_datalength, realtype *interp_x, realtype *interp_y, int out_datalength){
int i;
realtype temp;
for(i=0;i<out_datalength;i++){
  splint(x,y,coef,in_datalength,interp_x[i],&temp);
  interp_y[i] = temp;
  }

return 0;
}

realtype min(realtype a, realtype b){
  if(a>b){
    return b;
    }
  else
    return a;
    }
