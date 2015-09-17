#include "functions.h"

#define NR_END 1
#define FREE_ARG char*

/* Functions for printing error message if wrongs arguments are paased */
void WrongArgs(char *name){

  printf("\nUsage: %s -inputdata <data filename> -inputhbtimes <heartbeattime data filename> -outputfile <output filename> -pars <parameter file> -initalvalues <initial values file> -settings <settings file> \n", name);
  printf("\n\n Inputdata file should contain two columns of data, time in column one, and blood pressure in column two.\n Inputhbtimes file should contain one column of times corresponding to the times of heart contraction.\n");
  exit(0);
  };
/* This function will open all the files required for the program */
void ProcessArgs(int argc, char *argv[], FILE **input, FILE **input2, FILE **output, FILE **parfile, FILE **settings, FILE **initialvaluesfile){
    /* If there is not enough arguments, make an error */
    if(argc != 13){ 
      WrongArgs(argv[0]);
    };
  
    int i;

  for(i=1; i <= argc-2; i=i+2)
  {
      /* If argument i = "-inputdata: */
      if(strcmp(argv[i],"-inputdata") == 0)
      {
        *input = fopen(argv[i+1],"r"); 
        if(*input==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input data file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as inputdata.\n", argv[i+1]);
      }
      
      /* If argument i = "-inputhbtimes: */
      else if(strcmp(argv[i],"-inputhbtimes") == 0)
      {
        *input2 = fopen(argv[i+1],"r"); 
        if(*input2==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input heartbeat times data file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
        else{
//        printf("Opened file %s as hearbeattimes data file.\n", argv[i+1]);
        }
      }

      /* If argument i = "-outputfile: */
      else if(strcmp(argv[i],"-outputfile") == 0)
      {
        *output = fopen(argv[i+1],"w");
//        printf("Opened file %s as outputfile.\n", argv[i+1]);
      }
      
      /* If argument i = "-inputdata: */
      else if(strcmp(argv[i],"-pars") == 0)
      {
        *parfile = fopen(argv[i+1],"r");
        if(*parfile==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input parameter file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as parameter file.\n", argv[i+1]);
      }

      /* If argument i = "-inputdata: */
      else if(strcmp(argv[i],"-settings") ==0)
      {
        *settings = fopen(argv[i+1],"r");
        if(*settings==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input settings file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as settings file.\n", argv[i+1]);
      }

      /* If argument i = "-initialvalues: */
      else if(strcmp(argv[i],"-initialvalues") ==0)
      {
        *initialvaluesfile = fopen(argv[i+1],"r");
        if(*initialvaluesfile==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Initialvalues file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as initial values file.\n", argv[i+1]);
      }

      /* If the argument is something else, end the program */
      else
      {
        WrongArgs(argv[0]);
      }
  }
 // printf("Finished loading input files\n");
  }; 
/* Counts the number of lines in input file*/
int LineCount(FILE *f){

        /*  */
        char c;
        int lines = 0;

        /* As long a the char we are reading is the end of the file */
        while((c = fgetc(f)) != EOF)

        if(c == '\n')
                lines++;

        if(c != '\n')
                lines++;

        return lines;
        }
/* Extract data */
//int DataExtract(FILE **f, int datalength, realtype *time, realtype *data1, realtype *data2){
//
//  int i;
//   
//  /* Loop through all lines and use scanf to extract data */
//  for(i=0; i<datalength; i++){
//    fscanf(f, "%e %e %e\n", time[i], data1[i], data2[i]); 
//  }
//return 0;
//}


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/* Fetches the column number spec_col, from the file f which contains data formatted in no_cols columns, and saves it in the array pdata */
int ExtractColumn(FILE *f, int number_of_columns, int specified_column, int datalength, realtype *data){

        /* We need a few variables */
        realtype number;
        int i = 0,  anothernumber;

        /* To use float instead of realtype to read the numbers */
        // float number;


        /* Read from the start of the file */
        rewind (f);

        /* As long as i isnt larger than the number of data entries */
        while(i<number_of_columns*datalength){ 
                
                /* Scan the file for numbers of %le format, and save to &number */
                fscanf(f,"%le",&number);

                /* To use float instead of realtype uncomment this line below */
                // fscanf(f,"%le",&number);

                /* if the modulo of i with regards to the number of columns in the file matches the column we are interessted in */
                if(i%number_of_columns==specified_column){ 

                        /* We find what line number we are at, both i and no_cols are integers */
                        anothernumber = i / number_of_columns;

                        /* We save what we read */
                        data[anothernumber] = number;
                        }

                i++;
                }
        return 0;
        }

/* Function taken from the CVODE examples. It tests i memory is allocated correctly for the buildin stuff */
realtype  LinearInterpolation(realtype *X, realtype *Y, realtype T, int length){
        
        /* We need some variables */
        int j;
        realtype data_step, dt, dp, p;

        /* First we check if the time is too large to be interpolated to */
        if(X[length-1]<T){
                printf("# ERROR: The value of T=%f is outside the range of the data (function interp) \n",T);
                return 0;
                }

        /* If T equals the last element of X, return last element of Y */
        if(T==X[length-1]){
                return Y[length-1];
                }

        /* If T equals the first element of X, return first element of Y */
        if(T==X[0]){
                return Y[0];
                }


        /* Find the step size - Assume similar distribution of data throughought the data */
        data_step = X[1]-X[0];

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

void splint(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *y){
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

int splintdotequalspace(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *ydot){
  /* Implemented from numerical recipes */
  
  int klo, khi;
  realtype h, b, a;

  /* We will find the right place in the table by means of bisection.  This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous values of klo and khi and test if they remain appropriate on the next call. */
  klo = 0;
  khi = datalength-1;
 
  /* First we check if the time is too large to be interpolated to */
  if(xa[datalength-1]<x){
          //printf("# ERROR: The value of T=%f is outside the range of the data (function interp) \n",T);
          return ya[datalength-1];
          }

  /* If T equals the last element of X, return last element of Y */
  if(x==xa[datalength-1]){
          return ya[datalength-1];
          }

  /* If T equals the first element of X, return first element of Y */
  if(x==xa[0]){
          return ya[0];
          }


  /* Find the step size - Assume similar distribution of data throughought the data */
  h = (xa[datalength-1]-xa[0])/datalength; 
  //h = X[1]-X[0];

  /* Find the largest index where the X value is smaller than T */
  klo = floorf((x-xa[0])/h);
 khi = klo+1;
 
  //printf("klo = %i and khi = %i\n", klo, khi); 
  /* klo and khi now bracket the input value of x */

  /* The xa's must be distinct */
  if(h == RCONST(0.0)){
    printf("Bad XA input to routine splint");
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;

  /* Cubic spline polynomial is now evaluated */
  //y[0] = a*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
  //*y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
  //*y = (ya[khi]-ya[klo])/(xa[khi]-x) - ((RCONST(3.0)*a*a-RCONST(1.0))/RCONST(6.0))*(xa[khi]-xa[klo])*y2a[klo] +
  *ydot = (ya[khi]-ya[klo])/(xa[khi]-xa[klo]) - (h/RCONST(6.0))*((RCONST(3.0)*a*a-RCONST(1.0))*y2a[klo] - (RCONST(3.0)*b*b-RCONST(1.0))*y2a[khi]);
  //*ydot = (ya[khi]-ya[klo])/(a*h) - (h/RCONST(6.0))*((RCONST(3.0)*a*a-RCONST(1.0))*y2a[klo] - (RCONST(3.0)*b*b-RCONST(1.0))*y2a[khi]);
        
 // printf("\n\n%5.5e\n\n",a*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0));
 return(0);
}

int splintdot(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *ydot){
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
  //*y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
  //*y = (ya[khi]-ya[klo])/(xa[khi]-x) - ((RCONST(3.0)*a*a-RCONST(1.0))/RCONST(6.0))*(xa[khi]-xa[klo])*y2a[klo] +
  *ydot = (ya[khi]-ya[klo])/(xa[khi]-xa[klo]) - (h/RCONST(6.0))*((RCONST(3.0)*a*a-RCONST(1.0))*y2a[klo] - (RCONST(3.0)*b*b-RCONST(1.0))*y2a[khi]);
  //*ydot = (ya[khi]-ya[klo])/(a*h) - (h/RCONST(6.0))*((RCONST(3.0)*a*a-RCONST(1.0))*y2a[klo] - (RCONST(3.0)*b*b-RCONST(1.0))*y2a[khi]);
        
 // printf("\n\n%5.5e\n\n",a*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0));
 return(0);
}

void splintequal(realtype xa[], realtype ya[], realtype y2a[], int datalength, realtype x, realtype *y){
  /* Implemented from numerical recipes */
  
  int klo, khi;
  realtype h, b, a;

  /* We will find the right place in the table by means of bisection.  This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous values of klo and khi and test if they remain appropriate on the next call. */
  klo = 0;
  khi = datalength-1;
 
  h = xa[1]-xa[0]; 
  /* Find the largest index where the X value is smaller than T */
  klo = floorf(x/h);
  khi = klo+1;
 
  //printf("kha= %i and khi=%i\n",klo,khi); 
  /* klo and khi now bracket the input value of x */
 // h=xa[khi]-xa[klo];

  /* The xa's must be distinct */
  if(h == RCONST(0.0)){
    printf("Bad XA input to routine splint");
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  //printf("a = %5.5e and b=%5.5e\n",a,b);
  /* Cubic spline polynomial is now evaluated */
  *y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/RCONST(6.0);
}

int printoutput(FILE *outputfile, realtype *time, int datalength, realtype *output){
  int i;
#if defined(SUNDIALS_EXTENDED_PRECISION)
    for(i=0;i<datalength;i++){
      fprintf(outputfile, "%0.4Le  %14.6Le  %14.6Le  %14.6Le\n", time[i], IJth(output,datalength,i,0), IJth(output,datalength,i,1), IJth(output,datalength,i,2));
    }
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    for(i=0;i<datalength;i++){
      fprintf(outputfile, "%0.4le  %14.6le  %14.6le  %14.6le\n", time[i], IJth(output,datalength,i,0), IJth(output,datalength,i,1), IJth(output,datalength,i,2));
    }
#else
    for(i=0;i<datalength;i++){
      fprintf(outputfile, "%0.4le  %14.6e  %14.6e  %14.6e\n", time[i], IJth(output,datalength,i,0), IJth(output,datalength,i,1), IJth(output,datalength,i,2));
    }
#endif

  return 0;
}

/* After the differential equation has been solved, the solution has to be processed to compare to data. This code is copied directly from the Matlab function wrapper. */
int AfterSolve(realtype *output, UserData data){
 
  int i; 

  realtype *V, *C, *pic, *Ga, *pc, *qsim, fact, kR, Gv;

  fact = RCONST(0.25);
  kR = data->p[8];
  Gv = data->p[5];

  C = malloc(sizeof(C)*data->datalength);
  pic = malloc(sizeof(pic)*data->datalength);
  V = malloc(sizeof(V)*data->datalength);
  Ga = malloc(sizeof(Ga)*data->datalength);
  pc = malloc(sizeof(pc)*data->datalength);
  qsim = malloc(sizeof(qsim)*data->datalength);

  /* For each datapoint, calculate the output */ 
  for(i=0;i<data->datalength;i++){
    pic[i] = DoubleArrayAccess(output,data->no_equations+1,i,1);
    C[i] = DoubleArrayAccess(output,data->no_equations+1,i,2);
    V[i] = C[i]*(data->data[i]-pic[i]);
    Ga[i] = V[i]*V[i]/kR;
    pc[i] = (data->data[i]*Ga[i]+pic[i]*Gv)/(Ga[i]+Gv);
    qsim[i] = Ga[i]*(data->data[i]-pc[i]);

    /* Write the output to the array we store the results in */
    DoubleArrayAccess(output,data->no_equations+1,i,data->no_equations+1) = qsim[i]/fact;
  }   

  /* Free all the allocated memory */
  free(C);
  free(pic);
  free(V);
  free(Ga);
  free(pc);
  free(qsim);

  return(0);
}

/* This function will print the results to the outputfile */
int PrintSensToFile(FILE *outputfile, realtype **current_solution, UserData data){
  int i,j;

  /* For each timepoint */
  for(i=0;i<data->hbdatalength-1;i++){

    /* For each parameter */
    for(j=0;j<data->no_parameters+1;j++){

      /* Print the result */
      fprintf(outputfile,"%8.8e  ", current_solution[j][i]);
         // DoubleArrayAccess(current_solution,(data->no_parameters+1),i,j));
    }

    /* Add a new line before next time point */
    fprintf(outputfile,"\n");
  }
  return(0);
}

int PrintToFile(FILE *outputfile, realtype *current_solution, UserData data){
  int i,j;

  /* For each timepoint */
  for(i=0;i<data->hbdatalength;i++){

    /* For each equation */
    for(j=0;j<data->no_equations+1;j++){

      /* Print the result */
      fprintf(outputfile,"%8.8e  ",DoubleArrayAccess(current_solution,(data->no_equations+1),i,j));
    }

    /* Add a new line before next time point */
    fprintf(outputfile,"\n");
  }
  return(0);
}

/* Will return the likelihood using exp(-error) as likelihood function*/
//int likelihood(realtype *array1, realtype *array2, int data_length, realtype *result, realtype sigma, int stepno10s){
int likelihood(realtype *array1, realtype *array2, realtype error_value, realtype *result, realtype sigma){
  realtype temp_result, smallest_number;


        /* The smallest number double float can use - http://www.dummies.com/how-to/content/determining-types-of-numbers-in-c.html */
        smallest_number = RCONST(1.7e-308);

        temp_result = exp(-error_value/(2*sigma*sigma));

        /* Check if the results will be rounded to zero */
        if(temp_result<smallest_number){
                temp_result = smallest_number;
                printf("#The likelyhood was tiny, so it was rounded to zero.\n");
        }

        *result = temp_result;

        return 0;
        }

/* Will return the sum of the square of the difference between the data points in the two arrays */
realtype error(realtype *array1, realtype *array2, int data_length, int stepno10s){
        realtype sum;
        int i, diff;
        sum = RCONST(0.0);
	diff = data_length-stepno10s;
        for(i=0; i<diff; i++){
                // if(i%500==0){
                        // printf("Test. Value of proposal is: %12.4le. Value of data is: %12.4le. Total error is before this point: %12.4le\n", array1[i], array2[i], sum);
                // }
                sum = sum + (array1[i+stepno10s]- array2[i+stepno10s])*(array1[i+stepno10s]- array2[i+stepno10s]);
        }

        // printf("SUM: %12.4le", sum);
        return sum;
        }

int generate_proposed_pars(realtype *current_par, realtype proposed_par[], realtype *sigma, realtype *minimum_value, realtype *maximum_value, int no_of_pars, gsl_rng *r){
        int i;
        realtype temp;

        /* set up GSL RNG */
        // gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
        /* end of GSL setup */

            for(i=0;i<no_of_pars;i++){
                temp=RCONST(-1.0);

                /* Only update the parameter if the value is larger than one */
               while(temp<minimum_value[i] || temp>maximum_value[i]){
                       temp = current_par[i] + gsl_ran_gaussian(r, sigma[i]);
               }
                proposed_par[i] = temp;
                }

        return 0;
        }

int next_mc_link(realtype *proposed, realtype *current_parameter, realtype *proposed_error, realtype *current_error, int no_of_pars, gsl_rng *r, realtype sigma){
        realtype alpha, likelihood_ratio, random_number;
        int j;

        /* Calculate the ratio of the likelihood of the data given the new parameters, to the likelihood of the data given the old parameters. */
        likelihood_ratio = exp(-RCONST(0.5)*(proposed_error-current_error)/(2*sigma*sigma));

        /* choose alpha as the minimum of this ratio and one */
        alpha = fmin(RCONST(1.0),likelihood_ratio);
        random_number = gsl_rng_uniform(r);

        // printf("# ---------------------- \n# Alpha: %12.4le                random_numer: %12.4le\n# ----------------------\n", alpha, random_number);
        if(alpha>random_number){
                for(j=0;j<no_of_pars;j++){
                        current_parameter[j]=proposed[j];
                }
                current_error[0]=proposed_error[0];
        }
        else{
        }

        /* if a uniformly random drafted number [0-1] is smaller than alpha, accept the proposed parameters */
        /* else keep current_parameter as the parameter value */

        return 0;
        }

int PrintParametersToFile(FILE *outputfile, realtype *current_pars, int number_of_parameters, realtype error_value){
    int i;
    for(i=0;i<number_of_parameters;i++){
      fprintf(outputfile,"%5.5e ", current_pars[i]);
    }
    fprintf(outputfile,"%5.5e ", error_value);
    fprintf(outputfile,"\n");
    return 0;
}

int calculate_hr_from_potential(realtype *result, realtype *times, realtype *potential, realtype datalength){
  int i;
  //result[0] = RCONST(0.0);
  for(i=0;i<datalength-1; i++){
    result[i] = (potential[i+1]-potential[i])/(times[i+1]-times[i]);
    //result[i] = (potential[i+1]-potential[i]);
  };
  return 0;
}

int calculate_hr_from_times(realtype *result, realtype *times, realtype datalength){
  int i;
  //result[0] = RCONST(0.0);
  for(i=0; i<datalength-1; i++){
    result[i] = RCONST(1.0)/(times[i+1]-times[i]);
    //result[i] = (times[i+1]-times[i]);
  };
  return 0;
}

realtype  fast_interp(realtype *X, realtype *Y, realtype T, int length){
        
        /* We need some variables */
        int j;
        realtype data_step, dt, dp, p;

        /* First we check if the time is too large to be interpolated to */
        if(X[length-1]<T){
                printf("# ERROR: The value of T=%f is outside the range of the data (function interp) \n",T);
                return Y[length-1];
                }

        /* If T equals the last element of X, return last element of Y */
        if(T==X[length-1]){
                return Y[length-1];
                }

        /* If T equals the first element of X, return first element of Y */
        if(T==X[0]){
                return Y[0];
                }


        /* Find the step size - Assume similar distribution of data throughought the data */
        data_step = X[1]-X[0];

        /* Find the largest index where the X value is smaller than T */
        j = floorf(T/data_step);

        /* Calculate the distance from the floor value */
        dt = T-j*data_step;

        dp = (Y[j+1]-Y[j])/(X[j+1]-X[j]);

        /* calculate the interpolated value of p */
        p = Y[j] + dt*dp;
        return p;
}

int Sensitivity(void *cvode_mem, UserData data, realtype y0[], realtype **output, realtype **outputsol, realtype Rtol, realtype Atol, realtype T0, realtype T1){

        int i,j,*plist,is;

        plist = NULL;

        realtype *pbar;
        pbar = NULL;

        realtype *sens_tol;
        sens_tol = malloc(sizeof(realtype)*data->no_parameters);
        for(i=0;i<data->no_parameters;i++){
          sens_tol[i]=Atol;
        }

        /* The sundials solver needs the y-values in an N_Vector */
        N_Vector y,*yS;
        y = NULL;
        y = N_VNew_Serial(data->no_equations);
        if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

        yS = NULL;
        yS = N_VCloneVectorArray_Serial(data->no_parameters, y);


        int sensi_meth;
        sensi_meth = CV_SIMULTANEOUS;

        booleantype err_con;
        err_con = TRUE;

        realtype *sdata;

        /* We need the initial values read into the y-vector at the beginning */
        for(i=0;i<data->no_equations;i++){
          Ith(y,i) = y0[i];
        }

        /* flag is to store the results of cvode_mem. iout is to loop over all timesteps */
        int flag;

        realtype t,tout;
        tout = T1;

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
        flag = CVDense(cvode_mem, data->no_equations);
        if (check_flag(&flag, "CVDense", 1)) return(1);


        plist = (int *) malloc(data->no_parameters * sizeof(int));
        if(check_flag((void *)plist, "malloc", 2)) return(1);
        for(is=0; is<data->no_parameters; is++) plist[is] = is;

        pbar  = (realtype *) malloc(data->no_parameters * sizeof(realtype));
        if(check_flag((void *)pbar, "malloc", 2)) return(1);
        for(is=0; is<data->no_parameters; is++) pbar[is] = data->p[plist[is]];

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

        
        /* We save the initial values to the solution array, and the initial timepoint as well. The +2 is there because we want to store time as well as the state solution, and also make room for one more time series calculated post ode. */
        for(j=0;j<data->no_equations;j++){
          outputsol[j+1][0]= y0[j];
        }

        /* Initial time point */
        outputsol[0][0]= data->t[0];
        output[0][0] = data->t[0];

        /* Save initial values for sensitivities */
        for(j=0;j<data->no_parameters;j++){
          output[j+1][0] = RCONST(0.0);
        }

        /* Loop over output points, call CVode, save results, test for error */
        for (i=1; i<data->hbdatalength; i++) {
                /* Call the integrator, and check the flag */
                flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
                //if (check_flag(&flag, "CVode", 1)) break;

                // We want to stop this iterations if the solver breaks.
                if (check_flag(&flag, "CVode", 1)==1){
                 return 1;
                }

                /* We save the time, and the results for each equation to the output array */

                flag = CVodeGetSens(cvode_mem, &t, yS);
                if (check_flag(&flag, "CVodeGetSens", 1)) break;

                /* Save time */
                output[0][i] = tout;
                outputsol[0][i] = tout;

                /* Save sensitivities */
                for(j=0;j<data->no_parameters;j++){
                  sdata = NV_DATA_S(yS[j]);
                  output[j+1][i]=sdata[0];
                }

                /* Save the solution */
                for(j=0;j<data->no_equations;j++){
                  outputsol[j+1][i] = Ith(y,j);
                }

                /* Update time for next iteration */
                tout = data->hbtime[i];
        }
        

        //AfterSolve(output, data);
        free(sens_tol);
        return 0;
}

void PrintOutputS(N_Vector *uS, UserData data){
  realtype *sdata;
  int i;
  for(i=0;i<data->no_parameters;i++){
    sdata = NV_DATA_S(uS[i]);
    printf("%12.4le   ",sdata[0]);
      }
  printf("\n ");
}

realtype **matrix(long nrl, long nrh, long ncl, long nch){
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	realtype **m;

	/* allocate pointers to rows */
	m=(realtype **) malloc((size_t)((nrow+NR_END)*sizeof(realtype*)));
	if (!m) printerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(realtype *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(realtype)));
	if (!m[nrl]) printerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void printerror(char error_text[]){
/* Error handler */
	fprintf(stderr,"Run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void free_matrix(realtype **m, long nrl, long nrh, long ncl, long nch){
/* free a realtype matrix allocated by matrix() */
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

int CalcRelSens(realtype **relsens, realtype **solution, realtype **sens, UserData data){
//int AddRelSensToSum(realtype **senssum, realtype **solution, realtype **sens, UserData data){
/* Will calculate the relative sensitivity at the current parameter values and add to the senssum matrix */
int j,k;
    for(j=1;j<data->no_parameters+1;j++)
    {
      relsens[j][0] = sens[j][0]/solution[1][0];
      for(k=1;k<data->hbdatalength;k++)
      {
        ///* This is the intended equation. It is the relative sensitivity */
        relsens[j][k] = (sens[j][k] - sens[j][k-1])/(solution[1][k]-solution[1][k-1])*data->p[j-1];

        /*Trying to figure out what is wrong. This line will just add solutions.*/
        //senssum[j][k] = senssum[j][k] + solution[1][k];
        
        /*Trying to figure out what is wrong. This line will just add sensitivity.*/
        //senssum[j][k] = senssum[j][k] + sens[j][k];
      }
    }

    //for(j=1;j<data->no_parameters+1;j++)
    //{
    //  for(k=0;k<data->hbdatalength-1;k++)
    //  {
    //    ///* This is the intended equation. It is the relative sensitivity */
    //    senssum[j][k] = senssum[j][k] + (sens[j][k+1] - sens[j][k])/(solution[1][k+1]-solution[1][k])*data->p[j-1];

    //    /*Trying to figure out what is wrong. This line will just add solutions.*/
    //    //senssum[j][k] = senssum[j][k] + solution[1][k];
    //    
    //    /*Trying to figure out what is wrong. This line will just add sensitivity.*/
    //    //senssum[j][k] = senssum[j][k] + sens[j][k];
    //  }
    //}
      //printf("\n");
 return 0;
}

int AddSens(realtype **senssum, realtype **sens, UserData data){
/* Will calculate the relative sensitivity at the current parameter values and add to the senssum matrix */
int j,k;
    for(j=1;j<data->no_parameters+1;j++)
    {
      for(k=0;k<data->hbdatalength;k++)
      {
        /* For all sensitivities - add */
        senssum[j][k] = senssum[j][k] + sens[j][k];
      }
    }
      //printf("\n");
 return 0;
}

int PotToHr(realtype **pot, realtype **solution, realtype **HR, UserData data){
int j,k;
    for(j=1;j<data->no_parameters+1;j++)
    {
      /* Calculate initial sensitivity. Divid by time passed since 0 */
      HR[j][0] = pot[j][0]/solution[1][0];

      for(k=1;k<data->hbdatalength;k++)
      {
        /* Calculate remaining */
        HR[j][k] = (pot[j][k]-pot[j][k-1])/(solution[1][k]-solution[1][k-1]);
      }
    }
      //printf("\n");
 return 0;
}

void ProcessArgsML(int argc, char *argv[], FILE **input, FILE **input2, FILE **output, FILE **parfile, FILE **settings, FILE **initialvaluesfile){
    /* If there is not enough arguments, make an error */
    if(argc != 13){ 
      WrongArgs(argv[0]);
    };
  
    int i;

  for(i=1; i <= argc-2; i=i+2)
  {
      /* If argument i = "-inputdata: */
      if(strcmp(argv[i],"-inputdata") == 0)
      {
        *input = fopen(argv[i+1],"r"); 
        if(*input==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input data file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as inputdata.\n", argv[i+1]);
      }
      
      /* If argument i = "-inputhbtimes: */
      else if(strcmp(argv[i],"-inputhbtimes") == 0)
      {
        *input2 = fopen(argv[i+1],"r"); 
        if(*input2==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input heartbeat times data file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
        else{
//        printf("Opened file %s as hearbeattimes data file.\n", argv[i+1]);
        }
      }

      /* If argument i = "-outputfile: */
      else if(strcmp(argv[i],"-outputfile") == 0)
      {
        *output = fopen(argv[i+1],"w");
//        printf("Opened file %s as outputfile.\n", argv[i+1]);
      }
      
      /* If argument i = "-inputdata: */
      else if(strcmp(argv[i],"-pars") == 0)
      {
        *parfile = fopen(argv[i+1],"r");
        if(*parfile==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input parameter file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as parameter file.\n", argv[i+1]);
      }

      /* If argument i = "-inputdata: */
      else if(strcmp(argv[i],"-settings") ==0)
      {
        *settings = fopen(argv[i+1],"r");
        if(*settings==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Input settings file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as settings file.\n", argv[i+1]);
      }

      /* If argument i = "-initialvalues: */
      else if(strcmp(argv[i],"-initialvalues") ==0)
      {
        *initialvaluesfile = fopen(argv[i+1],"r");
        if(*initialvaluesfile==NULL)
        {
          /* If the file doesn't exist, produce an error */
          printf("ERROR: Initialvalues file \"%s\" doesn't exist.\n", argv[i+1]);
          exit(0);
        }
//        printf("Opened file %s as initial values file.\n", argv[i+1]);
      }

      /* If the argument is something else, end the program */
      else
      {
        WrongArgs(argv[0]);
      }
  }
  //printf("Finished loading input files\n");
  }; 

int Scale(realtype **sens, realtype scale, UserData data){
/* Will calculate the relative sensitivity at the current parameter values and add to the senssum matrix */
int j,k;
//realtype scale2;
//scale2 = (realtype) scale;
    for(j=1;j<data->no_parameters+1;j++)
    {
      for(k=0;k<data->hbdatalength;k++)
      {
        /* For all sensitivities - add */
        sens[j][k] = sens[j][k]/ scale;
      }
    }
      //printf("\n");
 return 0;
}
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
