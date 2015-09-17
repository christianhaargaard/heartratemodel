from __future__ import division
import scipy.integrate
import numpy
#from sobol_lib import *

#from c_functions import quasi_random
#from joblib import Parallel, delayed
#from multiprocessing import Pool
from multiprocessing import Pool, cpu_count
from tabulate import tabulate



def sobol_saltelli(M,generators,response):
    """
    [Si, STi] = sobol_saltelli(M,generators,response)
    
    Calculates first order and total Sobol Indices of a scalar function, using the method produced by Weirs et al., 2010.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the original Saltelli algorithm.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"
    p = numpy.size(generators)
    A = numpy.zeros([M,p])
    B = numpy.zeros([M,p])

    for i in xrange(p):
        A[:,i] = generators[i](M)
        B[:,i] = generators[i](M)

    #     print 'Here comes A,B\n'
    #     print A,B
    #     return A,B

    yA = numpy.zeros(M)
    yB = numpy.zeros(M)
    yC = numpy.zeros(M)

    for j in xrange(M):
        yA[j] = response(numpy.transpose(numpy.matrix(A[j,:])))
        yB[j] = response(numpy.transpose(numpy.matrix(B[j,:])))

    #     print 'Here comes yA,yB\n'
    #     print yA, yB

    f02 = (numpy.sum(yA)*numpy.sum(yB))/(M**2)
    #     print f02

    Si = numpy.zeros(p)
    STi = numpy.zeros(p)

    C = numpy.zeros([M,p])
    for i in xrange(p):
        Ci = numpy.copy(B)
        Ci[:,i] = numpy.copy(A[:,i])
        #         print 'Difference A-C:' + str(A-Ci) + 'B-Ci:' + str(B-Ci)
        for j in xrange(M):
            yC[j] = response(numpy.transpose(numpy.matrix(Ci[j,:])))

        #print 'sum(yA*yC)/M =' + str(numpy.sum(yA*yC)/M) + ', sum(yA*yA)/M =' + str(numpy.sum(yA*yA)/M) + '.'
        Si[i] = (numpy.sum(yA*yC)/M-f02)/(numpy.sum(yA*yA)/M-f02)
        STi[i] = 1 - (numpy.sum(yB*yC)/M-f02)/(numpy.sum(yA*yA)/M-f02)

        print 'Done for parameter:' + str(i+1) + '/' + str(p)

    return [Si, STi]

def sobol_weirs(M,no_p,response):
    """ 
    [Si, STi] = sobol_weirs(M,generators,response) 
    
    Calculates first order and total Sobol Indices of a scalar function, using the method produced by Weirs et al., 2010.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the method presented in the paper "Sensitivity analysis techniques applied to a system of hyperbolic conservation laws" by V. Gregory Weirs, James R. Kamm, Laura P. Swiler, Stefano Tarantola, Marco Ratto, Brian M. Adams, William J. Rider, Michael S. Eldred in Reliability Engineering and System Safety, vol. 107, 2012, pages 157-170.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"
    #A = numpy.zeros([M,no_p])
    #B = numpy.zeros([M,no_p])


    # Draw random parameter values ~ U(0,1)
    A = numpy.random.uniform(size=[M,no_p])
    B = numpy.random.uniform(size=[M,no_p])
    #print 'Here comes A,B\n'
    #print A,B
    #return A,B

    yA = numpy.zeros(M)
    yB = numpy.zeros(M)
    yAB = numpy.zeros(M)
    yBA = numpy.zeros(M)

    for j in xrange(M):
        yA[j] = response(A[j])
        yB[j] = response(B[j])

        # If the system cannot be solved, try a new parameter set
        while(yA[j] == 0):
            temp_p = numpy.zeros(no_p)
            for i in xrange(no_p):
                #temp_p[i] = generators[i](generator_input[i],1)
                temp_p[i] = numpy.random.uniform(3)
            yA[j] = response(temp_p)
            del(temp_p)

        while(yB[j] == 0):
            temp_p = numpy.zeros(no_p)
            for i in xrange(no_p):
                #temp_p[i] = generators[i](generator_input[i],1)
                temp_p[i] = numpy.random.uniform(3)
            yB[j] = response(temp_p)
            del(temp_p)
        #yA[j] = response(numpy.transpose(numpy.matrix(A[j,:])))
        #yB[j] = response(numpy.transpose(numpy.matrix(B[j,:])))
    
    yC = numpy.append(yA,yB)
    yC2 = numpy.core.mean(yC)**2
    denominator = numpy.sum(yC*yC- yC2)

    Si = numpy.zeros(no_p)
    STi = numpy.zeros(no_p)

    for i in xrange(no_p):
        Bi = numpy.copy(B)
        Bi[:,i] = numpy.copy(A[:,i])
        Ai = numpy.copy(A)
        Ai[:,i] = numpy.copy(B[:,i])

        for j in xrange(M):
            yAB[j] = response(Ai[j,:])
            yBA[j] = response(Bi[j,:])

            #yC[j] = response(numpy.transpose(numpy.matrix(Ci[j,:])))

        #Si[i] = 2*numpy.sum(yA*(yBA -yB))/denominator
        Si[i] = 1.0/M*sum(yA*(yBA-yB))/(1.0/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))
        #STi[i] = numpy.sum((yA-yAB)*(yA-yAB))/denominator
        STi[i] = 1.0/(2*M)*sum((yA-yAB)**2)/(1.0/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))

        print 'Done for parameter:' + str(i+1) + '/' + str(no_p)

    return [Si, STi]

def sobol_weirs_series(M,generators,generator_input,response):
    """ 
    [Si, STi] = sobol_weirs(M,generators,response) 
    
    Calculates first order and total Sobol Indices of a function returning an array, using the method produced by Weirs et al., 2010.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the method presented in the paper "Sensitivity analysis techniques applied to a system of hyperbolic conservation laws" by V. Gregory Weirs, James R. Kamm, Laura P. Swiler, Stefano Tarantola, Marco Ratto, Brian M. Adams, William J. Rider, Michael S. Eldred in Reliability Engineering and System Safety, vol. 107, 2012, pages 157-170.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"
    no_p = numpy.size(generators)


    A = numpy.zeros([M,no_p])
    B = numpy.zeros([M,no_p])

        
    #for i in xrange(M):
        #A[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*i4_sobol(no_p,i+2*M)[0]
        #B[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*i4_sobol(no_p,i+4*M)[0]

    for i in xrange(no_p):
        A[:,i] = generators[i](generator_input[i],M)
        B[:,i] = generators[i](generator_input[i],M)

    #print 'Here comes A,B\n'
    #print A,B
    #return A,B

    # Asses the size of the output.
    resp = response(A[0])
    n_t = numpy.size(resp)

    # Container for model output
    yA = numpy.zeros([M,n_t])
    yB = numpy.zeros([M,n_t])
    yAB = numpy.zeros([M,n_t])
    yBA = numpy.zeros([M,n_t])


    for j in xrange(M):
        yA[j,:] = response(A[j])
        yB[j,:] = response(B[j])

    yC = numpy.vstack([yA,yB])
    yC2 = numpy.core.mean(yC,0)**2
    denominator = numpy.sum(yC*yC- yC2,0)

    Si = numpy.zeros([no_p,n_t])
    STi = numpy.zeros([no_p,n_t])
    

    # For each parameter
    for i in xrange(no_p):
        Bi = numpy.copy(B)
        Bi[:,i] = numpy.copy(A[:,i])
        Ai = numpy.copy(A)
        Ai[:,i] = numpy.copy(B[:,i])

        # For parameter set
        for j in xrange(M):
            yAB[j,:] = response(Ai[j,:])
            yBA[j,:] = response(Bi[j,:])

            #yC[j] = response(numpy.transpose(numpy.matrix(Ci[j,:])))

        Si[i] = 2*numpy.sum(yA*(yBA -yB),0)/denominator
        #Si[i] = 1/M*sum(yA*(yBA-yB))/(1/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))
        STi[i] = numpy.sum((yA-yAB)*(yA-yAB),0)/denominator
        #STi[i] = 1/(2*M)*sum((yA-yAB)**2)/(1/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))

        print 'Done for parameter:' + str(i+1) + '/' + str(no_p)

    return [Si, STi]

def sobol_weirs_series_sobol_seq(M,generators,generator_input,response):
    """ 
    [Si, STi] = sobol_weirs(M,generators,response) 
    
    Calculates first order and total Sobol Indices of a function returning an array, using the method produced by Weirs et al., 2010.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the method presented in the paper "Sensitivity analysis techniques applied to a system of hyperbolic conservation laws" by V. Gregory Weirs, James R. Kamm, Laura P. Swiler, Stefano Tarantola, Marco Ratto, Brian M. Adams, William J. Rider, Michael S. Eldred in Reliability Engineering and System Safety, vol. 107, 2012, pages 157-170.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"
    no_p = numpy.size(generators)


    A = numpy.zeros([M,no_p])
    B = numpy.zeros([M,no_p])

    sobol_seq(A,0,0)
    sobol_seq(B,no_p,0)

        
    for i in xrange(M):
        A[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*A[i]
        B[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*B[i]

    # Asses the size of the output.
    A0 = numpy.copy(A[0])
    resp = response(A0)
    n_t = numpy.size(resp)

    # Container for model output
    yA = numpy.zeros([M,n_t])
    yB = numpy.zeros([M,n_t])
    yAB = numpy.zeros([M,n_t])
    yBA = numpy.zeros([M,n_t])


    for j in xrange(M):
        yA[j,:] = response(A[j])
        yB[j,:] = response(B[j])

    yC = numpy.vstack([yA,yB])
    yC2 = numpy.core.mean(yC,0)**2
    denominator = numpy.sum(yC*yC- yC2,0)

    Si = numpy.zeros([no_p,n_t])
    STi = numpy.zeros([no_p,n_t])
    

    # For each parameter
    for i in xrange(no_p):
        Bi = numpy.copy(B)
        Bi[:,i] = numpy.copy(A[:,i])
        Ai = numpy.copy(A)
        Ai[:,i] = numpy.copy(B[:,i])

        # For parameter set
        for j in xrange(M):
            yAB[j,:] = response(Ai[j,:])
            yBA[j,:] = response(Bi[j,:])

            #yC[j] = response(numpy.transpose(numpy.matrix(Ci[j,:])))

        Si[i] = 2*numpy.sum(yA*(yBA -yB),0)/denominator
        #Si[i] = 1/M*sum(yA*(yBA-yB))/(1/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))
        STi[i] = numpy.sum((yA-yAB)*(yA-yAB),0)/denominator
        #STi[i] = 1/(2*M)*sum((yA-yAB)**2)/(1/(2*M)*sum(yC*yC-numpy.core.mean(yC)**2))

        print 'Done for parameter:' + str(i+1) + '/' + str(no_p)
    return [Si, STi]

def sobol_weirs_series_sobol_seq_parallel(M,generators,generator_input,response):
    """ 
    DO NOT USE IN CURRENT VERSION
    [Si, STi] = sobol_weirs_series_sobol_seq_parallel(M,generators,generator_input,response)
    
    Calculates first order and total Sobol Indices of a function returning an array, using the method produced by Weirs et al., 2010.
    The code will spawn workers corresponding to the number of kernels available to python, and run model evaluations in parallel.
    Parameter sampling is done using a scrambled Halton sequence by XXXXXX et. al, using the GNU Scientific Library implementation.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the method presented in the paper "Sensitivity analysis techniques applied to a system of hyperbolic conservation laws" by V. Gregory Weirs, James R. Kamm, Laura P. Swiler, Stefano Tarantola, Marco Ratto, Brian M. Adams, William J. Rider, Michael S. Eldred in Reliability Engineering and System Safety, vol. 107, 2012, pages 157-170.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"
    no_p = numpy.size(generators)


    A = numpy.zeros([M,no_p])
    B = numpy.zeros([M,no_p])

    #halton_seq(A,0,0)
    #halton_seq(B,no_p,0)
    rng = quasi_random()

    #rng.sobol_seq(A,0,0)
    #rng.sobol_seq(B,no_p,0)
    skip = int(-numpy.floor(-numpy.log2(M)))
    rng.reversehalton_seq(A,0,skip)
    rng.reversehalton_seq(B,no_p,skip)
        
    for i in xrange(M):
        A[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*A[i]
        B[i] = generator_input[:,0] + (generator_input[:,1] - generator_input[:,0])*B[i]

    # Asses the size of the output.
    A0 = numpy.copy(A)
    resp = response(A0[0])
    n_t = numpy.size(resp)

    # Container for model output
    yA = numpy.zeros([M,n_t])
    yB = numpy.zeros([M,n_t])
    yAB = numpy.zeros([M,n_t])
    yBA = numpy.zeros([M,n_t])

    # Some parallel magic
    # from http://stackoverflow.com/questions/9742739/how-do-i-make-processes-able-to-write-in-an-array-of-the-main-program/9849971
    pool = Pool(processes=cpu_count())

    yA = parallel(A,response,pool)
    yB = parallel(B,response,pool)

    yC = numpy.vstack([yA,yB])
    yC2 = numpy.core.mean(yC,0)**2
    denominator = numpy.sum(yC*yC- yC2,0)

    Si = numpy.zeros([no_p,n_t])
    STi = numpy.zeros([no_p,n_t])
    

    # For each parameter
    for i in xrange(no_p):
        Bi = numpy.copy(B)
        Bi[:,i] = numpy.copy(A[:,i])
        Ai = numpy.copy(A)
        Ai[:,i] = numpy.copy(B[:,i])

        # For parameter set
        yAB = parallel(Ai,response,pool)
        yBA = parallel(Bi,response,pool)

        Si[i] = 2*numpy.sum(yA*(yBA -yB),0)/denominator
        STi[i] = numpy.sum((yA-yAB)*(yA-yAB),0)/denominator

        print 'Done for parameter:' + str(i+1) + '/' + str(no_p)

    pool.close()
    pool.join()
    return [Si, STi]

def sobol_weirs_series_uniform_parallel(M,no_p,response):
    """ 
    [Si, STi] = sobol_weirs(M,generators,response) 
    
    Calculates first order and total Sobol Indices of a function returning an array, using the method produced by Weirs et al., 2010.
    
    Parameters
    ----------
    M : The number of points in parameter space used.
    generators : An array containing generators for the parameters of the array required as input to the scalar function.
    response : A function accepting as input a numpy.array of parameter values.

    Returns
    -------
    Si : First order indices
    STi : Total indices
    
    References
    ----------
    The code is based on the method presented in the paper "Sensitivity analysis techniques applied to a system of hyperbolic conservation laws" by V. Gregory Weirs, James R. Kamm, Laura P. Swiler, Stefano Tarantola, Marco Ratto, Brian M. Adams, William J. Rider, Michael S. Eldred in Reliability Engineering and System Safety, vol. 107, 2012, pages 157-170.
    """

    print "Calculating first order Sobol Indices Si and STi using M=" + str(M) + " points.\n"

    # Draw random parameter values ~ U(0,1)
    A = numpy.random.uniform(size=[M,no_p])
    B = numpy.random.uniform(size=[M,no_p])

    # Asses the size of the output.
    A0 = numpy.copy(A)
    resp = response(A0[0])
    n_t = numpy.size(resp)

    # Container for model output
    yA = numpy.zeros([M,n_t])
    yB = numpy.zeros([M,n_t])
    yAB = numpy.zeros([M,n_t])
    yBA = numpy.zeros([M,n_t])

    # Some parallel magic
    # from http://stackoverflow.com/questions/9742739/how-do-i-make-processes-able-to-write-in-an-array-of-the-main-program/9849971
    #pool = Pool(processes=cpu_count())
    pool = Pool(processes=2)

    yA = parallel(A,response,pool)
    yB = parallel(B,response,pool)

    yC = numpy.hstack([yA,yB])
    yC2 = numpy.core.mean(yC,0)**2
    denominator = numpy.sum(yC*yC- yC2,0)

    Si = numpy.zeros([no_p,n_t])
    STi = numpy.zeros([no_p,n_t])
    

    # For each parameter
    for i in xrange(no_p):
        Bi = numpy.copy(B)
        Bi[:,i] = numpy.copy(A[:,i])
        Ai = numpy.copy(A)
        Ai[:,i] = numpy.copy(B[:,i])

        # For parameter set
        yAB = parallel(Ai,response,pool)
        yBA = parallel(Bi,response,pool)

        Si[i] = 2*numpy.sum(yA*(yBA -yB),0)/denominator
        STi[i] = numpy.sum((yA-yAB)*(yA-yAB),0)/denominator

        print 'Done for parameter:' + str(i+1) + '/' + str(no_p)

        # Clean up
        del Bi, Ai, yAB, yBA

    pool.close()
    pool.join()

    # Clean up
    del A, B, A0, resp, n_t, yA, yB, yC, denominator

    return [Si, STi]

def morris2(l,r,p,distribution_rescalers,rhs,yinit,t,response):
    print "Calculating Morris measures mean and variance for l=" + str(l) + ", r=" + str(r) + ".\n"

    # For l even, this delta prevent Type2 errors
    delta = l/(2*(l-1))

    #Sample base parameter values 
    #Q: Should these all be uniform, or should I be using another distribution for k?
    q = qvalue(l,[r,p]) # Need to be on the interval [0,1-delta]


    ###############################################################################
    #                       Do steps and calculate indices                        #
    ###############################################################################

    #Containers for responses and d values
    y = numpy.zeros([r,p])
    d = numpy.zeros([r,p])

    # Run through each parameter
    for i in xrange(p):

        # Container for the base solution
        y0 = numpy.zeros(r)

        # Sample parameter values from grid
        q = numpy.zeros([r,p])
        q = qvalue(l,[r,p]) 

        # Calculate base solution
        for j in xrange(r): 
            tq = numpy.copy(q[j])

            #rescale parameter to follow correct distribution - given by rescalers 
            for k in xrange(p):
                tq[k] = distribution_rescalers[k](tq[k])

            y0[j] = response(rhs,yinit,t,tq)

        # Reset step vector

        # Make sure the value+step is within boundaries
        for j in xrange(r):
            vec = numpy.zeros(p)
            if q[j,i]+delta>1:
                vec[i] = -delta
            else:
                vec[i] = delta

        # Run through each of the r realizations
        #for j in xrange(r):

            tq = numpy.copy(q[j] + vec)
            if numpy.abs(tq[1]) > 1:
                print 'error'

            #rescale parameter no follow correct distribution - given by rescalers 
            for k in xrange(p):
                tq[k] = distribution_rescalers[k](tq[k])

            y[j,i] = response(rhs,yinit,t,tq)
            d[j,i] = (y[j,i] - y0[j])/vec[i]

    #Containers for mean and variance of d values 
    uis = numpy.zeros(p)
    ui  = numpy.zeros(p)
    s2i = numpy.zeros(p)
    
    #Calculate mean indices and variances
    for i in xrange(p):
        uis[i] = 1/r*numpy.sum(numpy.abs(d[:,i]))
        ui[i] = 1/r*numpy.sum(d[:,i])
        s2i[i] = 1/(r-1)*numpy.sum((d[:,i]-ui[i])**2)

    return  [uis, numpy.sqrt(s2i)]

def morris_series_parallel_old(l,r,p,response):
    print "Calculating Morris measures mean and variance for l=" + str(l) + ", r=" + str(r) + ".\n"

    # For l even, this delta prevent Type2 errors
    delta = l/(2.0*(l-1))

    #skip = 0

    #Sample base parameter values 
    #q = qvalue_halton(l,[r,p],0) # Need to be on the interval [0,1-delta]
    q = qvalue(l,[r,p])
    #skip = r


    ###############################################################################
    #                       Do steps and calculate indices                        #
    ###############################################################################

    q0 = numpy.copy(q)
    resp = response(q0[0])
    n_t = numpy.size(resp)

    #Containers for responses and d values
    y = numpy.zeros([r,p,n_t])
    d = numpy.zeros([r,p,n_t])

    pool = Pool(processes=cpu_count())

    # Run through each parameter
    for i in xrange(p):

        # Container for the base solution
        y0 = numpy.zeros([r,n_t])

        # Sample parameter values from grid
        #q = qvalue_halton(l,[r,p],(i+1)*skip) 
        q = qvalue(l,[r,p])

        #for j in xrange(p):
            #numpy.random.shuffle(q[:,j])

        # Calculate base solution
        #for j in xrange(r): 
            #tq = numpy.copy(q[j])
            #y0[j] = response(tq)

        result = []
        for j in xrange(r):
            result.append(pool.apply_async(response,([q[j]])))

        y0 = []
        for j in xrange(r):
            y0.append(result[j].get())
            #y0[j] 

        # Reset step vector

        # Make sure the value+step is within boundaries
        result = []
        for j in xrange(r):
            vec = numpy.zeros(p)
            if q[j,i]+delta>1:
                vec[i] = -delta
            else:
                vec[i] = delta

            # Run through each of the r realizations

            tq = numpy.copy(q[j] + vec)
            if numpy.abs(tq[1]) > 1:
                print 'error'

            result.append(pool.apply_async(response,([tq])))

        for j in xrange(r):
            #y[j,i] = response(tq)
            y[j,i] = result[j].get()
            d[j,i] = (y[j,i] - y0[j])/vec[i]

    #Containers for mean and variance of d values 
    uis = numpy.zeros([p,n_t])
    ui  = numpy.zeros([p,n_t])
    s2i = numpy.zeros([p,n_t])

    pool.close()
    pool.join()

    #Calculate mean indices and variances
    for i in xrange(p):
        for j in xrange(n_t):
            uis[i,j] = numpy.sum(numpy.abs(d[:,i,j]))/r
            ui[i,j] = numpy.sum(d[:,i,j])/r
            s2i[i,j] = numpy.sum((d[:,i,j]-ui[i,j])**2)/(r-1)

    return  [uis, numpy.sqrt(s2i)]

def morris_series_parallel(l,r,p,response):
    """ 
    [Mus, S2] = morris_series_parallel(l,r,p,response) 
    
    Calculates elementary effects using Morris method.
    
    Parameters
    ----------
    l : The number of points per dimension on the parameter grid.
    r : The number of evaluations of elementary effects for each parameter.
    p : The number of parameters.
    response : a function that takes a numpy.array of length p with parameter-values uniformly distributed between 0 and 1.

    Returns
    -------
    Mus : Sampling mean (the mean of the ABSOLUTE value of the elementary effects).
    S2 : Sampling variance.
    
    References
    ----------
    The code is presented on the algorithm presented by Smith in the book "Uncertainty Quantification: Theory, Implementation, and Applications" page 331-334.

    Comments
    --------
    The code could probably be improved by implementing methods for selecting good initial points on the parameter grid.

    Author
    ------
    Christian Haargaard Olsen <chaarga@ncsu.edu>

    """
    print "Calculating Morris measures mean and variance for l=" + str(l) + ", r=" + str(r) + ".\n"

    # For l even, this delta prevent Type2 errors
    delta = l/(2.0*(l-1))

    #Sample base parameter values 
    q = qvalue(l,[r,p])

    ###############################################################################
    #                       Do steps and calculate indices                        #
    ###############################################################################

    #Determine size of output
    resp = response(q[0])
    n_t = numpy.size(resp)

    #Clean up
    del q, resp

    #Containers for responses and d values
    d = numpy.zeros([r,p,n_t])
    uis = numpy.zeros([p,n_t])
    ui = numpy.zeros([p,n_t])
    s2i = numpy.zeros([p,n_t])

    pool = Pool(processes=cpu_count())

    for i in xrange(r):

        # generate initial parameter configuration - the point to start in the grid
        q = qvalue(l,p)

        results = numpy.zeros([p+1,n_t])

        # Generate the permutation matrix - See Smith. 2013
        A = numpy.zeros([p+1,p])
        J = numpy.ones([p+1,p])
        for j in xrange(numpy.size(A,0)):
            for k in xrange(numpy.size(A,1)):
                if j>k:
                    A[j,k] = 1

        D = numpy.diag(numpy.floor(2*numpy.random.uniform(size=p)-1)*2+1)
        P = numpy.identity(p)
        numpy.random.shuffle(P)

        B = numpy.dot(J*q + delta/2*(numpy.dot((2*A-J),D) + J),P)

        # In case the permutation takes any values to be larger than 1, find the indexes that overshoots, and subtract 2*delta instead
        over = -numpy.floor(-numpy.mod(B,1)*(B//1))
        B = B-2*delta*over

        # Find the order of steps, so the elementary effects can be assigned to the correct parameter
        shift = numpy.transpose(numpy.diff(B.T,1))
        shiftindex = numpy.argmax(numpy.abs(shift),1)
        deltas = numpy.sum(shift,1)
        #print 'Here are the deltas:' + str(deltas)
        #print 'And the shifts:' + str(shiftindex)

        #parallel code
        result = []
        for j in xrange(p+1):
            #print B[j]
            result.append(pool.apply_async(response,([B[j]])))

        y0 = []
        for j in xrange(p+1):
            y0.append(result[j].get())
        

        #Serial code
        #y0 = []
        #for j in xrange(p+1):
            #y0.append(response(B[j]))
        
        # Calculate the elementary effects, and assign them to the correct parameter
        for j in xrange(p):
            d[i,shiftindex[j]] = (y0[j+1]-y0[j])/deltas[j]
        #print 'Done for r=' + str(i)

        # Clean up
        del q, results, A, J, D, P, B, over, shift, shiftindex, deltas, result, y0
        #del q, A, J, D, P, B, over, shift, shiftindex, deltas, result
         
    pool.close()
    pool.join()

    for i in xrange(p):
        for j in xrange(n_t):
            uis[i,j] = numpy.sum(numpy.abs(d[:,i,j]))/r
            ui[i,j] = numpy.sum(d[:,i,j])/r
            s2i[i,j] = numpy.sum((d[:,i,j]-ui[i,j])**2)/(r-1)

    del ui
    return  [uis, numpy.sqrt(s2i)]
    #return  [y0, y0]

def osm(basis = numpy.matrix([]), tolerance = 1e-10):
    p = numpy.size(basis,1)
    STS =  numpy.dot(numpy.transpose(basis),basis)
    [l, Q] = numpy.linalg.eig(STS);

    # ii = argsort(l)
    # ii = ii[::-1]
    # l = l[ii]
    # Q = Q[:,ii]
    E = numpy.array(numpy.dot(numpy.abs(Q),numpy.abs(l)))/numpy.sum(numpy.abs(l))
    E = E[0]

    # Select first item
    element = numpy.argmax(E)
    selected_index = [element]
    not_selected = range(p)
    not_selected.remove(element)

    # Now for angle
    angle = [0]     # First element has no angle
    maxI = [0]      # First element has no angle and thus no I-value

    if(1): # We do want to find angles, etc
        for i in range(1,p):        # We have already selected one parameter
            d = numpy.zeros(p)      # Container for angles
            for j in range(p):
                if j in selected_index:
    #                 print 'already selected:' + str(j)
                    d[j] = 0        # If parameter has already been added angle is zer
                else:
                    # find a first
                    a = numpy.linalg.lstsq(basis[:,selected_index],basis[:,j])
                    s = basis[:,selected_index]*a[0]
                    denom = numpy.linalg.norm(s)*numpy.linalg.norm(basis[:,j])

                    # If the column has near 0 2-norm, set d[j] = 0
                    if numpy.linalg.norm(basis[:,j]) < tolerance:
                        d[j] = 0.0

                    else:
                        if numpy.linalg.norm(a[0]) < tolerance:
                            d[j] = 1.0
                        else:

                            # Secure against numerical errors causing a projection longer than what is projected onto
                            proj_length = numpy.dot(numpy.transpose(s),basis[:,j])[0,0]

                            if numpy.abs(proj_length/denom) > 1:
                                d[j] = 0
                            else:
                                d[j] = numpy.sin(numpy.arccos(proj_length/denom))

            I = E*d;
            
            # If all values of I is 0, select the one with the largest E-value
            if numpy.linalg.norm(I)==0:
                print "All I-values are 0\n"

                # If all E-values are 0, pick first occuring number
                if numpy.linalg.norm(E[not_selected]) == 0:
                    element = numpy.min(not_selected)
                else:
                    print "Im here"
                    element = not_selected[numpy.argmax(E[not_selected])]
                    #no = numpy.argmax(E[not_selected])
                    #element = not_selected[no]
            else: 
                element = numpy.argmax(I)

            #print "element = " + str(element)
            #print "At round %(round)i of %(totalround)i." % {"round": i, "totalround": range(p)-1}

            maxI.append(numpy.max(I));
            angle.append(d[element]);
            selected_index.append(element)  # Should select the largest that has not been selected yet
            not_selected.remove(element)


    names = ['Rank', 'Par' , 'd', 'E', 'I']
    data = numpy.zeros([p,5])
    data[:,0] = range(1,p+1)
    # data[:,1] = numpy.mat(selected_index)+1
    data[:,1] = selected_index
    data[:,2] = angle
    data[:,3] = E[selected_index]
    data[:,4] = maxI
    return data, names

def scm(basis = numpy.matrix([]), condition_limit= 1e20):
    # Should return ordering of the columns that is in 'basis'
    p = numpy.size(basis,1)
    STS =  numpy.dot(numpy.transpose(basis),basis)
    
    K = numpy.linalg.cond(STS)
    if (K < condition_limit):
       covariance_matrix = numpy.linalg.pinv(STS)
       correlation_matrix = numpy.zeros(numpy.shape(covariance_matrix))
       for i in xrange(p):
            for j in xrange(p):
                correlation_matrix[i,j] = covariance_matrix[i,j]/numpy.sqrt(covariance_matrix[i,i]*covariance_matrix[j,j])
       return correlation_matrix, K
    else:
        print 'condition number exceeds limit' + '{0:2.2e}'.format(K) + '>' + '{0:2.2e}'.format(condition_limit) 
        return [], K


def qvalue(l,n=1):
    # Draw parameter values from the proper grid points
    #return numpy.floor(numpy.random.uniform(size=n)*(l/2))/(l-1)
    return numpy.floor(numpy.random.uniform(size=n)*(l))/(l-1)

def qvalue_halton(l,dim=numpy.array([1,1]),skip = 0):
    # Draw parameter values from the proper grid points, using reverse halton
    r = dim[0]
    p = dim[1]

    rng = quasi_random()
    A = numpy.zeros(dim)

    skip = skip + int(-numpy.floor(-numpy.log2(r)))

    rng.reversehalton_seq(A,0,skip)
    return numpy.floor(A*(l))/(l-1)

def parallel(A,response,pool):
    ## Some parallel magic
    ## http://stackoverflow.com/questions/9742739/how-do-i-make-processes-able-to-write-in-an-array-of-the-main-program/9849971#9849971
    #pool = Pool(16)

    result = []
    for i in range(len(A)):
        result.append(pool.apply_async(response,([A[i]])))

    m = []
    for i in range(len(A)):
        m.append(result[i].get())

    return numpy.array(m)
