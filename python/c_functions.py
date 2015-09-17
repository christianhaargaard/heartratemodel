import numpy
import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import pointer 
from ctypes import CDLL
import matplotlib.pyplot 

# must be a double array, with single dimension that is contiguous
array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=numpy.double, ndim=2, flags='CONTIGUOUS')

        
# load the library, using numpy mechanisms
libsol = npct.load_library("../c/libsol", ".")

#######################################################################
#                           Model solution                            #
#######################################################################
# setup the return typs and argument types
libsol.model.restype = c_int
libsol.model.argtypes = [c_int, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, c_int, array_1d_double, array_1d_double, c_int, array_1d_double]

#@profile
def solve(input_time, input_bp, input_bp_spline_coef, input_angle, input_angle_spline_coef, input_respiration, input_respiration_spline_coef, output_time, output_hr, parameters):
    return libsol.model(len(input_time), input_time, input_bp, input_bp_spline_coef, input_angle, input_angle_spline_coef, input_respiration, input_respiration_spline_coef, len(output_time), output_time, output_hr, len(parameters), parameters)


#######################################################################
#                   Derivative based sensitivities                    #
#######################################################################

#######################################################################
#                               Spline                                #
#######################################################################
#libsol.spline.restype = c_int
#libsol.spline.argtypes = [array_1d_double, array_1d_double, c_int, c_double, c_double, array_1d_double]
#def spline(input_time, input_data, coef): # input_angle, input_respiration, output_time, output_hr, parameters):
    #return libsol.spline(input_time, input_data, len(input_time), 1e30, 1e30, coef)

#######################################################################
#                     Sampling of Sobol Sequences                     #
#######################################################################
# This library allows to samle a MxN matrix where each column is a sobol sequence. Only good for N<40.
# One should supply a seed - an integer that determines how many steps of the chain to skip.
#lib_custom_gsl = npct.load_library("/Users/christian/Documents/PhD/research/code/c/lib_custom_gsl", ".")
#lib_custom_gsl.sobol_seq.restype = c_int
#lib_custom_gsl.sobol_seq.argtypes = [c_int, c_int, c_int, c_int, array_2d_double]
#def sobol_seq(output_array, start_dim, seed):



class wrapper(object):
    def __init__(self):
        path = "/Users/christian/Documents/PhD/research/code/articles/new_model/c/"
        self.so = CDLL('%s/libsol.so' %path)

    def model(self,input_time, input_bp, input_angle, input_respiration, output_time, output, parameters, flags=numpy.array([0.0])):
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')

        self.so.model.argtypes = [c_int, array_1d_double, array_1d_double, array_1d_double, array_1d_double, c_int, array_1d_double, array_1d_double, c_int, array_1d_double, array_1d_double]
        self.so.model.restype = c_int

        #return self.so.model(N_in, input_time, input_bp, input_bp_spline_coef, input_angle, input_angle_spline_coef, input_respiration, input_respiration_spline_coef, N_out, output_time, output_hr, N_parameters, parameters)
        return self.so.model(numpy.size(input_time), input_time, input_bp, input_angle, input_respiration, numpy.size(output_time), output_time, output, numpy.size(parameters), parameters, flags)

    def sens(self, input_time, input_bp, input_angle, input_respiration, output_time, output_hr, parameters, output_sensitivity, flags=numpy.array([0.0])):   
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')
        array_3d_double = npct.ndpointer(dtype=numpy.double, ndim=3, flags='CONTIGUOUS')

        self.so.sens.argtypes = [c_int, array_1d_double, array_1d_double, array_1d_double, array_1d_double, c_int, array_1d_double, array_1d_double, c_int, array_1d_double, array_3d_double, array_1d_double]
        self.so.sens.restype = c_int

        return self.so.sens(len(input_time), input_time, input_bp, input_angle, input_respiration, len(output_time), output_time, output_hr, len(parameters), parameters, output_sensitivity, flags)

    def spline(self,input_time, input_data, coef): 
        # input_angle, input_respiration, output_time, output_hr, parameters):
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')

        self.so.spline.restype = c_int
        self.so.spline.argtypes = [array_1d_double, array_1d_double, c_int, c_double, c_double, array_1d_double]
        return self.so.spline(input_time, input_data, len(input_time), 1e30, 1e30, coef)

    def spline_evaluation(self, input_time, input_data, coef, output_time, output_data):
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')
        n = numpy.size(input_time)
        nout = numpy.size(output_time)

        self.so.spline_evaluation.restype = c_int
        self.so.spline_evaluation.argtypes = [array_1d_double, array_1d_double, array_1d_double, c_int, array_1d_double, array_1d_double, c_int]

        return self.so.spline_evaluation(input_time, input_data, coef, n, output_time, output_data, nout)

    def plot_model_output(self, output_time, output_vector, input_time, input_bp, input_respiration, parameters,flags=numpy.array([0.0])):
        [p0,k,AM0,a1,a2,b1,b2,s1,s2,Tpm,TpM,xi,fp,Tsm,TsM,eta,fs,height_threshold,am,aM,ak,bm,bM,bk,td,tA,qp,kiN,tN,qs,tAF,mu,KA,tAS,tNS, KN, h0, hm, hM, hrM, fresp_thresshold, fresp_k ] = parameters
        n = numpy.size(output_time)

        print flags[0]
        if flags[0]> 0:
            print "using interpolated pressure!"
            bp_out = output_vector[:n]

        else:
            bp_coef = 0*numpy.empty_like(input_time)
            self.spline(input_time,input_bp,bp_coef)

            bp_out = 0*numpy.empty_like(output_time)
            self.spline_evaluation(input_time,input_bp,bp_coef,output_time,bp_out)

        
        #######################################################################
        #                                 BP                                  #
        #######################################################################
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time, bp_out,'b')
        matplotlib.pyplot.legend(['Blood pressure'])
        matplotlib.pyplot.show()

        #######################################################################
        #                            Wall Strains                             #
        #######################################################################

        wall_strain = (1.0-numpy.sqrt((p0**k +bp_out**k)/(p0**k + AM0*bp_out**k)))
        e1 = output_vector[n:2*n]              # Strain of both components in visco wall]
        e2 = output_vector[2*n:3*n]           # Strain of both components in visco wall
        ebr = wall_strain-e1
        
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time,wall_strain,'k',output_time,e1,'g',output_time,e2,'c',output_time,ebr,'r')
        matplotlib.pyplot.legend(['$\epsilon_w$','$\epsilon_1$', '$\epsilon_2$', '$\epsilon_{br}$'])
        matplotlib.pyplot.show()

        #######################################################################
        #                            Neuron firing                            #
        #######################################################################

        firing = s1*ebr + s2

        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time,firing,'k')
        matplotlib.pyplot.legend(['$f$'])
        matplotlib.pyplot.show()

        #######################################################################
        #                                 ANS                                 #
        #######################################################################

        parasymp_baro = Tpm + (1 -Tpm)*firing**xi/(firing**xi + fp**xi)
        symp_baro = 1 - (1 -Tsm)*firing**eta/(firing**eta + fs**eta)

        rho = 5
        factorial = 4*3*2
        delay_alpha = rho/td
        scale_constant = delay_alpha**rho/factorial
        delayed_symp = scale_constant*output_vector[7*n:8*n]

        respiration_coef = 0*numpy.empty_like(input_time)
        respiration_out = 0*numpy.empty_like(output_time)
        self.spline(input_time,input_respiration,respiration_coef)
        self.spline_evaluation(input_time,input_respiration,respiration_coef,output_time,respiration_out)
        parasymp_resp = 1-respiration_out**5/(respiration_out**5 + .3**5)

        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time,parasymp_baro,'k', output_time, symp_baro, 'r', output_time, delayed_symp, 'g', output_time, parasymp_resp, 'c')
        matplotlib.pyplot.legend(['$T_{p,br}$','$T_s$', '$T_s^d$', '$T_{r}$'])
        matplotlib.pyplot.show()

        #######################################################################
        #                          Respiratory input                          #
        #######################################################################

        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(input_time,input_respiration,'k')
        matplotlib.pyplot.legend(['Resp'])
        a = numpy.array(matplotlib.pyplot.axis())
        a[2] = 0
        a[3] = 0.5
        matplotlib.pyplot.axis(a)
        matplotlib.pyplot.show()

        #######################################################################
        #                           Concentrations                            #
        #######################################################################

        CA = output_vector[8*n:9*n]
        CN = output_vector[9*n:10*n]
        
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time,CA, 'k', output_time, CN,'r')
        matplotlib.pyplot.legend(['$C_A$','$C_N$'])
        matplotlib.pyplot.show()

        #######################################################################
        #                            Cell Pathways                            #
        #######################################################################
        
        CAS = output_vector[10*n:11*n]
        CNS = output_vector[11*n:12*n]
        CAF = mu*CA**2/(CA**2+KA**2)

        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time,CAS, 'k', output_time, CAF, 'r', output_time, CNS,'g', output_time, CAS+CAF, 'b')
        matplotlib.pyplot.legend(['$C_{AS}$','$C_{AF}$','$C_{NS}$', '$C_{AT}$'])
        matplotlib.pyplot.show()

        #######################################################################
        #                             Heart Rate                              #
        #######################################################################
        
        HR = output_vector[12*n:13*n]

        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(output_time[1:],numpy.diff(HR)/numpy.diff(output_time), 'k')
        matplotlib.pyplot.legend(['Heart rate'])
        matplotlib.pyplot.show()

        return 0

    def HRsens(self,y,S,p,numbers):
        npars = len(p)
        nout = len(y)/13
        SHR = numpy.zeros([npars,nout])
        hm = p[37]
        hM = p[38]
        h0 = p[39]
        mu = p[31]
        KA = p[32]

        CA = y[8*nout:9*nout]
        CAF = mu*CA**2/(CA**2+KA**2)
        CN = y[9*nout:10*nout]
        CAS = y[10*nout:11*nout]
        CNS = y[11*nout:12*nout]

        for i in numbers:
            dCA = S[8,i]
            dCAF = 2*CA*dCA*mu/(CA**2+KA**2)*(1-CA**2/(CA**2+KA**2))
            dCAS = S[10,i]
            dCNS = S[11,i]

            if i==32:
                dCAF = -mu*CA**2/((CA**2+KA**2)**2)*2*KA

            if i==31:
                CAF = CA**2/(CA**2+KA**2)

            SHR[i] = -hm*(dCAS + dCAF) + hM*dCNS - 1/h0*hm*hM*(dCAF*CN + CAF*dCNS)

            if i==36: # h0
                SHR[i] = numpy.ones(numpy.shape(dCA)) + 1/(p[36]**2)*hm*hM*CAF*CNS
                
            if i==37:
                SHR[i] = -(CAS + CAF) - 1/h0*hM*CAF*CN

            if i==38:
                SHR[i] = CNS - 1/h0*hm*CAF*CN
        return SHR

    def all_equations(self,input_time, input_bp, input_angle, input_respiration, output_time, output, parameters, flags=numpy.array([1.0])):
        self.model(input_time, input_bp, input_angle, input_respiration, output_time, output, parameters,flags)
        [p0,k,AM0,a1,a2,b1,b2,s1,s2,Tpm,TpM,xi,fp,Tsm,TsM,eta,fs,height_threshold,am,aM,ak,bm,bM,bk,td,tA,qp,kiN,tN,qs,tAF,mu,KA,tAS,tNS, KN, h0, hm, hM, hrM, fresp_thresshold, fresp_k] = parameters
        n = len(output_time)

        bp_coef = 0.0*numpy.empty_like(input_time)
        self.spline(input_time,input_bp,bp_coef)

        bp_out = 0.0*numpy.empty_like(output_time)
        self.spline_evaluation(input_time,input_bp,bp_coef,output_time,bp_out)

        bp_mean = output[0:n]
    
        #p0 = parameters[0]
        #k =  parameters[1]
        #AM0 = parameters[2]
        wall_strain = (1.0-numpy.sqrt((p0**k +bp_mean**k)/(p0**k + AM0*bp_mean**k)))
        #wall_strain = output[13*n:14*n]
        e1 = output[n:2*n]              # Strain of both components in visco wall]
        e2 = output[2*n:3*n]           # Strain of both components in visco wall
        ebr = wall_strain-e1
        
        #s1 = parameters[7]
        #s2 = parameters[8]
        firing = s1*ebr + s2

        #Tpm = parameters[
        parasymp_baro = Tpm + (1 -Tpm)*firing**xi/(firing**xi + fp**xi)
        symp_baro = 1 - (1 -Tsm)*firing**eta/(firing**eta + fs**eta)

        rho = 5
        factorial = 4*3*2
        if(td>0):
            delay_alpha = rho/td
            scale_constant = (delay_alpha**rho)/factorial
            delayed_symp = scale_constant*output[7*n:8*n]
        else:
            delayed_symp = symp_baro

        respiration_coef = 0*numpy.empty_like(input_time)
        respiration_out = 0*numpy.empty_like(output_time)
        self.spline(input_time,input_respiration,respiration_coef)
        self.spline_evaluation(input_time,input_respiration,respiration_coef,output_time,respiration_out)
        parasymp_resp = 1-respiration_out**fresp_k/(respiration_out**fresp_k + fresp_thresshold**fresp_k)


        angle_coef = 0*numpy.empty_like(input_time)
        self.spline(input_time,input_angle,angle_coef)

        angle_out = 0*numpy.empty_like(output_time)
        self.spline_evaluation(input_time,input_angle,angle_coef,output_time,angle_out)
        
        baro_weight = am + (aM-am)*(numpy.sin(angle_out)**ak)/((numpy.sin(angle_out)**ak) + (numpy.sin(height_threshold)**ak));
        resp_weight = bM - (bM-bm)*(numpy.sin(angle_out)**bk)/((numpy.sin(angle_out)**bk) + (numpy.sin(height_threshold)**bk));

        parasymp_total = baro_weight*parasymp_baro + resp_weight*parasymp_resp

        CA = output[8*n:9*n]
        CN = output[9*n:10*n]
        
        CAS = output[10*n:11*n]
        CNS = output[11*n:12*n]
        CAF = mu*CA**2/(CA**2+KA**2)

        HR = output[12*n:13*n]
        
        result = numpy.array([output_time, bp_out, bp_mean, wall_strain, e1, e2, ebr, firing, parasymp_baro, symp_baro, delayed_symp, parasymp_resp, parasymp_total, CA, CN, CAS, CAF, CNS, HR])
        names = numpy.array(['output_time', 'bp_smoothed', 'bp_mean', 'wall_strain', 'e1', 'e2', 'ebr', 'firing', 'parasymp_baro', 'symp_baro', 'delayed_symp', 'parasymp_resp', 'parasymp_total', 'CA', 'CN', 'CAS', 'CAF', 'CNS', 'HR'])
        return [result, names]


class quasi_random(object):

    def __init__(self):
        path = "/Users/christian/Documents/PhD/research/code/c/"
        self.so = CDLL('%s/lib_custom_gsl.so' %path)

    def sobol_seq(self,output_array, start_dim, seed):
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')
        array_2d_double = npct.ndpointer(dtype=numpy.double, ndim=2, flags='CONTIGUOUS')

        self.so.sobol_seq.argtypes = [c_int, c_int, c_int, c_int, array_2d_double]
        self.so.sobol_seq.restype = c_int
        return self.so.sobol_seq(len(output_array), len(output_array[0]), start_dim, seed, output_array)

    def halton_seq(self,output_array, start_dim, seed):
        array_1d_double = npct.ndpointer(dtype=numpy.double, ndim=1, flags='CONTIGUOUS')
        array_2d_double = npct.ndpointer(dtype=numpy.double, ndim=2, flags='CONTIGUOUS')

        self.so.halton_seq.argtypes = [c_int, c_int, c_int, c_int, array_2d_double]
        self.so.halton_seq.restype = c_int
        return self.so.halton_seq(len(output_array), len(output_array[0]), start_dim, seed, output_array)

    def reversehalton_seq(self,output_array, start_dim, seed):
        array_2d_double = npct.ndpointer(dtype=numpy.double, ndim=2, flags='CONTIGUOUS')

        self.so.reversehalton_seq.argtypes = [c_int, c_int, c_int, c_int, array_2d_double]
        self.so.reversehalton_seq.restype = c_int
        return self.so.reversehalton_seq(len(output_array), len(output_array[0]), start_dim, seed, output_array)
