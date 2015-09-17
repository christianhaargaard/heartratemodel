""" Example of wrapping a C library function that accepts a C double array as
input using the numpy.ctypeslib. """

import numpy 
import numpy.ctypeslib as npct
import numpy.linalg
import matplotlib.pyplot as pyplot
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import pointer 

# Import the python SA tools
import sys
#sys.path.append('../python')
from sensitivity_analysis import sobol_weirs
from c_functions import solve, sens, spline; 



#######################################################################
#             Some functions for manipulating input data              #
#######################################################################

def add_sinangle(time = numpy.array([]), current=numpy.array([]), tilttime = 0, angle_dif=numpy.pi/3.0, duration = 14):
    print 'tilttime ='+str(tilttime)
    output = numpy.empty_like(time)
    output[:] = 0.0
    for i in xrange(numpy.size(time)):
        if time[i]<tilttime:
            output[i] = current[i]
        elif tilttime<time[i]<tilttime+duration:
            output[i] = current[i] + numpy.sin( angle_dif * (time[i]-tilttime)/duration)
        else:
            output[i] = current[i] + numpy.sin(angle_dif)
    return output

#######################################################################
#             Functions for generating parameter samples              #
#######################################################################

#def unif(n=1):
    #return numpy.random.uniform(size=n)

def unif(config = numpy.array([]),n=1):
    return config[0] + (config[1]-config[0])*numpy.random.uniform(size=n)
    #return numpy.random.uniform(low=config[0],high=config[1],size=n)

#def scaled_unif(min=0, max=1):
    #return min + (max-min)*unif

#######################################################################
#                    The code that runs the model                     #
#######################################################################

# Read data, and arange it into numpy arrays
data = numpy.genfromtxt('/Users/christian/Documents/PhD/research/data/standard_tilt/bp_martin_thinned.csv',delimiter=',')
t = numpy.array(data[:,0])
bp = numpy.array(data[:,1])

# and coefficients for interpolation
bp_coef = 0*numpy.empty_like(t)
spline(t,bp,bp_coef)

# Timeseries for sin to angle - and coefficients for interpolation
tiltangle = numpy.ones(numpy.shape(t))*0.0
tiltangle = add_sinangle(t,tiltangle,tilttime=181.0)
tiltangle_coef = 0*numpy.empty_like(t)
spline(t,bp,tiltangle_coef)

# Timeseries for respiration - and coefficients for interpolation
respiration = .25 + .1*numpy.sin(2*numpy.pi*numpy.copy(t)/10.0)
respiration_coef = 0*numpy.empty_like(t)
spline(t,bp,tiltangle_coef)

# Timepoints for which we wants output
n = 100000
step = (t[-1]-t[0])/n
tout = numpy.arange(t[0],t[-1],step)


# Number of equahiotn
NEQ = 11

# Memory for solution
y = numpy.zeros((NEQ+1)*numpy.size(tout))

# Set parameters
p = numpy.array([
70,         # p0
6,          # k
1.83,       # Am0
.201,       # a1
.075,       # a2
2.99,       # b1
.313,       # b2
3.5,        # s1
.05,        # s2
0,          # Tpm
1,          # TpM
5,          # xi
.5,         # fp
0,          # Tsm
1,          # TsM
5,          # eta
.5,         # fs
numpy.pi/6, # height_threshold
.44,         # am
1.4,        # aM
30,         # ak
.5,         # bm
3,          # bM
30,         # bk
.405,       # td
6.25,       # tA
15,         # qp
.3,         # kiN
33.3,       # tN
1.5,        # qs
.01,        # tAF
.35,        # mu
10,         # KA
2.5,        # tAS
2.22,       # tNS
20,         # KN
75,         # h0
33.5,         # hm
115         # hM
])

def metric(parameters = numpy.array([])):
    n = numpy.size(tout)
    y = numpy.zeros((NEQ+1)*n)
    result = solve(t,bp,bp_coef,tiltangle,tiltangle_coef,respiration,respiration_coef,tout,y,parameters)
    if result==0:
        hr = y[11*n:12*n]
        return numpy.linalg.norm(hr,2)
    else:
        print 'ignoring this result'
        return 0

# Set parameter intervals
p_interval = numpy.array([
[50, 120],         # p0
[3, 8],          # k
[1, 2.5],       # Am0
[.1, 1],       # a1
[.01, .1],       # a2
[2.5, 10],       # b1
[.25, 1],       # b2
[1, 10],        # s1
[0.2, .5],        # s2
[0, .5],          # Tpm
[.5, 1],          # TpM
[3,10],          # xi
[.2,.8],         # fp
[0,.5],          # Tsm
[.5,1],          # TsM
[3,10],          # eta
[.2,.8],         # fs
[.523, .523], # height_threshold
[0,.6],         # am
[.6,3],        # aM
[10,50],         # ak
[0,1],         # bm
[.8,5],          # bM
[10,50],         # bk
[0.1,6],       # td
[1, 40],       # tA
[2.5, 90],         # qp
[.05,1.8],         # kiN
[5, 200],       # tN
[.25, 9.0],        # qs
[.001, .1],        # tAF
[.15,75],        # mu
[2, 50],         # KA
[.5, 12.5],        # tAS
[.4, 10],       # tNS
[4, 100],         # KN
[60, 110],         # h0
[10, 60],         # hm
[90, 200]]         # hM
)

#p_interval = numpy.empty_like(p_interval)
#for i in xrange(len(p)):
    #p_interval[i][0] = .99*p[i]
    #p_interval[i][1] = p[i]

npars  = len(p)
sobolgenerators = numpy.array([])
for i in xrange(npars):
    sobolgenerators=numpy.append(sobolgenerators,unif)

# Solve model
solve(t,bp,bp_coef,tiltangle,tiltangle_coef,respiration,respiration_coef,tout,y,p)

#[Si, STi] = sobol_weirs(5000,sobolgenerators,p_interval,metric)


p0 = p[0]
k = p[1]
Am0 = p[2]
wall_strain = (1.0-numpy.sqrt((p0**k +bp**k)/(p0**k + Am0*bp**k)))

y0 = y[:n]              # Strain of both components in visco wall]
y1 = y[n:2*n]           # Strain of both components in visco wall

s1 = p[7]
s2 = p[8]
firing = s1*(numpy.lib.function_base.interp(tout,t,wall_strain)-y0)+s2;
Tpm = p[9]
TpM = p[10]
xi = p[11]
fP = p[12]
parasymp_baro = Tpm + (TpM -Tpm)*firing**xi/(firing**xi + fP**xi);

y2 = y[2*n:3*n]
y3 = y[3*n:4*n]
y4 = y[4*n:5*n]
y5 = y[5*n:6*n]
y6 = y[6*n:7*n]         # Unscaled value of delayed sympathetic signal.

td = p[24]
rho = 5;
factorial = 4*3*2;
delay_alpha = rho/td;
scale_constant = delay_alpha**rho/factorial;
symp_delayed = scale_constant*y6;# Scaled value of delayed sympathetic signal.

# Looks good this far :D

y7 = y[7*n:8*n]
y8 = y[8*n:9*n]
y9 = y[9*n:10*n]
y10 = y[10*n:11*n]
y11 = y[11*n:12*n]

#Stuff for testing sens
p[13] = .01
p[9] = .01
