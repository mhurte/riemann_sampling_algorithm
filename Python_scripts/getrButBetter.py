import numpy as np
from scipy.optimize import newton
from scipy.stats import norm
import matplotlib.pyplot as plt
import time 
from mpmath import mp
from mpmath import findroot
mp.prec = 50

def lng(r,n,sg):
    d = int(n*(n+1)/2)
    #loggr = (-r**2/(2*sg**2))+np.log(np.sinh(r/2))*(d-1)
    loggr = (-r**2/(2*sg**2))+(r/2)*(d-1)-np.log(2)*(d-1)
    return loggr
def g(r, n, sg, dr):


    d = int(n*(n+1)/2)
    #gr = np.exp(-r**2/(2*sg**2)) * np.sinh(r/2)**(d-1)
    #loggr = (-r**2/(2*sg**2))+np.log(np.sinh(r/2))*(d-1)
    #gr = (np.exp(loggr))
    #gr = gr / (np.sum(gr) * dr)
    
    ##MP version
    g=[]
    sum = mp.mpf(0)
    for i in range(np.size(r)):
        temp = mp.fmul(mp.exp(mp.fdiv(-mp.power(r[i],2),2*mp.power(sg,2))),mp.power(mp.sinh(r[i]/2),d-1))
        sum = mp.fadd(sum,temp)
        g.append(temp)
    
    #print(g[30])
    #print(sum)
    for i in range(np.size(g)):
        g[i] = mp.fdiv(g[i],sum*dr)
    
    return np.array(g)
    #return gr

    
def g1(r, n, sg, dr): # first derivative
    d = int(n*(n+1)/2)
    gr = g(r, n, sg, dr)
    gr1 = np.multiply(gr, -r/sg**2 + 1/2*(d-1)*np.cosh(r/2)/np.sinh(r/2))   
    return gr1

def g2(r, n, sg, dr): # second derivative
    d = int(n*(n+1)/2)
    gr = g(r, n, sg, dr)
    gr1 = g1(r, n, sg, dr)
    gr2 = np.multiply(gr1, -r/sg**2 + 1/2*(d-1)*np.cosh(r/2)/np.sinh(r/2))
    gr2 = gr2 - np.multiply(gr, 1/sg**2 + 1/4*(d-1)*1.0/np.sinh(r/2)**2)
    return gr2



def generate_unit_norm_symmetric_matrices(n_dim):
    n_matrices = int(n_dim * (n_dim + 1) / 2)
    S = np.zeros((n_matrices, n_dim, n_dim))
    k = 0
    for i in range(n_dim):
        for j in range(i, n_dim):
            if i == j:
                S[k,i,i] = 1
            else:
                S[k,i,j] = 1/np.sqrt(2)
                S[k,j,i] = 1/np.sqrt(2)
            k = k + 1
    return S


def s_sampling(n_dim):

    # dimensionality of the space of SPD matrices
    d = int(n_dim*(n_dim+1)/2)
    
    # - (1) generate random unit norm vector (i.e. a direction)
    z = np.random.randn(d)
    z = z / np.linalg.norm(z)
    
    # - (2) list a set of base vectors for the space of symmetric matrices
    base = generate_unit_norm_symmetric_matrices(n_dim=n_dim)
    
    # - (3) create the symmetric matrix coming from uniform distribution
    s = sum([zi * bi for zi, bi in zip(z, base)])
    
    return s



##################################################
##################################################
##################################################
######### FIRST VERSION, REQUIRES MPMATH #########
##################################################
##################################################
##################################################


def r_sampling(dim, sigma = 2, iter = 1) :
    d = dim*(dim+1)/2
    cut = 100000

        #Finding the zero of this function is equivalent to finding the zero of the first derivative 
    zero_eq = lambda r :((d-1)*(sigma**2)*mp.cosh(r/2)/2) - r*mp.sinh(r/2) 
    #We use scipy Newton method to get an approximation of the zero
    #We obtain approx_x0 by using a Taylor approximaiton on zero_eq which gives a rough idea of the location of the solution
    if (sigma >= 0.1):
        approx_x0 = sigma**2* (d-1)/2

    else: 
        approx_x0 = sigma * np.sqrt(d-1)

    #g_zero = newton(zero_eq, approx_x0, maxiter = 1000)
    g_zero = mp.findroot(zero_eq,approx_x0)
    g_zero = float(g_zero)

    if g_zero > 10*sigma :
        rr = np.linspace(g_zero - 10*sigma, g_zero + 10*sigma, cut)

    else :
        rr = np.linspace(0,g_zero+10*sigma,cut)
    dr = rr[1]-rr[0]
    gr = g(rr, dim, sigma, dr)


    imax = int((-rr[0]+g_zero)/dr)
    gmax = gr[imax]
    rmax = rr[imax]
    cmax = 1/sigma**2
    sc = 1.0 / np.sqrt(cmax)


    #Envelope to be used in the rejection sampling
    er = norm(loc=rmax, scale=sc).pdf(rr)
    M = gmax * np.sqrt(2*np.pi*sc**2)
    
    R_sample = []
    counter_points = 0
    accepted_points = 0
    while (accepted_points < iter):
        uniform_sample = np.random.uniform(0.,1.)
        er_sample = norm.rvs(loc = rmax, scale = sc)
        index = int((-rr[0] + er_sample)/dr)
        counter_points+=1
        if(uniform_sample<=((gr[index])/(M*er[index]))):
           R_sample.append(er_sample)
           accepted_points += 1
    return R_sample, (accepted_points/counter_points)

##################################################
##################################################
##################################################
#########        Upgraded version        #########
##################################################
##################################################
##################################################


def r_sampling2(dim, sigma = 2, iter = 1) :
    d = dim*(dim+1)/2
    cut = 100000

        #Finding the zero of this function is equivalent to finding the zero of the first derivative 
    zero_eq = lambda r :((d-1)*(sigma**2)*mp.cosh(r/2)/2) - r*mp.sinh(r/2) 

    #We will only use the classic zero_eq if n<3 and sigma <3 else the following equation. 
    #In the following we approximate both sinh and cosh for exp as the values of r are great enough which greatly simplifies the global expression

    lnzero_eq = lambda r : np.log(d-1)+2*np.log(sigma)-np.log(r)+np.log(1/2)
    #We use scipy Newton method to get an approximation of the zero
    #We obtain approx_x0 by using a Taylor approximaiton on zero_eq which gives a rough idea of the location of the solution
    if (sigma <= 0.01):
        approx_x0 = sigma**2* (d-1)/2

    else: 
        approx_x0 = sigma * np.sqrt(d-1)

    #g_zero = newton(zero_eq, approx_x0, maxiter = 1000)
    lng_zero = newton(lnzero_eq,approx_x0)

    if lng_zero > 5*sigma :
        rr = np.linspace(lng_zero - 5*sigma, lng_zero + 5*sigma, cut)

    else :
        rr = np.linspace(0,lng_zero+5*sigma,cut)
    dr = rr[1]-rr[0]
    lngr = lng(rr, dim, sigma)


    imax = int((-rr[0]+lng_zero)/dr)
    lngrmax = lngr[imax]
    rmax = rr[imax]
    cmax = 1/sigma**2
    sc = 1.0 / np.sqrt(cmax)


    #Envelope to be used in the rejection sampling
    er = norm(loc=rmax, scale=sc).pdf(rr)
    erlog = np.log(er)
    M = lngrmax * np.sqrt(2*np.pi*sc**2)
    
    R_sample = []
    counter_points = 0
    accepted_points = 0
    while (accepted_points < iter):
        uniform_sample = np.random.uniform(0.,1.)
        er_sample = norm.rvs(loc = rmax, scale = sc)
        index = int((-rr[0] + er_sample)/dr)
        counter_points+=1
        if(np.log(uniform_sample)<=((lngr[index])-(M*erlog[index]))):
           R_sample.append(er_sample)
           accepted_points += 1
    return R_sample, (accepted_points/counter_points)



a=time.time()
sigma = 0.2
n = 10000
d = n*(n+1)/2
lower = -10*sigma
upper = 10*sigma
if (sigma <= 0.01)  : 
    approx = sigma * np.sqrt(d-1)
else :
    approx = sigma**2* (d-1)/2
    print("approx")
    print(approx)
if approx < upper : 
    lower = -approx

#Test with a given sigma and n
if(1 == 1):
    rr = np.linspace(approx+lower, approx+upper, 10000)
    dr = rr[1]-rr[0]
    test = r_sampling2(n,sigma, iter = 10000)
    gr = g(rr,n,sigma,dr)
    print(test[1])
    plt.hist(test[0], bins = 20, density = True)
    plt.plot(rr,gr)
    plt.show()