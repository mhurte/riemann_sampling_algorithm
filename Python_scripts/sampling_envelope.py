import numpy as np
from scipy.optimize import newton
from scipy.stats import norm
import matplotlib.pyplot as plt

def g(r, n, sg, dr):
    d = int(n*(n+1)/2)
    gr = np.exp(-r**2/(2*sg**2)) * np.sinh(r/2)**(d-1)
    gr = gr / (np.sum(gr) * dr)    
    return gr

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

def envelope_sampling(dim, sigma = 2, iter = 1, newton_x0 = 10 ) :
    d = dim*(dim-1)/2

    #Finding the zero of this function is equivalent to finding the zero of the first derivative 
    zero_eq = lambda r :((d-1)*(sigma**2)*np.cosh(r/2)/2) - r*np.sinh(r/2) 
    #We use scipy Newton method to get an approximation of the zero
    gr1_zero = newton(zero_eq, 50, maxiter = 100)



    rr = np.linspace(0.01, +100, 10000)
    dr = rr[1]-rr[0]
    gr = g(rr, dim, sigma, dr)


    imax = int(gr1_zero/0.01)
    gmax = gr[imax]
    rmax = rr[imax]
    cmax = 1/sigma**2
    sc = 1.0 / np.sqrt(cmax)


    #Envelope to be used in the rejection sampling
    er = norm(loc=rmax, scale=sc).pdf(rr)
    M = gmax * np.sqrt(2*np.pi*sc**2)
    
    R_sample = []
    counter_points = 0
    while (counter_points < iter):
        uniform_sample = np.random.uniform(0.,1.)
        er_sample = norm.rvs(loc = rmax, scale = sc)
        index = int(er_sample/0.01)
        if(uniform_sample<=((gr[index])/(M*er[index]))):
           R_sample.append(er_sample)
           counter_points += 1
    return R_sample
        

    
rr = np.linspace(0.01, +100, 1000)
gr = g(rr, 4, 1, 0.01)
test = envelope_sampling(4, iter = 1000)
plt.hist(test, bins = 20)
plt.show()
