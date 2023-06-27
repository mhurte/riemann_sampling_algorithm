import numpy as np
from scipy.optimize import newton
from scipy.stats import norm
import matplotlib.pyplot as plt
import time 
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

def envelope_sampling(dim, sigma = 2, iter = 1) :
    d = dim*(dim+1)/2
    cut = 10000
    #Finding the zero of this function is equivalent to finding the zero of the first derivative 
    zero_eq = lambda r :((d-1)*(sigma**2)*np.cosh(r/2)/2) - r*np.sinh(r/2) 
    #We use scipy Newton method to get an approximation of the zero
    #We obtain approx_x0 by using a Taylor approximaiton on zero_eq which gives a rough idea of the location of the solution
    if sigma >= 0.9 :
        approx_x0 = sigma**2* (d-1)/2

    else : 
        approx_x0 = sigma * np.sqrt(d-1)

    gr1_zero = newton(zero_eq, approx_x0, maxiter = 100)

    if gr1_zero > 10*sigma :
        rr = np.linspace(gr1_zero - 10*sigma, gr1_zero + 10*sigma, cut)

    else :
        rr = np.linspace(0,gr1_zero+10*sigma,cut)
    dr = rr[1]-rr[0]
    gr = g(rr, dim, sigma, dr)


    imax = int((-rr[0]+gr1_zero)/dr)
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

a = time.time()


sigma = 4.1
n = 3
d = n*(n+1)/2
lower = -10*sigma
upper = 10*sigma
if (sigma <= 1)  : 
    approx = sigma * np.sqrt(d-1)
else :
    approx = sigma**2* (d-1)/2

if approx < upper : 
    lower = -approx

#Test with a given sigma and n
if(1 == 0):
    rr = np.linspace(approx+lower, approx+upper, 10000)
    dr = rr[1]-rr[0]
    gr = g(rr, n, sigma,dr )
    test = envelope_sampling(n,sigma, iter = 10000)
    b = time.time()-a
    print(test[1])
    print(b)
    plt.hist(test[0], bins = 20, density = True)
    plt.plot(rr,gr, color = "red")
    plt.show()



#Test for execution speed
if(1 == 0):
    execution_time = []
    rr = np.linspace(0.1,3.9,39)
    print(rr)
    for i in range (1,40):
        a = time.time()
        test = envelope_sampling(5,sigma = 0.1*i, iter = 10000)
        b = time.time()-a
        execution_time.append(b)
        print(i)

    plt.plot(rr,execution_time)
    plt.show()


#Test for acceptance rate
if(1 == 1):
    acceptance_rate = []
    rr = np.linspace(0.1,3,30)
    for i in range (1,31):
        a = time.time()
        test = envelope_sampling(5,sigma = 0.1*i, iter = 1000)
        b = time.time()-a
        acceptance_rate.append(test[1])
        print(i)

    plt.plot(rr,acceptance_rate)
    plt.show()