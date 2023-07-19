import numpy as np
from scipy.optimize import newton
from scipy.stats import norm
import matplotlib.pyplot as plt
import time 
from mpmath import mp
from mpmath import findroot
import pyriemann as pr
mp.prec = 50



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


#####################################################################################################################################
#Sampling r



#Used in algorithm 3
def lng_small_sigma(r,n,sigma):
    d = n*(n+1)/2
    return (1/sigma)*((np.log(sigma)-r**2/(2*sigma**2)+r/2*sigma-np.log(2)))
    

#Used in algorithm 2
def lng(r,n,sg):
    d = int(n*(n+1)/2)
    
    return (-r**2/(2*sg**2))+(r/2)*(d-1)-np.log(2)*(d-1)

#Used in algorithm one 
def g(r, n, sg, dr):


    d = int(n*(n+1)/2)
    ##MP version
    g=[]
    sum = mp.mpf(0)
    for i in range(np.size(r)):
        temp = mp.fmul(mp.exp(mp.fdiv(-mp.power(r[i],2),2*mp.power(sg,2))),mp.power(mp.sinh(r[i]/2),d-1))
        sum = mp.fadd(sum,temp)
        g.append(temp)
    
    for i in range(np.size(g)):
        g[i] = mp.fdiv(g[i],sum*dr)
    
    return np.array(g)


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

    #Finding the zero of this function is equivalent to finding the zero of the first derivative of g
    zero_eq = lambda r :((d-1)*(sigma**2)*mp.cosh(r/2)/2) - r*mp.sinh(r/2) 
    #We use scipy Newton method to get an approximation of the zero
    #We obtain approx_x0 by using a Taylor approximation on zero_eq which gives a rough idea of the location of the solution
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
    return R_sample, (accepted_points/counter_points),g_zero


def r_sampling2(dim, sigma = 2, iter = 1):

    d = dim*(dim+1)/2
    cut = 10000
    #Basic case :
    approx_x0 = sigma**2*(d-1)/2

    if (n<=5 and sigma <= 7):

    
        zero_eq = lambda r : ((d-1)*(sigma**2)*np.cosh(r/2)/2) - r*np.sinh(r/2)


        g_zero = newton(zero_eq, approx_x0)
        if g_zero > 5*sigma :
            rr = np.linspace(g_zero - 5*sigma, g_zero + 5*sigma, cut)

        else :
            rr = np.linspace(0,g_zero+5*sigma,cut)
        dr = rr[1]-rr[0]

        gr = g(rr, dim, sigma,dr)

        imax = int((-rr[0]+g_zero)/dr)
        gmax = gr[imax]
        rmax = rr[imax]
        cmax = 1/sigma**2
        sc = 1.0 / np.sqrt(cmax)
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

    if (sigma >= 1 and n > 5) or (sigma <=1 and n >= 60):

        lnzero_eq = lambda r : np.log(d-1)+2*np.log(sigma)-np.log(r)+np.log(1/2)
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

    if ((n>=5 and n < 10 and sigma >= 1 and sigma <= 1) or (n>= 10 and n <= 50 and sigma <= 1 and sigma >= 0.1)) :
        zero_eq = lambda r :((d-1)*(sigma**2)*mp.cosh(r/2)/2) - r*mp.sinh(r/2) 
        #We use scipy Newton method to get an approximation of the zero
        #We obtain approx_x0 by using a Taylor approximaiton on zero_eq which gives a rough idea of the location of the solution

        approx_x0 = sigma**2* (d-1)/2



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
        return R_sample, (accepted_points/counter_points),g_zero
        #Last case for large n and small sigma, wip
    else :

        ln_zero_eq2_chg_vr = lambda r : np.log(sigma**2 * (r**2)/2)-np.log((1/2)*(d-1)*sigma**2)
        lnzero_eq = lambda r : np.log(d-1)+2*np.log(sigma)-np.log(r)+np.log(1/2)
        lng_zero = sigma * newton(ln_zero_eq2_chg_vr,approx_x0)
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





#########################################################################################################################################################
#########################################################################################################################################################
#CURS Algorithm


def riemannian_gaussian_matrix_sampling(dim, mean, sigma, amount = 1, accurate = False):
    nb_samples = 0
    samples = []
    mean_sqrt = pr.utils.base.sqrtm(mean)

    if (accurate):
        r = r_sampling(dim, sigma, amount)[0]
    else :
        r = r_sampling2(dim, sigma, amount)[0]
        
    for i in range(amount):
        S = s_sampling(dim)
        u = np.random.uniform(0.,1.)
        diag = np.linalg.eig(S)[0]
        detA = 1

        for k in range(dim):
            for j in range(dim):
                if j < k:
                    detA *= np.sinh(r[i] * np.abs(diag[k] - diag[j])/2)/((diag[k] - diag[j])/2)
        
        if (u <= np.abs(detA) / (np.sinh(r[i]))**(dim-1)):
            samples.append(mean_sqrt @ np.exp(r[i]*S) @ mean_sqrt)
            nb_samples+=1

    return samples







################TESTS###########################
##################################################
#############################################


sigma = 0.1
n = 6
d = n*(n+1)/2
lower = -10*sigma
upper = 10*sigma
approx = sigma**2* (d-1)/2

#Test with a given sigma and n

if approx < upper : 
    lower = -approx
if(1 == 0):
    #np.random.seed(50)
    rr = np.linspace(approx+lower, approx+upper, 10000)
    dr = rr[1]-rr[0]
    #test = r_sampling2(n,sigma, iter = 10000)
    test2 = r_sampling(n,sigma,iter = 10000)
    gr = g(rr,n,sigma,dr)
    #print(test[1])
    #plt.hist(test[0], bins = 20, color = "green", density = True)
    plt.hist(test2[0], bins = 20, color = "red", density = True)
    plt.plot(rr,gr)
    plt.show()

if(1 == 1):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    samples = pr.datasets.sample_gaussian_spd(n_matrices= 300,mean = np.array([[1,4,-5],[0,30,0],[2,0,10]]),sigma =0.5)
    print(samples)
    #3DPLOT for 2D case
    points = []
    points.append([])
    points.append([])
    points.append([])
    for i in range(0,300):
        diag = np.linalg.eig(samples[i])[1]
#        print("above")
        points[0].append(diag[0])
        points[1].append(diag[1])
        points[2].append(diag[2])

    ax.scatter(points[0],points[1],points[2], color = 'blue')

    samples = riemannian_gaussian_matrix_sampling(3,np.array([[1,4,-5],[0,30,0],[2,0,10]]),0.5,5000, True)
#    print("SECOND PART")
    #3DPLOT for 2D case
    points = []
    points.append([])
    points.append([])
    points.append([])
    for i in range(0,300):
        diag = np.linalg.eig(samples[i])[1]
#        print(diag)
        points[0].append(diag[0])
        points[1].append(diag[1])
        points[2].append(diag[2])

#    print(points[0])

    ax.scatter(points[0],points[1],points[2], color = 'red')
    ax.set_aspect('equal')
    plt.show()

