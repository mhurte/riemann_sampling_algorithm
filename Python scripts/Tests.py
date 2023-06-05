
import numpy as np
import matplotlib.pyplot as plt
import pyriemann.datasets as pr
import time
#import random
#import scipy as sp

################################################################################################################################################


##################
##################
#Replace the 0 by 1 in the if statement to test the different parts of the code 
#The "Score" used throughout the code is simply the execution time / acceptance rate and gives a rough idea of the efficiency of sampling algorithms for a similar number of points



#Stuff
if (1==0) : 
    t = np.arange(-5., 5., 0.01)
    plt.plot(t,np.exp(-(t*t)/2)*((np.sin(6*t)**2)+(3*(np.cos(t)**2)*(np.sin(4*t)**2))+1),'r',t,5*np.exp(-(t**2)/2),'b')
    plt.show()
    y1 = [np.exp(-(t**2)/2)*((np.sin(6*t)**2)+(3*(np.cos(t)**2)*(np.sin(4*t)**2))+1) for t in t]
    y2 = [5*np.exp(-(t**2)/2) for t in t]
    plt.plot(t,y1,t,y2)
    plt.show()






#Uniform values within a certain distribution naively by computing random values uniformly within a bigger rectangle first
if (1==0) : 
    nb_points = 10000
    leften_edge = -5.
    righten_edge = 5.
    time0 = time.time()
    rand = np.arange(leften_edge,righten_edge,(righten_edge-leften_edge)/nb_points)
    unif = [np.random.uniform(0,4) for x in rand]
    size = np.arange(0.,10000.,1.)
    sizeint = [int(x) for x in size]
    #plt.scatter(rand,unif)
    a=[(rand[i],unif[i]) for i in sizeint]

    #print(a)
    b=[]
    for i in range(len(a)) :
        if a[i][1]<=np.exp(-(a[i][0]*a[i][0])/2)*((np.sin(6*a[i][0])**2)+(3*(np.cos(a[i][0])**2)*(np.sin(4*a[i][0])**2))+1) :
            b.append(a[i])

    AcceptanceRate = len(b)/len(a)
    print(AcceptanceRate)
    execution_time = time.time() - time0
    print(execution_time)


    #This is an arbitrary score I use to compare my methods, basically, the lower the better
    print("executiontime*points/accuracy")
    print((execution_time*10000)/AcceptanceRate)
    plt.scatter(*zip(*b))
    plt.show()




############# Sampling with pyRiemann


#size0 = np.arange(0.,1000.,1.)
#size = [int(i) for i in size0]
#test = pr.generate_random_spd_matrix(n_dim = 10)
#print(test)






#################
if(1==0) :

    size0 = np.arange(0.,10000.,1.)
    size = [int(i) for i in size0]
    a = [np.random.normal(1,1) for i in size]
    plt.hist(a, bins = 100)
    plt.show()





#Normal values random
if(1==0) :

    time0 = time.time()
    nb_points = 10000
    uniform_simulated = []
    normal_simulated=[]

    for i in range(0, int(nb_points/2.)) :

        a = np.random.uniform(0.,1.)
        b = np.random.uniform(0.,1.)

        uniform_simulated.append(4.5*a)
        uniform_simulated.append(4.5*b)

        normal_simulated.append(np.sqrt(-2.*np.log(a))*np.sin(2.*np.pi*b))
        normal_simulated.append(np.sqrt(-2.*np.log(a))*np.cos(2.*np.pi*b))


    #plt.hist(normal_simulated,bins = 200,color='b')
    #plt.hist(uniform_simulated, bins = 200)
    x2=[]
    y2=[]
    for i in range(len(normal_simulated)) :
        if uniform_simulated[i]<=np.exp(-(normal_simulated[i]*normal_simulated[i])/2)*((np.sin(6.*normal_simulated[i])**2.)+(3.*(np.cos(normal_simulated[i])**2.)*(np.sin(4.*normal_simulated[i])**2.))+1.) :
            y2.append(uniform_simulated[i])
            x2.append(normal_simulated[i])

    AcceptanceRate = len(x2)/len(uniform_simulated)
    print(AcceptanceRate)

    execution_time = time.time() - time0
    print(execution_time)

    print("(executiontime*points)/accuracy")
    print((execution_time*nb_points)/AcceptanceRate)
    t = np.arange(-3.,3.,0.01)

    plt.plot(t,np.exp(-(t*t)/2.)*((np.sin(6.*t)**2.)+(3.*(np.cos(t)**2.)*(np.sin(4.*t)**2.))+1.),'r--')
    plt.scatter(x2,y2)
    plt.show()


    

###Performance test normal distribution numpy

if(1==0) :
    time0 = time.time()
    nb_points = 10000
    normal_simulated=[]
    uniform_simulated = []

    for i in range(0, int(nb_points)) :
        normal_simulated.append(np.random.normal(0.,1.))
        uniform_simulated.append(np.random.uniform(0.,1.))

    
    execution_time = time.time() - time0
    #print(execution_time)

    #print("(executiontime*points)/accuracy")
    #print((execution_time*nb_points)/AcceptanceRate)


    #plt.hist(normal_simulated,bins = 200,color='b')
    #plt.hist(uniform_simulated, bins = 200)
    #plt.show()



####Accept reject  (Most optimized so far )
##Given any distribution f, given g (easy to compute), such that we can find M>=1 st f(x)<= Mg(x) for any x ...
#ref : 2.3.2 
if(1==0) :
    time0 = time.time()
    nb_points = 10000
    leften_edge = -5.
    righten_edge = 5.
    M = 5.
    t = np.arange(leften_edge,righten_edge,(righten_edge-leften_edge)/nb_points)
    f = lambda x : np.exp(-(x*x)/2)*((np.sin(6*x)**2)+(3*(np.cos(x)**2)*(np.sin(4*x)**2))+1)
    g = lambda x : np.exp(-(x**2)/2)
    f_simulated1 = []
    for i in range(0, int(nb_points)) :

        uniform_sample = np.random.uniform(0.,1.)
        g_sample = np.random.normal(0,1)
        if (uniform_sample<=(f(g_sample)/(M*g(g_sample)))) :
            f_simulated1.append(g_sample)
    
    execution_time = time.time() - time0
    print(M,"execution_time",execution_time)

    AcceptanceRate = (len(f_simulated1)/nb_points)
    print(M,"AcceptanceRate",AcceptanceRate)

    score1 = (execution_time)/AcceptanceRate
    print(M,"(executiontime)/accuracy",score1)

    
    #plt.hist(f_simulated1,bins = 200)
    #plt.title("M=5")
    #plt.show()

    time0 = time.time()
    M=4.5
    f_simulated2 = []
    for i in range(0, int(nb_points)) :

        uniform_sample = np.random.uniform(0.,1.)
        g_sample = np.random.normal(0,1)
        if (uniform_sample<=(f(g_sample)/(M*g(g_sample)))) :
            f_simulated2.append(g_sample)

    execution_time = time.time() - time0
    print(M,"execution_time",execution_time)

    AcceptanceRate = (len(f_simulated2)/nb_points)
    print(M,"AcceptanceRate",AcceptanceRate)

    score2 = (execution_time)/AcceptanceRate
    print(M,"(executiontime)/accuracy",score2)

    #plt.hist(f_simulated2,bins = 200)
    #plt.title("M=4.5")
    #plt.show()

    time0 = time.time()
    M=4.
    f_simulated3 = []
    for i in range(0, int(nb_points)) :

        uniform_sample = np.random.uniform(0.,1.)
        g_sample = np.random.normal(0,1)
        if (uniform_sample<=(f(g_sample)/(M*g(g_sample)))) :
            f_simulated3.append(g_sample)
    


    execution_time = time.time() - time0
    print(M,"execution_time",execution_time)

    AcceptanceRate = (len(f_simulated3)/nb_points)
    print(M,"AcceptanceRate",AcceptanceRate)

    score3 = (execution_time)/AcceptanceRate
    print(M,"(executiontime)/accuracy",score3)

    #plt.hist(f_simulated3,bins = 200)
    #plt.title("M=4")
    #plt.show()



    time0 = time.time()
    M=2.
    f_simulated4 = []
    for i in range(0, int(nb_points)) :

        uniform_sample = np.random.uniform(0.,1.)
        g_sample = np.random.normal(0,1)
        if (uniform_sample<=(f(g_sample)/(M*g(g_sample)))) :
            f_simulated4.append(g_sample)
    


    execution_time = time.time() - time0
    print(M,"execution_time",execution_time)

    AcceptanceRate = (len(f_simulated4)/nb_points)
    print(M,"AcceptanceRate",AcceptanceRate)

    score4 = (execution_time)/AcceptanceRate
    print(M,"(executiontime)/accuracy",score4)

    #plt.hist(f_simulated4,bins = 200)
    #plt.title("M=2")
    #plt.show()



    plt.plot(t,f(t))
    plt.plot(t,5.*g(t),t,4.5*g(t),t,4.*g(t),t,3.*g(t))

    plt.show()

    ##Clean plots
    figure, axis = plt.subplots(2, 4)
    
    ##histo
    figure.set_size_inches((10,10))
    figure.suptitle("Comparison of the different returned distributions given choices of M")
    axis[0, 0].hist(f_simulated1,bins = 100,color = "r")
    axis[0, 0].set_title("M=5,s={:.4f}".format(score1))
  

    axis[0, 1].hist(f_simulated2,bins = 100,color = "g")
    axis[0, 1].set_title("M=4.5,s={:.4f}".format(score2))
  

    axis[0, 2].hist(f_simulated3,bins = 100, color = "b")
    axis[0, 2].set_title("M=4,s={:.4f}".format(score3))
  

    axis[0, 3].hist(f_simulated4,bins = 100, color = "y")
    axis[0, 3].set_title("M=2,s={:.4f}".format(score4))

    ##Functions
    axis[1, 0].plot(t,f(t),"black",t,5.*g(t), "r")
    axis[1, 0].set_title("M=5")

    axis[1, 1].plot(t,f(t),"black",t,4.5*g(t),"g")
    axis[1, 1].set_title("M=4.5")

    axis[1, 2].plot(t,f(t),"black",t,4.*g(t), "b")
    axis[1, 2].set_title("M=4")

    axis[1, 3].plot(t,f(t),"black",t,2.*g(t),"y")
    axis[1, 3].set_title("M=2")
    #The lower the choice of M, the higher the performance
    plt.show()

###Box Muller Normal distribution (2nd version)
if(1==0) :
    time0 = time.time()
    nb_points = 10000
    u1 = 2
    u2 = 0
    normal_simulated = []
    for i in range (nb_points) :
        u1 = 2
        while(u1**2+u2**2 > 1.) :

            u1 = np.random.uniform(-1.,1.)
            u2 = np.random.uniform(-1.,1.)

        s = u1**2 + u2**2
        z = np.sqrt((-2*np.log(s))/s)
        normal_simulated.append(z*u1)
        normal_simulated.append(z*u2)
        

    plt.hist(normal_simulated,bins = 100)
    plt.show()


################Box-Muller(3rd version) 
if(1==0) :
    time0 = time.time()
    nb_points = 10000

    normal_simulated = []
    for i in range (nb_points) :
        y1 = 0
        y2 = 0
        while(y2 <= ((1.-y1)**2)/2) :

            y1 = np.random.exponential(1.)
            y2 = np.random.exponential(1.)
        u = np.random.uniform(0.,1.)
        if (u < 0.5) :
            normal_simulated.append(y1)
        else : 
            normal_simulated.append(-y1)
        
    plt.hist(normal_simulated,bins = 100)
    plt.show()

####Student's t distribution algorithm

if(1==0) :
    time0 = time.time()
    nb_points = 100000
    nu = 10
    std_simulated = []
    counterpoints = 0
    while(counterpoints<nb_points) :
        u1 = np.random.uniform(0.,1.)
        u2 = np.random.uniform(0.,1.)

        if (u1 < 0.5) :
            x = 1/(4*u1-1)
            v = (x**-2)*u2
        else : 
            x = 4*u1 - 3
            v = u2

        if(v<(1 -(abs(x)/2)) or v < (1+(x**2/nu))**(-(nu+1)/2)) :
            std_simulated.append(x)
            counterpoints+=1


    plt.xlim([-5,5])
    plt.hist(std_simulated,bins = 100, range = [-5,5])
    plt.show()

#####Cheng and Feast's gamma 

if(1==1) :
    time0 = time.time()
    nb_points = 10000
    alpha = 10
    ga_simulated = []
    counterpoints = 0


    c1 = alpha -1
    c2 = (alpha - (1/(6*alpha)))/c1
    c3 = 2/c1
    c4 = 1 + c3
    c5 = 1/np.sqrt(alpha)
    while(counterpoints<nb_points) :
        u1=-1

        
        while(u1<=0 or u1 >= 1) :
            u1 = np.random.uniform(0.,1.)
            u2 = np.random.uniform(0.,1.)
            if(alpha>2.5) :
                u1 = u2 + c5*(1-1.86*u1)
        
        w =c2*u2/u1
        if(c3*u1+w+(w**-1) <= c4 or c3*np.log(u1)-np.log(w)+w <= 1):
            ga_simulated.append(c1*w)
            counterpoints+=1





    #plt.xlim([-8,8])
    plt.hist(ga_simulated,bins = 100)
    plt.show()

###########Ahren's and Dieter's gamma
#####This algorithm seems to be working for smaller values of alpha but doesn't seem like a gamma distribution at all for gamma >= 3 
if(1==1) :
    time0 = time.time()
    nb_points = 1000
    alpha = 10
    ga_simulated = []
    counterpoints = 0
    while(counterpoints<nb_points) :
        u1 = np.random.uniform(0.,1.)
        u2 = np.random.uniform(0.,1.)

        if (u1 > np.e/(np.e + alpha) ) :
            x = -np.log(((alpha + np.e)*(1-u1))/(alpha*np.e))
            y = x**(alpha - 1)
        else : 
            x = (((alpha + np.e)*u1)/np.e)**(1/alpha)
            y = np.exp(-x)

        if (u2 < y) :
            ga_simulated.append(x)
            counterpoints+=1






    plt.hist(ga_simulated,bins = 100)
    plt.show()