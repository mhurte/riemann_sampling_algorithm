import numpy as np
import matplotlib.pyplot as plt
import pyriemann.datasets as pr
import time
import accept_reject as ar

f = lambda x : np.exp(-(x*x)/2)*((np.sin(6*x)**2)+(3*(np.cos(x)**2)*(np.sin(4*x)**2))+1)

###Quick testing of the functions 
if(1==1) :
    sample = ar.accept_reject_algorithm(f,  M = 10,envelope_distribution ='exponential', nb_points=10000 ,lambda_exp= 5)
    print(sample[1])
    ar.sample_distribution_hist(sample[0])

    sample1 = ar.accept_reject_algorithm(f, M = 4.5, nb_points=10000 )
    print(sample1[1])
    ar.sample_distribution_hist(sample1[0])
