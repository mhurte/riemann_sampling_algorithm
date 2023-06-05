############Cleaner version of the rejection algorithm with regard to the contents of Tests.py with easy benchmarking possibilities
############For making this easier, I will use mostly numpy functions for generating the different distributions




import numpy as np
import matplotlib.pyplot as plt
import pyriemann.datasets as pr
import time
#import random
#import scipy as sp






"""Returns a sample from a given distribution using a rejection method.

    Parameters
    ----------
    distribution : function
        density function of the distribution we attempt to sample from
    envelope_distribution : string, default='normal'
        density function (supposed easy to compute) we will use to accelerate the computation of the sample 
    M : float, default=1.0
        Coefficient for adjusting the envelope distribution
    nb_points : int, default=1000
        Number of points that we want in the sample

    Returns
    -------
    Sample : ndarray, shape (1, nb_points)
        Numpy array of samples

"""


def accept_reject_algorithm (density,
                            M = 1.0,
                            envelope_distribution = 'normal',
                            nb_points = 1000,
                            efficiency_computation = True) :
    time0 = time.time()
    sample = []

    if(efficiency_computation == False ) :
        if (envelope_distribution == 'normal' ) :
            normal_density = lambda x : np.exp(-(x**2)/2)
            for i in range(0, nb_points) :

                uniform_sample = np.random.uniform(0.,1.)
                g_sample = np.random.normal(0,1)
                if (uniform_sample<=(density(g_sample)/(M*normal_density(g_sample)))) :
                    sample.append(g_sample)
    
        return sample
    
    else :



        if (envelope_distribution == 'normal' ) :
            normal_density = lambda x : np.exp(-(x**2.)/2.)
            counter_points = 0
            nb_points_generated = 0
            while (counter_points < nb_points) :
                nb_points_generated +=1
                uniform_sample = np.random.uniform(0.,1.)
                g_sample = np.random.normal(0.,1.)
                if (uniform_sample<=(density(g_sample)/(M*normal_density(g_sample)))) :
                    sample.append(g_sample)
                    counter_points+=1

        time_exec = time.time() - time0
        acceptance_rate = nb_points / nb_points_generated
        efficiency = time_exec / acceptance_rate
        return sample, efficiency


def sample_distribution(sample, bins = 200) :

    plt.hist(sample,bins)
    plt.title("Distribution")
    plt.show()

