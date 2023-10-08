# riemann_sampling_algorithm

Development of an algorithm which should be faster than the current Monte-Carlo metho for the sampling of matrices according to the Gaussian-Riemannian distribution using different means of computations.
Currently working for most cases but the edge cases are leading to approximations errors that we are currently under the process of solving in order to be able to implement this new method in the pyRiemann library.
It is still possible to compute the matrices according to this distribution using an external library that deals with really large numbers and let us avoid using approximations but this is much slower and doesn't let us compute faster, but shows us that our method gives proper results. If we can manage to find a better approximations depending on the parameters, then this method would allow for much faster computation in every scenario.
