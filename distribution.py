# -*- coding: utf-8 -*-
"""
@author: David Hu 
"""

print('Loaded distribution function G(x)...')

from math import exp, log
from scipy.stats import norm

# The distribution identical to that in the lecture notes
def G(x):
    '''
    It is the CDF
    Restrictions on parameters of the distribution:
        alpha > 0
        0 < gamma < beta
    '''
    alpha = 5.0     # bigger alpha will force k (S) to converge to 1 more quickely
    gamma = 1.0     # bigger gamma will extend the concave part of k(S)
    beta = 100.0    # bigger beta will extend the convex part of k(S)
    
    return (gamma+beta)/(1+beta) * (1+beta)*exp(-alpha/x) / (gamma+beta*exp(-alpha/x))


def invG(U):
    '''
    It is the inverse function of G(x) and will be used to draw samples from the distribution via inverse transform
    Parameters should be the same as in G(x)
    Restrictions on parameters of the distribution:
        alpha > 0
        0 < gamma < beta
    '''
    alpha = 5.0     # bigger alpha will force k (S) to converge to 1 more quickely
    gamma = 1.0     # bigger gamma will extend the concave part of k(S)
    beta = 100.0    # bigger beta will extend the convex part of k(S)
    
    return -alpha / log(gamma*U/(gamma+beta-U*beta))


"""
# An example of another distribution - normal distribution
def G(x):
    return norm.cdf(x, loc = 3.0, scale = 1.0) # CDF of normal distribution, loc is mean, scale is standard deviation

def invG(U):
    return norm.ppf(U, loc = 3.0, scale = 1.0) # inverse of CDF of normal distribution, loc is mean, scale is standard deviation
    
"""