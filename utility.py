# -*- coding: utf-8 -*-
"""
@author: David Hu 
"""

'''
Note: 
    Utility function is not explicitly used in the main file in this version. 
    Nonetheless, all the analyses were based on and derived from this function.
    If wanted, it can be used to calculate the agents' utility in Part 2-1 of the main file.
'''

print('Loaded utility function...')

from typeFunction import h, type_function
from distribution import G

def u(s, x, S_current, S_past, lumbda, c):
    return (type_function(x, lumbda, S_current, S_past) - c) * s

# This k function is the equilibrium condition derived from the utility function above. 
# If other utility functions are defined, this condition function should also be changed accordingly. 
def k(S_current, S_past, c, lumbda):
    if (lumbda*h(S_past) + (1-lumbda)*h(S_current)) == 0:
        return 0
    else:
        return 1 - G(c/(lumbda*h(S_past) + (1-lumbda)*h(S_current)))