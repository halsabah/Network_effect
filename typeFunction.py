# -*- coding: utf-8 -*-
"""
@author: David Hu 
"""

print('Loaded type function x*h(x)...')

# h(s) can be any increasing function of S on the interval [0,1]
def h(S):
    return S

def type_function(x, lumbda, S_current, S_past):
    return x * (lumbda*h(S_past) + (1-lumbda)*h(S_current))