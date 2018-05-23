# -*- coding: utf-8 -*-
"""
@author: David Hu 
"""
import numpy as np
import matplotlib.pyplot as plt
from math import exp, log, ceil
from scipy.optimize import fsolve

from utility import u, k
from typeFunction import h, type_function
from distribution import G, invG
from parameters import lumbda, c


'''
Part 1: Solving and plotting the equilibria given the functions and parameters defined by user
'''

# f(S) defined for solving the equilibria
def f(S):
    return k(S, S, c, lumbda) - S

# Solver for equilibria
tolerance = 1e-10 # tolerance for checking whether the solution provided by the solver is valid
test = 0.001 # small perturbation used for determine whether the equilibrium is stable using the rule of thumb
guess = list(np.arange(0.0, 1.0, 0.02))
S_equilibrium = []
for x in guess:
    S_equilibrium.append(float(fsolve(f, x)))
S_l = min(S_equilibrium)
S_h = max(S_equilibrium)
for i in range (1, len(S_equilibrium)):
    if S_equilibrium[i] > S_l + tolerance and S_equilibrium[i] < S_h - tolerance:
        S_m = S_equilibrium[i]
        break
    else:
        S_m = 1.1 #some non-sense number that will be recognized as "not exist" by the following program
print(' ')
print('The Equilibria under the functions (utility, distribution, type) you defined is: ')
if abs(f(S_l)) > tolerance:
    print('S_low: does not exist')
elif f(S_l + test) < 0 and f(S_l - test) > 0:
    print('S_low: ' + str(S_l) + ', and it is STABLE')
else:
    print('S_low: ' + str(S_l) + ', and it is UNSTABLE')
    
if abs(f(S_m)) > tolerance:
    print('S_medium: does not exist')
elif f(S_m + test) < 0 and f(S_m - test) > 0:
    print('S_medium: ' + str(S_m) + ', and it is STABLE')
else:
    print('S_medium: ' + str(S_m) + ', and it is UNSTABLE')
    
if abs(f(S_h)) > tolerance:
    print('S_high: does not exist')
elif f(S_h + test) < 0 and f(S_h - test) > 0:
    print('S_high: ' + str(S_h) + ', and it is STABLE')
else:
     print('S_high: ' + str(S_h) + ', and it is UNSTABLE')


# Plot solver result
S_plot = np.arange(0.001, 1, 0.001)
n = len(S_plot)
KS_plot = []
for i in range (n):
    KS_plot.append(k(S_plot[i], S_plot[i], c, lumbda))
plt.plot(S_plot, KS_plot, 'b', label = 'k(S)')
plt.plot(S_plot, S_plot, 'r', label = 'S')
if abs(f(S_l)) <= tolerance:
    plt.plot(S_l, S_l, 'go')
    plt.annotate('S_l', (S_l, k(S_l, S_l, c, lumbda)), xytext = (S_l-0.02, k(S_l, S_l, c, lumbda)+0.07))
if abs(f(S_m)) <= tolerance:
    plt.plot(S_m, k(S_m, S_m, c, lumbda), 'go')
    plt.annotate('S_m', (S_m, k(S_m, S_m, c, lumbda)), xytext = (S_m-0.06, k(S_m, S_m, c, lumbda)+0.07))
if abs(f(S_h)) <= tolerance: 
    plt.plot(S_h, k(S_h, S_h, c, lumbda), 'go')
    plt.annotate('S_h', (S_h, k(S_h, S_h, c, lumbda)), xytext = (S_h-0.02, k(S_h, S_h, c, lumbda)-0.07))
plt.title('The Equilibria')
plt.legend()
plt.show(block=False)

