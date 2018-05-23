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

'''
Part 2: User to choose which function to use
'''
choice = input('Please select which tool you want to use (For finite agents simulation, please enter "agent"; For analytical form analysis, please enter "ana"; If you want to do both, enter "both"): ')



'''
Part 2-1: Finite agents simulation
'''

if (choice == 'agent') or (choice == 'both'):
    print(' ')
    print('-----------------------------------')
    print('      Finite Agents Simulation     ')
    print('-----------------------------------')
    
    
    num_agents = int(input('Please enter the number of agents in this simulation (suggested: 100000): '))
    S0 = input('Please enter an initial state of the fraction that applies strategy s=1 (if you want to test the equilibria solved above, enter "S_l", "S_m", or "S_h"): ')
    epsilon = input('Please enter a small perturbation from the current state (can be either positive of negative): ')
    
    if S0 == 'S_l':
        S0 = S_l
    elif S0 == 'S_m':
        S0 = S_m
    elif S0 == 'S_h':
        S0 = S_h
    
    S0 = float(S0)
    epsilon = float(epsilon)
    
    # save the input values in case the user wants to test the same numbers in part 2-2
    S0_save = S0
    epsilon_save = epsilon
    
    print(' ')
    print('The number of agents you entered is: ' + str(num_agents))
    print('The initial state you entered is: ' + str(S0))
    print('The perturbation you entered is: ' + str(epsilon))
    
    maxIter = int(input('Now enter the maximum iterations that you wish the program to perform to get started (suggested: 100): '))
    print(' ')
    
    S = []
    S.append(S0)
    S.append(S0 + epsilon)
    
    # sample the types of the agents from G(x) using inverse transform
    x = []
    for j in range (num_agents):
        U = np.random.uniform()
        x.append(invG(U))
    
    # set up the initial state
    s = [1]*ceil(num_agents*(S0+epsilon)) + [0]*(num_agents - ceil(num_agents*(S0+epsilon)))
    
    # main iteration
    for i in range (2, maxIter+2):
        for j in range (num_agents):
            type_j = type_function(x[j], lumbda, S[i-1], S[i-2]) - c
            if type_j > 0:
                s[j] = 1
            else:
                s[j] = 0
        
        S.append(sum(s) / num_agents)
        print('Number of agents choosing s=1 at time ' + str(i-1) + ' is: ' + str(S[i]))
        if (S[i] == S[i-1]) and (S[i-1] == S[i-2]):
            break
        
    # print the result
    print(' ')
    print('The equilibrium state is: ' + str(S[-1]))
    print('The time spent to achieve the equilibirum is: ' + str(len(S)-4))


'''
Part 2-2: Analytical form analysis
'''
if (choice == 'ana') or (choice == 'both'):

    print(' ')
    print('-----------------------------------')
    print('      Analytical Form Analysis     ')
    print('-----------------------------------')
    
    if choice == 'ana':
        S0 = input('Please enter an initial state of the fraction that applies strategy s=1 (if you want to test the equilibria solved above, enter "S_l", "S_m", or "S_h"): ')
        epsilon = input('Please enter a small perturbation from the current state (can be either positive of negative): ')
    
    if choice == 'both':
        S0 = input('Please enter an initial state of the fraction that applies strategy s=1 (if you want to test the equilibria solved above, enter "S_l", "S_m", or "S_h"; if you want to use the same number as before, please enter "same"): ')
        epsilon = input('Please enter a small perturbation from the current state (can be either positive of negative; if you want to use the same number as before, please enter "same"): ')
    
    if S0 == 'S_l':
        S0 = S_l
    elif S0 == 'S_m':
        S0 = S_m
    elif S0 == 'S_h':
        S0 = S_h
    elif S0 == 'same':
        S0 = S0_save
    
    if epsilon == 'same':
        epsilon = epsilon_save
    
    S0 = float(S0)
    epsilon = float(epsilon)
    
    print(' ')
    print('The initial state you entered is: ' + str(S0))
    print('The perturbation you entered is: ' + str(epsilon))
    
    maxIter = int(input('Now enter the maximum iterations that you wish the program to perform to get started (suggested: 100): '))
    
    
    S = [];
    S.append(S0)
    S.append(S0 + epsilon)
    
    
    # plot initial state
    plt.plot(S_plot, KS_plot, 'b', label = 'k(S)')
    plt.plot(S_plot, S_plot, 'r', label = 'S')
    plt.plot(S[0], S[0], 'ko')
    plt.title('Initial State (t = 0)')
    plt.legend()
    plt.show(block=False)
    
    # plot perturbed initial state
    plt.plot(S_plot, KS_plot, 'b', label = 'k(S)')
    plt.plot(S_plot, S_plot, 'r', label = 'S')
    plt.plot(S[0]+epsilon, S[0]+epsilon, 'yo')
    plt.title('Perturbed Initial State (t = 0)')
    plt.legend()
    plt.show(block=False)
    
    # main iteration
    for i in range (2, maxIter+2):
        best_response = k(S[i-1], S[i-2], c, lumbda)
        S.append(best_response)
        
        # plotting
        plt.plot(S_plot, KS_plot, 'b', label = 'k(S)')
        plt.plot(S_plot, S_plot, 'r', label = 'S')
        plt.plot(S[i], S[i], 'yo')
        plt.title('Intermediate State (t = ' + str(i-1) + ')')
        plt.legend()
        plt.show(block=False)
        
        # check whether already in equilibrium
        if S[i] == S[i-1]:
            plt.plot(S_plot, KS_plot, 'b', label = 'k(S)')
            plt.plot(S_plot, S_plot, 'r', label = 'S')
            plt.plot(S[i], S[i], 'go')
            plt.title('Equilibrium State (t = ' + str(i-2) + ')')
            plt.show(block=False)
            break
    
    # print the result
    print('The equilibrium state is: ' + str(S[-1]))
    print('The time spent to achieve the equilibirum is: ' + str(len(S)-3))