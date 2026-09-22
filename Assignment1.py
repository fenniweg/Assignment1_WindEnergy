# Assignment 1  
#import libraries
import numpy as np
import matplotlib.pyplot as plt
#import the BEM_algorithm function from bem.py
from bem import BEM_algorithm


# Question 1
# Compute the highest obtainable power coefficient, Cp,max (λmax, θp,max).
# Using BEM two methods: the polynomial method and the Madsen et al. method. 
# Compare the results of the two methods by plotting the power coefficient, Cp, as a function of the tip speed ratio, λ, for a range of pitch angles, θp. 

#Define range of tip speed ratios (lambda) and pitch angles (theta_p) to test
tip_speed_ratio = [5,6,7,8,9,10]
theta_p =[-3,-2,-1,0,1,2,3,4]

method = 'Polynomial'  # or 'Madsen' remov later and loop 
#Outer loop over methods, inner loop over tip speed ratios and pitch angles
#Loop over all combinations of tip speed ratios and pitch angles to compute thrust, torque, and power using BEM_algorithm
#for method in ['Polynomial','Madsen']:
thrust_list = []
torque_list = []
power_list = []
tip_speed_list = []
theta_list = []
for s in tip_speed_ratio:
        for theta in theta_p:
            print(f"Tip speed ratio: {s}, Pitch angle: {theta}")
            print(f"Method: {method}")
            thrust, torque,Power= BEM_algorithm(s,theta,method)
            print(f"Thrust: {thrust}, Torque: {torque}")
            print(f"Power: {Power}")
            #Save for plotting later
            tip_speed_list.append(s)
            theta_list.append(theta)
            # power_list.append(Power)
            thrust_list.append(thrust)
            torque_list.append(torque)
            




