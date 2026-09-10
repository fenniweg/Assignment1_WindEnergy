# Assignment 1  
# Question 1
# Compute the highest obtainable power coefficient, Cp,max (λmax, θp,max), and the corresponding values of the tip speed ratio and pitch angle, λmax and θmax.
# Using BEM two methods: the polynomial method and the Madsen et al. method. Compare the results of the two methods by plotting the power coefficient, Cp, as a function of the tip speed ratio, λ, for a range of pitch angles, θp. 
#Glauert's method is used to compute the power coefficient, Cp, as a function of the tip speed ratio, λ, for a range of pitch angles, θp. 

from bem import BEM_algorithm
from glauert_rotor import c_p_lambda

#Data for 10 MW Turbine

R = 89.17
n_blades = 3
P_rated = 10e6
v_min = 4 #cut in wind speed
v_max = 25 #cut out wind speed
rho = 1.225



tip_speed_ratio = [5,6,7,8,9,10]
theta_p =[-3,-2,-1,0,1,2,3,4]



for s in speed:
    for t in theta_p:
        print(f"Tip speed ratio: {s}, Pitch angle: {t}")
        omega = s*v
        C_l, C_d =  (R, n_blades,rho,v_min,omega,theta_p,beta,chord,thickness,method):
        #Call BEM_algorithm to compute Cp for given speed and pitch angle
        #Call c_p_lambda to compute Cp for given speed and pitch angle