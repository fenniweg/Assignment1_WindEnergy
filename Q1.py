
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
tip_speed_ratio= np.linspace(5.0, 10.0, 25)      # Tip speed ratio range 
theta_p = np.linspace(-4.0, 3.0, 25)       # Pitch angle range [deg] 


Cp = np.zeros((len(tip_speed_ratio), len(theta_p), 2)) #Initialize array to store Cp values for each method


#Outer loop over methods, inner loop over tip speed ratios and pitch angles
#Loop over all combinations of tip speed ratios and pitch angles to compute thrust, torque, and power using BEM_algorithm
for method in ['Polynomial','Madsen']:
    print(f"Running BEM_algorithm for method: {method}")

    for s in tip_speed_ratio:
        for theta in theta_p:
            Cp_value = BEM_algorithm(s,theta,method) #Call BEM_algorithm function to compute Cp for combinaion 
            Cp[np.where(tip_speed_ratio==s)[0][0], np.where(theta_p==theta)[0][0], np.where(np.array(['Polynomial','Madsen'])==method)[0][0]] = Cp_value #Store Cp value in array

    print(f"Completed BEM_algorithm for method: {method}")

#extract optimum
# Find the maximum Cp value and its corresponding tip speed ratio and pitch angle for each method
optimum_index_polynomial = np.unravel_index(np.argmax(Cp[:,:,0]), Cp[:,:,0].shape)
optimum_index_madsen = np.unravel_index(np.argmax(Cp[:,:,1]), Cp[:,:,1].shape)

cp_max_polynomial = Cp[optimum_index_polynomial[0], optimum_index_polynomial[1], 0]
cp_max_madsen = Cp[optimum_index_madsen[0], optimum_index_madsen[1], 1]

optimum_lambda_polynomial = tip_speed_ratio[optimum_index_polynomial[0]]
optimum_theta_polynomial = theta_p[optimum_index_polynomial[1]]

optimum_lambda_madsen = tip_speed_ratio[optimum_index_madsen[0]]
optimum_theta_madsen = theta_p[optimum_index_madsen[1]]
print("Search for maximum Cp values completed.")    
print("\n=================== RESULTS ===================")
#results
print(f"Polynomial Method: Cp_max = {cp_max_polynomial:.4f}, λ_max = {optimum_lambda_polynomial:.2f}, θp_max = {optimum_theta_polynomial:.2f}")
print(f"Madsen Method: Cp_max = {cp_max_madsen:.4f}, λ_max = {optimum_lambda_madsen:.2f}, θp_max = {optimum_theta_madsen:.2f}") 


print("\n=================== PLOTTING ===================")
#plotting
#contour plot 
fig, ax = plt.subplots(1, 2, figsize=(15, 10))


            





