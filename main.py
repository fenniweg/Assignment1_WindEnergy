'''Calculate DTU Wind Energy 10MW wind turbine performance using BEM theory.
This code calculates the performance of the DTU Wind Energy 10MW wind turbine using Blade Element Momentum (BEM) theory. 
The code computes the power coefficient (Cp) as a function of tip speed ratio (lambda) and pitch angle (theta_p) using two methods: Polynomial and Interpolation. 
The rated wind speed and maximum rotational speed at rated wind speed are also calculated.
Stall ad feather pitching are implemented to limit the mechanical power to P_rated and the rotational speed to omega_max.
The results are plotted to visualize the performance characteristics of the wind turbine.'''

### IMPORT LIBRARIES AND DATA ###

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import brentq
from bem import BEM_algorithm, solve_pitch
from load_data import v_min,v_max,P_rated,R,rho,tip_speed_ratio, theta_p

####QUESTION 1: Compare the results of the two methods by plotting the power coefficient, Cp, as a function of the tip speed ratio, λ and pitch angles, θp.
#Initialize arrays
Cp = np.zeros((len(tip_speed_ratio), len(theta_p), 2)) #Initialize array to store Cp values for each method

print("Starting BEM_algorithm computations for all combinations of tip speed ratios and pitch angles...")

#Outer loop over methods, inner loop over tip speed ratios and pitch angles
#Loop over all combinations of tip speed ratios and pitch angles to compute thrust, torque, and power using BEM_algorithm
for method in ['Polynomial','Madsen']:
    print(f"Running BEM_algorithm for method: {method}")

    for s in tip_speed_ratio:
        for theta in theta_p:
            Cp_value = BEM_algorithm(s,theta,method) #Call BEM_algorithm function to compute Cp for combinaion 
            Cp[np.where(tip_speed_ratio==s)[0][0], np.where(theta_p==theta)[0][0], 0 if method=='Polynomial' else 1] = Cp_value #Store Cp value in array for correct method
    print(f"Completed BEM_algorithm for method: {method}")

df = pd.DataFrame(Cp[:,:,0], index=tip_speed_ratio, columns=theta_p) #Create DataFrame for Polynomial method
df2 = pd.DataFrame(Cp[:,:,1], index=tip_speed_ratio, columns=theta_p) #Create DataFrame for Madsen method


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

#Save results to load_data.py to use in Q2.py
with open('load_data.py', 'r') as file:
    lines = file.readlines()
    V_0 = 11.19
    optimum_theta = optimum_theta_polynomial


# Contour plot of Cp as a function of tip speed ratio and pitch angle for both methods
THETA_GRID, LAMBDA_GRID = np.meshgrid(theta_p, tip_speed_ratio)

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, cp_data, title, theta_opt, lambda_opt in [
    (axs[0], Cp[:, :, 0], 'Polynomial Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_polynomial, optimum_lambda_polynomial),
    (axs[1], Cp[:, :, 1], 'Madsen Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_madsen, optimum_lambda_madsen),
]:
    contour = ax.contourf(THETA_GRID, LAMBDA_GRID, cp_data, levels=20, cmap='viridis')
    fig.colorbar(contour, ax=ax, label='$C_p$')
    ax.scatter(theta_opt, lambda_opt, color='r', s=80, marker='*', label='Max $C_p$')
    ax.set_title(title)
    ax.set_xlabel('Pitch Angle $\\theta_p$ [deg]')
    ax.set_ylabel('Tip Speed Ratio $\\lambda$')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.savefig('Figures/Cp_contour_comparison.png', dpi=300)
print("Contour plots saved as 'Figures/Cp_contour_comparison.png'.")

print("Compute rated wind speed and maximum rotational speed at rated wind speed")
C_p_max = cp_max_polynomial #from Q1 results
optimum_lambda = optimum_lambda_polynomial #from Q1 results
optimum_theta = optimum_theta_polynomial #from Q1 results


#Calculate rotor area
A = np.pi*R**2

#Calculate rated wind speed
v_rated = (P_rated/(0.5*rho*A*C_p_max))**(1/3)

#Calculate maximum rotational speed at rated wind speed
omega_max = (optimum_lambda*v_rated)/R #[rad/s]
rpm_max = omega_max*60/(2*np.pi) #[rpm]

#Output results
print(f"Rated wind speed: {v_rated:.2f} m/s")
print(f"Maximum rotational speed at rated wind speed: {omega_max:.2f} rad/s ({rpm_max:.2f} rpm)")

print('Compute omega as a function of wind speed from cut-in to max speed')
print('Compute power as a function of wind speed from cut-in to max speed ')

v_sweep = np.linspace(v_min, v_rated, 100) #Sweep wind speed from cut-in to rated
omega_sweep = (optimum_lambda*v_sweep)/R #Calculate omega for each wind speed
rpm_sweep = omega_sweep*60/(2*np.pi) #Convert omega to rpm
P_sweep = 0.5*rho*A*C_p_max*v_sweep**3 #Calculate power for each wind speed
v_sweep = np.append(v_sweep, np.linspace(v_rated, v_max, 100)) #Sweep wind speed from rated to max
rpm_sweep = np.append(rpm_sweep, np.full(100, rpm_max)) #Omega is constant at omega_max above rated wind speed
P_sweep = np.append(P_sweep, np.full(100, P_rated)) #Power is constant at rated power above rated wind speed

#Plot 1, omega vs wind speed from cut-in to max speed
plt.figure(figsize=(10,6))
plt.plot(v_sweep, rpm_sweep, 'b-', linewidth=2.5, label='Rotational Speed $\\omega(V_0)$')
# Add a vertical dashed line to highlight the rated wind speed point
plt.axvline(v_rated, color='red', linestyle='--', 
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')

plt.axhline(rpm_max, color='green', linestyle='--',
            label=f'Maximum Rotational Speed = {rpm_max:.2f} rpm')

# Format the plot with titles, labels, and grid
plt.title('DTU 10MW: Rotational Speed vs Wind Speed', fontsize=14)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Rotational Speed [RPM]', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig('Figures/omega_vs_wind_speed.png', dpi=300)


#Plot power against wind speed from cut-in to max speed
plt.figure(figsize=(10,6))
plt.plot(v_sweep, P_sweep, 'g-', linewidth=2.5, label='Power $P(V_0)$')
plt.axvline(v_rated, color='red', linestyle='--', 
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')
plt.axhline(P_rated, color='blue', linestyle='--',
            label=f'Rated Power = {P_rated:.2f} W')
plt.grid(True, linestyle=':', alpha=0.7)
plt.title('DTU 10MW: Power vs Wind Speed', fontsize=14)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Power [W]', fontsize=12)
plt.legend(fontsize=11)
plt.savefig('Figures/power_vs_wind_speed.png', dpi=300)

print('Plots for Q2 saved in Figures folder')

### QUESTION 3: Implement stall and feather pitching to limit the mechanical power to P_rated and the rotational speed to omega_max.

# ---------------------------------------------------------------------------
# Find the pitch setting that makes the power equal to the rated power.
# For feather pitching, a larger pitch angle reduces Cp.
# For stall pitching, a smaller pitch angle reduces Cp.
# --------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Sweep between cut-in and cut-out wind speed, leave out rated wind speed since it is already known.
# ---------------------------------------------------------------------------
v_sweep = np.linspace(v_rated+0.1, v_max, 10)
theta_p_sweep_feather = np.zeros_like(v_sweep)
theta_p_sweep_stall = np.zeros_like(v_sweep)
Cp_sweep_feather = np.zeros_like(v_sweep)
Ct_sweep_feather = np.zeros_like(v_sweep)
Cp_sweep_stall = np.zeros_like(v_sweep)
Ct_sweep_stall = np.zeros_like(v_sweep)
P_sweep_feather = np.zeros_like(v_sweep)
P_sweep_stall = np.zeros_like(v_sweep)
T_sweep_feather = np.zeros_like(v_sweep)
T_sweep_stall = np.zeros_like(v_sweep)

for i, v in enumerate(v_sweep):
    lambda_i = omega_max * R / v
    cp_target = P_rated / (0.5 * rho * A * v**3)

    # Feather: increase pitch angle above the optimal setting to reduce Cp.
    # Stall: decrease pitch angle below the optimal setting to reduce Cp.
    theta_f = solve_pitch(0.0, 40.0, lambda_i, cp_target,optimum_theta = optimum_theta)
    theta_s = solve_pitch(-40.0, 0.0, lambda_i, cp_target,optimum_theta = optimum_theta)
    theta_p_sweep_feather[i] = theta_f
    theta_p_sweep_stall[i] = theta_s

    #compute final Cp, Ct, T and P for feather and stall pitching
    Cp_sweep_feather[i], Ct_sweep_feather[i],T_sweep_feather[i],P_sweep_feather[i] = BEM_algorithm(lambda_i, theta_f,Loads = True)
    Cp_sweep_stall[i], Ct_sweep_stall[i],T_sweep_stall[i],P_sweep_stall[i] = BEM_algorithm(lambda_i, theta_s,Loads = True)

# ---------------------------------------------------------------------------
# Plot the required quantities.
# ---------------------------------------------------------------------------

#4 plots: 1) P vs V0, 2) theta_p vs V0, 3) T vs V0, 4) Cp and Ct vs V0 with each feather and stall pitch

fig,axs = plt.subplots(2,2,figsize=(12,10))
axs[0,0].plot(v_sweep, P_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[0,0].plot(v_sweep, P_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[0,0].axhline(P_rated, color='k', ls='--', lw=1.5, label='Rated power')
axs[0,0].set_xlabel('Wind speed $V_0$ [m/s]')
axs[0,0].set_ylabel('Mechanical power $P$ [W]')
axs[0,0].set_title('Power regulation by pitching')
axs[0,0].grid(True, alpha=0.3)


axs[0,1].plot(v_sweep, theta_p_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[0,1].plot(v_sweep, theta_p_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[0,1].set_xlabel('Wind speed $V_0$ [m/s]')
axs[0,1].set_ylabel(r'Pitch angle $\theta_p$ [deg]')
axs[0,1].set_title(r'Pitch angle required for rated power at $\omega_{max}$')
axs[0,1].grid(True, alpha=0.3)

axs[1,0].plot(v_sweep, T_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[1,0].plot(v_sweep, T_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[1,0].set_xlabel('Wind speed $V_0$ [m/s]')
axs[1,0].set_ylabel('Thrust $T$ [N]')
axs[1,0].set_title('Thrust at rated-speed power limit')
axs[1,0].grid(True, alpha=0.3)

axs[1,1].plot(v_sweep, Cp_sweep_feather, 'b-', lw=2.5, label=r'$C_p$ feather')
axs[1,1].plot(v_sweep, Cp_sweep_stall, 'r--', lw=2.5, label=r'$C_p$ stall')
axs[1,1].plot(v_sweep, Ct_sweep_feather, 'b:', lw=2.0, label=r'$C_T$ feather')
axs[1,1].plot(v_sweep, Ct_sweep_stall, 'r:', lw=2.0, label=r'$C_T$ stall')
axs[1,1].set_xlabel('Wind speed $V_0$ [m/s]')
axs[1,1].set_ylabel('Coefficient value')    
axs[1,1].set_title('Dimensionless coefficients for each pitch strategy')
axs[1,1].grid(True, alpha=0.3)
for ax in axs.flat:
    ax.legend()

plt.tight_layout()

plt.savefig('Figures/Q3_pitch_control.png', dpi=300)