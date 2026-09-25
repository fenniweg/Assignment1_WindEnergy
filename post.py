''' This file is used for post pprocessing of the results from main. It generates plots needed for the report and saves them in the Figures folder. It also saves the results in a text file for further analysis.
The plots generated are:
1. Cp as a function of tip speed ratio and pitch angle for both Polynomial and Madsen methods.
2. Power, pitch angle, thrust and Cp/Ct as a function of wind speed for feather and stall pitching.'''



import matplotlib.pyplot as plt
import numpy as np

from load_data import P_rated

#load npz files with results from main.py
npzfile = np.load('results/Cp_plotting_data.npz')
theta_p = npzfile['theta_p']
tip_speed_ratio = npzfile['tip_speed_ratio']
THETA_GRID = npzfile['THETA_GRID']
LAMBDA_GRID = npzfile['LAMBDA_GRID']
Cp_polynomial = npzfile['Cp_polynomial']
Cp_madsen = npzfile['Cp_madsen']
cp_max_polynomial = npzfile['cp_max_polynomial']
cp_max_madsen = npzfile['cp_max_madsen']
optimum_lambda_polynomial = npzfile['optimum_lambda_polynomial']
optimum_theta_polynomial = npzfile['optimum_theta_polynomial']
optimum_lambda_madsen = npzfile['optimum_lambda_madsen']
optimum_theta_madsen = npzfile['optimum_theta_madsen']

npzfile = np.load('results/sweep_plotting_data.npz')
v_sweep1 = npzfile['v_sweep']
omega_sweep = npzfile['omega_sweep']
P_sweep = npzfile['P_sweep']
omega_max = npzfile['omega_max']

npzfile = np.load('results/Q3_pitch_control_plotting_data.npz')
v_sweep = npzfile['v_sweep']
theta_p_sweep_feather = npzfile['theta_p_sweep_feather']
theta_p_sweep_stall = npzfile['theta_p_sweep_stall']
Cp_sweep_feather = npzfile['Cp_sweep_feather']
Ct_sweep_feather = npzfile['Ct_sweep_feather']
Cp_sweep_stall = npzfile['Cp_sweep_stall']
Ct_sweep_stall = npzfile['Ct_sweep_stall']
P_sweep_feather = npzfile['P_sweep_feather']
P_sweep_stall = npzfile['P_sweep_stall']
T_sweep_feather = npzfile['T_sweep_feather']
T_sweep_stall = npzfile['T_sweep_stall']
v_rated = npzfile['v_rated']

print("Imported plotting data from 'results' folder for Q1, Q2 and Q3.")

'Q1: Contour plot of Cp as a function of tip speed ratio and pitch angle for both Polynomial and Madsen methods'
fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, cp_data, title, theta_opt, lambda_opt in [
    (axs[0], Cp_polynomial, 'Polynomial Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_polynomial, optimum_lambda_polynomial),
    (axs[1], Cp_madsen, 'Madsen Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_madsen, optimum_lambda_madsen),
]:
    contour = ax.contourf(THETA_GRID, LAMBDA_GRID, cp_data, levels=20, cmap='viridis', rasterized=True)
    fig.colorbar(contour, ax=ax, label='$C_p$')
    ax.scatter(theta_opt, lambda_opt, color='r', s=80, marker='*', label='Max $C_p$')
    ax.set_title(title)
    ax.set_xlabel('Pitch Angle $\\theta_p$ [deg]')
    ax.set_ylabel('Tip Speed Ratio $\\lambda [-]$')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
fig.savefig('Figures/Cp_contour_comparison_raster.pdf', dpi=300)

'Q2: Plot omega and P vs wind speed up to max speed'
# Plot 1, omega vs wind speed from cut-in to max speed
plt.figure(figsize=(10, 6))
plt.plot(v_sweep1, omega_sweep, 'b-', linewidth=2.8, label='Rotational Speed $\\omega(V_0)$')
# Add a vertical dashed line to highlight the rated wind speed point
plt.axvline(v_rated, color='red', linestyle='--', linewidth=2.0,
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')

plt.axhline(omega_max, color='green', linestyle='--', linewidth=2.0,
            label=f'Maximum Rotational Speed = {omega_max:.2f} rad')

# Format the plot with titles, labels, and grid
plt.title('DTU 10MW: Rotational Speed vs Wind Speed', fontsize=12)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Rotational Speed [rad/s]', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7, linewidth=1.0)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig('Figures/omega_vs_wind_speed.pdf', dpi=300, format='pdf', bbox_inches='tight')


# Plot power against wind speed from cut-in to max speed
plt.figure(figsize=(10, 6))
plt.plot(v_sweep1, P_sweep, 'g-', linewidth=2.8, label='Power $P(V_0)$')
plt.axvline(v_rated, color='red', linestyle='--', linewidth=2.0,
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')
plt.axhline(P_rated, color='blue', linestyle='--', linewidth=2.0,
            label=f'Rated Power = {P_rated:.2f} W')
plt.grid(True, linestyle=':', alpha=0.7, linewidth=1.0)
plt.title('DTU 10MW: Power vs Wind Speed', fontsize=12)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Power [W]', fontsize=12)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig('Figures/power_vs_wind_speed.pdf', dpi=300, format='pdf', bbox_inches='tight')

print('Plots for Q2 saved in Figures folder')


# ---------------------------------------------------------------------------
# Plot the required quantities.
# ---------------------------------------------------------------------------

# 4 plots: 1) P vs V0, 2) theta_p vs V0, 3) T vs V0, 4) Cp and Ct vs V0 with each feather and stall pitch
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
for ax in axs.flat:
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.3, linewidth=1.0)

axs[0, 0].plot(v_sweep, P_sweep_feather, 'b-', lw=2.8, label='Feather')
axs[0, 0].plot(v_sweep, P_sweep_stall, 'r--', lw=2.8, label='Stall')
axs[0, 0].axhline(P_rated, color='k', ls='--', lw=1.8, label='Rated power')
axs[0, 0].set_xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
axs[0, 0].set_ylabel('Mechanical power $P$ [W]', fontsize=12)
axs[0, 0].set_title('Power regulation by pitching', fontsize=12)

axs[0, 1].plot(v_sweep, theta_p_sweep_feather, 'b-', lw=2.8, label='Feather')
axs[0, 1].plot(v_sweep, theta_p_sweep_stall, 'r--', lw=2.8, label='Stall')
axs[0, 1].set_xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
axs[0, 1].set_ylabel(r'Pitch angle $\theta_p$ [deg]', fontsize=12)
axs[0, 1].set_title(r'Pitch angle required for rated power at $\omega_{max}$', fontsize=12)

axs[1, 0].plot(v_sweep, T_sweep_feather, 'b-', lw=2.8, label='Feather')
axs[1, 0].plot(v_sweep, T_sweep_stall, 'r--', lw=2.8, label='Stall')
axs[1, 0].set_xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
axs[1, 0].set_ylabel('Thrust $T$ [N]', fontsize=12)
axs[1, 0].set_title('Thrust at rated-speed power limit', fontsize=12)

axs[1, 1].plot(v_sweep, Cp_sweep_feather, 'b-', lw=2.8, label=r'$C_p$ feather')
axs[1, 1].plot(v_sweep, Cp_sweep_stall, 'r--', lw=2.8, label=r'$C_p$ stall')
axs[1, 1].plot(v_sweep, Ct_sweep_feather, 'b:', lw=2.4, label=r'$C_T$ feather')
axs[1, 1].plot(v_sweep, Ct_sweep_stall, 'r:', lw=2.4, label=r'$C_T$ stall')
axs[1, 1].set_xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
axs[1, 1].set_ylabel('Coefficient value [-]', fontsize=12)
axs[1, 1].set_title('Dimensionless coefficients for each pitch strategy', fontsize=12)
for ax in axs.flat:
    ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig('Figures/Q3_pitch_control.pdf', dpi=300, format='pdf', bbox_inches='tight')
print("Plots for Q3 saved in Figures folder")