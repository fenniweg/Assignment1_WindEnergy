#BEM algorithm steady

import numpy as np
from load_data import blade_dat, airfoil_data, R,n_blades,rho,v_min,v_max,V_0,F


def bem (r, R, B,rho,V_0,omega,theta_p,beta,chord,C_l,C_d,method):

    a = 0
    a_prime = 0
    epsilon = 1e-7
    

    psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r))

    
    

    C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
    C_t = C_l*np.sin(psi) - C_d*np.cos(psi)

    sigma = chord*B/(2*np.pi*r) 
    C_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
    
    F = 2/np.pi*np.arccos(np.exp(-B/2*(R-r)/(r*np.sin(abs(psi)))))

    i = 0

    while True:

        if a < 0.33: #no correction for a_star applied
            a_star = sigma*C_n/(4*F*np.sin(psi)**2)*(1-a)
        else: # a_star with correction
            if method == 'Polynomial': #third order polynomial
                dC_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
                a_star = dC_T/(4*F*(1-0.25*(5-3*a)*a))
                
            else: #Madsen et al
                dC_T= (1-a)**2*C_n*sigma/(np.sin(psi)**2)
                a_star = 0.246*(dC_T/F)+0.0586*(dC_T/F)**2+0.083*(dC_T/F)**3

        a_new = 0.1*a_star+(1-0.1)*a #Update a with relaxation factor 0.1

        #Update a_prime with relaxation factor 0.1, no correction for a_prime_star applied
        a_prime_star = sigma*C_t/(4*F*np.sin(psi)*np.cos(psi))*(1+a_prime)
        a_prime_new = 0.1*a_prime_star+(1-0.1)*a_prime

        #Check for convergence
        if abs(a_new - a) < epsilon and abs(a_prime_new - a_prime) < epsilon:
            print(f"Converged after {i} iterations")
            break
        else:
            a = a_new
            a_prime = a_prime_new

        if i > 1000:
            print("Warning: BEM did not converge after 1000 iterations")
            break

        #Calculate psi, C_n, C_t, sigma, C_T, F for the next iteration or return
        psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r))
        C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
        C_t = C_l*np.sin(psi) - C_d*np.cos(psi)
        
        sigma = chord*B/(2*np.pi*r) 
        C_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
            
        F = 2/np.pi*np.arccos(np.exp(-B/2*(R-r)/(r*np.sin(abs(psi)))))

        V_rel = V_0*(1-a)/np.sin(psi)

        p_n = 0.5*rho*(V_rel)**2*chord*C_n
        p_t = 0.5*rho*(V_rel)**2*chord*C_t

    return a, a_prime, p_t,p_n,F


def BEM_algorithm (s,theta_p,method):
    omega = s*V_0/R

    #Load blade data
    r_list = blade_dat['r'].values
    chord_list = blade_dat['c'].values
    beta_list = blade_dat['beta'].values
    t_over_c_list = blade_dat['t/c'].values

    p_n_list = []
    p_t_list = []

    #loop through each blade element
    for i in range(len(r_list)):
        r = r_list[i]
        chord = chord_list[i]
        beta = beta_list[i]
        t_over_c = t_over_c_list[i]

    #Initialze a and a_prime, convergence tolerance

        a = 0
        a_prime = 0
        epsilon = 1e-7
    
    #Calculate flowangle
        psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r))

    #Compute local angle of attack alpha
        alpha = psi-(beta+theta_p)

        #Lookup C_l and C_d from airfoil data based on alpha with double interpolation
        #Interpolate first the values for each of the  6 different airfoils based on alpha
        C_l_thickness = np.zeros((6))
        C_d_thickness = np.zeros((6))
        
        for k in range (6):
            C_l_thickness[k] = np.interp(alpha,airfoil_data[k][:,1],airfoil_data[k][:,0])
            C_d_thickness[k] = np.interp(alpha,airfoil_data[k][:,2],airfoil_data[k][:,0])

        print(C_l_thickness)    

        #Interpolate to actual thickness
        thickness_profile = [100,60,48,36,30.1,24.1]
        C_l = np.interp(t_over_c, C_l_thickness,thickness_profile)
        C_d = np.interp(t_over_c, C_d_thickness,thickness_profile)
            

        C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
        C_t = C_l*np.sin(psi) - C_d*np.cos(psi)

        sigma = chord*n_blades/(2*np.pi*r) 
        C_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
        
        F = 2/np.pi*np.arccos(np.exp(-n_blades/2*(R-r)/(r*np.sin(abs(psi)))))

        i = 0

        while True:

            if a < 0.33: #no correction for a_star applied
                a_star = sigma*C_n/(4*F*np.sin(psi)**2)*(1-a)
            else: # a_star with correction
                if method == 'Polynomial': #third order polynomial
                    dC_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
                    a_star = dC_T/(4*F*(1-0.25*(5-3*a)*a))
                    
                if method =='Madsen': #Madsen et al
                    dC_T= (1-a)**2*C_n*sigma/(np.sin(psi)**2)
                    a_star = 0.246*(dC_T/F)+0.0586*(dC_T/F)**2+0.083*(dC_T/F)**3

            a_new = 0.1*a_star+(1-0.1)*a #Update a with relaxation factor 0.1

            #Update a_prime with relaxation factor 0.1, no correction for a_prime_star applied
            a_prime_star = sigma*C_t/(4*F*np.sin(psi)*np.cos(psi))*(1+a_prime)
            a_prime_new = 0.1*a_prime_star+(1-0.1)*a_prime

            #Check for convergence
            if abs(a_new - a) < epsilon and abs(a_prime_new - a_prime) < epsilon:
                print(f"Converged after {i} iterations")
                break
            else:
                a = a_new
                a_prime = a_prime_new

            if i > 1000:
                print("Warning: BEM did not converge after 1000 iterations")
                break

            #Calculate psi, C_n, C_t, sigma, C_T, F for the next iteration or return
            psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r))
            C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
            C_t = C_l*np.sin(psi) - C_d*np.cos(psi)
            
            sigma = chord*n_blades/(2*np.pi*r) 
            C_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
                
            F = 2/np.pi*np.arccos(np.exp(-n_blades/2*(R-r)/(r*np.sin(abs(psi)))))

            V_rel = V_0*(1-a)/np.sin(psi)

            p_n = 0.5*rho*(V_rel)**2*chord*C_n
            p_t = 0.5*rho*(V_rel)**2*chord*C_t

        #Save p_n and p_t for each blade element
        p_n_list.append(p_n)
        p_t_list.append(p_t)

        #Integrate over the blade elements to get total thrust and torque
        thrust = np.trapz(p_n_list,r)
        torque = np.trapz(p_t_list,r)

        P = torque*omega


    return thrust, torque, P


# R = 31
# B = 3
# rho = 1.225
# V_0 = 8.0
# omega = 2.61
# theta_p = -3.0
# beta = 2.0
# chord = 1.7
# C_l = 0.5
# C_d= 0.01
# r = 24.5

# method = 'Polynomial' #third order polynomial


# a, a_prime ,p_t, p_n, F = bem (r, R, B,rho,V_0,omega,theta_p,beta,chord,C_l,C_d,method)
# print(f"Method {method}:")
# print(f"a: {a}, a': {a_prime}, p_t: {p_t}, p_n: {p_n}, F: {F}")

# method = 'Madsen'#Madsen et al


# a, a_prime ,p_t, p_n, F = bem (r, R, B,rho,V_0,omega,theta_p,beta,chord,C_l,C_d,method)
# print(f"Method {method}:")
# print(f"a: {a}, a': {a_prime}, p_t: {p_t}, p_n: {p_n}, F: {F}")




