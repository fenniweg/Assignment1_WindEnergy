#BEM algorithm steady

import numpy as np
from load_data import blade_dat, airfoil_data, R,n_blades,rho,v_min,v_max,V_0


def double_interpolation(alpha,t_over_c):
    C_l_thickness = np.zeros((6))
    C_d_thickness = np.zeros((6))
    
    for k in range (6):
        C_l_thickness[k] = np.interp(alpha,airfoil_data[k][:,0],airfoil_data[k][:,1])
        C_d_thickness[k] = np.interp(alpha,airfoil_data[k][:,0],airfoil_data[k][:,2])

    #Interpolate to actual thickness
    thickness_profile = [100,60,48,36,30.1,24.1] #thickness profile of the 6 airfoils in percentage
    #Use np.argsort to sort the thickness profile and corresponding C_l and C_d values
    sort_indices = np.argsort(thickness_profile)
    x_p_sorted = np.array(thickness_profile)[sort_indices]
    y_p_sorted_C_l = np.array(C_l_thickness)[sort_indices]
    y_p_sorted_C_d = np.array(C_d_thickness)[sort_indices]

    C_l = np.interp(t_over_c,x_p_sorted,y_p_sorted_C_l)
    C_d = np.interp(t_over_c,x_p_sorted,y_p_sorted_C_d)
    return C_l, C_d


def BEM_algorithm (s,theta_p,method = 'Polynomial', Loads = False):
    '''
    BEM_algorithm computes the power coefficient (Cp) and thrust coefficient (CT) for a given tip speed ratio (s), 
    pitch angle (theta_p), and method ('Polynomial' or 'Madsen'). Polynomial is standard.'''

    omega = s*V_0/R
    #Load blade data
    r_list = blade_dat['r'].values
    chord_list = blade_dat['c'].values
    beta_list = blade_dat['beta'].values
    t_over_c_list = blade_dat['t/c'].values

    p_n_list = np.zeros(len(r_list))
    p_t_list = np.zeros(len(r_list))

    #loop through each blade element
    for i in range(len(r_list)-1): #last element is tip, skip to avoid numerical issues
        r = r_list[i]
        chord = chord_list[i]
        beta = beta_list[i]
        t_over_c = t_over_c_list[i]

        sigma = chord*n_blades/(2*np.pi*r) 

    #Initialze a and a_prime, convergence tolerance

        a = 0
        a_prime = 0
        epsilon = 1e-6

        count = 0

        while True:
            count += 1
            #Calculate flowangle
            psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r)) #radians
    
            #Compute local angle of attack alpha
            alpha = np.degrees(psi)-(beta+theta_p) #degrees
    
            #Lookup C_l and C_d from airfoil data based on alpha with double interpolation
            C_l, C_d = double_interpolation(alpha,t_over_c)

            #Compute normal and tangential force coefficients
            C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
            C_t = C_l*np.sin(psi) - C_d*np.cos(psi)
    
            
            dC_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
    
            # Prandtl Tip Loss Correction F 
            if np.abs(np.sin(psi)) < 1e-6:
                F = 1.0
            else:
                exp_arg = -(n_blades / 2.0) * (R - r) / (r * np.abs(np.sin(psi)))
                F = (2.0 / np.pi) * np.arccos(np.exp(exp_arg))
                F = max(F, 1e-4)  # Avoid division by zero

            #Axial Induction Factor a and Tangential Induction Factor a_prime for chosen method
            if method  == 'Polynomial':
                if a <= 0.33: #no correction for a_star applied
                    a_star = sigma*C_n/(4*F*np.sin(psi)**2)*(1-a)
                else: # a_star with correction
                    a_star = dC_T / (4.0 * F * (1.0 - 0.25 * (5.0 - 3.0 * a) * a))
            elif method =='Madsen': #Madsen et al
                    CT_F = dC_T / F
                    a_star = 0.246 * CT_F + 0.0586 * (CT_F**2) + 0.0883 * (CT_F**3)

            else:
                raise ValueError("Invalid method. Choose 'Polynomial' or 'Madsen'.")

            #Update a and a_prime with relaxation factor 0.1
            a_new = 0.1*a_star+(1-0.1)*a #Update a with relaxation factor 0.1

            a_prime_star = sigma*C_t/(4*F*np.sin(psi)*np.cos(psi))*(1+a_prime)
            a_prime_new = 0.1*a_prime_star+(1-0.1)*a_prime

            #Check for convergence
            if abs(a_new - a) < epsilon and abs(a_prime_new - a_prime) < epsilon:
                # print(f"Converged after {i} iterations")
                a, a_prime = a_new, a_prime_new
                break
            else:
                a = a_new
                a_prime = a_prime_new

            if count > 1000:
                print(f"Warning: BEM did not converge after 1000 iterations in element {i}")
                #break and go to next element
        
                break
            


        #if converged, calculate p_n and p_t for this blade element
            

        V_rel = np.sqrt((V_0 * (1.0 - a))**2 + (omega * r * (1.0 + a_prime))**2)

        p_n = 0.5*rho*(V_rel)**2*chord*C_n
        p_t = 0.5*rho*(V_rel)**2*chord*C_t

        #Save p_n and p_t for each blade element
        p_n_list[i] = (p_n)
        p_t_list[i]=(p_t)
        

    #Force last element to be zero to avoid numerical issues at the tip
    p_n_list[-1] = 0.0
    p_t_list[-1] = 0.0
    # Integrate over the blade elements to get total thrust and torque
    thrust = n_blades*np.trapz(p_n_list,r_list)
    torque = n_blades*np.trapz(r_list*p_t_list,r_list)

    P = torque*omega

    #Dimensionless Coefficients
    A = np.pi*R**2
    Cp = P/(0.5*rho*A*V_0**3)
    CT = thrust/(0.5*rho*A*V_0**2)

    if Loads == True:
        return Cp,CT,thrust,P


    return Cp


'TESTS'   

# V_0 = 10
# theta_p = 0
# s = 5
# method = 'Madsen'
# #run test of BEM_algorithm with these values
# p_n, p_t= BEM_algorithm(s,theta_p,method)

# print(f"p_n: {p_n}, p_t: {p_t}")
# # R = 31
# # B = 3
# # rho = 1.225
# # V_0 = 8.0
# # omega = 2.61
# # theta_p = -3.0
# # beta = 2.0
# # chord = 1.7
# # C_l = 0.5
# # C_d= 0.01
# # r = 24.5

# # method = 'Polynomial' #third order polynomial


# # a, a_prime ,p_t, p_n, F = bem (r, R, B,rho,V_0,omega,theta_p,beta,chord,C_l,C_d,method)
# # print(f"Method {method}:")
# # print(f"a: {a}, a': {a_prime}, p_t: {p_t}, p_n: {p_n}, F: {F}")

# # method = 'Madsen'#Madsen et al


# # a, a_prime ,p_t, p_n, F = bem (r, R, B,rho,V_0,omega,theta_p,beta,chord,C_l,C_d,method)
# # print(f"Method {method}:")
# # print(f"a: {a}, a': {a_prime}, p_t: {p_t}, p_n: {p_n}, F: {F}")




