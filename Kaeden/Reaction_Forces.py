import math
import numpy as np
import sympy as sym


# Reaction Forces
def reaction_forces():

    ###########################################################################
    # Variables
    ###########################################################################
    pi = math.pi
    H = 150                 # hp
    n_input = 6700          # rpm
    d_input = 5.629165      # in
    d_output = 1.8763833    # in

    # Eq. 13-19: Calculate transverse angle using helix & normal angles
    phi_t = math.degrees(math.atan(math.tan(20*pi/180) / math.cos(30*pi/180)))
    ###########################################################################



    ###########################################################################
    # Calculate W force (lbf)
    ###########################################################################
    # Eq. 13-35: Calculate transmitted load
    W_t = (33000 * H * 12)/(math.pi * d_input * n_input)

    # Eq. 13-40: Calculate radial and axial loads
    W_r = W_t * math.tan(phi_t*pi/180)
    W_a = W_t * math.tan(30*pi/180)

    # Define W force vector
    W_input = np.array([W_a, W_r, -W_t])
    W_output = np.array([-W_a, -W_r, W_t])
    ###########################################################################



    ###########################################################################
    # Define vectors
    ###########################################################################
    # Create symbolic variables
    F_Ax, F_Ay, F_Az = sym.symbols('F_Ax F_Ay F_Az')
    F_Bx, F_By, F_Bz = sym.symbols('F_Bx F_By F_Bz')
    F_Cx, F_Cy, F_Cz = sym.symbols('F_Cx F_Cy F_Cz')
    F_Dx, F_Dy, F_Dz = sym.symbols('F_Dx F_Dy F_Dz')
    T = sym.Symbol('T')

    # Force vectors (lbf)
    F_A = np.array([0, F_Ay, F_Az])
    F_B = np.array([F_Bx, F_By, F_Bz])
    F_C = np.array([F_Cx, F_Cy, F_Cz])
    F_D = np.array([0, F_Dy, F_Dz])

    # Position Vectors (in)
    R_AB = np.array([3.75, 0, 0])
    R_AW = np.array([1.875, -d_input/2, 0])
    R_DC = np.array([-3.75, 0, 0])
    R_DW = np.array([-1.875, d_output/2, 0])

    # Define i, j, and k to access elements in vectors
    i = 0
    j = 1
    k = 2
    ###########################################################################



    ###########################################################################
    # INPUT 
    # Sum of moments [(R_AW x W) + (R_AB x F_B) + T = 0]
    ###########################################################################
    # Moment at A from W
    M_W = np.cross(R_AW, W_input)
    #print(M_AW)

    # Moment at A from F_B
    M_FB = np.cross(R_AB, F_B)
    #print(M_FB)

    # Set i components equal to 0 and solve
    sum_M_i = T + M_W[i] + M_FB[i]
    T_sol = sym.solve(sum_M_i, T)
    T_input = T_sol[0]
    #print(T_input)

    # Set j components equal to 0 and solve
    sum_M_j = M_W[j] + M_FB[j]
    F_Bz_sol = sym.solve(sum_M_j, F_Bz)
    F_B[k] = F_Bz_sol[0]
    #print(F_Bz_sol)

    # Set k components equal to 0 and solve
    sum_M_k = M_W[k] + M_FB[k]
    F_By_sol = sym.solve(sum_M_k, F_By)
    F_B[j] = F_By_sol[0]
    #print(F_By_sol)
    ###########################################################################



    ###########################################################################
    # INPUT
    # Sum of forces (F_A + F_B + W = 0)
    ###########################################################################
    # Set i components equal to 0 and solve
    sum_F_i = F_A[i] + F_B[i] + W_input[i]
    F_Bx_sol = sym.solve(sum_F_i, F_Bx)
    F_B[i] = F_Bx_sol[0]
    #print(F_Bx_sol)

    # Set j components equal to 0 and solve
    sum_F_j = F_A[j] + F_B[j] + W_input[j]
    F_Ay_sol = sym.solve(sum_F_j, F_Ay)
    F_A[j] = F_Ay_sol[0]
    #print(F_Ay_sol)

    # Set k components equal to 0 and solve
    sum_F_k = F_A[k] + F_B[k] + W_input[k]
    F_Az_sol = sym.solve(sum_F_k, F_Az)
    F_A[k] = F_Az_sol[0]
    #print(F_Az_sol)
    ###########################################################################



    ###########################################################################
    # OUTPUT
    # Sum of moments [(R_DW x W) + (R_DC x F_C) + T = 0]
    ###########################################################################
    # Moment at D from W
    M_W_output = np.cross(R_DW, W_output)
    #print(M_W_output)

    # Moment at D from F_C
    M_FC = np.cross(R_DC, F_C)
    #print(M_FC)

    # Set i components equal to 0 and solve
    sum_M_out_i = T + M_W_output[i] + M_FC[i]
    T_sol = sym.solve(sum_M_out_i, T)
    T_output = T_sol[0]
    #print(T_output)
    
    # Set j components equal to 0 and solve
    sum_M_out_j = M_W_output[j] + M_FC[j]
    F_Cz_sol = sym.solve(sum_M_out_j, F_Cz)
    F_C[k] = F_Cz_sol[0]
    #print(F_Cz_sol)
    
    # Set k components equal to 0 and solve
    sum_M_out_k = M_W_output[k] + M_FC[k]
    F_Cy_sol = sym.solve(sum_M_out_k, F_Cy)
    F_C[j] = F_Cy_sol[0]
    #print(F_Cy_sol)

    ###########################################################################



    ###########################################################################
    # OUTPUT
    # Sum of forces (F_C + F_D + W = 0)
    ###########################################################################
    # Set i components equal to 0 and solve
    sum_F_out_i = F_C[i] + F_D[i] + W_output[i]
    F_Cx_sol = sym.solve(sum_F_out_i, F_Cx)
    F_C[i] = F_Cx_sol[0]
    #print(F_Cx_sol)
    
    # Set j components equal to 0 and solve
    sum_F_out_j = F_C[j] + F_D[j] + W_output[j]
    F_Dy_sol = sym.solve(sum_F_out_j, F_Dy)
    F_D[j] = F_Dy_sol[0]
    #print(F_Dy_sol)
    
    # Set k components equal to 0 and solve
    sum_F_out_k = F_C[k] + F_D[k] + W_output[k]
    F_Dz_sol = sym.solve(sum_F_out_k, F_Dz)
    F_D[k] = F_Dz_sol[0]
    #print(F_Dz_sol)
    ###########################################################################

    
    ###########################################################################
    # Print Forces
    ###########################################################################
    print("\nInput Shaft:")
    print("T = %.4f" % T_input)
    print("W   = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(W_input[0], 4), 
                                                      round(W_input[1], 4), 
                                                      round(W_input[2], 4)))
    print("F_A = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(F_A[0], 4), 
                                                      round(F_A[1], 4), 
                                                      round(F_A[2], 4)))
    print("F_B = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(F_B[0], 4), 
                                                      round(F_B[1], 4), 
                                                      round(F_B[2], 4)))
    print("\nOutput Shaft:")
    print("T = %.4f" % T_output)
    print("W   = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(W_output[0], 4), 
                                                      round(W_output[1], 4), 
                                                      round(W_output[2], 4)))
    print("F_C = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(F_C[0], 4), 
                                                      round(F_C[1], 4), 
                                                      round(F_C[2], 4)))
    print("F_D = [{:<15.4f} {:<15.4f} {:.4f}]".format(round(F_D[0], 4), 
                                                      round(F_D[1], 4), 
                                                      round(F_D[2], 4)))
    print("\n\n\n")
    ###########################################################################


    # Sanity Check
    R_CD = np.array([3.75, 0, 0])
    R_CW = np.array([1.875, d_output/2, 0])
    test1 = np.cross(R_CW, W_output) + np.cross(R_CD, F_D)
    #print(test1)
    test2 = np.cross(R_DW, W_output) + np.cross(R_DC, F_C)
    #print(test2)
    
    
    return
###############################################################################

reaction_forces()




