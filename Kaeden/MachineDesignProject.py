import math


###############################################################################
# Reaction forces at bearings
###############################################################################
def reaction_forces(d_in, d_out):

 
 # INPUT SHAFT 
    W_a_in = 1629.3/d_in
    W_r_in = 1186.03/d_in
    W_t_in = -2822.03/d_in
    
    print("\nGear Forces (Input):")
    print("W_a = %.4f lbf" % W_a_in)
    print("W_r = %.4f lbf" % W_r_in)
    print("W_t = %.4f lbf" % W_t_in)

    # Reaction forces at A
    F_Ay = (-593.015/d_in) + 217.24
    F_Az = 1411.015/d_in

    # Reaction forces at B
    F_Bx = -1629.3/d_in
    F_By = (-593.015/d_in) - 217.24
    F_Bz = 1411.015/d_in

    print("\nInput Shaft:")
    print("F_Ay = %.4f lbf" % F_Ay)
    print("F_Az = %.4f lbf" % F_Az)
    print("F_Bx = %.4f lbf" % F_Bx)
    print("F_By = %.4f lbf" % F_By)
    print("F_Bz = %.4f lbf" % F_Bz)

 # OUTPUT SHAFT
    W_a_out = -1629.3/d_in
    W_r_out = -1186.03/d_in
    W_t_out = 2822.03/d_in
    
    print("\nGear Forces (Output):")
    print("W_a = %.4f lbf" % W_a_out)
    print("W_r = %.4f lbf" % W_r_out)
    print("W_t = %.4f lbf" % W_t_out)

    # Reaction forces at C
    F_Cx = 1629.3/d_in
    F_Cy = (593.015/d_in) + 217.24
    F_Cz = -1411.015/d_in

    # Reaction forces at D
    F_Dy = (593.015/d_in) - 217.24
    F_Dz = -1411.015/d_in

    print("\nOutput Shaft:")
    print("F_Cx = %.4f lbf" % F_Cx)
    print("F_Cy = %.4f lbf" % F_Cy)
    print("F_Cz = %.4f lbf" % F_Cz)
    print("F_Dy = %.4f lbf" % F_Dy)
    print("F_Dz = %.4f lbf\n" % F_Dz)

    return
###############################################################################



d_in = 9.526279442
d_out = 3.175426481

reaction_forces(d_in, d_out)