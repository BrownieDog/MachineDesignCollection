import math

# comprehensive list of all the AGMA coefficients
def AGMA_coefficients(W_t, Q_v, V, P_d, N_cycle_P, N_cycle_G, F, p_x, pt_angle, N_G, N_P, d_P, d_G, P_n, S):
    # W_t is tangential transmitted load (lbf)
    # Q_v is quality number of gears
    # V is inline pitch velocity (ft/min)
    # P_d is transverse diameteral pitch
    # N_cycle_P is number of cycles for the pinion
    # N_cycle_G is the number of cycles for the gear
    # F is face width of narrow member
    # p_x is axial pitch
    # pt_angle is transverse pressure angle
    # N_G is number of teeth of gear
    # N_P is number of teeth of pinion
    # d_P is pitch diameter of the pinion (in)
    # d_G is gear pitch diameter
    # P_n is normal diametrical pitch
    # S is distance between center of bearings

    K_o = overload_factor()                         # overload factor
    print(f"K_o = {K_o}")

    K_v = dynamic_factor(V, Q_v)                    # dynamic factor
    print(f"K_v = {K_v}, where V = {V} and Q_v = {Q_v}")

    K_s = size_factor()                             # size factor
    print(f"K_s = {K_s}")

    K_m = load_distribution_factor(d_P, F, S)       # load-distribution factor
    print(f"K_m = {K_m}, where d_P = {d_P}, F = {F}, and S = {S}")

    K_B = rim_thickness_factor()                    # rim-thickness factor
    print(f"K_B = {K_B}")

    J_P = bending_geometry_factor(p_x, F, N_P, N_G)      # J is geometry factor for bending stress including root fillet stress concentration factor - Fig 14-6 for the pinion
    print(f"J_P = {J_P}, where p_x = {p_x}, F = {F}, and N_P = {N_P}")

    J_G = bending_geometry_factor(p_x, F, N_G, N_P)      # J is geometry factor for bending stress including root fillet stress concentration factor - Fig 14-6 for the gear
    print(f"J_G = {J_G}, where p_x = {p_x}, F = {F}, and N_G = {N_G}")

    S_t = 45000     # guess using Table 14-3 flame or induction hardened            # bending strength (lbf/in^2)
    print(f"S_t = {S_t}, which is a guess")

    Y_N_P = bending_stress_cycle_factor(N_cycle_P)  # Y_N is bending stress cycle life factor for the pinion
    print(f"Y_N_P = {Y_N_P}, where N_cycle_P = {N_cycle_P}")

    Y_N_G = bending_stress_cycle_factor(N_cycle_G)  # Y_N is bending stress cycle life factor for the gear
    print(f"Y_N_G = {Y_N_G}, where N_cycle_G = {N_cycle_G}")

    K_T = temperature_factor()                      # temperature factor
    print(f"K_T = {K_T}")

    K_R = reliability_factor()                      # reliability factor
    print(f"K_R = {K_R}")

    C_p = elastic_coefficient()                     # elastic coefficient (sqrt(lbf/in^2))
    print(f"C_p = {C_p}")

    C_f = surface_condition_factor()                # surface condition factor
    print(f"C_f = {C_f}")

    I = contact_geometry_factor(pt_angle, N_G, N_P, d_G, d_P, P_n)  # I is contact geometry factor for pitting resistance
    print(f"I = {I}, where pt_angle = {pt_angle}, N_G = {N_G}, N_P = {N_P}, d_G = {d_G}, d_P = {d_P}, P_n = {P_n}")

    S_c = 170000    # guess using Table 14-6 for flame or induction hardened  # allowable contact stress (lbf/in^2)
    print(f"S_c = {S_c}, which is a guess")

    Z_N_P = contact_stress_cycle_factor(N_cycle_P)            # Z_N is wear/contact stress cycle life factor
    print(f"Z_N_P = {Z_N_P}, where N_cycle_P = {N_cycle_P}")

    Z_N_G = contact_stress_cycle_factor(N_cycle_G)            # Z_N is wear/contact stress cycle life factor
    print(f"Z_N_G = {Z_N_G}, where N_cycle_G = {N_cycle_G}")

    C_H_G = gear_hardness_ratio_factor(N_G, N_P, d_G, d_P)    # gear hardness ratio factors for pitting resistance
    print(f"C_H_G = {C_H_G}, where N_G = {N_G}, N_P = {N_P}, d_G = {d_G}, d_P = {d_P}")

    C_H_P = pinion_hardness_ratio_factor()          # pinion hardness ratio factors for pitting resistance
    print(f"C_H_P = {C_H_P}")

    # AGMA bending stress factor of safety, a stress ratio
    S_F_G = bending_safety_factor_AGMA(S_t, Y_N_G, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J_G)     # gear bending factor of safety
    print(f"S_F_G = {S_F_G}, where S_t = {S_t}, Y_N_G = {Y_N_G}, K_T = {K_T}, K_R = {K_R}, W_t = {W_t}, K_o = {K_o}, K_v = {K_v}, K_s = {K_s}, P_d = {P_d}, F = {F}, K_m = {K_m}, K_B = {K_B}, J_G = {J_G}")

    S_F_P = bending_safety_factor_AGMA(S_t, Y_N_P, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J_P)     # pinion bending factor of safety
    print(f"S_F_P = {S_F_P}, where S_t = {S_t}, Y_N_P = {Y_N_P}, K_T = {K_T}, K_R = {K_R}, W_t = {W_t}, K_o = {K_o}, K_v = {K_v}, K_s = {K_s}, P_d = {P_d}, F = {F}, K_m = {K_m}, K_B = {K_B}, J_P = {J_P}")

    # AGMA wear/contact factor of safety, a stress ratio
    S_H_G = contact_safety_factor_AGMA(S_c, Z_N_G, C_H_G, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I)
    print(f"S_H_G = {S_H_G}, where S_c = {S_c}, Z_N_G = {Z_N_G}, C_H_G = {C_H_G}, K_T = {K_T}, K_R = {K_R}, C_p = {C_p}, W_t = {W_t}, K_o = {K_o}, K_v = {K_v}, K_s = {K_s}, K_m = {K_m}, d_P = {d_P}, F = {F}, C_f = {C_f}, and I = {I}")

    S_H_P = contact_safety_factor_AGMA(S_c, Z_N_P, C_H_P, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I)
    print(f"S_H_G = {S_H_G}, where S_c = {S_c}, Z_N_G = {Z_N_P}, C_H_G = {C_H_P}, K_T = {K_T}, K_R = {K_R}, C_p = {C_p}, W_t = {W_t}, K_o = {K_o}, K_v = {K_v}, K_s = {K_s}, K_m = {K_m}, d_P = {d_P}, F = {F}, C_f = {C_f}, and I = {I}")

    return K_o, K_v, K_s, K_m, K_B, S_t, Y_N_P, Y_N_G, K_T, K_R, C_p, C_f, I, S_c, Z_N_P, Z_N_G, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P
    # SCOTT, I replaced Y_N with Y_N_P and then added Y_N_G after that.
    # The new Y_N_P and Y_N_G take into account that the gear and pinion have different number of cycles.
    # Consequently, I have added the variables into the AGMA coefficient arguments as N_cycle_P and N_cycle_G.
    # Before the change, the arguments are W_t, Q_v, V, P_d, d_P, N, F, p_x, pt_angle, N_G, N_P, d_G, P_n, S.
    # Now, they are W_t, Q_v, V, P_d, d_P, N_cycle_P, N_cycle_G, F, p_x, pt_angle, N_G, N_P, d_G, P_n, S.
    # The list defining the arguments has been updated.

    # Now, with the new cycle inputs, I have also updated the Z_N to accommodate the different number of cycles for the gear and pinion.
    # I have created Z_N_P and Z_N_G.
    # The return list has been updated from K_o, K_v, K_s, K_m, K_B, S_t, Y_N_P, Y_N_G, K_T, K_R, C_p, C_f, I, S_c, Z_N, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P
    # to K_o, K_v, K_s, K_m, K_B, S_t, Y_N_P, Y_N_G, K_T, K_R, C_p, C_f, I, S_c, Z_N_P, Z_N_G, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P
    # You may have to update your arguments in your main function.

# equation 14-15
def calc_bending_stress_AGMA(W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J):
    # W_t tangential transmitted load (lbf)
    # K_o is overload factor
    # K_v is dynamic factor
    # K_s is size factor
    # P_d is transverse diameteral pitch
    # F is face width of narrow member
    # K_m is load-distribution factor
    # K_B is rim-thickness factor
    # J is geometry factor for bending stress

    calc_bend_AGMA = W_t * K_o * K_v * K_s * (P_d / F) * (K_m * K_B / J)
    return calc_bend_AGMA

# equation 14-17
def calc_allowable_bending_stress_AGMA(S_t, S_F, Y_N, K_T, K_R):
    # S_t is bending strength (lbf/in^2)
    # S_F is AGMA bending factor of safety
    # Y_N is bending stress cycle life factor
    # K_T is temperature factor
    # K_R is reliability factor

    calc_allow_bend_AGMA = (S_t / S_F) * (Y_N / (K_T * K_R))
    return calc_allow_bend_AGMA

# equation 14-16
def calc_contact_stress_AGMA(C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I):
    # C_p is elastic coefficient (sqrt(lbf/in^2))
    # W_t tangential transmitted load (lbf)
    # K_o is overload factor
    # K_v is dynamic factor
    # K_s is size factor
    # K_m is load-distribution factor
    # d_P pitch diameter of the pinion (in)
    # F is face width of narrow member
    # C_f is surface condition factor
    # I is contact geometry factor
    calc_cont_stress_AGMA = C_p * math.sqrt(W_t * K_o * K_v * K_s * (K_m / (d_P * F)) * (C_f / I))
    return calc_cont_stress_AGMA

# equation 14-18
def calc_allowable_contact_stress_AGMA(S_c, Z_N, C_H, K_T, K_R, S_H):
    # S_c is allowable contact stress (lbf/in^2)
    # Z_N is wear/ contact stress cycle life factor
    # C_H is hardness ratio factors for pitting resistance
    # K_T is temperature factor
    # K_R is reliability factor
    # S_H is AGMA wear/contact factor of safety, a stress ratio

    calc_allow_cont_AGMA = (S_c / S_H) * ((Z_N * C_H) / K_T * K_R)
    return calc_allow_cont_AGMA

# equation 14-22
def speed_ratio(N_G, N_P, d_G, d_P):
    # N_G is number of teeth in the gear
    # N_P is number of teeth in the pinion
    # d_G is gear pitch diameter
    # d_P is pinion pitch diameter

    m_G = N_G / N_P
    # OR #
    m_G = d_G / d_P
    return m_G

# equation xxxxxxxxx
def addendum(P_n):
    # P_n is normal diametrical pitch

    a = 1 / P_n  # addendum
    return a

def contact_geometry_factor(pt_angle, N_G, N_P, d_G, d_P, P_n):
    # pt_angle transverse pressure angle
    # N_G is number of teeth in the gear
    # N_P is number of teeth in the pinion
    # d_G is gear pitch diameter
    # d_P is pinion pitch diameter
    # P_n is normal diametrical pitch

    a = addendum(P_n)   # calling the addendum function


    r_P = d_P / 2     # transverse pitch radius of pinion

    # equation 14-26
    r_bP = r_P * math.cos(math.radians(pt_angle))   # base-circle radii for pinion

    r_G = d_G / 2     # transverse pitch radius of gear
    r_bG = r_G * math.cos(math.radians(pt_angle))   # base-circle radii for gear

    # equation 14-25
    Z_1 = math.sqrt((r_P + a) ** 2 - r_bP ** 2)
    Z_2 = math.sqrt((r_G + a) ** 2 - r_bG ** 2)
    Z_3 = (r_P + r_G) * math.sin(math.radians(pt_angle))

    # Z condition statements.
    # if either of the first two term groups are larger than the last, replace it with the last term
    if Z_1 > Z_3:
        Z = Z_3 + Z_2 - Z_1
    elif Z_2 > Z_3:
        Z = Z_1 + Z_3 - Z_2
    else:
        Z = Z_1 + Z_2 - Z_3

    p_n = math.pi / P_n                 # normal circular pitch
    # round up to a standard

    # equation 14-24
    p_N = p_n * math.cos(math.radians(pt_angle))

    # equation 14-21
    m_N = p_N / (0.95 * Z)          # load sharing ratio

    # equation 14-22
    m_G = speed_ratio(N_G, N_P, d_G, d_P)

    # equation 14-23
    I_ext = (math.cos(math.radians(pt_angle)) * math.sin(math.radians(pt_angle))) / (2 * m_N) * (m_G / (m_G + 1))    # geometry factor for external gear

    return I_ext

def bending_geometry_factor(p_x, F, N_teeth_1, N_teeth_2):
    # p_x is axial pitch
    # F is narrow face width
    # N_teeth_1 is number of teeth of the gear being considered
    # N_teeth_2 is number of teeth of the other gear
    # if N_teeth_1 is the pinion, then N_teeth_2 i the gear

    #m_F = F / p_x       # unused now?

    # J modifiers to find J
    if N_teeth_1 <= 25:
        J_mod = 0.465           # Fig 14-7
    elif N_teeth_1 <= 45:
        J_mod = 0.5             # Fig 14-7
    elif N_teeth_1 <= 105:
        J_mod = 0.54            # Fig 14-7
    elif N_teeth_1 <= 325:
        J_mod = 0.565           # Fig 14-7
    else:
        J_mod = 0.58            # Fig 14-7

    # J factors to find J
    if N_teeth_2 <= 25:
        J_factor = 0.945        # Fig 14-8
    elif N_teeth_2 <= 40:
        J_factor = 0.965        # Fig 14-8
    elif N_teeth_2 <= 62:
        J_factor = 0.99         # Fig 14-8
    elif N_teeth_2 <= 112:
        J_factor = 1.0          # Fig 14-8
    elif N_teeth_2 <= 325:
        J_factor = 1.015        # Fig 14-8
    else:
        J_factor = 1.0275       # Fig 14-8

    J_P = J_mod * J_factor      # geometry bending factor
    return J_P

# equation 14-13 or Table 14-8 can be used
# equation 14-13 is used
def elastic_coefficient():  # Eq. 14-12 or Table 14-8
    v_P = 0.3               # pinion Poisson's ratio
    E_P = 30 * 10 ** 6      # pinion Modulus of Elasticity

    v_G = 0.3               # gear Poisson's ratio
    E_G = 30 * 10 ** 6      # gear Modulus of Elasticity

    # arbitrary variables to simplify the function
    C_p_1 = (1 - v_P ** 2) / E_P
    C_p_2 = (1 - v_G ** 2) / E_G

    C_p = math.sqrt(1 / (math.pi * (C_p_1 + C_p_2)))        # elastic coefficient
    return C_p

# equations 14-27 and 14-28
def dynamic_factor(V, Q_v):
    # V is inline pitch velocity (ft/min)
    # Q_v is quality number of gears

    # equation 14-28
    B_v = 0.25 * (12 - Q_v) ** (2 / 3)
    A_v = 50 + 56 * (1 - B_v)

    # equation 14-27
    K_v = ((A_v + math.sqrt(V)) / A_v) ** B_v       # dynamic factor
    return K_v

# uses Table of Overload Factors
def overload_factor():
    K_o = 1.25                       # overload factor for light shock
    return K_o

def load_distribution_factor(d_P, F, S):
    # d_P is pinion pitch diameter
    # F is face width of the narrowest member
    # S is distance between center of bearings

    # equation 14-31
    C_mc = 1            # uncrowned

    # equation 14-32
    if F <= 1:          # in
        C_pf = (F / (10 * d_P)) - 0.025
    elif 1 < F <= 17:   # in
        C_pf = (F / (10 * d_P)) - 0.0375 + 0.0125 * F
    else:
        print(f"face width (F) is out of range")

    # equation 14-33
    S_1 = 0             # centered
    # since our straddle mounted pinion configuration is centered, S_1/S = 0, which is < 0.175
    if S_1 / S < 0.175:
        C_pm = 1
    elif S_1 / S >= 0.175:
        C_pm = 1.1

    # Table 14-9 for commercial, enclosed units
    A_m = 0.127
    B_m = 0.0158
    C_m = -0.930 * 10 ** -4

    # equation 14-34
    C_ma = A_m + B_m * F + C_m * F ** 2

    # equation 14-35, because not adjusted at assembly
    C_e = 1

    # equation 14-30
    K_m = 1 + C_mc * (C_pf * C_pm + C_ma * C_e)     # load distribution factor
    return K_m

# the gear undergoes fewer cycles and is not hardened as much as the pinion
def gear_hardness_ratio_factor(N_G, N_P, d_G, d_P):
    # N_G is number of teeth in the gear
    # N_P is number of teeth in the pinion
    # d_G is gear pitch diameter
    # d_P is pinion pitch diameter

    # equation 14-22
    m_G = speed_ratio(N_G, N_P, d_G, d_P)

    H_BP = 1 #should be changed
    H_HG = 1 #should be changed
    BH = H_BP / H_HG    # Brinell hardness
    BH = 1        # guess
    if BH < 1.2:
        A_H = 0
    elif 1.2 <= BH <= 1.7:
        A_H = (8.98 * 10 ** -3) * BH - (8.29 * 10 ** -3)
    elif BH > 1.7:
        A_H = 0.00698

    # equation 14-36 - use the gear is through hardened
    C_H = 1.0 + A_H * (m_G - 1.0)   # gear hardness factor
    # if surface hardened, it is something else
    return C_H

# C_H only applies to the gear, so C_H = 1
def pinion_hardness_ratio_factor():
    C_H = 1     # pinion hardness factor
    return C_H

# we need to adjust for whether its pinion or gear here and in J
# Fig. 14-14
def bending_stress_cycle_factor(N):
    # N is number of cycles

    Y_N = 1.6831 * (N ** -0.0323)      # Fig 14-14      # bending cycle factor
    return Y_N

# Fig. 14-15
def contact_stress_cycle_factor(N):
    # N is number of cycles

    Z_N = 2.466 * (N ** -0.056)  # Fig 14-15             # contact cycle factor
    return Z_N

# equation 14-38 or Table 14-10
# we are using equation 14-38, because we used the equation to match the decision for the elastic coefficient
def reliability_factor():
    R = 0.90        # reliability
    K_R = 0.658 - 0.0759 * math.log(1 - R, math.e)      # reliability factor
    return K_R

def temperature_factor():
    # project instructions state this is 1
    K_T = 1     # for temp up to 250 F                  # temperature factor
    return K_T

def size_factor():
    # project instructions state this is 1
    K_s = 1                                             # size factor
    return K_s

def surface_condition_factor():
    # project instructions state this is 1
    C_f = 1                                             # surface condition factor
    return C_f

def rim_thickness_factor():
    t_R = None          # look at fig
    h_t = None          # look at fig
    #m_B = t_R / h_t     # backup ratio

    # equation 14-39
    m_B = 1.2           # assuming this

    # equation 14-40
    if m_B < 1.2:
        K_B = 1.6 * math.log((2.242 / m_B), math.e)
    else:
        K_B = 1
    return K_B

# equation 14-41
def bending_safety_factor_AGMA(S_t, Y_N, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J):
    # S_t is bending strength (lbf/in^2)
    # Y_N is bending stress cycle life factor --- this will change with pinion and gear
    # K_T is temperature factor
    # K_R is reliability factor
    # W_t tangential transmitted load (lbf)
    # K_o is overload factor
    # K_v is dynamic factor
    # K_s is size factor
    # K_m is load-distribution factor
    # P_d is transverse diameteral pitch
    # F is face width of narrow member
    # K_m is load-distribution factor
    # K_B is rim-thickness factor
    # J is geometry factor for bending stress

    calc_bend_AGMA = calc_bending_stress_AGMA(W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J)

    fully_corrected_bending_strength = (S_t * Y_N) / (K_T * K_R)

    S_F = fully_corrected_bending_strength / calc_bend_AGMA     # bending safety factor
    return S_F

# equation 14-42
def contact_safety_factor_AGMA(S_c, Z_N, C_H, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I):
    # S_c is allowable contact stress (lbf/in^2)
    # Z_N is wear/ contact stress cycle life factor --- this will change with pinion and gear
    # C_H is hardness ratio factors for pitting resistance --- this will change with pinion and gear
    # K_T is temperature factor
    # K_R is reliability factor
    # C_p is elastic coefficient (sqrt(lbf/in^2))
    # W_t tangential transmitted load (lbf)
    # K_o is overload factor
    # K_v is dynamic factor
    # K_s is size factor
    # K_m is load-distribution factor
    # d_P pitch diameter of the pinion (in)
    # F is face width of narrow member
    # C_f is surface condition factor
    # I is contact geometry factor

    calc_cont_AGMA = calc_contact_stress_AGMA(C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I)

    fully_corrected_contact_strength = (S_c * Z_N * C_H) / (K_T * K_R)

    # must be greater than 1.2
    S_H = fully_corrected_contact_strength / calc_cont_AGMA     # contact safety factor

    return S_H


# end of AGMA Equations and start of gear size Equations
def gear_ratio(rpm_in, rpm_out):
    # rpm_in is the number of rpms of the shaft in or number of teeth in and out...
    # rpm_out is the number of rpms of the shaft out

    #equation 13-5
    m = (rpm_out / rpm_in)
    return m

def pinion_min_teeth(rpm_in, rpm_out, k, helix_angle, pt_angle):   # number of minimum pinion teeth
    # rpm_in is the number of rpms of the shaft in or number of teeth in and out...
    # rpm_out is the number of rpms of the shaft out
    # k is uncrowned teeth factor
    # helix_angle_deg is the helix angle in degrees
    # pt_angle_deg is the transverse pressure angle in degrees

    m = gear_ratio(rpm_in, rpm_out)     # gear ratio

    #convert angles to radians since python wants it radians for its trig functions
    #equation 13-22 to find minimum number of pinion teeth for specified gear ration

    N_teeth = ((2 * k * math.cos(math.radians(helix_angle)))/((1 + 2 * m) * (math.sin(math.radians(pt_angle)) ** 2)) * (m + math.sqrt(m ** 2 + (1 + 2 * m) * math.sin(math.radians(pt_angle)) ** 2)))

    # round up to the nearest tooth
    N_teeth_round = round(N_teeth)
    if N_teeth_round < N_teeth:
        N_teeth = N_teeth_round + 1
    else:
        N_teeth = N_teeth_round
    return N_teeth

def gear_max_teeth(N_P, k, helix_angle, pt_angle):     # maximum number of teeth on the gear
    # N_P is the number of teeth in the pinion
    # k is uncrowned teeth factor
    # helix_angle_deg is the helix angle in degrees
    # pt_angle_deg is the transverse pressure angle in degrees

    #convert angles to radians since python wants it radians for its trig functions
    #equation 13-23 to find maximum number of gear teeth for specified gear ratio
    N_teeth = ((N_P ** 2) * (math.sin(math.radians(pt_angle)) ** 2) - 4 * (k ** 2) * (math.cos(math.radians(helix_angle)) ** 2)) / (4 * k * math.cos(math.radians(helix_angle)) - 2 * N_P * math.sin(math.radians(pt_angle)) ** 2)

    # round up to the nearest tooth
    N_teeth_round = round(N_teeth)
    if N_teeth_round < N_teeth:
        N_teeth = N_teeth_round + 1
    else:
        N_teeth = N_teeth_round
    return N_teeth

def gear_sizes(N_P, N_G, rpm_in, rpm_out):
    # N_P is the minimum number of teeth of the pinion
    # N_G is the maximum number of teeth of the gear
    # rpm_in is the number of rpms of the shaft in
    # rpm_out is the number of rpms of the shaft out

    m = gear_ratio(rpm_in, rpm_out)     # gear ratio

    # define the smallest number of teeth on the pinion and gear using gear ratio equation
    P_size = N_P                # pinion size
    G_size = round(N_P * m)     # gear size
    if G_size > N_G:
        print("the gear size is invalid, try again")
    return P_size, G_size

def check_gear_ratio(input_size, output_size, rpm_in, rpm_out, e):
    # input_size is the input gear size (number of teeth, diameter)
    # output_size is the output gear size (number of teeth, diameter)
    # rpm_in is the number of rpms of the shaft in
    # rpm_out is the number of rpms of the shaft out
    # e is

    m = gear_ratio(rpm_in, rpm_out)  # gear ratio

    actual_ratio = input_size / output_size
    percent_error = abs((actual_ratio - m) / m) * 100
    if percent_error < e:
        print("The specified Gear ratio is within allowable bounds. "
              "The percent difference between the actual gear ratio and ideal gear ratio is " + str(round(percent_error, 4)) + "%")
    else:
        print("the ratio is unnacceptable")
        print(percent_error)
def pitchline_velocity(d_gear, G_rpm):
    # d_gear is the diameter of the gear
    # G_rpm is the rpm of the gear

    #using equation 13-34 the pitchline velocity is calculated in feet/second
    V = math.pi * d_gear * G_rpm / 12    # pitchline velocity (ft/s)
    return V

def transmitted_load_tangential(H, d_gear, G_rpm):
    # H is horsepower of input
    # d_gear is the diameter of the gear
    # G_rpm is the rpm of the gear

    V = pitchline_velocity(d_gear, G_rpm)       # pitchline velocity

    #equation 13-35 is used to calculate the Wt value
    W_t = 33000 * H / V     # tangential transmitted load (lbf)
    return W_t
def normal_diametral_pitch(G_1, G_2, y, c, helix_angle):
    # G_1 is gear 1
    # G_2 is gear 2
    # y is
    # c is
    # helix_angle_deg is the helix angle in degrees

    # this equation was derived from the gearbox diagram for the total height
    P_n = ((G_1 / math.cos(math.radians(helix_angle))) + (G_2 / math.cos(math.radians(helix_angle))) + 2) / (y - c)
    P_n = round(P_n)          # P_n is normal diametrical pitch
    return P_n

def cycles_lifetime(rpm):
    # rpm is the rotations per second of either of the pinion or gear

    #find number of cycles from hour life
    cycles = rpm * 60 * 1000    # converting units from xxx to yyy
    return cycles

def transverse_pressure_angle(pn_angle, helix_angle):
    # pn_angle is the normal pressure angle
    # helix_angle is the helix angle in degrees

    pt_angle = math.degrees(math.atan((math.tan(math.radians(pn_angle)) / math.cos(math.radians(helix_angle)))))    # transverse pressure angle
    return pt_angle

def gear_design(rpm_in, rpm_out, k, helix_angle, N_P, N_G, G_1, G_2, e, d_gear, G_rpm, H, y, c, rpm, pn_angle):
    pt_angle = transverse_pressure_angle(pn_angle, helix_angle)  # transverse pressure angle
    print(f"pt_angle = {pt_angle}, where pn_angle = {pn_angle} and helix_angle = {helix_angle}")

    m = gear_ratio(rpm_in, rpm_out)     # gear ratio
    print(f"m = {m}, where rpm_in = {rpm_in} and rpm_out = {rpm_out}")

    N_P_min = pinion_min_teeth(rpm_in, rpm_out, k, helix_angle, pt_angle)   # minimum number of teeth for the pinion
    print(f"N_P_min = {N_P_min}, where rpm_in = {rpm_in} and rpm_out = {rpm_out}, k = {k}, pn_angle = {pn_angle} and helix_angle = {helix_angle}")

    N_G_max = gear_max_teeth(N_P, k, helix_angle, pt_angle)     # maximum number of teeth for the gear
    print(f"N_G_max = {N_G_max}, where N_P = {N_P}, k = {k}, pn_angle = {pn_angle} and helix_angle = {helix_angle}")

    P_size, G_size = gear_sizes(N_P, N_G, rpm_in, rpm_out)      # pinion and gear sizes
    print(f"P_size = {P_size} and G_size = {G_size}, where rpm_in = {rpm_in} and rpm_out = {rpm_out}")

    check_gear_ratio(G_1, G_2, rpm_in, rpm_out, e)  # ????????????
    print(f"")

    V = pitchline_velocity(d_gear, G_rpm)       # pitchline velocity
    print(f"V = {V}, where d_gear = {d_gear} and G_rpm = {G_rpm}")

    W_t = transmitted_load_tangential(H, d_gear, G_rpm)     # tangential transmitted load
    print(f"W_t = {W_t}, where H = {H}, d_gear = {d_gear} and G_rpm = {G_rpm}")

    P_n = normal_diametral_pitch(G_1, G_2, y, c, helix_angle)       # normal diametrical pitch
    print(f"P_n = {P_n}, where ")

    cycles = cycles_lifetime(rpm)       # lifetime cycles of part
    print(f"cycles = {cycles}, where rpm = {rpm}")


def main():
    allowableWidth = 15
    clearanceAndWallThickness = 1.5
    hp = 150


    #define rpm then gear ratio
    rpmIn = 6700
    rpmOutIdeal = 20000
    m = gear_ratio(rpmIn, rpmOutIdeal)
    print("The gear ratio required is " + str(m))


    #define k value and helix and normal pressure angles for equations 13-22 and 13-23
    k = 1
    helixAngle = 30
    pn_angle = 20


    #equation 13-19 is used to find the transverse pressure angle
    pt_angle = transverse_pressure_angle(pn_angle, helixAngle)
    print("The transverse pressure angle is " + str(pt_angle))


    #minimum number of teeth on the pinion is calculated from equation 13-22
    minPinionTeeth = pinion_min_teeth(rpmIn,rpmOutIdeal, k, helixAngle, pt_angle)
    print("The minimum number of teeth allowable on the pinion is " + str(minPinionTeeth))


    #maximum number of teeth on the gear is calculated using equation 13-23
    maxGearTeeth = gear_max_teeth(minPinionTeeth, k, helixAngle, pt_angle)
    print("The maximum number of teeth allowable on the gear is " + str(maxGearTeeth))


    #gear and pinion sizes are calculated using specified ratio
    minPinionTeeth += 2     # minimum number of pinion teeth
    N_P, N_G = gear_sizes(minPinionTeeth, maxGearTeeth, rpmIn, rpmOutIdeal)



    #allowable tolerance is defined from the gear box specifications
    print("The number of teeth on the pinion and gear are " + str(N_P) + " and " + str(N_G))
    allowablePercentError = 1


    #chosen gear and pinion ratio is checked to ensure it is within tolerance
    print("The allowable percent error between the ideal and actual gear ratios is " + str(allowablePercentError) + "%")
    check_gear_ratio(N_G, N_P, rpmIn, rpmOutIdeal, allowablePercentError)
    actualGearRatio = gear_ratio(N_P, N_G)


    #calculate the actual rpm output of the gear ratio
    rpmOutActual = rpmIn * (actualGearRatio)
    print("The actual gear ratio is " + str(actualGearRatio) + " and the actual output speed is " + str(rpmOutActual) + " rpm")


    #the normal diametral pitch is then selected from table 13-2. Based off the design requirements, height is not a primary concern, but should be limited
    #gear and pinion diameter affect forces on bearings
    #units in teeth/inch
    normalDiametralPitch = normal_diametral_pitch(N_P, N_G, allowableWidth, clearanceAndWallThickness, helixAngle)
    #transverse diametral pitch is calculated using equation 13-18
    P_n = 8
    P_d = P_n * math.cos(math.radians(helixAngle))
    print("Using a Normal pitch of " + str(P_n) + " the Transverse diametral pitch is " + str(P_d) + " teeth per inch")


    #gear and pinion diameters are calculated using equation 13-1
    d_P = N_P / P_d
    d_G = N_G / P_d
    print("The pinion diameter is " + str(d_P) + " inches and the gear diameter is " + str(d_G) + " inches")


    #Pitchline Velocities are calculated
    V = pitchline_velocity(d_G, rpmIn)
    print("The pitchline velocity is " + str(V) + " feet/second")



    #transmitted load is calculated
    W_t = transmitted_load_tangential(hp, d_G, rpmIn)
    print("The transmitted load is " + str(W_t))
    N_cycle_P = cycles_lifetime(rpmOutActual)
    N_cycle_G = cycles_lifetime(rpmIn)
    p_x = ((math.pi/ P_d) / math.tan(math.radians(helixAngle)))
    normalCircularPitch = math.pi / normalDiametralPitch
    Q_v = 11
    F = 1.5
    S = 3.375
    ### OLD ###
    # W_t is tangential transmitted load (lbf)
    # Q_v is quality number of gears
    # V is inline pitch velocity (ft/min)
    # P_d is transverse diameteral pitch
    # N is number of cycles
    # F is face width of narrow member
    # p_x is axial pitch
    # d_P is pitch diameter of the pinion (in)
    # pt_angle is transverse pressure angle
    # N_G is number of teeth of gear
    # N_P is number of teeth of pinion
    # d_G is gear pitch diameter
    # P_n is normal diametrical pitch
    # S is distance between center of bearings
    # J is geometry factor for bending stress including root fillet stress concentration factor

    ### NEW ###
    # W_t is tangential transmitted load (lbf)
    # Q_v is quality number of gears
    # V is inline pitch velocity (ft/min)
    # P_d is transverse diameteral pitch
    # N_cycle_P is number of cycles for the pinion
    # N_cycle_G is the number of cycles for the gear
    # F is face width of narrow member
    # p_x is axial pitch
    # pt_angle is transverse pressure angle
    # N_G is number of teeth of gear
    # N_P is number of teeth of pinion
    # d_P is pitch diameter of the pinion (in)
    # d_G is gear pitch diameter
    # P_n is normal diametrical pitch
    # S is distance between center of bearings

    K_o, K_v, K_s, K_m, K_B, S_t, Y_N_P, Y_N_G, K_T, K_R, C_p, C_f, I, S_c, Z_N_P, Z_N_G, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P = AGMA_coefficients(W_t, Q_v, V, P_d, N_cycle_P, N_cycle_G, F, p_x, pt_angle, N_G, N_P, d_P, d_G, P_n, S)
    #coefficients for gear
    print("Both the gear and pinion face width are " + str(F))
    print("Both the gear and pinion are made from steel and are grade 1 flame or induction hardened as per table 14-3 and 14-6")
    print(S_F_G)
    print(S_H_G)
    print(S_F_P)
    print(S_H_P)
    print(S_H_G)
#def AGMA_coefficients(W_t, Q_v, V, P_d, d_P, N, F, p_x, pt_angle, N_G, N_P, d_G, P_n, S):
main()