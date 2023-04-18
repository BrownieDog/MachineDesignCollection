import math

def AGMA_coefficients(W_t, Q_v, V, P_d, d_P, N, F, p_x, pt_angle, N_G, N_P, d_G, P_n, S):
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


    K_o = overload_factor()                         # overload factor
    K_v = dynamic_factor(V, Q_v)                    # dynamic factor
    K_s = 1                                         # size factor
    K_m = load_distribution_factor(d_P, F, S)       # load-distribution factor
    K_B = rim_thickness_factor()                    # rim-thickness factor
    J_P = pinion_bending_geometry_factor(p_x, F, N_P)             # J is geometry factor for bending stress including root fillet stress concentration factor - Fig 14-6
    J_G = gear_bending_geometry_factor(p_x, F, N_G)

    S_t = 45000     # guess                         # bending strength (lbf/in^2) - Table 14-3 or 14-4 and Fig 14-2, 14-3, and 14-4
    Y_N = bending_stress_cycle_factor(N)            # Y_N is bending stress cycle life factor
    K_T = temperature_factor()                      # temperature factor
    K_R = reliability_factor()                      # reliability factor

    C_p = elastic_coefficient()                     # elastic coefficient (sqrt(lbf/in^2))
    C_f = 1                      # surface condition factor
    I = contact_geometry_factor(pt_angle, N_G, N_P, d_G, d_P, P_n)  # I is contact geometry factor for pitting resistance

    S_c = 170000    # guess                         # allowable contact stress (lbf/in^2) - Table 14-6, 14-7, and Fig 14-5
    Z_N = contact_stress_cycle_factor(N)            # Z_N is wear/ contact stress cycle life factor
    C_H_G = gear_hardness_ratio_factor(N_G, N_P, d_G, d_P)       # gear hardness ratio factors for pitting resistance
    C_H_P = pinion_hardness_ratio_factor()          # pinion hardness ratio factors for pitting resistance

    # AGMA bending stress factor of safety, a stress ratio
    S_F_G = bending_safety_factor_AGMA(S_t, Y_N, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J_G)
    S_F_P = bending_safety_factor_AGMA(S_t, Y_N, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J_P)

    # AGMA wear/contact factor of safety, a stress ratio
    S_H_G = contact_safety_factor_AGMA(S_c, Z_N, C_H_G, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I)
    S_H_P = contact_safety_factor_AGMA(S_c, Z_N, C_H_P, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I)

    return K_o, K_v, K_s, K_m, K_B, S_t, Y_N, K_T, K_R, C_p, C_f, I, S_c, Z_N, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P

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

def calc_allowable_bending_stress_AGMA(S_t, S_F, Y_N, K_T, K_R):
    # S_t is bending strength (lbf/in^2)
    # S_F is AGMA bending factor of safety
    # Y_N is bending stress cycle life factor
    # K_T is temperature factor
    # K_R is reliability factor

    calc_allow_bend_AGMA = (S_t / S_F) * (Y_N / (K_T * K_R))
    return calc_allow_bend_AGMA

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

def calc_allowable_contact_stress_AGMA(S_c, Z_N, C_H, K_T, K_R, S_H):
    # S_c is allowable contact stress (lbf/in^2)
    # Z_N is wear/ contact stress cycle life factor
    # C_H is hardness ratio factors for pitting resistance
    # K_T is temperature factor
    # K_R is reliability factor
    # S_H is AGMA wear/contact factor of safety, a stress ratio

    calc_allow_cont_AGMA = (S_c / S_H) * ((Z_N * C_H) / K_T * K_R)
    return calc_allow_cont_AGMA

def speed_ratio(N_G, N_P, d_G, d_P):
    # N_G is number of teeth in the gear
    # N_P is number of teeth in the pinion
    # d_G is gear pitch diameter
    # d_P is pinion pitch diameter

    m_G = N_G / N_P
    # OR #
    m_G = d_G / d_P
    return m_G

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

    a = addendum(P_n)

    r_P = d_P / 2     # transverse pitch radius of pinion
    r_bP = r_P * math.cos(pt_angle) # base-circle radii for pinion

    r_G = d_G / 2     # transverse pitch radius of gear
    r_bG = r_G * math.cos(pt_angle) # base-circle radii for gear

    Z_1 = (r_P + a) ** 2 - r_bP ** 2
    Z_2 = (r_G + a) ** 2 - r_bG ** 2
    Z_3 = (r_P + r_G) * math.sin(pt_angle)

    if Z_1 > Z_3:
        Z = math.sqrt(Z_3) + math.sqrt(Z_2) - Z_1
    elif Z_2 > Z_3:
        Z = math.sqrt(Z_1) + math.sqrt(Z_3) - Z_2
    else:
        Z = math.sqrt(Z_1) + math.sqrt(Z_2) - Z_3

    p_n = N_G / d_G                 # normal circular pitch
    # round up to a standard

    p_N = p_n * math.cos(pt_angle)
    m_N = p_N / (0.95 * Z)          # load sharing ratio

    m_G = speed_ratio(N_G, N_P, d_G, d_P)

    I_ext = (math.cos(pt_angle) * math.sin(pt_angle)) / (2 * m_N) * (m_G / (m_G + 1))    # external gear
    I_int = (math.cos(pt_angle) * math.sin(pt_angle)) / (2 * m_N) * (m_G / (m_G - 1))    # internal gear

    return I_ext                    # external gear

def pinion_bending_geometry_factor(p_x, F, N_P):
    # p_x is axial pitch
    # F is narrow face width
    # N_P is number of teeth in the pinion

    m_F = F / p_x

    if m_F >= 2:
        # J_mod = None            # Fig 14-7
        # J_factor = None       # Fig 14-8

        if N_P <= 25:
            J_mod = 0.465
        elif N_P <= 45:
            J_mod = 0.5
        elif N_P <= 105:
            J_mod = 0.54
        elif N_P <= 325:
            J_mod = 0.565
        else:
            J_mod = 0.58

        if N_P <= 25:
            J_factor = 0.945
        elif N_P <= 40:
            J_factor = 0.965
        elif N_P <= 62:
            J_factor = 0.99
        elif N_P <= 112:
            J_factor = 1.0
        elif N_P <= 325:
            J_factor = 1.015
        else:
            J_factor = 1.0275

    else:
        print(":(")

    J_P = J_mod * J_factor
    return J_P

def gear_bending_geometry_factor(p_x, F, N_G):
    # p_x is axial pitch
    # F is narrow face width
    # N_G is number of teeth in the gear

    m_F = F / p_x

    if m_F >= 2:
        # J_mod = None            # Fig 14-7
        # J_factor = None       # Fig 14-8
        if N_G <= 25:
            J_mod = 0.465
        elif N_G <= 45:
            J_mod = 0.5
        elif N_G <= 105:
            J_mod = 0.54
        elif N_G <= 325:
            J_mod = 0.565
        else:
            J_mod = 0.58

        if N_G <= 25:
            J_factor = 0.945
        elif N_G <= 40:
            J_factor = 0.965
        elif N_G <= 62:
            J_factor = 0.99
        elif N_G <= 112:
            J_factor = 1.0
        elif N_G <= 325:
            J_factor = 1.015
        else:
            J_factor = 1.0275

    else:
        print(":(")

    J_G = J_mod * J_factor
    return J_G

def elastic_coefficient():  # Eq. 14-12 or Table 14-8
    v_P = 0.3               # pinion Poisson's ratio
    E_P = 30 * 10 ** 6      # pinion Modulus of Elasticity

    v_G = 0.3               # gear Poisson's ratio
    E_G = 30 * 10 ** 6      # gear Modulus of Elasticity

    C_p_1 = (1 - v_P ** 2) / E_P
    C_p_2 = (1 - v_G ** 2) / E_G

    C_p = math.sqrt(1 / (math.pi * (C_p_1 + C_p_2)))        # elastic coefficient
    return C_p

def dynamic_factor(V, Q_v):
    # V is inline pitch velocity (ft/min)
    # Q_v is quality number of gears

    B_v = 0.25 * (12 - Q_v) ** (2 / 3)
    A_v = 50 + 56 * (1 - B_v)

    K_v = ((A_v + math.sqrt(V)) / A_v) ** B_v       # dynamic factor
    return K_v

def overload_factor():
    K_o = 1.25     # check page 882-758             # overload factor
    return K_o

def load_distribution_factor(d_P, F, S):
    # d_P is pinion pitch diameter
    # F is face width of the narrowest member
    # S is distance between center of bearings

    C_mc = 1            # uncrowned

    if F <= 1:          # in
        C_pf = (F / (10 * d_P)) - 0.025
    elif 1 < F <= 17:   # in
        C_pf = (F / (10 * d_P)) - 0.0375 + 0.0125 * F
    else:
        print(":(")     # I am sad

    S_1 = 0             # centered

    if S_1 / S < 0.175:
        C_pm = 1
    elif S_1 / S >= 0.175:
        C_pm = 1.1

    # commercial, enclosed units
    A_m = 0.127
    B_m = 0.0158
    C_m = -0.930 * 10 ** -4
    C_ma = A_m + B_m * F + C_m * F ** 2

    C_e = 1     # not adjusted at assembly

    K_m = 1 + C_mc * (C_pf * C_pm + C_ma * C_e)     # load distribution factor
    return K_m

def gear_hardness_ratio_factor(N_G, N_P, d_G, d_P):
    # N_G is number of teeth in the gear
    # N_P is number of teeth in the pinion
    # d_G is gear pitch diameter
    # d_P is pinion pitch diameter

    m_G = speed_ratio(N_G, N_P, d_G, d_P)

    H_BP = None
    H_HG = None
    BH = H_BP / H_HG    # Brinell hardness
    BH = 1        # guess
    if BH < 1.2:
        A_H = 0
    elif 1.2 <= BH <= 1.7:
        A_H = (8.98 * 10 ** -3) * BH - (8.29 * 10 ** -3)
    elif BH > 1.7:
        A_H = 0.00698

    C_H = 1.0 + A_H * (m_G - 1.0)   # gear hardness factor
    # if surface hardened, it is something else
    return C_H

def pinion_hardness_ratio_factor():
    C_H = 1     # pinion hardness factor
    return C_H

def bending_stress_cycle_factor(N):
    # N is number of cycles

    Y_N = 1.6831 * (N ** -0.0323)      # Fig 14-14      # bending cycle factor
    return Y_N

def contact_stress_cycle_factor(N):
    # N is number of cycles

    Z_N = 2.66 * (N ** -0.056)  # Fig 14-15             # contact cycle factor
    return Z_N

def reliability_factor():
    R = 0.90        # reliability
    K_R = 0.658 - 0.0759 * math.log(1 - R, math.e)      # reliability factor
    return K_R

def temperature_factor():
    K_T = 1     # for temp up to 250 F                  # temperature factor
    return K_T

def rim_thickness_factor():
    t_R = None          # look at fig
    h_t = None          # look at fig
    m_B = t_R / h_t     # backup ratio
    m_B = 1.2           # assuming this
    if m_B < 1.2:
        K_B = 1.6 * math.log((2.242 / m_B), math.e)
    else:
        K_B = 1
    return K_B

def bending_safety_factor_AGMA(S_t, Y_N, K_T, K_R, W_t, K_o, K_v, K_s, P_d, F, K_m, K_B, J):
    # S_t is bending strength (lbf/in^2)
    # Y_N is bending stress cycle life factor
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

    S_F = fully_corrected_bending_strength / calc_bend_AGMA
    return S_F

def contact_safety_factor_AGMA(S_c, Z_N, C_H, K_T, K_R, C_p, W_t, K_o, K_v, K_s, K_m, d_P, F, C_f, I):
    # S_c is allowable contact stress (lbf/in^2)
    # Z_N is wear/ contact stress cycle life factor
    # C_H is hardness ratio factors for pitting resistance
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
    S_H = fully_corrected_contact_strength / calc_cont_AGMA

    return S_H
    # must be greater than 1.2

def gearratio(In, Out):
    #equation 13-5
    ratio = (Out/In)
    return ratio

def minimumNumberOfTeethOnPinion(m, k, helixAngle, transversePressureAngle):
    #convert angles to radians since python wants it radians for its trig functions
    helixAngleRads = math.radians(helixAngle)
    transversePressureAngleRads = math.radians(transversePressureAngle)
    #equation 13-22 to find minimum number of pinion teeth for specified gear ration
    teethNumber = ((2*k*math.cos(helixAngleRads))/((1+2*m)*(math.sin(transversePressureAngleRads)**2))*(m+math.sqrt(m**2+(1+2*m)*math.sin(transversePressureAngleRads)**2)))
    #code to round up to nearest tooth
    roundedTeethNumber = round(teethNumber)
    if roundedTeethNumber < teethNumber:
        teethNumber = roundedTeethNumber + 1
    else:
        teethNumber = roundedTeethNumber
    return teethNumber

def maximumNumberOfTeethOnGear(Np, k, helixAngle, transversePressureAngle):
    #convert angles to radians since python wants it radians for its trig functions
    helixAngleRads = math.radians(helixAngle)
    transversePressureAngleRads = math.radians(transversePressureAngle)
    #equation 13-23 to find maximum number of gear teeth for specified gear ratio
    teethNumber = ((Np**2)*(math.sin(transversePressureAngleRads)**2)-4*(k**2)*(math.cos(helixAngleRads)**2))/(4*k*math.cos(helixAngleRads)-2*Np*math.sin(transversePressureAngleRads)**2)
    #code to round up to nearest tooth
    roundedTeethNumber = round(teethNumber)
    if roundedTeethNumber < teethNumber:
        teethNumber = roundedTeethNumber + 1
    else:
        teethNumber = roundedTeethNumber
    return teethNumber

def gearSizes(Np, Ng, gearRatio):
    #define the smallest number of teeth on the pinion and gear using gear ratio equation
    pinionSize = Np
    gearSize = round(Np * gearRatio)
    if gearSize > Ng:
        print("the gear size is invalid, try again")
    return pinionSize, gearSize

def checkGearRatio(N1, N2, m, e):
    actualRatio = N2/N1
    percentError = abs((actualRatio-m)/m)*100
    if percentError < e:
        print("The specified Gear ratio is within allowable bounds. "
              "The percent difference between the actual gear ratio and ideal gear ratio is " + str(round(percentError, 4)) + "%")
    else:
        print("BAD RATIO DUMBASS")
        print(percentError)
def pitchlineVelocity(d, n):
    #using equation 13-34 the pitchline velocity is calculated in feet/second
    V = math.pi * d * n / 12
    return V

def transmittedLoad(H, V):
    #equation 13-35 is used to calculate the Wt value
    Wt = 33000 * H / V
    return Wt
def findNormalDiametralPitch(N2, N3, y, c, helixAngle):
    #this equation was derived from the gearbox diagram for the total height
    Pn=((N2/math.cos(math.radians(helixAngle)))+(N3/math.cos(math.radians(helixAngle)))+2)/(y-c)
    Pn = round(Pn)
    return Pn

def cyclesLifetime(rpm):
    #find number of cycles from hour life
    cycles = rpm * 60 * 1000
    return cycles

def main():
    W_t = 1
    Q_v =1
    V = 1
    P_d = 1
    N = 1
    F = 2
    p_x = 1
    d_P = 1
    pt_angle = 1
    N_G = 33
    N_P = 11
    d_G = 1
    P_n  = 1
    S = 1
    allowableWidth = 15
    clearanceAndWallThickness = 1.5
    hp = 150


    #define rpm then gear ratio
    rpmIn = 6700
    rpmOutIdeal = 20000
    m = gearratio(rpmIn, rpmOutIdeal)
    print("The gear ratio required is " + str(m))


    #define k value and helix and normal pressure angles for equations 13-22 and 13-23
    k=1
    helixAngle = 30
    normalPressureAngle = 20


    #equation 13-19 is used to find the transverse pressure angle
    transversePressureAngle =math.degrees(math.atan((math.tan(math.radians(normalPressureAngle))/math.cos(math.radians(helixAngle)))))
    print("The transverse pressure angle is " +str(transversePressureAngle))


    #minimum number of teeth on the pinion is calculated from equation 13-22
    minPinionTeeth = minimumNumberOfTeethOnPinion(m,k,helixAngle,transversePressureAngle)
    print("The minimum number of teeth allowable on the pinion is " + str(minPinionTeeth))


    #maximum number of teeth on the gear is calculated using equation 13-23
    maxGearTeeth = maximumNumberOfTeethOnGear(minPinionTeeth, k, helixAngle, transversePressureAngle)
    print("The maximum number of teeth allowable on the gear is " + str(maxGearTeeth))


    #gear and pinion sizes are calculated using specified ratio
    minPinionTeeth +=2
    numberOfPinionTeeth, numberOfGearTeeth = gearSizes(minPinionTeeth, maxGearTeeth, m)
    numberOfGearTeeth -= 0


    #allowable tolerance is defined from the gear box specifications
    print("The number of teeth on the pinion and gear are " + str(numberOfPinionTeeth) + " and " + str(numberOfGearTeeth))
    allowablePercentError = 1


    #chosen gear and pinion ratio is checked to ensure it is within tolerance
    print("The allowable percent error between the ideal and actual gear ratios is " + str(allowablePercentError) + "%")
    checkGearRatio(numberOfPinionTeeth, numberOfGearTeeth, m, allowablePercentError)
    actualGearRatio = gearratio(numberOfPinionTeeth, numberOfGearTeeth)


    #calculate the actual rpm output of the gear ratio
    rpmOutActual = rpmIn * (actualGearRatio)
    print("The actual gear ratio is " + str(actualGearRatio) + " and the actual output speed is " + str(rpmOutActual) + " rpm")


    #the normal diametral pitch is then selected from table 13-2. Based off the design requirements, height is not a primary concern, but should be limited
    #gear and pinion diameter affect forces on bearings
    #units in teeth/inch
    normalDiametralPitch = findNormalDiametralPitch(numberOfPinionTeeth, numberOfGearTeeth, allowableWidth, clearanceAndWallThickness, helixAngle)
    #transverse diametral pitch is calculated using equation 13-18
    normalDiametralPitch = 8
    transverseDiametralPitch = normalDiametralPitch * math.cos(math.radians(helixAngle))
    print("Using a Normal pitch of " + str(normalDiametralPitch) +  " the Transverse diametral pitch is " + str(transverseDiametralPitch) + " teeth per inch")


    #gear and pinion diameters are calculated using equation 13-1
    pinionDiameter = numberOfPinionTeeth / transverseDiametralPitch
    gearDiameter = numberOfGearTeeth / transverseDiametralPitch
    print("The pinion diameter is " + str(pinionDiameter) + " inches and the gear diameter is " + str(gearDiameter) + " inches")


    #Pitchline Velocities are calculated
    pitchline = pitchlineVelocity(gearDiameter, rpmIn)
    print("The pitchline velocity is " + str(pitchline) + " feet/second")



    #transmitted load is calculated
    Wt = transmittedLoad(hp, pitchline)
    print("The transmitted load is " + str(Wt))
    cyclesPinion = cyclesLifetime(rpmOutActual)
    cyclesGear = cyclesLifetime(rpmIn)
    axialPitch = transverseDiametralPitch/ (math.pi * math.tan(math.radians(helixAngle)))
    Q_V = 11
    F = 2
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

    K_o, K_v, K_s, K_m, K_B, S_t, Y_N, K_T, K_R, C_p, C_f, I, S_c, Z_N, C_H_G, C_H_P, S_F_G, S_F_P, S_H_G, S_H_P = AGMA_coefficients(Wt, Q_V,
            pitchline, transverseDiametralPitch, cyclesGear, F, axialPitch, pinionDiameter, transversePressureAngle, numberOfGearTeeth,
            numberOfPinionTeeth, gearDiameter, normalDiametralPitch, S)