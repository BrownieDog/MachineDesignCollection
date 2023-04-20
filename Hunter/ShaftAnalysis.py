import math

# stress concentration factors for different geometries
def shoulder_fillet_sharp(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.02    # notch radius
    K_t = 2.7       # K_t for bending
    K_ts = 2.2      # K_ts for shear
    return K_t, K_ts

def shoulder_fillet_round(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.02    # notch radius
    K_t = 1.7       # K_t for bending
    K_ts = 1.5      # K_ts for shear
    return K_t, K_ts

def retaining_ring_groove():   # Table 7-1
    r_notch = 0.01        # notch radius
    K_t = 5.0       # K_t for bending
    K_ts = 3.0      # K_ts for shear
    return K_t, K_ts

def end_mill_keyseat(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.02    # notch radius
    K_t = 2.14      # K_t for bending
    K_ts = 3.0      # K_ts for shear
    return K_t, K_ts

def specimen_endurance_limit(ult):      # Eqn 6-10      # rotary-beam test specimen endurance limit
    # ult is ultimate tensile strength

    if ult <= 200:          # ksi
        s_ei = 0.5 * ult    # ksi
    elif ult > 200:         # ksi
        s_ei = 100          # ksi
    else:
        print(f"ult is out of range")
    return s_ei

# find coefficients of endurance limit
def surface_factor(ult):                # Table 6-2 for machined surface in ksi
    # ult is ultimate tensile strength

    a = 2.00                # from Table 6-2 for machined surface
    b = -0.217              # from Table 6-2 for machined surface

    k_a = a * (ult ** b)    # surface factor
    return k_a

def size_factor(d):         # Eqn 6-19
    # d is diameter of shaft section

    # for rotating round specimens in bending
    if 0.11 <= d <= 2:  # inch
        k_b = 0.879 * (d ** -0.107)  # inch
    elif 2 < d < 10:    # inch
        k_b = 0.91 * (d ** -0.157)  # inch
    else:
        print("The shaft diameter is not in range")
    return k_b

def loading_factor():       # k_c - load factor
    k_c = 1  # combined loading with von mises stress calculations
    return k_c

def modifying_factors():
    k_a = surface_factor()      # surface factor
    k_b = size_factor()         # size factor
    k_c = loading_factor()      # loading factor
    k_d = 1                     # temperature factor
    k_e = 1                     # reliability factor
    k_f = 1                     # miscellaneous factor
    return k_a, k_b, k_c, k_d, k_e, k_f

def endurance_limit(ult):
    # ult is ultimate tensile strength

    k_a, k_b, k_c, k_d, k_e, k_f = modifying_factors()
    s_ei = specimen_endurance_limit(ult)

    s_e = k_a * k_b * k_c * k_d * k_e * k_f * s_ei      # endurance limit
    return s_e

def neuber_constant(ult):
    # ult is ultimate tensile strength

    if 50 <= ult <= 250:    # ksi
        neuber_bend_1 = (3.08 * 10 ** -3) * ult
        neuber_bend_2 = (1.51 * 10 ** -5) * (ult ** 2)
        neuber_bend_3 = (2.67 * 10 ** -8) * (ult ** 3)
        neuber_bend = 0.246 - neuber_bend_1 + neuber_bend_2 - neuber_bend_3

    else:
        print("Out of range :(")

    if 50 <= ult <= 220:    # ksi
        neuber_tor_1 = (2.51 * 10 ** -3) * ult
        neuber_tor_2 = (1.35 * 10 ** -5) * (ult ** 2)
        neuber_tor_3 = (2.67 * 10 ** -8) * (ult ** 3)
        neuber_tor = 0.190 - neuber_tor_1 + neuber_tor_2 - neuber_tor_3

    else:
        print("Out of range :(")

    return neuber_bend, neuber_tor

### fatigue factor
def fatigue_factor(K_t, K_ts, r_notch, ult):
    # K_t is bending stress concentration factor
    # K_ts is shear stress concentration factor
    # r_notch is the notch radius
    # ult is the ultimate tensile strength

    sqrt_r = math.sqrt(r_notch)  # r is notch radius
    neuber_bend, neuber_tor = neuber_constant(ult)

    # bending - use a for bending/ axial
    K_f = 1 + ((K_t - 1) / (1 + (neuber_bend / sqrt_r)))        # fatigue stress concentration factor
    q = (K_f - 1) / (K_t - 1)                               # q is notch sensitivity

    # shear - use a for torsion
    K_fs = 1 + ((K_ts - 1) / (1 + (neuber_tor / sqrt_r)))       # fatigue stress concentration factor
    q_shear = (K_fs - 1) / (K_ts - 1)                       # q is notch sensitivity

    return K_f, q, K_fs, q_shear

def material_properties():
    # HR 1030 Steel
    ult = 68        # ksi   ultimate strength
    s_y = 37.5      # ksi   yield strength
    return ult, s_y

def alternating_load():
    alternating_moment = 0
    alternating_torque = 0

def steady_load():
    steady_moment = 0
    steady_torque = 0

def alternating_stress(K_f, K_fs, a_m, a_t, d):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # d is diameter of shaft section

    a_norm = ((32 * K_f * a_m) / (math.pi * (d ** 3)))
    a_shear = ((16 * K_fs * a_t) / (math.pi * (d ** 3)))
    return a_norm, a_shear

def steady_stress(K_f, K_fs, s_m, s_t, d):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section

    s_norm = (32 * K_f * s_m) / (math.pi * (d ** 3))
    s_shear = (16 * K_fs * s_t) / (math.pi * (d ** 3))
    return s_norm, s_shear

def total_stress(K_f, K_fs, a_m, a_t, s_m, s_t, d):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section

    a_norm, a_shear = alternating_stress(K_f, K_fs, a_m, a_t, d)
    s_norm, s_shear = steady_stress(K_f, K_fs, s_m, s_t, d)
    return a_norm, a_shear, s_norm, s_shear

def goodman_fatigue_safety_factor(K_f, K_fs, a_m, a_t, s_m, s_t, d, ult):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section
    # ult is the ultimate tensile strength

    a_norm, a_shear, s_norm, s_shear = total_stress(K_f, K_fs, a_m, a_t, s_m, s_t, d)

    s_e = endurance_limit(ult)

    goodman_a_m = 4 * ((K_f * a_m) ** 2)        # alternating moment goodman
    goodman_a_t = 3 * ((K_fs * a_t) ** 2)       # alternating torque goodman
    goodman_a = (1 / (s_e * 1000)) * math.sqrt(goodman_a_m + goodman_a_t)   # the 1000 converts units from psi to ksi

    goodman_s_m = 4 * ((K_f * s_m) ** 2)        # steady moment goodman
    goodman_s_t = 3 * ((K_fs * s_t) ** 2)       # steady torque goodman
    goodman_s = (1 / (s_e * 1000)) * math.sqrt(goodman_s_m + goodman_s_t)   # the 1000 converts units from psi to ksi

    goodman_coeff = 16 / (math.pi * (d ** 3))   # a coefficient used in the final calculation. It simply saves space

    S_F_goodman = 1 / goodman_coeff * (goodman_a + goodman_s)
    return S_F_goodman

def first_cycle_yield_vm(K_f, K_fs, a_m, a_t, s_m, s_t, d, s_y):        # checking 1st cycle yield with Von Mises
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section
    # s_y is yield strength 

    a_norm, a_shear, s_norm, s_shear = total_stress(K_f, K_fs, a_m, a_t, s_m, s_t, d)

    tot_norm = a_norm + s_norm          # total normal stress
    tot_shear = a_shear + s_shear       # total shear stress
    max_vm = math.sqrt((tot_norm ** 2) + (3 * (tot_shear ** 2)))        # vm stress
    n_vm = (s_y * 1000) / max_vm        # safety factor against first cycle yield using vm
    return n_vm

def first_cycle_yield_conservative(K_f, K_fs, a_m, a_t, s_m, s_t, d, s_y):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section
    # s_y is yield strength

    a_norm, a_shear, s_norm, s_shear = total_stress(K_f, K_fs, a_m, a_t, s_m, s_t, d)

    alt = math.sqrt((a_norm ** 2) + (3 * (a_shear ** 2)))       # total alternating stress
    steady = math.sqrt((s_norm ** 2) + (3 * (s_shear ** 2)))    # total steady stress
    total = alt + steady                                        # total stress

    S_F_conservative = (s_y * 1000) / total       # safety factor against first cycle yield using the conservative approach
    return S_F_conservative

def shear_diagram(f):

    # analysis point: g, h, i, j, k, l, m, n, o, p
    return

class AnalysisPoint():
    def __init__(self, x, f, d, geometry):
        self.x = x
        self.f = f
        self.d = d
        self.geometry = geometry


def find_shear_and_moment_values():
    analysis_point_input_list = analysis_point_input()          # list of the analysis points for the input shaft
    analysis_point_output_list = analysis_point_output()        # list of the analysis points for the output shaft

    analysis_point = analysis_point_input_list
    analysis_point = analysis_point_output_list

    # initialize shear and moment lists
    shear_force = []
    bending_moment =[]

    for i in range(0, len(analysis_point)):                 # analysis point object for the input and output shafts are
                                                            # created below. Since they are the same length, we can use
                                                            # either; however, I don't know how to adjust if they were
                                                            # different lengths
        point = analysis_point[i]   # loop through the analysis point list

        if point.x == 0:
            shear_force = point.f   # appending the shear force list with the force value at the specific point of interest
            bending_moment = 0      # bending moment at x = 0

            new_shear = shear_force
            new_moment  = bending_moment
            new_location = point.x

        else:
            old_shear = new_shear
            old_moment = new_moment

            old_location = new_location
            new_location point.x
            moment_distance = new_location - old_location

            shear_force =old_shear + point.f
            bending_moment = old_moment + (point.f * moment_distance)

def input_shaft_diameters():
    di_1 =
    di_2 =
    di_3 =
    di_4 =
    di_5 =
    di_6 =
    di_7 =
    return di_1, di_2, di_3, di_4, di_5, di_6, di_7

def output_shaft_diameters():
    do_1 =
    do_2 =
    do_3 =
    do_4 =
    do_5 =
    do_6 =
    do_7 =
    return do_1, do_2, do_3, do_4, do_5, do_6, do_7

def analysis_point_input():
    di_1, di_2, di_3, di_4, di_5, di_6, di_7 = input_shaft_diameters()

    point_start_input = AnalysisPoint(x = 0, f = f, d = di_1,geometry = )
    point_g = AnalysisPoint(x = 2, f = , d = di_1, geometry = keyway)
    point_h = AnalysisPoint(x = 2.75, f = , d = , geometry = )
    point_i = AnalysisPoint(x = 3, f = , d = , geometry = clip)
    point_j = AnalysisPoint(x = 3.75, f = , d = , geometry = )
    point_k = AnalysisPoint(x = , f = , d = , geometry = )
    point_l = AnalysisPoint(x = , f = , d = , geometry = )
    point_m = AnalysisPoint(x = , f = , d = , geometry = )
    point_n = AnalysisPoint(x = , f = , d = , geometry = )
    point_o = AnalysisPoint(x = , f = , d = , geometry = )
    point_p = AnalysisPoint(x = 6.75, f = , d = , geometry = )
    point_end_input = AnalysisPoint(x = 8, f = , d = , geometry = )

    analysis_point_input = [point_start_input, point_g, point_h, point_i, point_j, point_k, point_l, point_m, point_n,
                            point_o, point_p, point_end_input]
    return analysis_point_input

def analysis_point_output():
    do_1, do_2, do_3, do_4, do_5, do_6, do_7 = output_shaft_diameters()

    point_start_output = AnalysisPoint(x = 0, f = , d = , geometry = )
    point_q = AnalysisPoint(x = 2, f = , d = , geometry = keyway)
    point_r = AnalysisPoint(x = 2.75, f = , d = , geometry = )
    point_s = AnalysisPoint(x = 3, f = , d = , geometry = clip)
    point_t = AnalysisPoint(x = 3.75, f = , d = , geometry = )
    point_u = AnalysisPoint(x = , f = , d = , geometry = )
    point_v = AnalysisPoint(x = , f = , d = , geometry = )
    point_w = AnalysisPoint(x = , f = , d = , geometry = )
    point_x = AnalysisPoint(x = , f = , d = , geometry = )
    point_y = AnalysisPoint(x = , f = , d = , geometry = )
    point_z = AnalysisPoint(x = 6.75, f = , d = do_7, geometry = )
    point_end_output = AnalysisPoint(x = 8, f = , d = do_7, geometry = )

    analysis_point_input = [point_start_output, point_q, point_r, point_s, point_t, point_u, point_v, point_w, point_x,
                            point_y, point_z, point_end_output]
    return analysis_point_input