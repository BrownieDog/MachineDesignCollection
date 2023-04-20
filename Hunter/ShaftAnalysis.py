import math

# stress concentration factors for different geometries
def shoulder_fillet_sharp(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.02    # notch radius
    K_t = 2.7       # K_t for bending
    K_ts = 2.2      # K_ts for shear
    return K_t, K_ts, r_notch

def shoulder_fillet_round(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.1  # notch radius
    K_t = 1.7           # K_t for bending
    K_ts = 1.5          # K_ts for shear
    return K_t, K_ts, r_notch

def end_mill_keyseat(d):   # Table 7-1
    # d is diameter of shaft section

    r_notch = d * 0.02  # notch radius
    K_t = 2.14          # K_t for bending
    K_ts = 3.0          # K_ts for shear
    return K_t, K_ts, r_notch

def retaining_ring_groove():   # Table 7-1
    r_notch = 0.01      # notch radius
    K_t = 5.0           # K_t for bending
    K_ts = 3.0          # K_ts for shear
    return K_t, K_ts, r_notch

def specimen_endurance_limit(ult):      # Eqn 6-10      # rotary-beam test specimen endurance limit
    # ult is ultimate tensile strength

    if ult <= 200:          # ksi
        s_ei = 0.5 * ult    # ksi
    elif ult > 200:         # ksi
        s_ei = 100          # ksi
    else:
        print("ult is out of range")
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
    print("diameter: ", d)
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

def temperature_factor():
    k_d = 1  # temperature factor
    return k_d

def reliability_factor():
    k_e = 1  # reliability factor
    return k_e

def miscellaneous_factor():
    k_f = 1  # miscellaneous factor
    return k_f

def modifying_factors(ult,d):
    k_a = surface_factor(ult)      # surface factor
    k_b = size_factor(d)         # size factor
    k_c = loading_factor()      # loading factor
    k_d = temperature_factor()                     # temperature factor
    k_e = reliability_factor()                     # reliability factor
    k_f = miscellaneous_factor()                     # miscellaneous factor
    return k_a, k_b, k_c, k_d, k_e, k_f

def endurance_limit(ult,d):
    # ult is ultimate tensile strength

    k_a, k_b, k_c, k_d, k_e, k_f = modifying_factors(ult,d)
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

    s_e = endurance_limit(ult,d)

    goodman_a_m = 4 * ((K_f * a_m) ** 2)        # alternating moment goodman
    goodman_a_t = 3 * ((K_fs * a_t) ** 2)       # alternating torque goodman
    goodman_a = (1 / (s_e * 1000)) * math.sqrt(goodman_a_m + goodman_a_t)   # the 1000 converts units from psi to ksi

    goodman_s_m = 4 * ((K_f * s_m) ** 2)        # steady moment goodman
    goodman_s_t = 3 * ((K_fs * s_t) ** 2)       # steady torque goodman
    goodman_s = (1 / (s_e * 1000)) * math.sqrt(goodman_s_m + goodman_s_t)   # the 1000 converts units from psi to ksi

    goodman_coeff = 16 / (math.pi * (d ** 3))   # a coefficient used in the final calculation. It simply saves space

    S_F_goodman = 1 / goodman_coeff * (goodman_a + goodman_s)
    return S_F_goodman

def first_cycle_yield_conservative(K_f, K_fs, a_m, a_t, s_m, s_t, d, s_y):
    # K_f is fatigue bending stress concentration factor
    # K_fs is fatigue shear stress concentration factor
    # a_m is alternating moment
    # a_t is alternating torque
    # s_m is steady moment
    # s_t is steady torque
    # d is diameter of shaft section
    # s_y is yield strength

    a_norm = ((32 * K_f * a_m) / (math.pi * (d ** 3)))
    print(f"a_norm = {a_norm}")
    a_shear = ((16 * K_fs * a_t) / (math.pi * (d ** 3)))
    print(f"a_shear = {a_shear}")

    s_norm = (32 * K_f * s_m) / (math.pi * (d ** 3))
    print(f"s_norm = {s_norm}")
    s_shear = (16 * K_fs * s_t) / (math.pi * (d ** 3))
    print(f"s_shear = {s_shear}")

    alt = math.sqrt((a_norm ** 2) + (3 * (a_shear ** 2)))       # total alternating stress


    steady = math.sqrt((s_norm ** 2) + (3 * (s_shear ** 2)))    # total steady stress
    total = alt + steady                                        # total stress

    S_F_conservative = (s_y * 1000) / total       # safety factor against first cycle yield using the conservative approach
    return S_F_conservative

def main(d, a_m, a_t, s_m, s_t, geometry="input geometry"):
    ult, s_y = material_properties()

    if geometry == "sharp":         # sharp shoulder fillet
        K_t, K_ts, r_notch = shoulder_fillet_sharp(d)
        print(f"For sharp shoulder, K_t = {K_t}, K_ts = {K_ts}, r_notch = {r_notch}, where d = {d}")

    elif geometry == "round":       # round shoulder fillet
        K_t, K_ts, r_notch = shoulder_fillet_round(d)
        print(f"For round shoulder, K_t = {K_t}, K_ts = {K_ts}, r_notch = {r_notch}, where d = {d}")

    elif geometry == "keyseat":     # end mill keyseat
        K_t, K_ts, r_notch = end_mill_keyseat(d)
        print(f"For end-mill keyseat, K_t = {K_t}, K_ts = {K_ts}, r_notch = {r_notch}, where d = {d}")

    elif geometry == "retaining":
        K_t, K_ts, r_notch = retaining_ring_groove()
        print(f"For retaining ring groove, K_t = {K_t}, K_ts = {K_ts}, r_notch = {r_notch}")

    else:
        print(f"Please input a valid geometry")

    # this will vary on geometry!
    K_f, q, K_fs, q_shear = fatigue_factor(K_t, K_ts, r_notch, ult)
    print(f"K_f = {K_f}, q = {q}, K_fs = {K_fs}, q_shear = {q_shear}, where K_t = {K_t}, K_ts = {K_ts}, r_notch = {r_notch}, and ult = {ult}")

    # this will vary on geometry!
    a_norm, a_shear = alternating_stress(K_f, K_fs, a_m, a_t, d)
    print(f"a_norm = {a_norm}, a_shear = {a_shear}, where K_f = {K_f}, K_fs = {K_fs}, a_m = {a_m}, a_t = {a_t}, d = {d}")

    # this will vary on geometry!
    s_norm, s_shear = steady_stress(K_f, K_fs, s_m, s_t, d)
    print(f"s_norm = {s_norm}, s_shear = {s_shear}, where K_f = {K_f}, K_fs = {K_fs}, s_m = {s_m}, s_t = {s_t}, d = {d}")

    # this will vary on geometry!
    S_F_goodman = goodman_fatigue_safety_factor(K_f, K_fs, a_m, a_t, s_m, s_t, d, ult)
    print(f"S_F_goodman = {S_F_goodman}, where K_f = {K_f}, K_fs = {K_fs}, a_m = {a_m}, a_t = {a_t}, s_m = {s_m}, s_t = {s_t} d = {d}, ult = {ult}")

    # this will vary on geometry!
    S_F_conservative = first_cycle_yield_conservative(K_f, K_fs, a_m, a_t, s_m, s_t, d, s_y)
    print(f"S_F_conservative = {S_F_conservative}, where K_f = {K_f}, K_fs = {K_fs}, a_m = {a_m}, a_t = {a_t}, s_m = {s_m}, s_t = {s_t} d = {d}, s_y = {s_y}")

    return S_F_goodman, S_F_conservative

