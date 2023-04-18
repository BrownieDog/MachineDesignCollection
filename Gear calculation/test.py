import math
v_P = 0.3  # pinion Poisson's ratio
E_P = 30 * 10 ** 6  # pinion Modulus of Elasticity

v_G = 0.3  # gear Poisson's ratio
E_G = 30 * 10 ** 6  # gear Modulus of Elasticity

C_p_1 = (1 - v_P ** 2) / E_P
C_p_2 = (1 - v_G ** 2) / E_G

C_p = math.sqrt(1 / (math.pi * (C_p_1 + C_p_2)))
print(C_p)