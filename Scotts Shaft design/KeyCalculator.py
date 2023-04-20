import math

#find transmitted torque on the input shaft using equation 3-42 T = (63025 * H) /n
#T is torque in lbf*in
#H is power/Hp
#n is shaft speed in revolutions per minute

#H and n are given in the design requirements
H = 150
n= 6700
T_input = 63025 * H / n
print(T_input)

#the torque on the output shaft is W_t * R_P where W_t is the transmitted loaded and R_P is the pinion radius
W_t = 501.3231747118936
R_P = 1.8763883748662835/2
T_output = W_t * R_P
print(T_output)

#the input and output shaft diameters at the gears are listed below
#once the shaft sizes are determined update these
D_input = 1.5
D_output = 1.25

#from table 7-6 the key sizes are determined and listed below. the keys are square keys
#once the shaft sizes are determined update these
t_Input = 5/16
t_Output = 1/4

#the key material is 1006 HR steel, with Sy = 24 ksi
Sy = 24 * 10**3
Ssy = Sy /math.sqrt(3) #equation from lecture 13, slide 4

#the required safety factor is 1.5
S_f = 1.5

#calculate key length for input shaft
#direct shear failure is checked first using equations from lecture 13, slide 4
#crushing is checked using equations from lecture 13
Tau_Max = Ssy
F = T_input / (D_input/2)

L_Input_Shear = (F * S_f) / (t_Input * Ssy)
L_Input_Crushing = (2 * F * S_f) / (t_Input * Sy)

#the minimum key length is determined
if L_Input_Crushing > L_Input_Shear:
    L_Input = L_Input_Crushing
else:
    L_Input = L_Input_Shear

print("The input shaft key length must be " + str(L_Input) + " inches long")

#calculate Key length for output shaft
#direct shear failure is checked first using equations from lecture 13, slide 4
#crushing is checked using equations from lecture 13
Tau_Max = Ssy
F = T_output / (D_output/2)
L_Output_Shear = (F * S_f) / (t_Output * Ssy)
L_Output_Crushing = (2 * F * S_f) / (t_Output * Sy)

#the minimum key length is determined
if L_Output_Crushing > L_Output_Shear:
    L_Output = L_Output_Crushing
else:
    L_Output = L_Output_Shear

print("The Output shaft key length must be " + str(L_Output) + " inches long")
