import math

#find transmitted torque on the input shaft using equation 3-42 T = (63025 * H) /n
#T is torque in lbf*in
#H is power/Hp
#n is shaft speed in revolutions per minute

#H and n are given in the design requirements
#units in Horsepower
H = 150
#units in rpm
n= 6700

#equation 3-42 is used to calculate the input torque from the input horsepower and input rpm
#Torque units are in-lbf
T_input = 63025 * H / n
print(T_input)


#W_t is the transmitted load brought in the AGMA and gear calculations file
#R_P is the pinion diameter from the AGMA and Gear Calculations file. The pinion is the gear on the output shaft

#units in lbf
W_t = 501.3231747118936

#units in inches
R_P = 1.8763883748662835/2

#the torque on the output shaft is W_t * R_P where W_t is the transmitted loaded and R_P is the pinion radius
#this relation comes from a modified form of equation 13-33 derived from figure 13-29 on page 711 in the book
#units in in-lbf
T_output = W_t * R_P
print(T_output)

#the input and output shaft diameters at the gears are listed below, they are from the forces and diameter calculation file
#units in inches
D_input = 1.417
#units in inches
D_output = 0.944

#from table 7-6 the key sizes are determined and listed below. Note that the keys are square keys
#units in inches
t_Input = 3/8
#units in inches
t_Output = 1/4

#the key material is 1006 HR steel, with Sy = 24 ksi
#unit conversion from Ksi to psi
Sy = 24 * 10**3

Ssy = Sy /math.sqrt(3) #equation from lecture 13, slide 4, also on equation sheet for exam 1, but they are not labeled with numbers

#the required safety factor is 1.5
S_f = 1.5

#calculate key length for input shaft
#Tau_Max is set to the Ssy value as that is the yield strength for the key equations
Tau_Max = Ssy

#The F is calculated using the equation on the equation sheet and in lecture 13, slide 4, which is used in later equations
F = T_input / (D_input/2)

#equation for key shear from equation sheet for exam 1, also given in lecture 13, slide 4
L_Input_Shear = (F * S_f) / (t_Input * Ssy)
#equation for key crushing from equation sheet for exam 1, also given in lecture 13, slide 4
L_Input_Crushing = (2 * F * S_f) / (t_Input * Sy)

#the minimum key length is determined using logic to compare the two checked values. the larger length is the required length
if L_Input_Crushing > L_Input_Shear:
    L_Input = L_Input_Crushing
else:
    L_Input = L_Input_Shear

print("The input shaft key length must be " + str(L_Input) + " inches long")

#calculate Key length for output shaft
#direct shear failure is checked first using equations from lecture 13, slide 4
#crushing is checked using equations from lecture 13
#Tau_Max is set to the Ssy value as that is the yield strength for the key equations
Tau_Max = Ssy

#The F is calculated using the equation on the equation sheet and in lecture 13, slide 4, which is used in later equations
F = T_output / (D_output/2)

#equation for key shear from equation sheet for exam 1, also given in lecture 13, slide 4
L_Output_Shear = (F * S_f) / (t_Output * Ssy)
#equation for key crushing from equation sheet for exam 1, also given in lecture 13, slide 4
L_Output_Crushing = (2 * F * S_f) / (t_Output * Sy)

#the minimum key length is determined by finding the larger of the two lengths
if L_Output_Crushing > L_Output_Shear:
    L_Output = L_Output_Crushing
else:
    L_Output = L_Output_Shear

print("The Output shaft key length must be " + str(L_Output) + " inches long")
