from ShaftAnalysis import *
from DiameterFinder import *

import math
sin = math.sin
cos = math.cos
tan = math.tan

hp = 150                                    ## Input power (HP)
Ohi = 6700                                  ## input rotation speed (RPM)
m = 3                                       ## gear ration
Oho = m*Ohi                                 ## output rotation speed
k = 1                                       ## full depth teeth
psi = math.radians(30)                      ## helix angle
phin = math.radians(20)                     ## Normal pressure angle
phit = (math.atan(tan(phin)/cos(psi)))      ## equation 13-19 Transverse pressure angle

face = 1.5                                  ## gear face width

Hlb = 550*hp*60*12/(2*math.pi)              ## input power in lbf
Ti = Hlb/Ohi                                ## input torque
To = Hlb/Oho                                ## output torque

dg = 5.629165124598851                      ## gear diameter
dp = 1.8763883748662835                     ## pinion diameter
V = math.pi*dp*Oho/12                       ## EQ. 13-34 Pitch line velocity

## Gear forces

Wt = 33000*hp/V                             ## eq. 13-35 Transverse load on gear
Wr = Wt*tan(phit)                           ## eq. 13-40 Radial load on gear
Wa = Wt*tan(psi)                            ## eq. 13-40 Radial load on gear
W = Wt/(cos(phin)*cos(psi))                 ## eq. 13-40 overall load on gear

## Reaction forces at bearings
## forces found using sum of moments and sum of forces

Fbx = Wa
Fby = (1.875*Wr + (dg/2) * Wa)/3.75
Fbz = 1.875*Wt/3.75

Fay = -(Wr-Fby)
Faz = Wt-Fbz

Fdy = (1.875*Wr-(dp/2)*Wa)/3.75
Fdz = 1.875*Wt/3.75

Fcx = Wa
Fcy = (Wr-Fdy)
Fcz = Wt-Fdz

print("\nGear forces")
print("Wt: ", Wt)
print("Wr: ", Wr)
print("Wa: ", Wa)

print("\nInput shaft")
print("Fay: ", Fay)
print("Faz: ", Faz)
print("Fbx: ", Fbx)
print("Fby: ", Fby)
print("Fbz: ", Fbz)

print("\nOutput shaft")
print("Fcx: ", Fcx)
print("Fcy: ", Fcy)
print("Fcz: ", Fcz)
print("Fdy: ", Fdy)
print("Fdz: ", Fdz)


## Moment analysis at each point
def InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg):
    A = 3.375
    E = 5.25
    if (x - A) > 0:
        if (x - E) > 0:
            My = (x - A) * Fay + (x - E) * Wr - (dg / 2) * Wa
            Mz = (x - A) * Faz - (x - E) * Wt
            a_m = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (x - A) * Fay
            Mz = (x - A) * Faz
            a_m = math.sqrt(My ** 2 + Mz ** 2)
    else: a_m = 0
    return a_m

def OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp):
    D = .5+.75+3+.325
    F = 3.25
    if (D - x) > 0:
        if (F - x) > 0:
            My = (D - x) * Fdy - (F - x) * Wr - (dp / 2) * Wa
            Mz = (D - x) * Fdz - (F - x) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (D - x) * Fdy
            Mz = (D - x) * Fdz
            M = math.sqrt(My ** 2 + Mz ** 2)
    else: M = 0
    return M

## Torque analysis at each point
def InputPointTorque(x,Hlb,Ohi):
    if (x > 1) and (x < 5.25):
        T = Hlb/Ohi
    else:
        T = 0
    return T

def OutputPointTorque(x,Hlb,Oho):
    if x > 2.75 and x < 7:
        T = Hlb/Oho
    else:
        T = 0
    return T

## input shaft point distance from left side
G = 2
H = 2.75
I = 3
J = 3.75
L = 3.75+((3-face)/2)
K = L-.25
M = L+.25
N = L + face
O = N + .25
P = 6.75

## output shaft point distance from left side
Q = 1.25
S = 1.25+((3-face)/2)
R = S-.25
U = S + face
T = U - .25
V = U + .25
W = 3+.75+.5
X = W + .75
Y = X + .25
Z = X + 1

## mean bending moment is 0 because it is a rotating circular shaft
s_m = 0
## alternating torque is 0 because there is a constant speed and constant torque load applied
a_t = 0

## point analysis at each point on input shaft
## d at each point is diameter based on bearing bore and the correlating shoulder
x = G
print("\nPoint G")
scf = 3
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.051
SafetyFactors(scf, s_t, a_m, d)

x = H
print("\nPoint H")
scf = 2
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.051
SafetyFactors(scf, s_t, a_m, d)

x = I
print("\nPoint I")
scf = 4
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.181
SafetyFactors(scf, s_t, a_m, d)

x = J
print("\nPoint J")
scf = 1
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.181
SafetyFactors(scf, s_t, a_m, d)

x = K
print("\nPoint K")
scf =0
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.417
SafetyFactors(scf, s_t, a_m, d)

x = L
print("\nPoint L")
scf = 4
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.417
SafetyFactors(scf, s_t, a_m, d)

x = M
print("\nPoint M")
scf = 3
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.417
SafetyFactors(scf, s_t, a_m, d)

x = N
print("\nPoint N")
scf = 2
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.417
SafetyFactors(scf, s_t, a_m, d)

x = O
print("\nPoint O")
scf = 2
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.417
SafetyFactors(scf, s_t, a_m, d)

x = P
print("\nPoint P")
scf = 1
a_m = InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg)
s_t = InputPointTorque(x, Hlb, Ohi)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.1811
SafetyFactors(scf, s_t, a_m, d)


## point analysis at each point on the output shaft

x = Q
print("\nPoint Q")
scf = 1
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.575
SafetyFactors(scf, s_t, a_m, d)

x = R
print("\nPoint R")
scf = 2
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 1.133
SafetyFactors(scf, s_t, a_m, d)

x = S
print("\nPoint S")
scf = 2
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.944
SafetyFactors(scf, s_t, a_m, d)

x = T
print("\nPoint T")
scf = 3
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.944
SafetyFactors(scf, s_t, a_m, d)

x = U
print("\nPoint U")
scf = 4
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.944
SafetyFactors(scf, s_t, a_m, d)

x = V
print("\nPoint V")
scf = 0
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.944
SafetyFactors(scf, s_t, a_m, d)

x = W
print("\nPoint W")
scf = 1
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.787
SafetyFactors(scf, s_t, a_m, d)

x = X
print("\nPoint X")
scf = 4
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.787
SafetyFactors(scf, s_t, a_m, d)

x = Y
print("\nPoint Y")
scf = 2
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.72
SafetyFactors(scf, s_t, a_m, d)

x = Z
print("\nPoint Z")
scf = 3
a_m = OutputPointMoments(x, Fdy, Wr, Wa, Fdz, Wt, dp)
s_t = OutputPointTorque(x, Hlb, Oho)
print("Moment: ", a_m, "lbf")
print("Torque: ", s_t, "in-lb")
findDiameter(scf, s_t, a_m)
d = 0.72
SafetyFactors(scf, s_t, a_m, d)







