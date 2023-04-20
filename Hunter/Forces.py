from ShaftAnalysis import *

import math
sin = math.sin
cos = math.cos
tan = math.tan

hp = 150
Ohi = 6700
m = 3
Oho = m*Ohi
k = 1 ## full depth teeth
psi = math.radians(30) ## helix angle
phin = math.radians(20) ## Normal pressure angle
phit = (math.atan(tan(phin)/cos(psi)))

face = 1

H = 550*hp*60*12/(2*math.pi)
Ti = H/Ohi
To = H/Oho

dg = 5.629165124598851
dp = 1.8763883748662835
V = math.pi*dp*Oho/12

## Gear forces

Wt = 33000*hp/V
Wr = Wt*tan(phit)
Wa = Wt*tan(psi)
W = Wt/(cos(phin)*cos(psi))

## Reaction forces at bearing

Fbx = Wa
Fby = (1.875*Wr + (dg/2) * Wa)/3.75
Fbz = 1.875*Wt/3.75

Fay = -(Wr-Fby)
Faz = Wt-Fbz

Fdy = (1.875*Wr+(dp/2)*Wa)/3.75
Fdz = 1.875*Wt/3.75

Fcx = Wa
Fcy = (Wr-Fdy)
Fcz = Wt-Fdz

print("gear forces")
print("Wt: ", Wt)
print("Wr: ", Wr)
print("Wa: ", Wa)

print("Input shaft")
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



def InputPointMoments(xi, Fay, Wr, Wa, Faz, Wt, dg):
    A = 3.375
    E = 5.25
    if (xi-A) > 0:
        if (xi-E) > 0:
            My = (xi-A) * Fay + (xi-E) * Wr - (dg/2) * Wa
            Mz = (xi-A) * Faz - (xi-E) * Wt
            a_m = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (xi-A) * Fay
            Mz = (xi-A) * Faz
            a_m = math.sqrt(My ** 2 + Mz ** 2)
    else: a_m = 0
    return a_m

def OutputPointMoments(xo,Fdy,Wr,Wa,Fdz,Wt, dp):
    D = .5+.75+3+.325
    F = 3.25
    if (D-xo) > 0:
        if (F-xo) > 0:
            My = (D-xo) * Fdy - (F-xo) * Wr - (dp/2) * Wa
            Mz = (D-xo) * Fdz - (F-xo) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (D-xo) * Fdy
            Mz = (D-xo) * Fdz
            M = math.sqrt(My ** 2 + Mz ** 2)
    else: M = 0
    return M
def InputPointTorque(xi,H,Ohi):
    if (xi > 1) and (xi < 4.5):
        T = H/Ohi
    else:
        T = 0
    return T

def OutputPointTorque(xo,H,Oho):
    if xo > 2.75 and xo < 7:
        T = H/Oho
    else:
        T = 0
    return T

## input shaft points
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

## output shaft points
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

s_m = 0
a_t = 0

## Point G

xi = J
S_F_goodman = 0
S_F_conservative = 0
d = .2
while S_F_goodman < 1.5 or S_F_conservative < 1.5:
    a_m = InputPointMoments(xi,Fay, Wr, Wa, Faz, Wt, dg)
    s_t = InputPointTorque(xi,H,Ohi)
    print("\nTi: ", Ti)
    print("a_m", a_m)
    print("diameter: ", d)
    S_F_goodman, S_F_conservative = main(d,a_m, a_t, s_m, s_t, geometry = "round")
    d += .1




