

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

H = 550*hp*60/(2*math.pi)
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



def InputPointMoments(x, Fay, Wr, Wa, Faz, Wt, dg):
    A = 3.325
    E = 5.25
    if (x-A) > 0:
        if (x-E) > 0:
            My = (x-A) * Fay + (x-E) * Wr - (dg/2) * Wa
            Mz = (x-A) * Faz - (x-E) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (x-A) * Fay
            Mz = (x-A) * Faz
            M = math.sqrt(My ** 2 + Mz ** 2)
    return M

def OutputPointMoments(x,Fdy,Wr,Wa,Fdz,Wt, dp):
    D = .5+.75+3+.325
    F = 3.25
    if (D-x) > 0:
        if (F-x) > 0:
            My = (D-x) * Fdy - (F-x) * Wr - (dp/2) * Wa
            Mz = (D-x) * Fdz - (F-x) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (D-x) * Fdy
            Mz = (D-x) * Fdz
            M = math.sqrt(My ** 2 + Mz ** 2)
    return M
def InputPointTorque(x,H,Ohi):
    if (x > 1) and (x < 4.5):
        T = H/Ohi
    else:
        T = 0
    return T

def OuputPointTorque(x,H,Oho):
    if x > 2.75 and x <
        T = H/Oho
    else:
        T = 0
    return T


