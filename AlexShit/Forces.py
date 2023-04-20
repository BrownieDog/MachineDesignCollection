

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



def InputPointMoments:
    A = 3
    E = 4.5
    if (x-A) > 0:
        if (x-E) > 0:
            My = (x-A) * Fay + (x-E) * Wr - (d/2) * Wa
            Mz = (x-A) * Faz - (X-E) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (x-A) * Fay
            Mz = (x-A) * Faz
            M = math.sqrt(My ** 2 + Mz ** 2)

def OutputPointMoments:
    if (x-C) > 0:
        if (x-F) > 0:
            My = (x-C) * Fcy - (x-F) * Wr - (d/2) * Wa
            Mz = (x-C) * Fcz - (X-F) * Wt
            M = math.sqrt(My ** 2 + Mz ** 2)
        else:
            My = (x-C) * Fcy
            Mz = (x-C) * Fcz
            M = math.sqrt(My ** 2 + Mz ** 2)

def InputPointTorque:
    if x > 1 && x < 4.5:
        T = H/Ohi

def OuputPointTorque:
    if x > 4.5 && x < end
        T = H/Oho


