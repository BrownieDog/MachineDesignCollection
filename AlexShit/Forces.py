

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

dg = 9.526
dp = 3.175
V = math.pi*dp*Oho/12
Wt = 33000*hp/V

Wr = Wt*tan(phit)
Wa = Wt*tan(psi)
W = Wt/(cos(phin)*cos(psi))



Fbx = Wa
Fby = (1.875*Wr + (dg/2) * Wa)/3.75
Fbz = 1.875*Wt/3.75

Fay = -(Wr-Fby)
Faz = Wt-Fbz

Fdy = (1.875*Wr+(dp/2)*Wa)/3.75
Fdz = 1.875*Wt/3.75

Fcx = Wa
Fcy = -(Wr-Fdy)
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


