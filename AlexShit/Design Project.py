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

H = 550*hp*60/(2*math.pi)
Ti = H/Ohi
To = H/Oho
print("Input torque: ", Ti)
print("Output torque: ", To)

phit = (math.atan(tan(phin)/cos(psi)))
print("Transverse pressure angle: ", math.degrees(phit))

## Gear and pinion teeth
Np = ((2*k*cos(psi))/((1+2*m)*((sin(phit))**2)))*(m+((m**2+(1+2*m)*sin(phit)**2)**(0.5)))
Np = int(Np)+1
print("Number of teeth on pinion: ", Np)
Ng = 3 * Np
print("Number of teeth on gear: ", Ng)


## Gear and pinion diameter

dp = 3.175
print("Pinion diameter: ", dp)

dg = 9.526
print("Gear diameter: ", dg)

P = Ng/dg
print("diametral pitch: ", P)
if P > 4:
    P = int(P) + 1
else:
    if P > 3:
        P = 4
    if P > 2.5:
        P = 3
    if P > 2.25:
        P = 2.5
    if P > 2:
        P = 2.25
    else:
        P = 2

V = math.pi*dp*Oho/12
print("Linear pitch velocity: ", V)

## Gear forces

Wt = 33000*hp/V
print("Transmitted force: ", Wt)

Wr = Wt*tan(phit)
print("Radial force: ", Wr)

Wa = Wt*tan(psi)
print("Axial force: ", Wa)

W = Wt/(cos(phin)*cos(psi))
print("Total force magnitude: ", W)

## Bearing forces

## input shaft
print("Fbx: ", Wa)

Fby = (1.5*Wr+(dg/2)*Wa)/3
print("Fby: ", Fby)

Fbz = 1.5*Wt/3
print("Fbz: ", Fbz)

Faz = Wt-Fbz
print("Faz: ", Faz)

Fay = Wr-Fby
print("Fay: ", Fay)

## Output Shaft

Fdy = (1.5*Wr-(dp/2)*Wa)/3
print("Fdy: ", Fdy)

Fdz = 1.5*Wt/3
print("Fdz: ", Fdz)

print("Fcx: ", Wa)

Fcy = Wr-Fdy
print("Fcy: ", Fcy)

print("Fcz: ", Fdz)



## Input shaft forces

A = 1.375
B = 4.75 +.375
C = 3.25

## G

## Y bending moment

Gy = A * Fay + C * Wr - B * Fby
print("G Bending moment y: ", Gy)
Gz = A * Faz - C * Wt +B * Fbz
print("G Bending moment z: ", Gz)
Gm = math.sqrt(Gy ** 2 + Gz ** 2)
print("G combined bending moment: ", Gm)








