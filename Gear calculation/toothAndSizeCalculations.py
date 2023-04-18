import math
#functions and their equations are defined below
def gearratio(In, Out):
    #equation 13-5
    ratio = (Out/In)
    return ratio

def minimumNumberOfTeethOnPinion(m, k, helixAngle, transversePressureAngle):
    #convert angles to radians since python wants it radians for its trig functions
    helixAngleRads = math.radians(helixAngle)
    transversePressureAngleRads = math.radians(transversePressureAngle)
    #equation 13-22 to find minimum number of pinion teeth for specified gear ration
    teethNumber = ((2*k*math.cos(helixAngleRads))/((1+2*m)*(math.sin(transversePressureAngleRads)**2))*(m+math.sqrt(m**2+(1+2*m)*math.sin(transversePressureAngleRads)**2)))
    #code to round up to nearest tooth
    roundedTeethNumber = round(teethNumber)
    if roundedTeethNumber < teethNumber:
        teethNumber = roundedTeethNumber + 1
    else:
        teethNumber = roundedTeethNumber
    return teethNumber

def maximumNumberOfTeethOnGear(Np, k, helixAngle, transversePressureAngle):
    #convert angles to radians since python wants it radians for its trig functions
    helixAngleRads = math.radians(helixAngle)
    transversePressureAngleRads = math.radians(transversePressureAngle)
    #equation 13-23 to find maximum number of gear teeth for specified gear ratio
    teethNumber = ((Np**2)*(math.sin(transversePressureAngleRads)**2)-4*(k**2)*(math.cos(helixAngleRads)**2))/(4*k*math.cos(helixAngleRads)-2*Np*math.sin(transversePressureAngleRads)**2)
    #code to round up to nearest tooth
    roundedTeethNumber = round(teethNumber)
    if roundedTeethNumber < teethNumber:
        teethNumber = roundedTeethNumber + 1
    else:
        teethNumber = roundedTeethNumber
    return teethNumber

def gearSizes(Np, Ng, gearRatio):
    #define the smallest number of teeth on the pinion and gear using gear ratio equation
    pinionSize = Np
    gearSize = round(Np * gearRatio)
    if gearSize > Ng:
        print("the gear size is invalid, try again")
    return pinionSize, gearSize

def checkGearRatio(N1, N2, m, e):
    actualRatio = N2/N1
    percentError = abs((actualRatio-m)/m)*100
    if percentError < e:
        print("The specified Gear ratio is within allowable bounds. "
              "The percent difference between the actual gear ratio and ideal gear ratio is " + str(round(percentError, 4)) + "%")
    else:
        print("BAD RATIO DUMBASS")
        print(percentError)
def pitchlineVelocity(d, n):
    #using equation 13-34 the pitchline velocity is calculated in feet/second
    V = math.pi * d * n / 12
    return V

def transmittedLoad(H, V):
    #equation 13-35 is used to calculate the Wt value
    Wt = 33000 * H / V
    return Wt
def findNormalDiametralPitch(N2, N3, y, c, helixAngle):
    #this equation was derived from the gearbox diagram for the total height
    Pn=((N2/math.cos(math.radians(helixAngle)))+(N3/math.cos(math.radians(helixAngle)))+2)/(y-c)
    Pn = round(Pn)
    return Pn

def cyclesLifetime(rpm):
    #find number of cycles from hour life
    cycles = rpm * 60 * 1000
    return cycles
def main():
    allowableWidth = 15
    clearanceAndWallThickness = 1.5
    hp = 150


    #define rpm then gear ratio
    rpmIn = 6700
    rpmOutIdeal = 20000
    m = gearratio(rpmIn, rpmOutIdeal)
    print("The gear ratio required is " + str(m))


    #define k value and helix and normal pressure angles for equations 13-22 and 13-23
    k=1
    helixAngle = 30
    normalPressureAngle = 20


    #equation 13-19 is used to find the transverse pressure angle
    transversePressureAngle =math.degrees(math.atan((math.tan(math.radians(normalPressureAngle))/math.cos(math.radians(helixAngle)))))
    print("The transverse pressure angle is " +str(transversePressureAngle))


    #minimum number of teeth on the pinion is calculated from equation 13-22
    minPinionTeeth = minimumNumberOfTeethOnPinion(m,k,helixAngle,transversePressureAngle)
    print("The minimum number of teeth allowable on the pinion is " + str(minPinionTeeth))


    #maximum number of teeth on the gear is calculated using equation 13-23
    maxGearTeeth = maximumNumberOfTeethOnGear(minPinionTeeth, k, helixAngle, transversePressureAngle)
    print("The maximum number of teeth allowable on the gear is " + str(maxGearTeeth))


    #gear and pinion sizes are calculated using specified ratio
    minPinionTeeth +=2
    numberOfPinionTeeth, numberOfGearTeeth = gearSizes(minPinionTeeth, maxGearTeeth, m)
    numberOfGearTeeth -= 0


    #allowable tolerance is defined from the gear box specifications
    print("The number of teeth on the pinion and gear are " + str(numberOfPinionTeeth) + " and " + str(numberOfGearTeeth))
    allowablePercentError = 1


    #chosen gear and pinion ratio is checked to ensure it is within tolerance
    print("The allowable percent error between the ideal and actual gear ratios is " + str(allowablePercentError) + "%")
    checkGearRatio(numberOfPinionTeeth, numberOfGearTeeth, m, allowablePercentError)
    actualGearRatio = gearratio(numberOfPinionTeeth, numberOfGearTeeth)


    #calculate the actual rpm output of the gear ratio
    rpmOutActual = rpmIn * (actualGearRatio)
    print("The actual gear ratio is " + str(actualGearRatio) + " and the actual output speed is " + str(rpmOutActual) + " rpm")


    #the normal diametral pitch is then selected from table 13-2. Based off the design requirements, height is not a primary concern, but should be limited
    #gear and pinion diameter affect forces on bearings
    #units in teeth/inch
    normalDiametralPitch = findNormalDiametralPitch(numberOfPinionTeeth, numberOfGearTeeth, allowableWidth, clearanceAndWallThickness, helixAngle)
    #transverse diametral pitch is calculated using equation 13-18
    normalDiametralPitch = 8
    transverseDiametralPitch = normalDiametralPitch * math.cos(math.radians(helixAngle))
    print("Using a Normal pitch of " + str(normalDiametralPitch) +  " the Transverse diametral pitch is " + str(transverseDiametralPitch) + " teeth per inch")


    #gear and pinion diameters are calculated using equation 13-1
    pinionDiameter = numberOfPinionTeeth / transverseDiametralPitch
    gearDiameter = numberOfGearTeeth / transverseDiametralPitch
    print("The pinion diameter is " + str(pinionDiameter) + " inches and the gear diameter is " + str(gearDiameter) + " inches")


    #Pitchline Velocities are calculated
    pitchlinePinion = pitchlineVelocity(gearDiameter, rpmIn)
    pitchlineGear = pitchlineVelocity(gearDiameter, rpmOutActual)
    print("The pitchline velocity is " + str(pitchlinePinion) + " feet/second")



    #transmitted load is calculated
    Wt = transmittedLoad(hp, pitchlinePinion)
    print("The transmitted load is " + str(Wt))
    cyclesPinion = cyclesLifetime(rpmOutActual)
    cyclesGear = cyclesLifetime(rpmIn)


main()