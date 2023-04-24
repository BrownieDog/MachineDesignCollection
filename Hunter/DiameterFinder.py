# Scott Michelsen
#this file finds the minimum diameter using a sefety factor of 1.5 and checking the conservative first cycle yield and Goodman Criteria
#The user inputs the fully reversed bending moment, the steady torque, Kf, Kfs,
# the stress concentration factor, and the relevant diameter. the shaft was hot-rolled, than machined to size

import math

#stress concentrator is 1-4, 1 = Sharp radius, 2 = Wide Radius, 3 = Keyway, 4 = retaining ring
#S_T = Steady Torque
#A_M = Alternating Bending Moment
def findDiameter(Stress_Concentrator, S_T, A_M ):
    #print("\nCalculated Diameter")
    #this section of code request all the necessary data to perform the calculations
    #print("Please enter all requested values in English/Imperial Units")
    FullyReversedBendingMoment = A_M
    SteadyTorque = S_T
    RelevantShaftDiameter = .11

    #this section of code allows the user to pick what the stress concentrator is
    # and then assigns values to Kt, Kts, and R for the Kf and Kfs calculations
    loopvar = True
    while loopvar == True:

        if Stress_Concentrator == 0:
            Kt = 1
            Kts = 1
            StressConcentrationRadius = 0
            loopvar = False
        if Stress_Concentrator ==1:
            Kt=2.7
            Kts=2.2
            StressConcentrationRadius = RelevantShaftDiameter * .02
            loopvar = False
        elif Stress_Concentrator ==2:
            Kt = 1.7
            Kts = 1.5
            StressConcentrationRadius = RelevantShaftDiameter * .1
            loopvar = False
        elif Stress_Concentrator ==3:
            Kt = 2.14
            Kts = 3.0
            StressConcentrationRadius = RelevantShaftDiameter * .02
            loopvar = False
        elif Stress_Concentrator ==4:
            Kt = 5.0
            Kts = 3.0
            StressConcentrationRadius = 0.01
            loopvar = False
        else:
            print("Please enter a number between 1 and 4")


    #variable for loop condition
    SafetyFactorMet = False
    GoodmanSafetyFactorMet = False
    ConservativeSafetyFactorMet = False
    SafetyFactorGoodman = 0
    SafetyFactorConservative = 0
    #while loop to iterate until the diameter meets both safety factor criteria
    while not SafetyFactorMet:
        #the stress concentration radius must be calculated for each iteration, which this sort does. menuSelect is defined by user input before the loop
        #this ensures that the stress concentration factors are updated each iteration
        if Stress_Concentrator == 0:
            Kt = 1
            Kts = 1
            StressConcentrationRadius = 0
            loopvar = False
        elif Stress_Concentrator ==1:
            Kt=2.7
            Kts=2.2
            StressConcentrationRadius = RelevantShaftDiameter * .02
            loopvar = False
        elif Stress_Concentrator ==2:
            Kt = 1.7
            Kts = 1.5
            StressConcentrationRadius = RelevantShaftDiameter * .1
            loopvar = False
        elif Stress_Concentrator ==3:
            Kt = 2.14
            Kts = 3.0
            StressConcentrationRadius = RelevantShaftDiameter * .02
            loopvar = False
        elif Stress_Concentrator ==4:
            Kt = 5.0
            Kts = 3.0
            StressConcentrationRadius = 0.01
            loopvar = False
        else:
            print("Please enter a number between 1 and 4")

        # Material properties given in kpsi from table A-20 for HR-1030 steel
        UltimateTensileStrength = 68
        YieldStrength = 37.5
        # Find the endurance limit
        # Se' is calculated by equation 6-10
        if UltimateTensileStrength <= 200:
            EnduranceStrengthPrime = 0.5 * UltimateTensileStrength
        else:
            EnduranceStrengthPrime = 100
        # define a and b value from table 6.2 for machine finish, then find Ka

        a = 2.00
        b = -0.217

        ka = a * UltimateTensileStrength ** b

        # calculate Kb using equation 6-19. the if else sorts for the different cases based off diamter size. all given diameters are inches
        if RelevantShaftDiameter >= 0.11 and RelevantShaftDiameter <= 2:
            kb = (RelevantShaftDiameter / 0.3) ** -0.107
        elif RelevantShaftDiameter > 2 and RelevantShaftDiameter <= 10:
            kb = (RelevantShaftDiameter * 0.91) ** -0.157
        else:
            print("The given shaft diameter is invalid")
            kb = 1
        #set other k values for finding the endurance limit. These are based off the problem statement.
        kc = 1
        kd = 1
        ke = .897
        kf = 1
        #calculate the Endurance Limit in Ksi using equation 6-17
        EnduranceLimit = EnduranceStrengthPrime * ka * kb * kc * kd * ke * kf

        if Kt == 1:
            Kf = 1
            Kfs = 1
        else:

            # find Kf, which is the bending stress concentration factor
            # find square root of a for bending using equation 6-35
            SquareRootofaforBending = 0.246 - 3.08 * (10 ** -3) * UltimateTensileStrength + 1.51 * (10 ** -5) * (
                        UltimateTensileStrength ** 2) - 2.67 * (10 ** -8) * (UltimateTensileStrength ** 3)
            # calculate Kf using equation 6-34
            Kf = 1 + (Kt - 1) / (1 + (SquareRootofaforBending / math.sqrt(StressConcentrationRadius)))

            # find Kfs, which is the torsional stress concentration factor
            # find square root of a for torsion using equation 6-36
            SquareRootofaforTorsion = 0.190 - 2.51 * (10 ** -3) * UltimateTensileStrength + 1.35 * (10 ** -5) * (
                        UltimateTensileStrength ** 2) - 2.67 * (10 ** -8) * (UltimateTensileStrength ** 3)
            # calculate Kfs using equation 6-34
            Kfs = 1 + (Kts - 1) / (1 + (SquareRootofaforTorsion / math.sqrt(StressConcentrationRadius)))


        # calculate out the Goodman Safety Factor using equation 7-7, the right hand side is calculated first, then 1 is divided by that side to solve for n
        # do not forget to convert the strengths from ksi to psi since the given forces are in pounds and the strengths are in Ksi
        onedividedbyn = (16 / (math.pi * (RelevantShaftDiameter ** 3))) * (
                    (1 / (EnduranceLimit * 10 ** 3)) * math.sqrt(4 * ((Kf * FullyReversedBendingMoment) ** 2)) + (
                        1 / (UltimateTensileStrength * 10 ** 3)) * math.sqrt(3 * ((Kfs * SteadyTorque) ** 2)))
        n = 1 / onedividedbyn
        SafetyFactorGoodman = n

        #calculate the Conservative First Cycle Yield Safety Factor
        #calculate Sigma a' and Sigma m' using equations 7-4 and 7-5

        SigmaAPrime = math.sqrt(((32 * Kf * FullyReversedBendingMoment) / (math.pi * RelevantShaftDiameter ** 3)) ** 2)
        SigmaMPrime = math.sqrt(3 * (((16 * Kfs * SteadyTorque) / (math.pi * RelevantShaftDiameter ** 3)) ** 2))

        # calculate the conservative Von Mises by adding SigmaA' and SigmaM'
        ConservativeVonMises = SigmaAPrime + SigmaMPrime

        #calculate the safety factor using equation 6-43 mod, Yield Strength is multiplied by 10^3 to account for it being in Ksi and the conservative von mises being in psi
        n = (YieldStrength * (10 ** 3)) / ConservativeVonMises
        SafetyFactorConservative = n

        #increase shaft diameter for the next iteration. this will increase once more before printing the final answer
        #which make the answer slightly more conservative
        RelevantShaftDiameter += 0.00001

        #if statements to check if the diameter has a safety factor of 1.5 for both Goodman and Conservative yielding
        if SafetyFactorGoodman > 1.5:
            GoodmanSafetyFactorMet = True
        if SafetyFactorConservative > 1.5:
            ConservativeSafetyFactorMet = True

        #check is both safety factors are at least 1.5 and if so, changes the SafetyFactorMet Variable to true to end the loop
        if GoodmanSafetyFactorMet and ConservativeSafetyFactorMet:
            SafetyFactorMet = True
    #print the safety factors and the diameter
    # print("The Goodman Safety Factor is " + str(SafetyFactorGoodman))
    # print("The Conservative Yield Safety Factor is " + str(SafetyFactorConservative))
    # print("The minimum diameter is " + str(RelevantShaftDiameter))
    # print("Kt value: ", Kt)
    # print("Kts value: ", Kts)
    # print("Kf value: ", Kf)
    # print("Kfs value: ", Kfs)
    return RelevantShaftDiameter

def SafetyFactors(Stress_Concentrator, S_T, A_M, RelevantShaftDiameter ):
    print("\nSafety Factors based on input diameter")
    FullyReversedBendingMoment = A_M
    SteadyTorque = S_T

#the stress concentration radius must be calculated for each iteration, which this sort does. menuSelect is defined by user input before the loop
        #this ensures that the stress concentration factors are updated each iteration
    if Stress_Concentrator == 0:
        Kt = 1
        Kts = 1
        StressConcentrationRadius = 0
        loopvar = False
    elif Stress_Concentrator ==1:
        Kt=2.7
        Kts=2.2
        StressConcentrationRadius = RelevantShaftDiameter * .02
        loopvar = False
    elif Stress_Concentrator ==2:
        Kt = 1.7
        Kts = 1.5
        StressConcentrationRadius = RelevantShaftDiameter * .1
        loopvar = False
    elif Stress_Concentrator ==3:
        Kt = 2.14
        Kts = 3.0
        StressConcentrationRadius = RelevantShaftDiameter * .02
        loopvar = False
    elif Stress_Concentrator ==4:
        Kt = 5.0
        Kts = 3.0
        StressConcentrationRadius = 0.01
        loopvar = False
    else:
        print("Please enter a number between 1 and 4")

    # Material properties given in kpsi from table A-20 for HR-1030 steel
    UltimateTensileStrength = 68
    YieldStrength = 37.5
    # Find the endurance limit
    # Se' is calculated by equation 6-10
    if UltimateTensileStrength <= 200:
        EnduranceStrengthPrime = 0.5 * UltimateTensileStrength
    else:
        EnduranceStrengthPrime = 100
    # define a and b value from table 6.2 for machine finish, then find Ka

    a = 2.00
    b = -0.217

    ka = a * UltimateTensileStrength ** b

    # calculate Kb using equation 6-19. the if else sorts for the different cases based off diamter size. all given diameters are inches
    if RelevantShaftDiameter >= 0.11 and RelevantShaftDiameter <= 2:
        kb = (RelevantShaftDiameter / 0.3) ** -0.107
    elif RelevantShaftDiameter > 2 and RelevantShaftDiameter <= 10:
        kb = (RelevantShaftDiameter * 0.91) ** -0.157
    else:
        print("The given shaft diameter is invalid")
        kb = 1
    #set other k values for finding the endurance limit. These are based off the problem statement.
    kc = 1
    kd = 1
    ke = .897
    kf = 1
    #calculate the Endurance Limit in Ksi using equation 6-17
    EnduranceLimit = EnduranceStrengthPrime * ka * kb * kc * kd * ke * kf

    if Kt == 1:
        Kf = 1
        Kfs = 1
    else:

        # find Kf, which is the bending stress concentration factor
        # find square root of a for bending using equation 6-35
        SquareRootofaforBending = 0.246 - 3.08 * (10 ** -3) * UltimateTensileStrength + 1.51 * (10 ** -5) * (
                    UltimateTensileStrength ** 2) - 2.67 * (10 ** -8) * (UltimateTensileStrength ** 3)
        # calculate Kf using equation 6-34
        Kf = 1 + (Kt - 1) / (1 + (SquareRootofaforBending / math.sqrt(StressConcentrationRadius)))

        # find Kfs, which is the torsional stress concentration factor
        # find square root of a for torsion using equation 6-36
        SquareRootofaforTorsion = 0.190 - 2.51 * (10 ** -3) * UltimateTensileStrength + 1.35 * (10 ** -5) * (
                    UltimateTensileStrength ** 2) - 2.67 * (10 ** -8) * (UltimateTensileStrength ** 3)
        # calculate Kfs using equation 6-34
        Kfs = 1 + (Kts - 1) / (1 + (SquareRootofaforTorsion / math.sqrt(StressConcentrationRadius)))


    # calculate out the Goodman Safety Factor using equation 7-7, the right hand side is calculated first, then 1 is divided by that side to solve for n
    # do not forget to convert the strengths from ksi to psi since the given forces are in pounds and the strengths are in Ksi
    onedividedbyn = (16 / (math.pi * (RelevantShaftDiameter ** 3))) * (
                (1 / (EnduranceLimit * 10 ** 3)) * math.sqrt(4 * ((Kf * FullyReversedBendingMoment) ** 2)) + (
                    1 / (UltimateTensileStrength * 10 ** 3)) * math.sqrt(3 * ((Kfs * SteadyTorque) ** 2)))
    n = 1 / onedividedbyn
    SafetyFactorGoodman = n

    #calculate the Conservative First Cycle Yield Safety Factor
    #calculate Sigma a' and Sigma m' using equations 7-4 and 7-5

    SigmaAPrime = math.sqrt(((32 * Kf * FullyReversedBendingMoment) / (math.pi * RelevantShaftDiameter ** 3)) ** 2)
    SigmaMPrime = math.sqrt(3 * (((16 * Kfs * SteadyTorque) / (math.pi * RelevantShaftDiameter ** 3)) ** 2))

    # calculate the conservative Von Mises by adding SigmaA' and SigmaM'
    ConservativeVonMises = SigmaAPrime + SigmaMPrime

    #calculate the safety factor using equation 6-43 mod, Yield Strength is multiplied by 10^3 to account for it being in Ksi and the conservative von mises being in psi
    n = (YieldStrength * (10 ** 3)) / ConservativeVonMises
    SafetyFactorConservative = n

    print("Diameter: ", RelevantShaftDiameter)
    print("Conservative safety factor: ", SafetyFactorConservative)
    print("Goodman safety factor: ", SafetyFactorGoodman)
