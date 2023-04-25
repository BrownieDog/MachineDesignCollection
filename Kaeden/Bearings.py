import math
import sympy as sym

###############################################################################
# Dimensionless life, x_D
###############################################################################
def life (bearing):
    
# Shaft speed
    # If the bearing is on the input shaft (A and B)
    if (bearing == 'A') or (bearing == 'B'):
        # The shaft speed is 6700 rpm
        n = 6700
    # Else the bearing is on the output shaft (C and D)
    elif (bearing == 'C') or (bearing == 'D'):
        # The shaft speed is 20100 rpm
        n = 20100
    # Terminate if typo
    else:
        exit()

# Convert lifetime from hours to cycles
    # Lifetime in hours
    L = 1000      
    # Lifetime in cycles: multiply life by (60 min / 1 hour) and (1 rev / 1 min)
    L_D = L * 60 * n

# Minimum life rating (cycles)
    L10 = 1000000       

# Calculate dimensionless life
    x_D = L_D / L10

    # Return the dimensionless life
    return x_D
###############################################################################


###############################################################################
# Axial and radial forces
###############################################################################
def forces(bearing):
    
# Reaction forces (magnitudes) from forces.py
    # Reaction forces on bearing A
    F_A = [0,           111.8928,       250.6616]
    # Reaction forces on bearing B
    F_B = [289.4391,    322.5873,       250.6616]
    # Reaction forces on bearing C
    F_C = [289.4391,    177.7604,       250.6616]
    # Reaction forces on bearing D
    F_D = [0,           32.9341,        250.6616]

    # Variables to access i, j, and k force components
    i = 0       # i component
    j = 1       # j component
    k = 2       # k component

# Calculate axial and radial loads
# Equation used is vector addition:
# F_a = F_x
# F_r = sqrt(F_y^2 + F_z^2)
    # Bearing A
    if bearing == 'A':
        # Axial force on A
        F_a = F_A[i]
        # Radial force on A
        F_r = math.sqrt((F_A[j])**2 + (F_A[k])**2)
    # Bearing B
    elif bearing == 'B':
        # Axial force on B
        F_a = F_B[i]
        # Radial force on B
        F_r = math.sqrt((F_B[j])**2 + (F_B[k])**2)
    # Bearing C
    elif bearing == 'C':
        # Axial force on C
        F_a = F_C[i]
        # Radial force on C
        F_r = math.sqrt((F_C[j])**2 + (F_C[k])**2)
    # Bearing D
    elif bearing == 'D':
        # Axial force on D
        F_a = F_D[i]
        # Radial force on D
        F_r = math.sqrt((F_D[j])**2 + (F_D[k])**2)

    # Return axial and radial forces of the bearing in question
    return F_a, F_r
###############################################################################


###############################################################################
# Table 11-2 Interpolation:
# This section of code finds the values of e and Y2 by interpolating the values
# found in Table 11-2 based on the value F_a / C0.
###############################################################################
def interpolation (Fa_c0):
    
    # e is the abscissa and Y2 is the slope used in Eq. 11-9

    if Fa_c0 <= 0.014:
        e = 0.19
        Y2 = 2.30
    
    elif (Fa_c0 > 0.014) and (Fa_c0 <= 0.021):
        e = ((Fa_c0 - 0.014)/(0.021 - 0.014) * (0.21 - 0.19)) + 0.19
        Y2 = ((Fa_c0 - 0.014)/(0.021 - 0.014) * (2.15 - 2.30)) + 2.30

    elif (Fa_c0 > 0.021) and (Fa_c0 <= 0.028):
        e = ((Fa_c0 - 0.021)/(0.028 - 0.021) * (0.22 - 0.21)) + 0.21
        Y2 = ((Fa_c0 - 0.021)/(0.028 - 0.021) * (1.99 - 2.15)) + 2.15

    elif (Fa_c0 > 0.028) and (Fa_c0 <= 0.042):
        e = ((Fa_c0 - 0.028)/(0.042 - 0.028) * (0.24 - 0.22)) + 0.22
        Y2 = ((Fa_c0 - 0.028)/(0.042 - 0.028) * (1.85 - 1.99)) + 1.99

    elif (Fa_c0 > 0.042) and (Fa_c0 <= 0.056):
        e = ((Fa_c0 - 0.042)/(0.056 - 0.042) * (0.26 - 0.24)) + 0.24
        Y2 = ((Fa_c0 - 0.042)/(0.056 - 0.042) * (1.71 - 1.85)) + 1.85

    elif (Fa_c0 > 0.056) and (Fa_c0 <= 0.070):
        e = ((Fa_c0 - 0.056)/(0.070 - 0.056) * (0.27 - 0.26)) + 0.26
        Y2 = ((Fa_c0 - 0.056)/(0.070 - 0.056) * (1.63 - 1.71)) + 1.71

    elif (Fa_c0 > 0.070) and (Fa_c0 <= 0.084):
        e = ((Fa_c0 - 0.070)/(0.084 - 0.070) * (0.28 - 0.27)) + 0.27
        Y2 = ((Fa_c0 - 0.070)/(0.084 - 0.070) * (1.55 - 1.63)) + 1.63

    elif (Fa_c0 > 0.084) and (Fa_c0 <= 0.110):
        e = ((Fa_c0 - 0.084)/(0.110 - 0.084) * (0.30 - 0.28)) + 0.28
        Y2 = ((Fa_c0 - 0.084)/(0.110 - 0.084) * (1.45 - 1.55)) + 1.55

    elif (Fa_c0 > 0.11) and (Fa_c0 <= 0.17):
        e = ((Fa_c0 - 0.11)/(0.17 - 0.11) * (0.34 - 0.30)) + 0.30
        Y2 = ((Fa_c0 - 0.11)/(0.17 - 0.11) * (1.31 - 1.45)) + 1.45

    elif (Fa_c0 > 0.17) and (Fa_c0 <= 0.28):
        e = ((Fa_c0 - 0.17)/(0.28 - 0.17) * (0.38 - 0.34)) + 0.34
        Y2 = ((Fa_c0 - 0.17)/(0.28 - 0.17) * (1.15 - 1.31)) + 1.31

    elif (Fa_c0 > 0.28) and (Fa_c0 <= 0.42):
        e = ((Fa_c0 - 0.28)/(0.42 - 0.28) * (0.42 - 0.38)) + 0.38
        Y2 = ((Fa_c0 - 0.28)/(0.42 - 0.28) * (1.04 - 1.15)) + 1.15

    elif (Fa_c0 > 0.42) and (Fa_c0 <= 0.56):
        e = ((Fa_c0 - 0.42)/(0.56 - 0.42) * (0.44 - 0.42)) + 0.42
        Y2 = ((Fa_c0 - 0.42)/(0.56 - 0.42) * (1.00 - 1.04)) + 1.04

    else:
        exit()

    # Return the interpolated values of e and Y2
    return e, Y2
###############################################################################


###############################################################################
# Select bearing
###############################################################################
def select_bearing(x_D, F_a, F_r):

    # Variables
    R_D = 0.975             # Reliability
    a_f = 1.4               # Application factor
    V = 1                   # Rotation factor: Inner ring rotates
    a = 3                   # Value for ball bearing       

    # Table 11-6: Weibull Parameters (Manufacturer 2 Values)      
    x0 = 0.02               # Guaranteed value of x
    theta = 4.459           # Characteristic parameter
    b = 1.483               # Shape parameter

    # Table 11-1 values (to be used later in Eq. 11-9)
    X1 = 1                  # Ordinate intercept 1
    Y1 = 0                  # Slope 1
    X2 = 0.56               # Ordinate intercept 2

    # If only a radial load exists
    if F_a == 0:

        # Eq. 11-6: Calculate C10 value using the radial load
        C10_real = a_f * F_r * (x_D/(x0+((theta-x0)*(math.log(1/R_D))**(1/b))))**(1/a)

        # Convert from lbf to kN (because the manufacturer specs are in kN)
        C10_real = (C10_real * 4.448) / 1000

        # Print load requirements and dynamic load rating
        print("\nAxial Force = %.4f lbf" % F_a)
        print("Radial Force = %.4f lbf" % F_r)
        print("Dynamic load rating needs to be greater than or equal to %.4f kN\n" % C10_real)

        # Iterate until a bearing meets specifications
        while True:
            # Get manufacturer static and dynmaic ratings
            C0 = float(input("C0 (static rating) from manufacturer (kN): "))
            C10 = float(input("C10 (dynamic rating) from manufacturer (kN): "))

            # Check that the calculated dynamic rating is less than the manufacturer's rating
            if C10_real >= C10:
                # If the bearing does not meet the required dynamic rating, choose a new bearing
                print("\nSelect a larger bearing")
            else:
                # Else bearing works- do not iterate again
                break

    # Else, the bearing experiences a combined load (radial and axial)    
    else:
        # Choose an initial guess of the equivalent radial load
        # The simple addition of the axial and radial forces provides a descent first guess
        F_e = F_a + F_r

        # Eq. 11-6: Calculate the C10 value based on the initial equivalent load guess
        C10_guess = a_f * F_e * (x_D/(x0+((theta-x0)*(math.log(1/R_D))**(1/b))))**(1/a)

        # Convert from lbf to kN (because the manufacturer specs are in kN)
        C10_guess = (C10_guess * 4.448) / 1000
        
        # Print load requirements and initial dynamic load rating guess
        print("\nAxial Force = %.4f lbf" % F_a)
        print("Radial Force = %.4f lbf" % F_r)
        print("Initial C10 guess: %.4f kN\n" % C10_guess)

        # Iterate until a bearing meets specifications
        while True:
            # C0 and C10 values from manufacturer (based on previous guess)
            C0 = float(input("C0 (static rating) from manufacturer (kN): "))
            C10 = float(input("C10 (dynamic rating) from manufacturer (kN): "))

            # Convert from kN to lbf (to perfrom calculations)
            C0_lbf = (C0 * 1000) / 4.48

            # Compute axial force divided by the static load rating (F_axial / C0)
            Fa_c0 = F_a / C0_lbf

            # Interpolate e and Y2 using Table 11-1 based on F_a/C0 value
            e, Y2 = interpolation(Fa_c0)

            # Compute the value of the axial force divided by the radial force and rotation factor
            # That is: F_axial / (V * F_radial)
            value = F_a / (V * F_r)

            # If the value is less than or greater to e (found from Table 11-1)
            if value <= e:
                # Use X1 and Y2 (later in Eq. 11-9)
                X = X1
                Y = Y1

            # Else
            else:
                # Use X2 and Y2 (later in Eq. 11-9)
                X = X2
                Y = Y2

            # Eq. 11-9: Calculate equivalent radial load
            F_e = (X * V * F_r) + (Y * F_a)

            # Eq. 11-6: Calculate C10 (dynamic load rating)
            C10_real = a_f * F_e * (x_D/(x0+((theta-x0)*(math.log(1/R_D))**(1/b))))**(1/a)

            # Convert to kN (because the manufacturer specs are in kN)
            C10_real = (C10_real * 4.48) / 1000
        
            # Check the actual C10 rating against the manufacturer C10 rating
            # If the real C10 is less than or equal to the manufacturer C10
            if C10_real <= C10:
                # The bearing will work- stop iterating
                break
            # If the real C10 is greater than the manufacturer C10
            else:
                # The bearing will not work- iterate again
                print("\nTry the next largest bearing.")

    # Return the manufacturer static and dynamic ratings and the actual dynamic rating
    return C0, C10, C10_real
###############################################################################


###############################################################################
# Main
###############################################################################
while True:
    # Bearing to be evaluated
    bearing = input("Which bearing is being evaluated (A, B, C, D): ")

    # Dimenionless (design) life of the bearing
    x_D = life(bearing)

    # Axial and radial forces on the bearing
    F_a, F_r = forces(bearing)

    # Select the bearing
    C0, C10, C10_real = select_bearing(x_D, F_a, F_r)

    # Print results
    print("\n\n\nBearing " + bearing)
    print("Axial Force = %.4f lbf" % F_a)
    print("Radial Force = %.4f lbf" % F_r)
    print("Bearing Type: Ball")
    print("Manufacturer C0: %.4f kN" % C0)
    print("Manufacturer C10: %.4f kN" % C10)
    print("Actual C10: %.4f kN\n\n\n" % C10_real)

    # Analyze another bearing?
    repeat = input("Select another bearing? (yes or no) ")
    if repeat == "no":
        break

# End
###############################################################################
