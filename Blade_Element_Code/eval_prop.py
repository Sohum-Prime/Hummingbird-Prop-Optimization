"""
Evaluates prop using Blade Element Method for static thrust

:AUTHOR: William Skelly
:VERSION: 1.0
"""

import numpy as np
import math
import find_zero as fz
import physics as phys
import prop_class
import matplotlib.pyplot as plt

def thrust_momentum(vel_ind, radius, dr):
    """
    Calculates thrust contribution from element according to momentum theory

    args:
        vel_ind:    induced velocity
        radius:     radius of element from hub
        dr:         thickness of element in radial dimension
    returns: thrust contribution from element
    """

    delta_T = 4 * math.pi * radius * phys.air_density * pow(vel_ind,2) * dr
    return delta_T

def torque_momentum(vel_ind, swirl_fact, angle_vel, radius, dr):
    """
    Calculates torque contribution from element according to momentum theory

    args:
        vel_ind:    induced velocity
        swirl_fact: swirl factor
        angle_vel:  angular velocity (rad/s)
        radius:     radius of element from hub
        dr:         thickness of element in radial dimension
    returns: torque contribution from element
    """

    delta_Q = 4 * math.pi * pow(radius,3) * phys.air_density *\
              vel_ind * swirl_fact * angle_vel * dr
    return delta_Q

def thrust_element(vel_tot, radius, dr, chord, num_blades, cl, cd, phi):
    """
    Calculates thrust contribution from element according to blade element
    theory

    args:
        vel_tot:    vector sum of induced velocity and radial velocity
        radius:     radius of element from hub
        dr:         thickness of element in radial dimension
        chord:      chord of one blade
        num_blades: number of blades
        cl:         coefficient of lift
        cd:         coefficient of drag
        phi:        geometric pitch angle minus angle of attack
    returns: thrust contribution from element
    """

    coeff_term = cl * math.cos(phi) - cd * math.sin(phi)
    const_term = 0.5 * phys.air_density * pow(vel_tot,2) * chord * num_blades
    delta_T = const_term * coeff_term * dr
    return delta_T

def torque_element(vel_tot, radius, dr, chord, num_blades, cl, cd, phi):
    """
    Calculates torque contribution from element according to blade element
    theory

    args:
        vel_tot:    vector sum of induced velocity and radial velocity
        radius:     radius of element from hub
        dr:         thickness of element in radial dimension
        chord:      chord of one blade
        num_blades: number of blades
        cl:         coefficient of lift
        cd:         coefficient of drag
        phi:        geometric pitch angle minus angle of attack
    returns: torque contribution from element
    """

    #Note addition of radius term in delta_Q, there is a typo in the
       #aerodynamics4students page where they forget the extra r term
    coeff_term = cd * math.cos(phi) + cl * math.sin(phi)
    const_term = 0.5 * phys.air_density * pow(vel_tot,2) * chord * num_blades
    delta_Q = const_term * coeff_term * radius * dr
    return delta_Q

def calc_vel_tan(radius, angle_vel, swirl_fact):
    """
    Calculates the velocity component tangential to the prop arc

    args:
        radius: radius at station (m)
        angle_vel: angular velocity (m/s)
        swirl_fact: swirl factor
    returns: the tangential velocity component (m/s)
    """

    v_tan_ideal = angle_vel * radius
    v_tan = v_tan_ideal * (1 - swirl_fact)
    return v_tan

def calc_vel_tot(vel_ind, vel_tan):
    """
    Calculates resutant velocity at station

    args:
        vel_ind: induced velocity (m/s)
        vel_tan: tangential velocity
    returns: total velocity at the station (m/s)
    """

    vel_tot = math.sqrt(pow(vel_ind,2) + pow(vel_tan,2))
    return vel_tot

def calc_phi(vel_ind, vel_tan):
    """
    Calculates the inflow angle, phi

    args:
        vel_ind: induced velocity (m/s)
        vel_tan: tangential velocity (m/s)
    returns: inflow angle in radians
    """

    return math.atan(vel_ind / vel_tan)

def thrust_error(vel_ind, swirl_fact, angle_vel, radius, dr, chord, num_blades,
                 cl_func, cd_func, pitch):
    """
    Calculates the error between momentum theory and element theory in thrust

    args:
        vel_ind: induced velocity (m/s)
        swirl_fact: swirl factor
        angle_vel: angular velocity (rad/s)
        radius: radius at station (m)
        dr: width of differential element in the r direction (m)
        chord: chord of blade element (m)
        num_blades: number of blades
        cl_func: function of AoA for coefficient of lift of foil
        cd_func: function of AoA for coefficient of drag of foil
        pitch: pitch of prop (m)
    returns: error in thrust prediction between the two methods
    """

    #calculate velocity components and inflo angle
    vel_tan = calc_vel_tan(radius, angle_vel, swirl_fact) #m/s
    vel_tot = calc_vel_tot(vel_ind, vel_tan) #m/s
    phi     = calc_phi(vel_ind, vel_tan) #rad
    theta   = np.arctan(pitch / (2 * radius * np.pi)) #rad
    AoA     = theta - phi #rad
    cl      = cl_func(AoA)
    cd      = cd_func(AoA)

    #calculate thrust by each method
    thrust_from_element  = thrust_element(vel_tot, radius, dr, chord,
                                          num_blades, cl, cd, phi)
    thrust_from_momentum = thrust_momentum(vel_ind, radius, dr)

    #calulate error and return
    error = thrust_from_momentum - thrust_from_element
    return error

def torque_error(vel_ind, swirl_fact, angle_vel, radius, dr, chord, num_blades,
                 cl_func, cd_func, pitch):
    """
    Calculates the error between momentum theory and element theory for torque

    args:
        vel_ind: induced velocity (m/s)
        swirl_fact: swirl factor
        angle_vel: angular velocity (rad/s)
        radius: radius at station (m)
        dr: width of differential element in the r direction (m)
        chord: chord of blade element (m)
        num_blades: number of blades
        cl_func: function of AoA for coefficient of lift of foil
        cd_func: function of AoA for coefficient of drag of foil
        pitch: pitch of prop (m)
    returns: error in torque prediction between the two methods
    """

    #calculate velocity components and inflow angle
    vel_tan = calc_vel_tan(radius, angle_vel, swirl_fact) #m/s
    vel_tot = calc_vel_tot(vel_ind, vel_tan) #m/s
    phi     = calc_phi(vel_ind, vel_tan) #rad
    theta   = np.arctan(pitch / (2 * radius * np.pi)) #rad
    AoA     = theta - phi #rad
    cl      = cl_func(AoA)
    cd      = cd_func(AoA)

    #calculate torque by each method
    torque_from_element  = torque_element(vel_tot, radius, dr, chord,
                                          num_blades, cl, cd, phi)
    torque_from_momentum = torque_momentum(vel_ind, swirl_fact, angle_vel,
                                           radius, dr)
    #print("Torque element: " + str(torque_from_element))
    #print("Torque from momentum: " + str(torque_from_momentum))

    #calulate error and return
    error = torque_from_momentum - torque_from_element
    return error

def create_thrust_error_func(angle_vel, radius, dr, chord, num_blades,
                             cl_func, cd_func, pitch):
    """
    Returns a function of the induced velocity and swirl factor for the error
    in the thrust with the angular velocity, radius, dr, chord, number of blades
    coefficient of lift, and coefficient of drag set as constant from the
    inputs to create_thrust_error_func

    args:
        angle_vel: angular velocity (rad/s)
        radius: radius at station (m)
        dr: width of differential element in the r direction (m)
        chord: chord of blade element (m)
        num_blades: number of blades
        cl_func: function of AoA for coefficient of lift of foil
        cd_func: function of AoA for coefficient of drag of foil
        pitch: pitch of prop (m)
    returns: a function of induced velocity and swirl factor for error in thrust
    """

    #Idea from Devlin Ih
    def error_of_2_vars(vel_ind, swirl_fact):
        return thrust_error(vel_ind, swirl_fact, angle_vel, radius, dr, chord,
                            num_blades, cl_func, cd_func, pitch)
    return error_of_2_vars

def create_torque_error_func(angle_vel, radius, dr, chord, num_blades,
                             cl_func, cd_func, pitch):
    """
    Returns a function of the induced velocity and swirl factor for the error
    in the torque with the angular velocity, radius, dr, chord, number of blades
    coefficient of lift, and coefficient of drag set as constant from the
    inputs to create_torque_error_func

    args:
        angle_vel: angular velocity (rad/s)
        radius: radius at station (m)
        dr: width of differential element in the r direction (m)
        chord: chord of blade element (m)
        num_blades: number of blades
        cl_func: function of AoA for coefficient of lift of foil
        cd_func: function of AoA for coefficient of drag of foil
        pitch: pitch of prop (m)
    returns: a function of induced velocity and swirl factor for error in torque
    """

    #Idea from Devlin Ih
    def error_of_2_vars(vel_ind, swirl_fact):
        return torque_error(vel_ind, swirl_fact, angle_vel, radius, dr, chord,
                            num_blades, cl_func, cd_func, pitch)
    return error_of_2_vars

def solve_element(angle_vel, radius, dr, chord, num_blades, cl_func, cd_func,
                  pitch, debug = False):
    """
    Solves for the thrust, torque, induced velocity, and swirl factor of a
    single differential element

    args:
        angle_vel:  angular velocity of prop (m/s)
        radius:     radius at station (m)
        dr:         radial thickness of differential element (m)
        chord:      chord of differential element (m)
        num_blades: number of blades
        cl_func:    function of AoA for coefficient of lift for foil
        cd_func:    function of AoA for coefficient of drag for foil
        pitch:      pitch of prop (m)
        debug:      whether to display intermediate output (default False)
    returns: a list of the form [thrust (N), torque (N*m),
             induced velocity (m/s), swirl factor]
    """

    # Create initial guess
    guess_init = get_initial_guess(angle_vel, pitch)

    # Create error functions
    thrust_error_simple = create_thrust_error_func(angle_vel, radius, dr,
                                    chord, num_blades, cl_func, cd_func, pitch)
    torque_error_simple = create_torque_error_func(angle_vel, radius, dr,
                                    chord, num_blades, cl_func, cd_func, pitch)

    # Solve for vel_ind and siwrl_fact with Newton's Method
    if debug:
        print("Section Radius: " + str(radius))
    guess_final = fz.find_zero2D(thrust_error_simple, torque_error_simple,
                                 guess_init, debug = debug)
    vel_ind     = guess_final[0]
    swirl_fact  = guess_final[1]

    # Calculate final thrust and torque
    """THIS IS THE LOCATION OF THE EMPERICAL CORRECTION FACTOR!"""
    emp_corr_fact = 1 #NOTE: Emperical Correction Factor! (1=no corrrection)
    thrust = thrust_momentum(vel_ind, radius, dr) / emp_corr_fact
    torque = torque_momentum(vel_ind, swirl_fact, angle_vel, radius, dr) /\
             emp_corr_fact
    if debug:
        vel_tan = calc_vel_tan(radius, angle_vel, swirl_fact) #m/s
        phi     = calc_phi(vel_ind, vel_tan) #rad
        theta   = np.arctan(pitch / (2 * radius * np.pi)) #rad
        AoA     = theta - phi #rad
        cl      = cl_func(AoA)
        cd      = cd_func(AoA)
        print("Thrust: " + str(thrust))
        print("Torque: " + str(torque))
        print("Induced Velocity: " + str(vel_ind))
        print("Swirl Factor: " + str(swirl_fact))
        print("AoA: " + str(AoA * (180/np.pi)) + "deg")
        print("Phi: " + str(phi * (180/np.pi)) + "deg")
        print("Theta: " + str(theta * (180/np.pi)) + "deg")
        print("Cl: " + str(cl))
        print("Cd: " + str(cd))
        print("Chord: " + str(chord))

    # Pack into list and return
    result = [thrust, torque, vel_ind, swirl_fact]
    return result

def get_initial_guess(angle_vel, pitch):
    """
    Returns initial guess for differential element

    args:
        angle_vel: angular velocity (rad/s)
        pitch:     pitch of prop (m)
    returns: initial guess as a numpy array with format [vel_ind, swirl_fact]
    """
    
    slip_factor = 0.01 #this is a complete guess
    vel_ind     = angle_vel * 2 * np.pi * pitch * slip_factor #from pitch speed 
    swirl_fact  = 0.4 #this is also a complete guess
    guess       = np.array([vel_ind, swirl_fact])
    return guess

def evaluate_prop_slow(In_prop, debug = False, num_stations = 100):
    """
    Calculates thrust, torque, induced velocity, and swirl factor profiles for
    a prop by saving them in lists. This method is slower than merely solving
    for the final values because it has to store the value at each station

    Total values for prop are calcualted via a Reimann sum over radius

    args:
        In_prop: a Prop object
        debug:   whether or not to display intermediate output (default False)
        num_stations: number of differential elements to use to approximate
            integral over radius of prop (default 100)
    returns: a list of lists of the form [thrust_list, torque_list,
             induced_velocity_list, swirl_factor_list]
    """

    # Read data from prop object
    angle_vel  = In_prop.angle_vel
    num_blades = In_prop.num_blades
    diameter   = In_prop.diameter
    cl_func    = In_prop.get_cl
    cd_func    = In_prop.get_cd
    pitch      = In_prop.pitch

    # Initialize loop parameters
    dr           = (diameter/2) / num_stations
    radius       = 0
    thrusts      = [0] * num_stations
    torques      = [0] * num_stations
    vel_inds     = [0] * num_stations
    swirl_facts  = [0] * num_stations
    
    # Loop through stations
    for i in range(num_stations):
        #iteration specific geometry
        radius += dr
        chord  = In_prop.get_chord(radius)

        # Calculate
        if debug:
            print("\n") #blank line for formatting debug output
        result_list = solve_element(angle_vel, radius, dr, chord, num_blades,
                                    cl_func, cd_func, pitch, debug = debug)
        thrusts[i]     = result_list[0]
        torques[i]     = result_list[1]
        vel_inds[i]    = result_list[2]
        swirl_facts[i] = result_list[3]

    # Pack into list of lists and return
    result = [thrusts, torques, vel_inds, swirl_facts]
    return result

def evaluate_prop_fast(In_prop, debug = False, num_stations = 100):
    """
    Calculates thrust, torque, induced velocity, and swirl factor profiles for
    a prop by saving them in lists. This method is faster because doesn't store
    the value at each station

    Total values for prop are calcualted via a Reimann sum over radius

    args:
        In_prop: a Prop object
        debug:   whether or not to display intermediate output (default False)
        num_stations: number of differential elements to use to approximate
            integral over radius of prop (default 100)
    returns: a list of the form [thrust, torque, induced_velocity, swirl_factor]
    """

    # Read data from prop object
    angle_vel  = In_prop.angle_vel
    num_blades = In_prop.num_blades
    diameter   = In_prop.diameter
    cl_func    = In_prop.get_cl
    cd_func    = In_prop.get_cd
    pitch      = In_prop.pitch

    # Initialize loop parameters
    dr           = (diameter/2) / num_stations
    radius       = 0
    thrust       = 0
    torque       = 0
    vel_ind      = 0
    swirl_fact   = 0
    
    # Loop through stations
    for i in range(num_stations):
        #iteration specific geometry
        radius += dr
        chord  = In_prop.get_chord(radius)

        # Calculate
        if debug:
            print("\n") #blank line for formatting debug output
        result_list = solve_element(angle_vel, radius, dr, chord, num_blades,
                                    cl_func, cd_func, pitch, debug = debug)
        thrust      += result_list[0]
        torque      += result_list[1]
        vel_ind     += result_list[2] #store sum (do division later for avg)
        swirl_fact  += result_list[3] #store sum (do division later for avg)

    # Average, Pack into list, and return
    vel_ind = vel_ind / num_stations
    result = [thrust, torque, vel_ind, swirl_fact]
    return result

def make_cropped_4blader():
    """
    Creates Prop object with attributes and methods correspoding to 4 bladed
    prop tested for Hummingbird during summer of 2022

    returns: Prop object corresponding to cropped 4-blader with 8in pitch
    """

    #basic properties
    num_blades  = 4
    diameter    = 0.254 #0.0889 m is 70% station
    pitch       = 0.21336 #m

    #foil properties correspond to NACA4412 (airfoil used for APC props)
    cl_slope    = 5.72957549575 #cl per rad
    cl_yint     = 0.5
    cd_shape    = 0.58361174676 #cd per rad^2

    #blade shape based on memory
    chord_table = [[0, 0.003175], [0.0254, 0.01905],
                   [0.1016, 0.0254], [0.127, 0.01905]] #all in m

    #assign default angular velocity for ease of testing
    angle_vel   = 539.62 #max tested angular velocity (rad/s)

    #assign to Prop object
    cropped_4_blade = prop_class.Prop(diameter, num_blades, pitch, cl_slope,
                                      cl_yint, cd_shape)
    cropped_4_blade.set_angle_vel(angle_vel)
    cropped_4_blade.set_chord_layout(chord_table)

    #return
    return cropped_4_blade

def rpm_sweep(rpm_max=5152, num_rpms=10):
    """
    Sweeps RPM of cropped 4-bladed prop
    """

    # Create Prop object to evaluate
    cropped_4_blade = make_cropped_4blader()

    # Calculate
    rad_to_rpm    = 60 / (2*np.pi)
    max_angle_vel = rpm_max / rad_to_rpm #rad/s
    num_rpms = num_rpms
    d_rpm    = max_angle_vel / num_rpms
    rpm_list = []
    thrusts  = []
    torques  = []
    for index in range(1,num_rpms):
        angle_vel = (index) * d_rpm
        rpm_list.append(angle_vel * rad_to_rpm)
        cropped_4_blade.set_angle_vel(angle_vel)
        result    = evaluate_prop_fast(cropped_4_blade, debug = False,
                                num_stations = 1000)
        thrusts.append(result[0])
        torques.append(result[1])

    plt.plot(rpm_list, thrusts)
    plt.ylabel("Thrust (N)")
    plt.xlabel("RPM")
    plt.show()

def eval_70_sta():
    """
    Evaluates based on fixed values for 70% station of cropped 4-blader this
    function is for debugging purposes
    """

    # Create Prop object to evaluate
    cropped_4_blade = make_cropped_4blader()

    # Assume good numbers for v_i, swirl factor
    vel_ind    = 5     #m/s based on simple momentum theory (avg. inflox speed)
    swirl_fact = 0.045 #value needed to make dQ_m = 0.0167N
    guess_init = [vel_ind, swirl_fact]

    # Define section params
    radius = 0.0889 #m this is 70% station
    dr     = 0.0127 #m this is 10% of total prop radius
    cl_func    = cropped_4_blade.get_cl
    cd_func    = cropped_4_blade.get_cd
    angle_vel  = cropped_4_blade.angle_vel
    chord      = cropped_4_blade.get_chord(radius)
    num_blades = cropped_4_blade.num_blades
    pitch      = cropped_4_blade.pitch

    # Calculated section params
    vel_tan = calc_vel_tan(radius, angle_vel, swirl_fact) #m/s
    phi     = calc_phi(vel_ind, vel_tan) #rad
    theta   = np.arctan(pitch / (2 * radius * np.pi)) #rad
    AoA     = theta - phi #rad
    cl      = cl_func(AoA)
    cd      = cd_func(AoA)

    # Calculations
    thrust_momentum = thrust_momentum(vel_ind, radius, dr)
    torque_momentum = torque_momentum(vel_ind, swirl_fact, angle_vel, radius,
                                      dr)

def eval_cropped_4blader():
    """
    Calculates and displays thrust for cropped 4 blader prop
    """

    # Create Prop object to evaluate
    cropped_4_blade = make_cropped_4blader()

    # Calculate
    result = evaluate_prop_slow(cropped_4_blade, debug = True,
                                num_stations = 20)
    thrust = sum(result[0])
    torque = sum(result[1])
    vel_inds    = result[2]
    swirl_facts = result[3]

    # Display
    print("\n\nTOTAL THRUST: " + str(thrust))
    print("TOTAL TORQUE: " + str(torque))
    plt.plot(vel_inds)
    plt.ylabel("Induced Velocity (m/s)")
    plt.title("Induced Velocity Distribution")
    plt.show()

def main():
    rpm_sweep()
    #eval_cropped_4blader()
    
if __name__ == "__main__":
    main()
