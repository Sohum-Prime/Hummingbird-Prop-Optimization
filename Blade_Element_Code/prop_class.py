"""
Datatype for propellers to go with eval_prop code

:AUTHOR: William Skelly
:VERSION: 1.0
"""

import numpy as np

class Prop():
    def __init__(self, diameter, num_blades, pitch, cl_slope, cl_yint,
                 cd_shape):
        """
        Initializes Prop object

        args:
            diameter:   diameter in meters
            num_blades: number of blades
            pitch:      pitch of prop (m)
            cl_slope: slope of the cl line (cl per radian)
            cl_yint:  y intercept of cl line (cl)
            cd_shape: shape parameter of cd curve (cd/rad^2)
        """
        self.num_blades = num_blades
        self.diameter   = diameter
        self.pitch      = pitch
        self.get_cl     = self.create_cl_func(cl_slope, cl_yint)
        self.get_cd     = self.create_cd_func(cd_shape)

    def set_chord_layout(self, chord_layout):
        """
        Sets list of lists representing chord layout of the form: [[r1,c1],
            [r2,c2], etc...] NOTE: chord_layout must be sorted by radius from
            low to high
        """
        self.chord_layout = chord_layout

    def set_angle_vel(self, angle_vel):
        """
        Set angular velocity in radians/sec
        """
        self.angle_vel = angle_vel

    def get_chord(self, radius):
        """
        Returns chord at radius

        args:
            radius: radius of station in meters
        returns: chord at that station
        """

        #find first station in table outboard of radius
        for index in range(len(self.chord_layout)):
            point = self.chord_layout[index]
            if point[0] > radius:
                #approximate between given points with a straight line
                point_prev = self.chord_layout[index - 1]
                slope  = (point_prev[1] - point[1]) / (point_prev[0] - point[0])
                r_delta = radius - point[0]
                chord = point[1] + (slope * r_delta)
                return chord
        # if radius outside defined points, return outermost given chord
        return self.chord_layout[-1][1]

    def create_cl_func(self, cl_slope, cl_yint):
        """
        Creates the function that returns the cl for a given AoA

        args:
            cl_slope: slope of the cl line (per radian)
            cl_yint:  y intercept of cl line (cl)
        """

        def cl_of_AoA(AoA):
            stall_AoA      = 0.174533 #rad
            if AoA < stall_AoA:
                return (cl_slope * AoA) + cl_yint
            else: #return constant Cl post-stall
                stall_cl = (cl_slope * stall_AoA) + cl_yint
                #((-cl_slope/40) * (AoA - stall_AoA)) + stall_cl
                return stall_cl
        return cl_of_AoA

    def create_cd_func(self, cd_shape):
        """
        Creates the function that returns the cd for a given AoA

        args:
            cd_shape: shape parameter for the cd parabola (cd/rad^2)
        """

        def cd_of_AoA(AoA):
            return pow(AoA,2) * cd_shape
        return cd_of_AoA

def main():
    cd_shape    = 0.58361174676
    cropped_4_blade = Prop(.25, 4, .2133, 5.729,
                                      .5, cd_shape)
    chord_table = [[0, 0.003175], [0.0254, 0.01905],
                   [0.1016, 0.0254], [0.127, 0.01905]]
    cropped_4_blade.set_chord_layout(chord_table)
    print("Chord: " + str(cropped_4_blade.get_chord(0.05)))
    print("Cd: " + str(cropped_4_blade.get_cd(0.174533)))

if __name__ == "__main__":
    main()
