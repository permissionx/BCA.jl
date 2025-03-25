#!/usr/bin/python3
from std_simulation import *
if __name__ == '__main__':
    with open("p_debug.csv", 'w') as file:
        file.write("ntargets,p0x,p0y,p0z,p1x,p1y,p1z\n")
    myprocess = MySimulation()