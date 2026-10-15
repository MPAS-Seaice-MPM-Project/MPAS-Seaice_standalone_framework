from netCDF4 import Dataset
import numpy as np
from math import sin, cos, pi, pow

import sys
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

import argparse

#-------------------------------------------------------------

def linear(x,y):

    Lx = 1.0
    Ly = 1.0

    a = (x - 0.5 * Lx)  - (y - 0.5 * Ly)

    return a

#-------------------------------------------------------------

def quadratic(x,y):

    Lx = 1.0
    Ly = 1.0

    a = pow(x - 0.5 * Lx, 2) + 2 * pow(y - 0.5 * Ly, 2)

    return a

#-------------------------------------------------------------

def sinusoid(x,y):

    A = 2.56
    B = 2.56
    C = 2.56
    D = 2.56

    Lx = 1.0
    Ly = 1.0

    a = sin((2.0 * pi * x * A) / Lx) * sin((2.0 * pi * y * B) / Ly) + \
        sin((2.0 * pi * x * C) / Lx) * sin((2.0 * pi * y * D) / Ly)

    return a

#-------------------------------------------------------------

def ice_area(x, y, functionType):

    if (functionType == "linear"):
            a = linear(x,y)
    elif (functionType == "quadratic"):
            a = quadratic(x,y)
    elif (functionType == "sinusoid"):
            a = sinusoid(x,y)

    return a

#-------------------------------------------------------------

def create_particle(gridfile, partfile, functionType):

    initial_particle_positions(gridfile,
                               partfile,
                               "number",
                                4,
                               "even",
                               "uniform",
                                1.0)

    # read position data and iceAreaCell from particle file
    fileOut = Dataset(partfile, "r+", format="NETCDF3_CLASSIC")
    nParticles = len(fileOut.dimensions["nParticles"])
    posnMP = fileOut.variables["posnMP"][:]
    iceAreaCellMP = fileOut.variables["iceAreaCellMP"][:]
    iceAreaCategoryMP = fileOut.variables["iceAreaCategoryMP"][:]

    # modify iceAreaCell at material points
    for iParticle in range(0, nParticles):
              xp = posnMP[iParticle, 0]
              yp = posnMP[iParticle, 1]

              icearea= ice_area(xp, yp, functionType)

              iceAreaCellMP[iParticle] = icearea
              iceAreaCategoryMP[iParticle, 0] = icearea

    fileOut.variables["iceAreaCellMP"][:] = iceAreaCellMP
    fileOut.variables["iceAreaCategoryMP"][:] = iceAreaCategoryMP

    fileOut.close()

#-------------------------------------------------------------

def create_particles(functionType):

    gridTypes = ["hex","quad"]

    grids = {"hex": ["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}


    for gridType in gridTypes:
        for grid in grids[gridType]:

            gridfile = "grid_%s_%s.nc" %(gridType,grid)
            partfile = "particles_%s_%s.nc" %(gridType,grid)

            create_particle(gridfile, partfile, functionType)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', dest='functionType')

    args = parser.parse_args()

    create_particles(args.functionType)
