from netCDF4 import Dataset
import numpy as np
from math import sin, cos, pi, pow

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

    a = pow(x - 0.5 * Lx, 2) + 0.5 * pow(y - 0.5 * Ly, 2)

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

def create_ic(gridfile, icfile, functionType):

    # load grid file
    grid = Dataset(gridfile, "r")

    nCells = len(grid.dimensions["nCells"])
    nVertices = len(grid.dimensions["nVertices"])
    vertexDegree = len(grid.dimensions["vertexDegree"])

    cellsOnVertex = grid.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    xCell = grid.variables["xCell"][:]
    yCell = grid.variables["yCell"][:]

    xVertex = grid.variables["xVertex"][:]
    yVertex = grid.variables["yVertex"][:]

    grid.close()

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    # calculate output variables
    iceAreaCell = np.empty(nCells)
    iceAreaVertex = np.empty(nVertices)

    for iVertex in range(0,nVertices):

        x = (xVertex[iVertex] - xMin)
        y = (yVertex[iVertex] - yMin)

        icearea= ice_area(x, y, functionType)

        iceAreaVertex[iVertex] = icearea

    for iCell in range(0, nCells):

        x = (xCell[iCell] - xMin)
        y = (yCell[iCell] - yMin)

        icearea = ice_area(x, y, functionType)

        iceAreaCell[iCell] = icearea

    # create output file
    fileOut = Dataset(icfile, "w", format="NETCDF3_CLASSIC")

    fileOut.functionType = functionType

    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("nCells", nCells)

    var = fileOut.createVariable("iceAreaVertex","d",dimensions=["nVertices"])
    var[:] = iceAreaVertex[:]

    var = fileOut.createVariable("iceAreaCell","d",dimensions=["nCells"])
    var[:] = iceAreaCell[:]

    fileOut.close()

#-------------------------------------------------------------

def create_ics(functionType):

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
            icfile = "ic_%s_%s.nc" %(gridType,grid)

            create_ic(gridfile, icfile, functionType)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', dest='functionType')

    args = parser.parse_args()

    create_ics(args.functionType)
