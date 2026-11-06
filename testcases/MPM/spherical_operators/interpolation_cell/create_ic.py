from netCDF4 import Dataset
import numpy as np
import math
from math import sin, cos, tan, pi, fabs, pow, sqrt, factorial

#-------------------------------------------------------------------------------

def icearea_analytical(x, y, z, lat, lon, c):

    icearea = c[0] + c[1]*x + c[2]*y + c[3]*z + c[4]*(lat+90) + c[5]*(lon+180) + c[6]*(2 + cos(lat)*cos(lat)*cos(2*lon))

    return icearea

#-------------------------------------------------------------------------------

def grid_rotation_forward(x, y, z, rotateCartesianGrid):

    # rotate xyz coordinates from geographical grid to rotated grid with poles on real equator

    if (rotateCartesianGrid):

       xp = -z
       yp = y
       zp = x

    else:

       xp = x
       yp = y
       zp = z

    return xp, yp, zp

#-------------------------------------------------------------------------------

def latlon_from_xyz(x, y, z, r):

    # given xyz coordinates determine the latitude and longitude

    lon = math.atan2(y, x)
    lat = math.asin(z/r)

    return lat, lon

#-------------------------------------------------------------------------------

def create_ic():

    tests = ["1", "1+x", "1+y", "1+z", "lat", "lon", "nonlin"]

    gridSizes = [2562, 10242, 40962, 163842]

    rotateCartesianGrid = True
    r = 1.0

    for test in tests:
        if (test == "1"):
           c = [1, 0, 0, 0, 0, 0, 0]
        elif (test == "1+x"):
           c = [1, 1, 0, 0, 0, 0, 0]
        elif (test == "1+y"):
           c = [1, 0, 1, 0, 0, 0, 0]
        elif (test == "1+z"):
           c = [1, 0, 0, 1, 0, 0, 0]
        elif (test == "lat"):
           c = [0, 0, 0, 0, 1, 0, 0]
        elif (test == "lon"):
           c = [0, 0, 0, 0, 0, 1, 0]
        elif (test == "nonlin"):
           c = [0, 0, 0, 0, 0, 0, 1]

        print("Mesh IC, testcase: ", test)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)

            # input
            filenameIn = "grid.%i.nc" %(gridSize)

            fileIn = Dataset(filenameIn,"r")

            nCells = len(fileIn.dimensions["nCells"])
            nVertices = len(fileIn.dimensions["nVertices"])

            xCell = fileIn.variables["xCell"][:]
            yCell = fileIn.variables["yCell"][:]
            zCell = fileIn.variables["zCell"][:]

            xVertex = fileIn.variables["xVertex"][:]
            yVertex = fileIn.variables["yVertex"][:]
            zVertex = fileIn.variables["zVertex"][:]

            latCell = fileIn.variables["latCell"][:]
            lonCell = fileIn.variables["lonCell"][:]

            latVertex = fileIn.variables["latVertex"][:]
            lonVertex = fileIn.variables["lonVertex"][:]

            fileIn.close()

            # ice area 
            iceAreaCell = np.zeros(nCells)
            uAirVelocity = np.zeros(nCells)
            vAirVelocity = np.zeros(nCells)

            for iCell in range(0, nCells):

                xp, yp, zp = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)
                lat, lon = latlon_from_xyz(xp, yp, zp, r)

                iceareacell = icearea_analytical(xp, yp, zp, lat, lon, c)

                iceAreaCell[iCell] = iceareacell
                uAirVelocity[iCell] = iceareacell
                vAirVelocity[iCell] = 0.0

            # output
            filenameOut = "ic_%s_%i.nc" %(test,gridSize)

            fileOut = Dataset(filenameOut, "w", format="NETCDF3_CLASSIC")

            fileOut.createDimension("nVertices", nVertices)
            fileOut.createDimension("nCells", nCells)
            fileOut.createDimension("TWO", 2)

            var = fileOut.createVariable("iceAreaCell","d",dimensions=["nCells"])
            var[:] = iceAreaCell[:]
            var = fileOut.createVariable("uAirVelocity","d",dimensions=["nCells"])
            var[:] = uAirVelocity[:]
            var = fileOut.createVariable("vAirVelocity","d",dimensions=["nCells"])
            var[:] = vAirVelocity[:]

            fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
