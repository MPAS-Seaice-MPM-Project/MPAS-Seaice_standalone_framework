from netCDF4 import Dataset
import numpy as np
import math
from math import sin, cos, tan, pi, fabs, pow, sqrt, factorial

#-------------------------------------------------------------------------------

def velocities_analytical(x, y, z, lat, lon, c):

    u = c[0] +   c[1]*x + c[2]*y + c[3]*z + c[4]*(x*x + y*y + x*y)     + c[5]*lat
    v = c[6]*x + c[7]*y + c[8]*z          + c[9]*(1/cos(lat)+tan(lat)) + c[10]*lon

    return u, v

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

    tests = ["1-x", "y-z", "latlon", "nonlin"]
    #tests = ["1-x"]

    gridSizes = [2562, 10242, 40962, 163842]

    rotateCartesianGrid = True
    r = 1.0

    for test in tests:
        if (test == "1-x"):
            c = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        elif (test == "y-z"):
            c = [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0]
        elif (test == "latlon"):
           c = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        elif (test == "nonlin"):
           c = [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]

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

            # velocities
            uVelocity = np.zeros(nVertices)
            vVelocity = np.zeros(nVertices)

            for iVertex in range(0, nVertices):

                xp, yp, zp = grid_rotation_forward(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], rotateCartesianGrid)
                lat, lon = latlon_from_xyz(xp, yp, zp, r)

                u, v = velocities_analytical(xp, yp, zp, lat, lon, c)

                uVelocity[iVertex] = u
                vVelocity[iVertex] = v

            solveVelocityPrevious = np.ones(nVertices,dtype="i")

            # output
            filenameOut = "ic_%s_%i.nc" %(test,gridSize)

            fileOut = Dataset(filenameOut, "w", format="NETCDF3_CLASSIC")

            fileOut.createDimension("nVertices", nVertices)
            fileOut.createDimension("nCells", nCells)
            fileOut.createDimension("TWO", 2)

            var = fileOut.createVariable("uVelocity","d",dimensions=["nVertices"])
            var[:] = uVelocity[:]

            var = fileOut.createVariable("vVelocity","d",dimensions=["nVertices"])
            var[:] = vVelocity[:]

            var = fileOut.createVariable("delvDyn","d",dimensions=["nVertices", "TWO"])
            var[:, 0] = uVelocity[:]
            var[:, 1] = vVelocity[:]

            var = fileOut.createVariable("solveVelocityPrevious","i",dimensions=["nVertices"])
            var[:] = solveVelocityPrevious[:]

            fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
