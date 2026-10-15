from netCDF4 import Dataset
import numpy as np
from scipy.special import sph_harm
import math
from math import sin, cos, tan, pi, fabs, pow, sqrt, factorial

#-------------------------------------------------------------------------------

def icearea_analytical(x, y, z, lat, lon, c):

    icearea = c[0] + c[1]*x + c[2]*y + c[3]*z + c[4]*lat + c[5]*lon + c[6]*(2 + cos(lat)*cos(lat)*cos(2*lon))

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

    tests = ["1", "x", "y", "z", "lat", "lon", "nonlin"]

    gridSizes = [2562, 10242, 40962, 163842]

    rotateCartesianGrid = True
    r = 1.0

    for test in tests:
        if (test == "1"):
           c = [1, 0, 0, 0, 0, 0, 0]
        elif (test == "x"):
           c = [0, 1, 0, 0, 0, 0, 0]
        elif (test == "y"):
           c = [0, 0, 1, 0, 0, 0, 0]
        elif (test == "z"):
           c = [0, 0, 0, 1, 0, 0, 0]
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

            nVertices = len(fileIn.dimensions["nVertices"])
            nCells = len(fileIn.dimensions["nCells"])

            xVertex = fileIn.variables["xVertex"][:]
            yVertex = fileIn.variables["yVertex"][:]
            zVertex = fileIn.variables["zVertex"][:]

            xCell = fileIn.variables["xCell"][:]
            yCell = fileIn.variables["yCell"][:]
            zCell = fileIn.variables["zCell"][:]

            latVertex = fileIn.variables["latVertex"][:]
            lonVertex = fileIn.variables["lonVertex"][:]

            fileIn.close()

            # ice area
            iceAreaVertex = np.zeros(nVertices)
            for iVertex in range(0, nVertices):

                xp, yp, zp = grid_rotation_forward(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], rotateCartesianGrid)
                lat, lon = latlon_from_xyz(xp, yp, zp, r)

                iceareavertex = icearea_analytical(xp, yp, zp, lat, lon, c)
                iceAreaVertex[iVertex] = iceareavertex

            iceAreaCell = np.zeros(nCells)
            for iCell in range(0, nCells):
                xp, yp, zp = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)
                lat, lon = latlon_from_xyz(xp, yp, zp, r)

                iceareacell = icearea_analytical(xp, yp, zp, lat, lon, c)
                iceAreaCell[iCell] = iceareacell

            # output
            filenameOut = "ic_%s_%i.nc" %(test,gridSize)

            fileOut = Dataset(filenameOut, "w", format="NETCDF3_CLASSIC")

            fileOut.createDimension("nVertices", nVertices)
            fileOut.createDimension("nCells", nCells)

            var = fileOut.createVariable("iceAreaVertex","d",dimensions=["nVertices"])
            var[:] = iceAreaVertex[:]
            var = fileOut.createVariable("iceAreaCell","d",dimensions=["nCells"])
            var[:] = iceAreaCell[:]

            fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
