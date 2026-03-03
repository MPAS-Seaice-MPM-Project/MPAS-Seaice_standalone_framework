from netCDF4 import Dataset
import numpy as np
import math
from math import sin, cos

#-------------------------------------------------------------------------------

def wind_stress_analytical(lat, earthRadius):

    stressx = 0.0
    stressy = cos(lat) * sin(lat) / earthRadius

    return stressx, stressy

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
    R = math.sqrt(x*x + y*y + z*z)
    lon = math.atan2(y, x)
    lat = math.asin(z/R)

    return lat, lon

#-------------------------------------------------------------------------------

def latlon_vector_to_xyz_vector(u, v, lat, lon):

# convert a latlon vector to a xyz vector

    vx = (-u) * math.sin(lon) - v * math.sin(lat) * math.cos(lon)
    vy =   u  * math.cos(lon) - v * math.sin(lat) * math.sin(lon)
    vz =                        v * math.cos(lat)

    return vx, vy, vz

#-------------------------------------------------------------------------------

def xyz_vector_to_latlon_vector(vx, vy, vz, lat, lon):

# convert a xyz vector vector to a latlon vector

    u = (-math.sin(lon)) * vx + math.cos(lon)  * vy

    v = (-math.sin(lat)) * math.cos(lon) * vx + (-math.sin(lat)) * math.sin(lon) * vy + math.cos(lat)  * vz

    return u, v

#-------------------------------------------------------------------------------

def xyz_vector_rotation_forward(vx, vy, vz, rotateCartesianGrid):

# rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    if (rotateCartesianGrid):

       vxp = -vz
       vyp =  vy
       vzp =  vx

    else:

       vxp = vx
       vyp = vy
       vzp = vz

    return vxp, vyp, vzp

#-------------------------------------------------------------------------------

def latlon_vector_rotation_forward(u, v, lat, lon, x, y, z, r, rotateCartesianGrid):

# rotate a latlon vector from geographical grid to rotated grid with poles on real equator

# u, v are components of velocity on the geographical grid
# lat, lon are lat/lon of point on the geographical grid
# x, y, z position of the point on the geographical grid
# r  is earth radius

    # perform rotation of the point from geographical grid to rotated grid
    xp, yp, zp = grid_rotation_forward(x, y, z, rotateCartesianGrid)

    # calculate latitude and longitude of the point in the rotated grid
    latp, lonp = latlon_from_xyz(xp, yp, zp, r)

    # convert lat lon vector to xyz vector on gegraphical grid
    vx, vy, vz = latlon_vector_to_xyz_vector(u, v, lat, lon)

    # perform rotation of geographical xyz vector to rotated grid
    vxp, vyp, vzp = xyz_vector_rotation_forward(vx, vy, vz, rotateCartesianGrid)

    # convert xyz vector to lat lon vector on rotated grid
    up, vp = xyz_vector_to_latlon_vector(vxp, vyp, vzp, latp, lonp)

    return up, vp

#-------------------------------------------------------------------------------

def create_ics(earthRadius, rotateCartesianGrid):

    forceTypes = ["cos_lat"]

    gridSizes = [2562, 10242, 40962, 163842]

    iceDensity = 917.0

    for forceType in forceTypes:

        print("Mesh IC, forceType for velocity: ", forceType)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)

            # input
            filenameIn = "grid.%i.nc" %(gridSize)

            fileIn = Dataset(filenameIn,"r")

            nCells = len(fileIn.dimensions["nCells"])
            nVertices = len(fileIn.dimensions["nVertices"])

            cellsOnVertex = fileIn.variables["cellsOnVertex"][:,:]

            xVertex = fileIn.variables["xVertex"][:]
            yVertex = fileIn.variables["yVertex"][:]
            zVertex = fileIn.variables["zVertex"][:]

            xCell = fileIn.variables["xCell"][:]
            yCell = fileIn.variables["yCell"][:]
            zCell = fileIn.variables["zCell"][:]

            latCell = fileIn.variables["latCell"][:]

            latVertex = fileIn.variables["latVertex"][:]
            lonVertex = fileIn.variables["lonVertex"][:]

            fileIn.close()

            # output
            filenameOut = "ic_%s_%i.nc" %(forceType, gridSize)

            fileOut = Dataset(filenameOut, "w", format="NETCDF3_CLASSIC")

            fileOut.createDimension("nCells", nCells)
            fileOut.createDimension("nVertices", nVertices)
            fileOut.createDimension("nCategories", size=1)
            fileOut.createDimension("ONE", size=1)

            # ice area and volume
            #iceAreaCell = np.zeros(nCells)
            #iceVolumeCell = np.zeros(nCells)
            #iceAreaCategory = np.zeros([nCells, 1])
            #iceVolumeCategory = np.zeros([nCells, 1])

            #for iCell in range (0, nCells):
            #    iceAreaCell[iCell] = 1.0
            #    iceVolumeCell[iCell] = 1.0

            #var = fileOut.createVariable("iceAreaCell",   'd', dimensions=("nCells"))
            #var[:] = iceAreaCell[:]
            #var = fileOut.createVariable("iceVolumeCell", 'd', dimensions=("nCells"))
            #var[:] = iceVolumeCell[:]

            #iceAreaCategory[:,0]   = iceAreaCell[:]
            #var = fileOut.createVariable("iceAreaCategory",   'd', dimensions=("nCells","nCategories","ONE"))
            #var[:] = iceAreaCategory

            #iceVolumeCategory[:,0] = iceVolumeCell[:]
            #var = fileOut.createVariable("iceVolumeCategory", 'd', dimensions=("nCells","nCategories","ONE"))
            #var[:] = iceVolumeCategory

            # wind_stress
            airStressVertexU = np.zeros(nVertices)
            airStressVertexV = np.zeros(nVertices)

            for iVertex in range(0, nVertices):

                lat, lon = latlon_from_xyz(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], earthRadius)

                stressx, stressy = wind_stress_analytical(lat, earthRadius)

                stressx, stressy = latlon_vector_rotation_forward(
                                       stressx, stressy,
                                       lat, lon,
                                       xVertex[iVertex], yVertex[iVertex], zVertex[iVertex],
                                       earthRadius, rotateCartesianGrid)

                airStressVertexU[iVertex] = (40.0 * 40.0) * iceDensity * stressx
                airStressVertexV[iVertex] = (40.0 * 40.0) * iceDensity * stressy


            var = fileOut.createVariable("airStressVertexU","d",dimensions=["nVertices"])
            var[:] = airStressVertexU[:]

            var = fileOut.createVariable("airStressVertexV","d",dimensions=["nVertices"])
            var[:] = airStressVertexV[:]

            fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics(earthRadius, rotateCartesianGrid)
