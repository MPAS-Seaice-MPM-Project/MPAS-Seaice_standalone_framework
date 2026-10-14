from netCDF4 import Dataset
import numpy as np
from math import atan2, sqrt, asin, radians, sin, cos
from tqdm import tqdm

#-------------------------------------------------------------------------------

def subtriangle_centers(a, b, c, n):

    def point(i,j):
        return a + (i/n) * (b-a) + (j/n) * (c-a)

    centers = []

    for i in range(n):
        for j in range(n-i):

            triangle = np.array([
                point(i, j),
                point(i+1, j),
                point(i, j+1)
            ])

            centers.append(triangle.mean(axis=0))

            if (i + j <= n - 2):
                triangle = np.array([
                    point(i+1, j),
                    point(i, j+1),
                    point(i+1, j+1)
                ])

                centers.append(triangle.mean(axis=0))

    return np.asarray(centers)

#-------------------------------------------------------------------------------

def xyz_to_latlon(x, y, z):

    lon = atan2(y, x)
    lat = asin(z/sqrt(x*x+y*y+z*z))

    return lat, lon

#-------------------------------------------------------------------------------

def latlon_to_xyz(lat, lon):

    r = 6371229.0
    x = r * cos(lat) * cos(lon)
    y = r * cos(lat) * sin(lon)
    z = r * sin(lat)

    return x, y, z

#-------------------------------------------------------------------------------

def create_interpolation_particles():

    icebergLength = 1000.0
    icebergHeight = 1000.0
    distLimit = 4e5

    lat0 = -68.633867
    lon0 = -69.291799

    lat0 = radians(lat0)
    lon0 = radians(lon0)
    x0, y0, z0 = latlon_to_xyz(lat0, lon0)

    nSubTriangle = 10

    fileMesh = Dataset("grid.nc","r")

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:] - 1

    xCell = fileMesh.variables["xCell"][:]
    yCell = fileMesh.variables["yCell"][:]
    zCell = fileMesh.variables["zCell"][:]

    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]
    zVertex = fileMesh.variables["zVertex"][:]

    fileMesh.close()

    cellIDCreationIB = []
    creationIndexIB = []
    latIceberg = []
    lonIceberg = []

    nCellsUsed = 0
    for iCell in tqdm(range(0,nCells)):

        dist2 = ((xCell[iCell]-x0)*(xCell[iCell]-x0)+
                 (yCell[iCell]-y0)*(yCell[iCell]-y0)+
                 (zCell[iCell]-z0)*(zCell[iCell]-z0))
        if (dist2 < distLimit*distLimit):
            nCellsUsed += 1

            creationIndex = 1

            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

                iVertexOnCell1 = iVertexOnCell
                iVertexOnCell2 = (iVertexOnCell + 1) % nEdgesOnCell[iCell]

                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]
                iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                positions = subtriangle_centers(np.array([xCell[iCell],yCell[iCell],zCell[iCell]]),
                                                np.array([xVertex[iVertex1],yVertex[iVertex1],zVertex[iVertex1]]),
                                                np.array([xVertex[iVertex2],yVertex[iVertex2],zVertex[iVertex2]]),
                                                nSubTriangle)

                #print(iCell, iVertexOnCell, iVertex1, iVertex2, positions)

                for i in range(0,positions.shape[0]):
                    x = positions[i,0]
                    y = positions[i,1]
                    z = positions[i,2]
                    lat, lon = xyz_to_latlon(x, y, z)
                    latIceberg.append(lat)
                    lonIceberg.append(lon)
                    cellIDCreationIB.append(iCell)
                    creationIndexIB.append(creationIndex)
                    creationIndex += 1

    nIcebergs = len(cellIDCreationIB)
    print("nIcebergs: ", nIcebergs, ", nCellsUsed: ", nCellsUsed)

    cellIDCreationIB = np.array(cellIDCreationIB)
    creationIndexIB = np.array(creationIndexIB)
    latIceberg = np.array(latIceberg)
    lonIceberg = np.array(lonIceberg)

    fileout = Dataset("icebergs_interpolation.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",nIcebergs)

    var = fileout.createVariable("cellIDCreationIB","d",dimensions=["nIcebergs"])
    var[:] = cellIDCreationIB[:]

    var = fileout.createVariable("creationIndexIB","d",dimensions=["nIcebergs"])
    var[:] = creationIndexIB[:]

    var = fileout.createVariable("latIceberg","d",dimensions=["nIcebergs"])
    var[:] = latIceberg[:]

    var = fileout.createVariable("lonIceberg","d",dimensions=["nIcebergs"])
    var[:] = lonIceberg[:]

    var = fileout.createVariable("icebergLength","d",dimensions=["nIcebergs"])
    var[:] = icebergLength

    var = fileout.createVariable("icebergHeight","d",dimensions=["nIcebergs"])
    var[:] = icebergHeight

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_interpolation_particles()
