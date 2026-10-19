from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def create_bathymetry_interpolation_file():

    earthRadius = 6371229.0

    filein = Dataset("grid.nc","r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    zCell = filein.variables["zCell"][:]
    zVertex = filein.variables["zVertex"][:]

    filein.close()

    bathymetry = np.zeros(nCells)
    uVelocity = np.zeros(nVertices)

    for iCell in range(0,nCells):
        bathymetry[iCell] = (zCell[iCell] + earthRadius) / earthRadius

    for iVertex in range(0,nVertices):
        uVelocity[iVertex] = (zVertex[iVertex] + earthRadius) / earthRadius


    fileout = Dataset("bathymetry_interpolation.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)

    var = fileout.createVariable("bathymetry","d",dimensions=["nCells"])
    var[:] = bathymetry[:]

    var = fileout.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = uVelocity[:]

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_bathymetry_interpolation_file()
