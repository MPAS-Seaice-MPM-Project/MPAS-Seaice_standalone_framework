from netCDF4 import Dataset
import numpy as np
import subprocess
import os
import argparse

#-------------------------------------------------------------------------------

def create_square_quad_mesh(nx, ny,
                            lx, ly,
                            x0, y0):

    filenameOut = "grid.quad.%ix%i.nc" %(nx,ny)

    dx = lx / float(nx)
    dy = ly / float(ny)

    nCells = nx * ny
    nVertices = (nx+1) * (ny+1)
    vertexDegree = 4

    xCell = np.zeros(nCells)
    yCell = np.zeros(nCells)
    zCell = np.zeros(nCells)

    xVertex = np.zeros(nVertices)
    yVertex = np.zeros(nVertices)
    zVertex = np.zeros(nVertices)

    cellsOnVertex = np.zeros((nVertices,vertexDegree))

    iCell = 0
    for ix in range(0, nx):
        for iy in range(0, ny):
            xCell[iCell] = x0 + dx*(float(ix)+0.5)
            yCell[iCell] = y0 + dy*(float(iy)+0.5)
            iCell += 1

    iVertex = 0
    for ix in range(0, nx+1):
        for iy in range(0, ny+1):
            xVertex[iVertex] = x0 + dx*float(ix)
            yVertex[iVertex] = y0 + dy*float(iy)
            iVertex += 1

    dxi = [0, 1, 1, 0]
    dyi = [0, 0, 1, 1]

    cellsOnVertex[:,:] = 0

    iVertex = 0
    for ix in range(0, nx+1):
        for iy in range(0, ny+1):

            ixc = ix-1
            iyc = iy-1

            for iCellOnVertex in range(0, vertexDegree):

                ixcc = ixc + dxi[iCellOnVertex]
                iycc = iyc + dyi[iCellOnVertex]

                if (ixcc <  0 or
                    ixcc >= nx or
                    iycc <  0 or
                    iycc >= ny):
                    iCell = -1
                else:
                    iCell = ixcc*ny + iycc

                cellsOnVertex[iVertex, iCellOnVertex] = iCell+1

            iVertex += 1

    # create mesh file description
    filenameTmp1Out = "grid.tmp1.nc"
    fileOut = Dataset(filenameTmp1Out, "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("vertexDegree", vertexDegree)

    var = fileOut.createVariable("xCell", "d", dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileOut.createVariable("yCell", "d", dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileOut.createVariable("zCell", "d", dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileOut.createVariable("xVertex", "d", dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileOut.createVariable("yVertex", "d", dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileOut.createVariable("zVertex", "d", dimensions=["nVertices"])
    var[:] = zVertex[:]

    var = fileOut.createVariable("cellsOnVertex", "d", dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileOut.on_a_sphere = "NO"

    fileOut.close()

    # create full mesh with the mesh converter
    MPAS_TOOLS_DIR = os.environ.get('MPAS_TOOLS_DIR')
    if (MPAS_TOOLS_DIR is None):
        raise Exception("MPAS_TOOLS_DIR environment variable must be set")
    MpasMeshConverter = MPAS_TOOLS_DIR+"/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x"
    if (not os.path.isfile(MpasMeshConverter)):
        raise Exception("MpasMeshConverter executable must be built")

    subprocess.run([MpasMeshConverter,filenameTmp1Out,filenameOut])

    # add attributes
    filein = Dataset(filenameOut,"a")

    filein.nx = nx
    filein.ny = ny
    filein.lx = lx
    filein.ly = ly
    filein.x0 = x0
    filein.y0 = y0

    filein.close()

    return filenameOut

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--nx', dest='nx', type=int, required=True)
    parser.add_argument('--ny', dest='ny', type=int, required=True)

    parser.add_argument('--lx', dest='lx', type=float, required=True)
    parser.add_argument('--ly', dest='ly', type=float, required=True)

    parser.add_argument('--x0', dest='x0', type=float, default=0.0)
    parser.add_argument('--y0', dest='y0', type=float, default=0.0)

    args = parser.parse_args()

    create_square_quad_mesh(args.nx, args.ny,
                            args.lx, args.ly,
                            args.x0, args.y0)
