from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#--------------------------------------------------------

def L2_norm_cell(numerical, analytical, nCells, areaCell):

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):
       norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)
       denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)
    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_cell(filenameGrid, filenameIC, filename):

    fileGrid = Dataset(filenameGrid, "r")
    nCells = len(fileGrid.dimensions["nCells"])
    areaCell = fileGrid.variables["areaCell"][:]
    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")
    iceAreaCellAnalytical = fileIC.variables["iceAreaCell"][:]
    fileIC.close()

    fileMPAS = Dataset(filename, "r")
    iceAreaCell = fileMPAS.variables["iceAreaCell"][0,:]
    fileMPAS.close()

    norm = L2_norm_cell(iceAreaCell, iceAreaCellAnalytical, nCells, areaCell)

    print('cell norm ',norm)
    return norm

#--------------------------------------------------------

def L2_norm_vertex(numerical, analytical, nVertices, areaTriangle):

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):
       norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)
       denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)
    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_vertex(filenameGrid, filenameIC, filename):

    fileGrid = Dataset(filenameGrid, "r")
    nVertices = len(fileGrid.dimensions["nVertices"])
    areaTriangle = fileGrid.variables["areaTriangle"][:]
    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")
    iceAreaVertexAnalytical = fileIC.variables["iceAreaVertex"][:]
    fileIC.close()

    fileMPAS = Dataset(filename, "r")
    iceAreaVertex = fileMPAS.variables["iceAreaVertex"][0,:]
    fileMPAS.close()

    norm = L2_norm_vertex(iceAreaVertex, iceAreaVertexAnalytical, nVertices, areaTriangle)

    print('vertex norm ',norm)
    return norm

#--------------------------------------------------------

def get_resolution(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    cellsOnEdge = fileMPAS.variables["cellsOnEdge"][:]
    cellsOnEdge[:] = cellsOnEdge[:] - 1

    dcEdge = fileMPAS.variables["dcEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        iCell1 = cellsOnEdge[iEdge,0]
        iCell2 = cellsOnEdge[iEdge,1]
        if (iCell1 != -1 and iCell2 != -1):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def reconstruction_scaling():

    mpl.rc('font', family='Times New Roman', size=6)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    gridTypes = ["hex","quad"]
    #gridTypes = ["quad"]

    grids = {"hex" :["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}
    #grids = {"hex" :["0082x0094"],
    #         "quad":["0080x0080"]}


    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    fig, axes = plt.subplots(1,1,figsize=(3,3))

    xMin = 2e-3
    xMax = 3.5e-3
    yMin = 4e-3

    # linear scaling
    scale = yMin / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    axes.loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
    axes.loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        x = []
        yVertex = []
        yCell = []

        for grid in grids[gridType]:

            print("    Grid: ", grid)

            filenameGrid = 'grid_%s_%s.nc' %(gridType,grid)
            filenameIC = "ic_%s_%s.nc" %(gridType,grid)
            filename = "output_%s_%s/output.2000.nc" %(gridType, grid)

            print("      ", filename, filenameIC, filenameGrid)

            normVertex = get_norm_vertex(filenameGrid, filenameIC, filename)
            normCell = get_norm_cell(filenameGrid, filenameIC, filename)

            x.append(get_resolution(filenameGrid))
            yVertex.append(normVertex)
            yCell.append(normCell)

        axes.loglog(x, yVertex, marker='o', color='blue', ls=lineStyles[gridType], markersize=5.0, label=('vertex, ' + gridType))
        axes.loglog(x, yCell, marker='^', color='red', ls=lineStyles[gridType], markersize=5.0, label=('cell, ' + gridType))

    axes.legend(frameon=False, loc=4, fontsize=6, handlelength=4)

    axes.set_xlabel("Grid resolution")
    axes.set_ylabel(r"$L_2$ error norm")

    axes.set_xticks([.01, .001])
    plt.tight_layout()
    plt.savefig("reconstruction_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_scaling()
