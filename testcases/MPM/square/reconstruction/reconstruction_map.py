from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sys
sys.path.append("../../spherical_operators/strain")
from strain_map import get_mpas_patch_collection
from strain_map import plot_subfigure

#---------------------------------------------------------------

def reconstruction_map():

    # hex grid
    fileGrid = Dataset("grid_hex_0082x0094.nc","r")
    # quad grid
    #fileGrid = Dataset("grid_quad_0080x0080.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    cellsOnCell = fileGrid.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]
    # overwrite lat arrays for plotting planar data
    latVertex = np.ones(nVertices)
    latCell = np.ones(nCells)

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    fileGrid.close()

    # ic hex
    fileIC = Dataset("ic_hex_0082x0094.nc","r")
    # ic quad
    #fileIC = Dataset("ic_quad_0080x0080.nc","r")

    iceAreaCellAnalytical = fileIC.variables["iceAreaCell"][:]
    iceAreaVertexAnalytical = fileIC.variables["iceAreaVertex"][:]

    fileIC.close()

    # computed fields hex
    fileMPM = Dataset("./output_hex_0082x0094/output.2000.nc","r")
    # computed fields hex
    #fileMPM = Dataset("./output_quad_0080x0080/output.2000.nc","r")

    iceAreaCellMPM = fileMPM.variables["iceAreaCell"][0,:]
    iceAreaVertexMPM  = fileMPM.variables["iceAreaVertex"][0,:]

    iceAreaCellDiff = iceAreaCellAnalytical - iceAreaCellMPM
    iceAreaVertexDiff = iceAreaVertexAnalytical - iceAreaVertexMPM

    print("Cell Diff:  ",
          np.amin(iceAreaCellDiff), np.amax(iceAreaCellDiff))
    print("Vertex Diff:  ",
          np.amin(iceAreaVertexDiff), np.amax(iceAreaVertexDiff))

    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minArea= -1.0
    maxArea=  1.0

    minAreaDiff = -1.e-2
    maxAreaDiff =  1.e-2

    fig, axes = plt.subplots(2, 2)

    fig.set_size_inches(6, 6)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexAnalytical, minArea, maxArea, xMin, xMax, yMin, yMax, \
                   False, False, r'iceAreaVertex', '(a)', False)
    plot_subfigure(axes[0,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellAnalytical, minArea, maxArea, xMin, xMax, yMin, yMax, \
                   False, True, r'iceAreaCell', '(b)', True)
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexDiff, minAreaDiff, maxAreaDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'iceAreaVertex Diff', '(c)', False)
    plot_subfigure(axes[1,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellDiff, minAreaDiff, maxAreaDiff, xMin, xMax, yMin, yMax, \
                   False, True, r'iceAreaCell Diff', '(d)', True)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("reconstruction_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_map()
