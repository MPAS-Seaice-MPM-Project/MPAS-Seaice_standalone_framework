from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import os

degreesToRadians = math.pi / 180.0

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > 20.0*degreesToRadians):

            polygonVertices = []

            useVertex = False
            for iCellOnVertex in range(0,vertexDegree[iVertex]):

                iCell = cellsOnVertex[iVertex,iCellOnVertex] - 1

                polygonVertices.append((xCell[iCell],yCell[iCell]))

                if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                    yCell[iCell] >= minY and yCell[iCell] <= maxY):
                    useVertex = True

            if (useVertex):
                polygon = Polygon(polygonVertices)
                patches.append(polygon)

                colours.append(mpasArray[iVertex])

                minval = min(minval,mpasArray[iVertex])
                maxval = max(maxval,mpasArray[iVertex])

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=True)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)

    patchCollection.set_clim(vmin=vmin,vmax=vmax)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
    axes.add_collection(patchCollection)
    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title, fontsize=8)

    if (subfigureLabel != None):
        axes.text(0.02, 0.89, subfigureLabel, verticalalignment='bottom', horizontalalignment='left',transform=axes.transAxes, fontsize=8)

    if (colorbar):
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(patchCollection,cax=cax)
        if (unityBar):
            cb.ax.set_yticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()

#---------------------------------------------------------------

def reconstruction_map():

    res = "40962" # specify the resolution to plot
    tests = ["lat","nonlin"] # specify two tests to plot

    # grid
    filename = "grid.%s.nc" %(res)
    fileGrid = Dataset(filename,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    # ic
    filename = "ic_%s_%s.nc" %(tests[0],res)
    fileIC = Dataset(filename,"r")
    iceAreaVertexAnalytical1 = fileIC.variables["iceAreaVertex"][:]
    print("iceArea ", tests[0], " vertex: ",
          np.amin(iceAreaVertexAnalytical1), np.amax(iceAreaVertexAnalytical1))
    iceAreaCellAnalytical1= fileIC.variables["iceAreaCell"][:]
    print("iceArea ", tests[0], "cell: ",
          np.amin(iceAreaCellAnalytical1), np.amax(iceAreaCellAnalytical1))
    fileIC.close()

    filename = "ic_%s_%s.nc" %(tests[1],res)
    fileIC = Dataset(filename,"r")
    iceAreaVertexAnalytical2 = fileIC.variables["iceAreaVertex"][:]
    print("iceArea ", tests[1], " vertex: ",
          np.amin(iceAreaVertexAnalytical2), np.amax(iceAreaVertexAnalytical2))
    iceAreaCellAnalytical2 = fileIC.variables["iceAreaCell"][:]
    print("iceArea ", tests[1], " cell: ",
          np.amin(iceAreaCellAnalytical2), np.amax(iceAreaCellAnalytical2))
    fileIC.close()

    # mpm
    filenameMPM = "./output_%s_%s/output.2000.nc" %(tests[0],res)
    if (not os.path.exists(filenameMPM)):
        raise Exception("Missing output file: %s" %(filenameMPM))
    fileMPM = Dataset(filenameMPM,"r")
    iceAreaVertexMPM1 = fileMPM.variables["iceAreaVertex"][0,:]
    iceAreaVertexDiff1 = (iceAreaVertexMPM1 - iceAreaVertexAnalytical1)
    print("MPM Reconstruction Error ", tests[0], " vertex:  ",
          np.amin(iceAreaVertexDiff1), np.amax(iceAreaVertexDiff1))
    iceAreaCellMPM1 = fileMPM.variables["iceAreaCell"][0,:]
    iceAreaCellDiff1 = (iceAreaCellMPM1- iceAreaCellAnalytical1)
    print("MPM Reconstruction Error ", tests[0], " cell:  ",
          np.amin(iceAreaCellDiff1), np.amax(iceAreaCellDiff1))
    fileMPM.close()

    filenameMPM = "./output_%s_%s/output.2000.nc" %(tests[1],res)
    if (not os.path.exists(filenameMPM)):
        raise Exception("Missing output file: %s" %(filenameMPM))
    fileMPM = Dataset(filenameMPM,"r")
    iceAreaVertexMPM2 = fileMPM.variables["iceAreaVertex"][0,:]
    iceAreaVertexDiff2 = (iceAreaVertexMPM2 - iceAreaVertexAnalytical2)
    print("MPM Reconstruction Error ", tests[1], "  vertex:  ",
          np.amin(iceAreaVertexDiff2), np.amax(iceAreaVertexDiff2))
    iceAreaCellMPM2 = fileMPM.variables["iceAreaCell"][0,:]
    iceAreaCellDiff2 = (iceAreaCellMPM2 - iceAreaCellAnalytical2)
    print("MPM Reconstruction Error ", tests[1], " cell:  ",
          np.amin(iceAreaCellDiff2), np.amax(iceAreaCellDiff2))
    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(4, 2)

    fig.set_size_inches(5, 8)

    minIceArea = -1.57
    maxIceArea =  1.57
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexAnalytical1, minIceArea, maxIceArea, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Vertex ' + tests[0], '(a)', True)
    minIceArea = 0.0
    maxIceArea = 3.0
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexAnalytical2, minIceArea, maxIceArea, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Vertex ' + tests[1], '(b)', True)
    minIceAreaDiff = -4.0e-6
    maxIceAreaDiff =  4.0e-6
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexDiff1, minIceAreaDiff, maxIceAreaDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM vertex diff', '(c)', True)
    minIceAreaDiff = -5.0e-5
    maxIceAreaDiff =  5.0e-5
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, iceAreaVertexDiff2, minIceAreaDiff, maxIceAreaDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM vertex diff', '(d)', True)

    minIceArea = -1.57
    maxIceArea =  1.57
    plot_subfigure(axes[2,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellAnalytical1, minIceArea, maxIceArea, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Cell ' + tests[0], '(e)', True)
    minIceArea = 0.0
    maxIceArea = 3.0
    plot_subfigure(axes[2,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellAnalytical2, minIceArea, maxIceArea, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Cell ' + tests[1], '(f)', True)
    minIceAreaDiff = -1.0e-3
    maxIceAreaDiff =  1.0e-3
    plot_subfigure(axes[3,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellDiff1, minIceAreaDiff, maxIceAreaDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM cell diff', '(g)', True)
    minIceAreaDiff = -6.0e-4
    maxIceAreaDiff =  6.0e-4
    plot_subfigure(axes[3,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, iceAreaCellDiff2, minIceAreaDiff, maxIceAreaDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM cell diff', '(h)', True)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("reconstruction_map.png",dpi=400)
    plt.savefig("reconstruction_map.pdf")

    plt.clf()
    plt.cla()
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_map()
