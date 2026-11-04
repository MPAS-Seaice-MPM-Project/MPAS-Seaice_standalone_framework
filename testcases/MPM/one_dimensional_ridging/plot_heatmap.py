import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import numpy as np
import matplotlib
import random
import argparse

#-------------------------------------------------------------------------------

def plot_axis(axis,
              fig,
              patches,
              array,
              cmin, cmax,
              tMin, tMax,
              xMin, xMax,
              title,
              units):

    pc = PatchCollection(patches, cmap="jet", edgecolor=None)
    pc.set_array(np.array(array))
    pc.set_clim((cmin, cmax))

    axis.add_collection(pc)

    axis.set_xlim((tMin, tMax))
    axis.set_ylim((xMin, xMax))

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Position (m)")

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(units)

    axis.set_title(title)

    return axis, fig

#-------------------------------------------------------------------------------

def plot_heatmap(gridFilename,
                 filenameIn,
                 filenameOut="ridging.pdf",
                 nSkip=10):

    maxThickness = 4.0

    yVal = 55000.0

    fileGrid = Dataset(gridFilename,"r")

    nCells = len(fileGrid.dimensions["nCells"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:] - 1

    dvEdge = fileGrid.variables["dvEdge"][:]
    cellWidth = dvEdge[0]

    fileGrid.close()


    filein = Dataset(filenameIn,"r")

    nTimes = len(filein.dimensions["Time"])

    secsSinceStartOfSimulation = filein.variables["daysSinceStartOfSim"][:] * 24*3600

    time0 = secsSinceStartOfSimulation[0]
    time1 = secsSinceStartOfSimulation[nSkip-1]
    writeInterval = secsSinceStartOfSimulation[1] - secsSinceStartOfSimulation[0]

    timeWidth = time1 - time0 + writeInterval

    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    tMin =  sys.float_info.max
    tMax = -sys.float_info.max

    iceAreaCell = filein.variables["iceAreaCell"][:]
    iceVolumeCell = filein.variables["iceVolumeCell"][:]

    uVelocity = filein.variables["uVelocity"][:]

    filein.close()


    maxVolume = 0.0
    maxThickness = 0.0

    # element
    elementPatches = []
    iceConcentrations = []
    iceVolumes = []
    iceThicknesses = []
    uVelocities = []

    for iTime in range(0,nTimes,nSkip):

        # time
        time = secsSinceStartOfSimulation[iTime]

        tMin = min(tMin,time-0.5*timeWidth)
        tMax = max(tMax,time+0.5*timeWidth)

        # elements
        for iCell in range(0,nCells):

            if (yCell[iCell] == yVal):

                elementPatches.append(Rectangle((time-0.5*timeWidth,xCell[iCell]-0.5*cellWidth),
                                                 timeWidth, cellWidth))

                iceThickness = 0.0
                if (iceAreaCell[iTime,iCell] > 0.0):
                    iceThickness = iceVolumeCell[iTime,iCell] / iceAreaCell[iTime,iCell]

                uVelocityCell = 0.0
                for iVertexOnCell in range(0,4):
                    iVertex = verticesOnCell[iCell,iVertexOnCell]
                    uVelocityCell += uVelocity[iTime,iVertex]
                uVelocityCell /= 4.0

                iceConcentrations.append(iceAreaCell[iTime,iCell])
                iceVolumes.append(iceVolumeCell[iTime,iCell])
                iceThicknesses.append(iceThickness)
                uVelocities.append(uVelocityCell)

                xMin = min(xMin,xCell[iCell]-0.5*cellWidth)
                xMax = max(xMax,xCell[iCell]+0.5*cellWidth)

                maxVolume    = max(maxVolume,    iceVolumeCell[iTime,iCell])
                maxThickness = max(maxThickness, iceThickness)


    fig, axes = plt.subplots(2,2,figsize=(15,20))

    # ice concentration
    axes[0,0], fig = plot_axis(axes[0,0],
                               fig,
                               elementPatches,
                               iceConcentrations,
                               0.0, 1.0,
                               tMin, tMax,
                               xMin, xMax,
                               "Ice concentration",
                               "(-)")

    # ice volume
    axes[0,1], fig = plot_axis(axes[0,1],
                               fig,
                               elementPatches,
                               iceVolumes,
                               0.0, maxVolume,
                               tMin, tMax,
                               xMin, xMax,
                               "Ice volume",
                               "(m)")

    # ice thickness
    axes[1,0], fig = plot_axis(axes[1,0],
                               fig,
                               elementPatches,
                               iceThicknesses,
                               0.0, maxThickness,
                               tMin, tMax,
                               xMin, xMax,
                               "Ice thickness",
                               "(m)")

    # velocity
    axes[1,1], fig = plot_axis(axes[1,1],
                               fig,
                               elementPatches,
                               uVelocities,
                               None, None,
                               tMin, tMax,
                               xMin, xMax,
                               "Velocity",
                               "(m/s)")

    plt.tight_layout()
    plt.savefig(filenameOut)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)
    parser.add_argument('-i', dest='filenameIn', required=True)
    parser.add_argument('-o', dest='filenameOut', default="ridging.pdf")
    parser.add_argument('-n', dest='nSkip', type=int, default=10)

    args = parser.parse_args()

    plot_heatmap(args.gridFilename,
                 args.filenameIn,
                 args.filenameOut,
                 args.nSkip)
