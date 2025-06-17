from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import numpy as np
import os
from average_particle_ice_area_to_cell import average_particle_ice_area_to_cell
import glob

#-------------------------------------------------------------

def plot_subfigure(axis, array, nCells, nEdgesOnCell, verticesOnCell, xCell, yCell, zCell, xVertex, yVertex, zVertex, cmin, cmax, cmap):

    xMin =  1.0e30
    xMax = -1.0e30
    yMin =  1.0e30
    yMax = -1.0e30

    cmap = plt.get_cmap(cmap)

    patches = []
    colors = []

    for iCell in range(0,nCells):

        if (yCell[iCell] > 0.0):

            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append((xVertex[iVertex],zVertex[iVertex]))

            colors.append(array[iCell])

            patches.append(Polygon(vertices))

            xMin = min(xMin,xVertex[iVertex])
            xMax = max(xMax,xVertex[iVertex])

            yMin = min(yMin,zVertex[iVertex])
            yMax = max(yMax,zVertex[iVertex])

    pc = PatchCollection(patches, cmap=cmap)
    pc.set_array(np.array(colors))
    pc.set_clim(cmin, cmax)

    axis.add_collection(pc)

    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)
    axis.set_aspect("equal")

    axis.ticklabel_format(style='plain')
    axis.tick_params(axis='x', \
                     which='both', \
                     bottom=False, \
                     top=False, \
                     labelbottom=False)
    axis.tick_params(axis='y', \
                     which='both', \
                     left=False, \
                     right=False, \
                     labelleft=False)

#-------------------------------------------------------------

def plot_testcase(runtype):

    #nGrids = [2562,10242,40962,163842]
    nGrids = [2562]
    tsteps = [3600, 1800, 900, 450]
    testTypes = ["cosine_bell","slotted_cylinder"]
    iTimes = [0,-1]

    for tstep in tsteps:

        print("tstep: ", tstep)
        filesExist = False
        for testType in testTypes:
            filenamein = "./output_%s_%i_%s/output.2000.nc" %(testType,tstep,runtype)
            if (os.path.exists(filenamein)):
                filesExist = True

        if (filesExist):

            fig, axes = plt.subplots(2,4)

            iTestType = -1
            for testType in testTypes:
                iTestType += 1

                print("  Test type: ", testType)

                # mesh
                filenamein = "./output_%s_%i_%s/output.2000.nc" %(testType,tstep,runtype)

                filein = Dataset(filenamein,"r")

                nCells = len(filein.dimensions["nCells"])

                nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
                verticesOnCell = filein.variables["verticesOnCell"][:]
                xCell = filein.variables["xCell"][:]
                yCell = filein.variables["yCell"][:]
                zCell = filein.variables["zCell"][:]
                xVertex = filein.variables["xVertex"][:]
                yVertex = filein.variables["yVertex"][:]
                zVertex = filein.variables["zVertex"][:]
                verticesOnCell[:] = verticesOnCell[:] - 1

                filein.close()

                # iceArea
                filenameParticlesTemplate = "./output_%s_%i_%s/particles_output*" %(testType,tstep,runtype)
                filenames = sorted(glob.glob(filenameParticlesTemplate))

                iceAreaCell0 = average_particle_ice_area_to_cell(filenames[0], nCells)

                iMethod = -1
                for iTime in iTimes:
                    iMethod +=1

                    print("iTime: ", iTime)

                    iceAreaCell = average_particle_ice_area_to_cell(filenames[iTime], nCells)

                    plot_subfigure(axes[iMethod,iTestType*2],
                                   iceAreaCell,
                                   nCells,
                                   nEdgesOnCell,
                                   verticesOnCell,
                                   xCell,
                                   yCell,
                                   zCell,
                                   xVertex,
                                   yVertex,
                                   zVertex,
                                   0.0,
                                   1.0,
                                   "viridis")

                    iceAreaCellDiff = iceAreaCell - iceAreaCell0

                    if (iMethod != 0):
                        plot_subfigure(axes[iMethod,iTestType*2+1],
                                       iceAreaCellDiff,
                                       nCells,
                                       nEdgesOnCell,
                                       verticesOnCell,
                                       xCell,
                                       yCell,
                                       zCell,
                                       xVertex,
                                       yVertex,
                                       zVertex,
                                       -1.0,
                                       1.0,
                                       "bwr")
                    else:
                        axes[iMethod, iTestType*2+1].axis('off')

            plt.savefig("advection_%6.6i.%s.png" %(tstep,runtype),dpi=300)
            plt.cla()
            plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase("")
