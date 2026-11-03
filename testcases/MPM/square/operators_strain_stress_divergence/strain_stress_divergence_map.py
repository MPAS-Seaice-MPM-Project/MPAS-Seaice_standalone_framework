from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

degreesToRadians = math.pi / 180.0

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertexx, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        polygonVertices = []

        useVertex = False
        for iCellOnVertex in range(0,vertexDegree):

            iCell = cellsOnVertex[iVertex,iCellOnVertex]

            polygonVertices.append((xCell[iCell],yCell[iCell]))

            if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                yCell[iCell] >= minY and yCell[iCell] <= maxY):
                useVertex = True

        if (useVertex and (useVertexx[iVertex] == 1)):
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

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
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

def strain_stress_divergence_map():

    # grid hex
    fileGrid = Dataset("grid_hex_0082x0094.nc","r")

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

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    fileGrid.close()

    # ic hex
    fileIC = Dataset("ic_hex_0082x0094.nc","r")

    uVelocityAnalytical = fileIC.variables["uVelocity"][:]
    vVelocityAnalytical = fileIC.variables["vVelocity"][:]

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Velocity: ",
          np.amin(uVelocityAnalytical), np.amax(uVelocityAnalytical),
          np.amin(vVelocityAnalytical), np.amax(vVelocityAnalytical))
    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_hex_wachspress_0082x0094/output.2000.nc","r")

    interiorCell = fileWach.variables["interiorCell"][0,:]

    uVelocityWachs = fileWach.variables["uVelocity"][0,:]
    vVelocityWachs = fileWach.variables["vVelocity"][0,:]

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    uVelocityWachsDiff = (uVelocityWachs - uVelocityAnalytical)
    vVelocityWachsDiff = (vVelocityWachs - vVelocityAnalytical)

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs velocity: ",
          np.amin(uVelocityWachsDiff), np.amax(uVelocityWachsDiff),
          np.amin(vVelocityWachsDiff), np.amax(vVelocityWachsDiff))
    print("Wachs stressDiv: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    useVertex = np.ones(nVertices,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    # Weak
    fileWeak = Dataset("./output_hex_weak_0082x0094/output.2000.nc","r")

    uVelocityWeak = fileWeak.variables["uVelocity"][0,:]
    vVelocityWeak = fileWeak.variables["vVelocity"][0,:]

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    uVelocityWeakDiff = (uVelocityWeak - uVelocityAnalytical)
    vVelocityWeakDiff = (vVelocityWeak - vVelocityAnalytical)

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak Velocity:  ",
          np.amin(uVelocityWeakDiff), np.amax(uVelocityWeakDiff),
          np.amin(vVelocityWeakDiff), np.amax(vVelocityWeakDiff))

    print("Weak stressDiv:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()

   # VAR MPM
    fileMPM = Dataset("./output_hex_mpmvar_0082x0094/output.2000.nc","r")

    uVelocityVarMPM = fileMPM.variables["uVelocity"][0,:]
    vVelocityVarMPM = fileMPM.variables["vVelocity"][0,:]

    stressDivergenceUVarMPM = fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVVarMPM = fileMPM.variables["stressDivergenceV"][0,:]

    uVelocityVarMPMDiff = (uVelocityVarMPM - uVelocityAnalytical)
    vVelocityVarMPMDiff = (vVelocityVarMPM - vVelocityAnalytical)

    stressDivergenceUVarMPMDiff = (stressDivergenceUVarMPM - stressDivergenceUAnalytical)
    stressDivergenceVVarMPMDiff = (stressDivergenceVVarMPM - stressDivergenceVAnalytical)

    print("MPM VAR velocity:   ",
          np.amin(uVelocityVarMPMDiff), np.amax(uVelocityVarMPMDiff),
          np.amin(vVelocityVarMPMDiff), np.amax(vVelocityVarMPMDiff))

    print("MPM VAR stressDiv:   ",
          np.amin(stressDivergenceUVarMPMDiff), np.amax(stressDivergenceUVarMPMDiff),
          np.amin(stressDivergenceVVarMPMDiff), np.amax(stressDivergenceVVarMPMDiff))

    fileMPM.close()

   # Weak MPM
    fileMPM = Dataset("./output_hex_mpmweak_0082x0094/output.2000.nc","r")

    uVelocityWeakMPM = fileMPM.variables["uVelocity"][0,:]
    vVelocityWeakMPM = fileMPM.variables["vVelocity"][0,:]

    stressDivergenceUWeakMPM = fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeakMPM = fileMPM.variables["stressDivergenceV"][0,:]

    uVelocityWeakMPMDiff = (uVelocityWeakMPM - uVelocityAnalytical)
    vVelocityWeakMPMDiff = (vVelocityWeakMPM - vVelocityAnalytical)

    stressDivergenceUWeakMPMDiff = (stressDivergenceUWeakMPM - stressDivergenceUAnalytical)
    stressDivergenceVWeakMPMDiff = (stressDivergenceVWeakMPM - stressDivergenceVAnalytical)

    print("MPM WEAK velocity:   ",
          np.amin(uVelocityWeakMPMDiff), np.amax(uVelocityWeakMPMDiff),
          np.amin(vVelocityWeakMPMDiff), np.amax(vVelocityWeakMPMDiff))

    print("MPM WEAK stressDiv:   ",
          np.amin(stressDivergenceUWeakMPMDiff), np.amax(stressDivergenceUWeakMPMDiff),
          np.amin(stressDivergenceVWeakMPMDiff), np.amax(stressDivergenceVWeakMPMDiff))

    fileMPM.close()

   # MPM MPM
    fileMPM = Dataset("./output_hex_mpmmpm_0082x0094/output.2000.nc","r")

    uVelocityMPMMPM= fileMPM.variables["uVelocity"][0,:]
    vVelocityMPMMPM = fileMPM.variables["vVelocity"][0,:]

    stressDivergenceUMPMMPM= fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVMPMMPM = fileMPM.variables["stressDivergenceV"][0,:]

    uVelocityMPMMPMDiff = (uVelocityMPMMPM - uVelocityAnalytical)
    vVelocityMPMMPMDiff = (vVelocityMPMMPM - vVelocityAnalytical)

    stressDivergenceUMPMMPMDiff = (stressDivergenceUMPMMPM - stressDivergenceUAnalytical)
    stressDivergenceVMPMMPMDiff = (stressDivergenceVMPMMPM - stressDivergenceVAnalytical)

    print("MPM MPM velocity:   ",
          np.amin(uVelocityMPMMPMDiff), np.amax(uVelocityMPMMPMDiff),
          np.amin(vVelocityMPMMPMDiff), np.amax(vVelocityMPMMPMDiff))

    print("MPM MPM stressDiv:   ",
          np.amin(stressDivergenceUMPMMPMDiff), np.amax(stressDivergenceUMPMMPMDiff),
          np.amin(stressDivergenceVMPMMPMDiff), np.amax(stressDivergenceVMPMMPMDiff))

    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    minVelDiff = -1.e-1
    maxVelDiff =  1.e-1

    minStressDiv = -750.0
    maxStressDiv =  750.0

    minStressDivDiff = -40.0
    maxStressDivDiff =  40.0


    fig, axes = plt.subplots(6, 4)

    fig.set_size_inches(9, 6.75)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityAnalytical, minVelocity, maxVelocity, xMin, xMax, yMin, yMax, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityAnalytical, minVelocity, maxVelocity, xMin, xMax, yMin, yMax, \
                   False, False, r'$v^\prime$', '(b)', False)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, xMin, xMax, yMin, yMax, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Analytical', r'(c)$\times20$', False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, xMin, xMax, yMin, yMax, \
                   False, False, r'$(\nabla \cdot \sigma)_{v^\prime}$ Analytical', r'(d)$\times20$', True)

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityWachsDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($u^\prime$ direction)', '(e)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityWachsDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($v^\prime$ direction)', '(f)', True)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($u^\prime$ direction)', '(g)', False)
    plot_subfigure(axes[1,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($v^\prime$ direction)', '(h)', True)

    plot_subfigure(axes[2,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityWeakDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($u^\prime$ direction)', '(i)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityWeakDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($v^\prime$ direction)', '(j)', True)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($u^\prime$ direction)', '(k)', False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($v^\prime$ direction)', '(l)', True)

    plot_subfigure(axes[3,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityMPMMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'MPM MPM ($u^\prime$ direction)', '(m)', False)
    plot_subfigure(axes[3,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityMPMMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'MPM MPM ($v^\prime$ direction)', '(n)', True)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUMPMMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'MPM MPM ($u^\prime$ direction)', '(o)', False)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVMPMMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'MPM MPM ($v^\prime$ direction)', '(p)', True)

    plot_subfigure(axes[4,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityWeakMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak MPM. ($u^\prime$ direction)', '(q)', False)
    plot_subfigure(axes[4,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityWeakMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak MPM. ($v^\prime$ direction)', '(r)', True)
    plot_subfigure(axes[4,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWeakMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak MPM. ($u^\prime$ direction)', '(s)', False)
    plot_subfigure(axes[4,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWeakMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak MPM. ($v^\prime$ direction)', '(t)', True)

    plot_subfigure(axes[5,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocityVarMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Var MPM. ($u^\prime$ direction)', '(u)', False)
    plot_subfigure(axes[5,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocityVarMPMDiff, minVelDiff, maxVelDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Var MPM. ($v^\prime$ direction)', '(v)', True)
    plot_subfigure(axes[5,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUVarMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Var MPM. ($u^\prime$ direction)', '(w)', False)
    plot_subfigure(axes[5,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVVarMPMDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Var MPM. ($v^\prime$ direction)', '(x)', True)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_stress_divergence_map.png",dpi=400)
    #plt.savefig("strain_stress_divergence_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_map()
