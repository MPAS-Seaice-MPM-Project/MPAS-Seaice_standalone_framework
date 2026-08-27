from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
from strain_map import get_mpas_patch_collection
from strain_map import plot_subfigure
from strain_map import scatter_plot

degreesToRadians = math.pi / 180.0

def stress_map():

    # grid
    fileGrid = Dataset("grid.40962.nc","r")

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
    fileIC = Dataset("ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    stress11VertexAnalytical = fileIC.variables["stress11VertexAnalytical"][:]
    stress22VertexAnalytical = fileIC.variables["stress22VertexAnalytical"][:]
    stress12VertexAnalytical = fileIC.variables["stress12VertexAnalytical"][:]

    stress11CellAnalytical = fileIC.variables["stress11CellAnalytical"][:]
    stress22CellAnalytical = fileIC.variables["stress22CellAnalytical"][:]
    stress12CellAnalytical = fileIC.variables["stress12CellAnalytical"][:]

    print("StressVertex Analytical: ",
          np.amin(stress11VertexAnalytical), np.amax(stress11VertexAnalytical),
          np.amin(stress22VertexAnalytical), np.amax(stress22VertexAnalytical),
          np.amin(stress12VertexAnalytical), np.amax(stress12VertexAnalytical))

    print("StressCell Analytical: ",
          np.amin(stress11CellAnalytical), np.amax(stress11CellAnalytical),
          np.amin(stress22CellAnalytical), np.amax(stress22CellAnalytical),
          np.amin(stress12CellAnalytical), np.amax(stress12CellAnalytical))

    fileIC.close()

    # mpm

    fileIC = Dataset("particles_40962.nc","r")

    stress11AnalyticalMP = fileIC.variables["stressAnalyticalMP"][:,0]
    stress22AnalyticalMP = fileIC.variables["stressAnalyticalMP"][:,1]
    stress12AnalyticalMP = fileIC.variables["stressAnalyticalMP"][:,2]
    posnMP = fileIC.variables["posnMP"][:,:]

    fileIC.close()

    fileMPM = Dataset("./output_mpm_40962/particles_output.2000-01-01_01.00.00.nc","r")

    stress11 = fileMPM.variables["stressTimesThicknessMP"][0,:,0]
    stress22 = fileMPM.variables["stressTimesThicknessMP"][0,:,1]
    stress12 = fileMPM.variables["stressTimesThicknessMP"][0,:,2]

    stress11Diff = stress11 - stress11AnalyticalMP
    stress22Diff = stress22 - stress22AnalyticalMP
    stress12Diff = stress12 - stress12AnalyticalMP

    print("MPM: ",
          np.amin(stress11Diff), np.amax(stress11Diff),
          np.amin(stress22Diff), np.amax(stress22Diff),
          np.amin(stress12Diff), np.amax(stress12Diff))

    fileMPM.close()

    # mpmvar

    fileMPMvar = Dataset("./output_mpmvar_40962/output.2000.nc")

    stress11Vertex = fileMPMvar.variables["stress11varAvgVertex"][0,:]
    stress22Vertex = fileMPMvar.variables["stress22varAvgVertex"][0,:]
    stress12Vertex = fileMPMvar.variables["stress12varAvgVertex"][0,:]

    fileMPMvar.close()

    stress11VertexDiff = stress11Vertex - stress11VertexAnalytical
    stress22VertexDiff = stress22Vertex - stress22VertexAnalytical
    stress12VertexDiff = stress12Vertex - stress12VertexAnalytical

    print("MPMvar: ",
          np.amin(stress11VertexDiff), np.amax(stress11VertexDiff),
          np.amin(stress22VertexDiff), np.amax(stress22VertexDiff),
          np.amin(stress12VertexDiff), np.amax(stress12VertexDiff))

    # mpmweak

    fileMPMweak = Dataset("./output_mpmweak_40962/output.2000.nc")

    stress11weak = fileMPMweak.variables["stress11weak"][0,:]
    stress22weak = fileMPMweak.variables["stress22weak"][0,:]
    stress12weak = fileMPMweak.variables["stress12weak"][0,:]

    fileMPMweak.close()

    stress11weakDiff = stress11weak - stress11CellAnalytical
    stress22weakDiff = stress22weak - stress22CellAnalytical
    stress12weakDiff = stress12weak - stress12CellAnalytical

    print("MPMeak: ",
          np.amin(stress11weakDiff), np.amax(stress11weakDiff),
          np.amin(stress22weakDiff), np.amax(stress22weakDiff),
          np.amin(stress12weakDiff), np.amax(stress12weakDiff))

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(6, 6)
    fig.set_size_inches(12, 8)

    minStress = -3.5
    maxStress =  3.5

    minDiff = -0.10
    maxDiff =  0.10

    # Velocities
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', True)
    axes[0,2].axis('off')
    axes[0,3].axis('off')
    axes[0,4].axis('off')
    axes[0,5].axis('off')

    # Analytical stresss
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress11VertexAnalytical, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                    False, False, r'$\sigma_{11}$', '(c)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress11VertexAnalytical, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{11}$', '(d)', False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress22VertexAnalytical, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\sigma_{22}$', '(e)', False)
    plot_subfigure(axes[1,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress22VertexAnalytical, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{22}$', '(f)', False)
    plot_subfigure(axes[1,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress12VertexAnalytical, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\sigma_{12}$', '(g)', False)
    plot_subfigure(axes[1,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress12VertexAnalytical, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{12}$', '(h)', True)

    # MPM
    scatter_plot(axes[2,0], fig, posnMP, stress11, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\sigma_{11}$ MPM.', '(i)', False)
    scatter_plot(axes[2,1], fig, posnMP, stress11, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{11}$ MPM.', '(j)', False)
    scatter_plot(axes[2,2], fig, posnMP, stress22, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\sigma_{22}$ MPM.', '(k)', False)
    scatter_plot(axes[2,3], fig, posnMP, stress22, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{22}$ MPM.', '(l)', False)
    scatter_plot(axes[2,4], fig, posnMP, stress12, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\sigma_{12}$ MPM.', '(m)', False)
    scatter_plot(axes[2,5], fig, posnMP, stress12, minStress, maxStress, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\sigma_{12}$ MPM.', '(n)', True)

    scatter_plot(axes[3,0], fig, posnMP, stress11Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{11}$ MPM.', '(i)', False)
    scatter_plot(axes[3,1], fig, posnMP, stress11Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{11}$ MPM.', '(j)', False)
    scatter_plot(axes[3,2], fig, posnMP, stress22Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{22}$ MPM.', '(k)', False)
    scatter_plot(axes[3,3], fig, posnMP, stress22Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{22}$ MPM.', '(l)', False)
    scatter_plot(axes[3,4], fig, posnMP, stress12Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{12}$ MPM.', '(m)', False)
    scatter_plot(axes[3,5], fig, posnMP, stress12Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{12}$ MPM.', '(n)', True)
    #MPM var
    plot_subfigure(axes[4,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress11VertexDiff, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{11}$ MPMvar', '(o)', False)
    plot_subfigure(axes[4,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress11VertexDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{11}$ MPMvar', '(p)', False)
    plot_subfigure(axes[4,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress22VertexDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{22}$ MPMvar', '(q)', False)
    plot_subfigure(axes[4,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress22VertexDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{22}$ MPMvar', '(r)', False)
    plot_subfigure(axes[4,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress12VertexDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{12}$ MPMvar', '(s)', False)
    plot_subfigure(axes[4,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stress12VertexDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{12}$ MPMvar', '(t)', True)

    # MPM weak
    plot_subfigure(axes[5,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress11weakDiff, minStress, maxStress, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{11}$ MPMweak', '(o)', False)
    plot_subfigure(axes[5,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress11weakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{11}$ MPMweak', '(p)', False)
    plot_subfigure(axes[5,2], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress22weakDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{22}$ MPMweak', '(q)', False)
    plot_subfigure(axes[5,3], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress22weakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{22}$ MPMweak', '(r)', False)
    plot_subfigure(axes[5,4], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress12weakDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\sigma_{12}$ MPMweak', '(s)', False)
    plot_subfigure(axes[5,5], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, stress12weakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\sigma_{12}$ MPMweak', '(t)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    stress_map()
