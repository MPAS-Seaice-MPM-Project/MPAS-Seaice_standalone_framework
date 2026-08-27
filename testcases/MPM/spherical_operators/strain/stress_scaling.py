from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np

#--------------------------------------------------------

def L2_norm_vertex(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_vertex(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    stress11VertexAnalytical = fileIC.variables["stress11VertexAnalytical"][:]
    stress22VertexAnalytical = fileIC.variables["stress22VertexAnalytical"][:]
    stress12VertexAnalytical = fileIC.variables["stress12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]
    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    stress11varAvg = fileMPAS.variables["stress11varAvgVertex"][0,:]
    stress22varAvg = fileMPAS.variables["stress22varAvgVertex"][0,:]
    stress12varAvg = fileMPAS.variables["stress12varAvgVertex"][0,:]

    fileMPAS.close()

    normE11 = L2_norm_vertex(stress11varAvg, stress11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE22 = L2_norm_vertex(stress22varAvg, stress22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE12 = L2_norm_vertex(stress12varAvg, stress12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (math.fabs(latCell[iCell]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_cell(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    stress11CellAnalytical = fileIC.variables["stress11CellAnalytical"][:]
    stress22CellAnalytical = fileIC.variables["stress22CellAnalytical"][:]
    stress12CellAnalytical = fileIC.variables["stress12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]

    areaCell = fileMPAS.variables["areaCell"][:]

    stress11weak = fileMPAS.variables["stress11weak"][0,:]
    stress22weak = fileMPAS.variables["stress22weak"][0,:]
    stress12weak = fileMPAS.variables["stress12weak"][0,:]

    normE11 = L2_norm_cell(stress11weak, stress11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE22 = L2_norm_cell(stress22weak, stress22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE12 = L2_norm_cell(stress12weak, stress12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_particle(numerical, analytical, nParticles, lat, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iParticle in range(0,nParticles):

        if (math.fabs(lat[iParticle]) > latitudeLimit * degreesToRadians):

            norm  = norm + math.pow(numerical[iParticle] - analytical[iParticle],2)

            denom = denom + math.pow(analytical[iParticle],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_particle(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    stressAnalyticalMP = fileIC.variables["stressAnalyticalMP"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    latParticle = fileMPAS.variables["latCellMP"][:]

    stressTimesThicknessMP = fileMPAS.variables["stressTimesThicknessMP"][:]

    normE11 = L2_norm_particle(stressTimesThicknessMP[0,:,0],
                               stressAnalyticalMP[:,0], nParticles, latParticle[0,:], latitudeLimit)
    normE22 = L2_norm_particle(stressTimesThicknessMP[0,:,1],
                               stressAnalyticalMP[:,1], nParticles, latParticle[0,:], latitudeLimit)
    normE12 = L2_norm_particle(stressTimesThicknessMP[0,:,2], 
                               stressAnalyticalMP[:,2], nParticles, latParticle[0,:], latitudeLimit)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def get_resolution(filename, latitudeLimit):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = math.pi / 180.0

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        if (math.fabs(latEdge[iEdge]) > latitudeLimit * degreesToRadians):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def scaling_lines(axis, xMin, xMax, yMin):

    # linear scaling
    scale = yMin / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    axis.loglog([xMin, xMax], [scaleMinLin,  scaleMaxLin],  linestyle=':', color='k')
    axis.loglog([xMin, xMax], [scaleMinQuad, scaleMaxQuad], linestyle=':', color='k')

#--------------------------------------------------------

def stress_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    stresses = ["stress11","stress22","stress12"]

    resolutions = [2562,10242,40962,163842]

    methods = ["mpmvar","mpmweak","mpm"]

    lineColours = ["black","grey","red","blue","green"]

    lineStyles = ["solid","dashed","dashdot"]

    markers = ["o","^","*"]

    labels = [r'$\sigma_{11}$', r'$\sigma_{22}$', r'$\sigma_{12}$']

    latitudeLimit = 20.0

    fig, axes = plt.subplots(figsize=(4,3))

    iStress = 0
    for stress in stresses:

        xMin = 4e-2
        xMax = 8e-2

        # linear scaling
        scale = 1e-3 / math.pow(xMin,1)
        scaleMinLin = math.pow(xMin,1) * scale
        scaleMaxLin = math.pow(xMax,1) * scale

        # quadratic scaling
        scale = 1e-3 / math.pow(xMin,2)
        scaleMinQuad = math.pow(xMin,2) * scale
        scaleMaxQuad = math.pow(xMax,2) * scale

        iMethod = 0
        for method in methods:

            axes.loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
            axes.loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

            x = []
            y = []

            for resolution in resolutions:

                if (method == "mpmvar"):
                   filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                   filenameIC = "./ic_%i.nc" %(resolution)
                   print(filename, filenameIC)
                   normE11, normE22, normE12 = get_norm_vertex(filenameIC, filename, latitudeLimit)
                elif (method == "mpmweak"):
                   filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                   filenameIC = "./ic_%i.nc" %(resolution)
                   print(filename, filenameIC)
                   normE11, normE22, normE12 = get_norm_cell(filenameIC, filename, latitudeLimit)
                else:
                   filenames = sorted(glob.glob("./output_%s_%i/particles_output*" %(method,resolution)))
                   filename = filenames[-1]
                   filenameIC = "particles_%s.nc" %resolution
                   print(filename, filenameIC)
                   normE11, normE22, normE12 = get_norm_particle(filenameIC, filename, latitudeLimit)
                   filename = "./output_%s_%i/output.2000.nc" %(method,resolution)

                x.append(get_resolution(filename, latitudeLimit))

                if (stress == "stress11"):
                    y.append(normE11)
                elif (stress == "stress22"):
                    y.append(normE22)
                elif (stress == "stress12"):
                    y.append(normE12)

            axes.loglog(x,y, marker=markers[iMethod], color=lineColours[iStress], ls=lineStyles[iMethod], markersize=5.0, label=labels[iStress]+" "+method)


            axes.legend(frameon=True, loc="upper left", fontsize=8, handlelength=4, bbox_to_anchor=(1, 1))

            axes.set_xlabel("Grid resolution")
            axes.set_ylabel(r"$L_2$ error norm")

            iMethod = iMethod + 1
        iStress = iStress + 1

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_scaling()
