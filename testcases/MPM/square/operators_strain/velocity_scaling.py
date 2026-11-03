from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

#--------------------------------------------------------

def L2_norm_vertex(numerical, analytical, nVertices, areaTriangle, useVertex):

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (useVertex[iVertex]):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_vertex(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    uVelocityAnalytical = fileIC.variables["uVelocity"][:]
    vVelocityAnalytical = fileIC.variables["vVelocity"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    uVelocity = fileMPAS.variables["uVelocity"][0,:]
    vVelocity = fileMPAS.variables["vVelocity"][0,:]

    useVertex = get_use_vertex(filename)

    normU = L2_norm_vertex(uVelocity, uVelocityAnalytical, nVertices, areaTriangle, useVertex)
    normV = L2_norm_vertex(vVelocity, vVelocityAnalytical, nVertices, areaTriangle, useVertex)

    fileMPAS.close()

    return normU, normV

#--------------------------------------------------------

def get_use_vertex(filenameIn):

    fileIn = Dataset(filenameIn,"r")
    nCells = len(fileIn.dimensions["nCells"])
    nVertices = len(fileIn.dimensions["nVertices"])
    nEdgesOnCell = fileIn.variables["nEdgesOnCell"][:]
    cellsOnCell = fileIn.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1
    verticesOnCell = fileIn.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1
    interiorCell = fileIn.variables["interiorCell"][0,:]
    fileIn.close()

    useVertex = np.ones(nVertices,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    return useVertex

#--------------------------------------------------------

def L2_norm_particle(numerical, analytical, nParticles):

    norm  = 0.0
    denom = 0.0

    for iParticle in range(0,nParticles):

        norm  = norm + math.pow(numerical[iParticle] - analytical[iParticle],2)

        denom = denom + math.pow(analytical[iParticle],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_particle(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    uvAnalyticalMP = fileIC.variables["uvVelMP"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    uvVelMP = fileMPAS.variables["uvVelMP"][:]

    normU = L2_norm_particle(uvVelMP[0,:,0], uvVelAnalyticalMP[:,0], nParticles, latParticle[0,:])
    normV = L2_norm_particle(uvVelMP[0,:,1], uvVelAnalyticalMP[:,1], nParticles, latParticle[0,:])

    fileMPAS.close()

    return normU, normV

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

def velocity_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    velocities = ["u","v"]

    operatorMethods = ["mpmvar","mpmweak","mpmmpm"]

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


    dirname = {"wachspress":"wachspress",
               "pwl":"pwl",
               "weak":"weak",
               "wachspress_avg":"wachspress",
               "pwl_avg":"pwl",
               "weak_avg":"weak",
               "mpmvar":"mpmvar",
               "mpmweak":"mpmweak",
               "mpmmpm":"mpmmpm"}

    legendLabels = {"mpmvar":"MPM Var",
                    "mpmweak":"MPM Weak",
                    "mpmmpm":"MPM MPM"}

    lineColours = {"mpmvar":"black",
                   "mpmweak":"red",
                   "mpmmpm":"blue"}

    markers = {"mpmvar":"o",
               "mpmweak":"x",
               "mpmmpm":"s"}

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    velocityLabels = {"u":r"(a) $u$",
                      "v":r"(b) $v$"}

    fig, axes = plt.subplots(2,1,figsize=(5,10))

    iVelo = 0
    for velocity in velocities:

        #xMin = 6e-3
        #xMax = 1e-2
        #yMin = 4e-4
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

        axes[iVelo].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iVelo].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        for gridType in gridTypes:

            print("Grid type: ", gridType)

            iPlot = 0
            for operatorMethod in operatorMethods:

                print("  Operator Method: ", operatorMethod)

                x = []
                y = []

                for grid in grids[gridType]:

                    print("    Grid: ", grid)

                    filenamePartIC = "particles_%s_%s.nc" %(gridType,grid)
                    filenamesPart = sorted(glob.glob("./output_%s_%s_%s/particles_output*" %(gridType, dirname[operatorMethod], grid)))
                    filenamePart = filenamesPart[-1]
                    filename = "output_%s_%s_%s/output.2000.nc" %(gridType, dirname[operatorMethod], grid)
                    filenameIC = "ic_%s_%s.nc" %(gridType,grid)

                    print("      ", filename, filenameIC)

                    normU, normV = get_norm_vertex(filenameIC, filename)

                    x.append(get_resolution(filename))

                    if (velocity == "u"):
                        y.append(normU)
                    elif (velocity == "v"):
                        y.append(normV)

                if (gridType == "hex"):
                    legendLabel = "%s" %(legendLabels[operatorMethod])
                else:
                    legendLabel = "_nolegend_"

                axes[iVelo].loglog(x, y, marker=markers[operatorMethod], color=lineColours[operatorMethod], ls=lineStyles[gridType], markersize=5.0, label=legendLabel)

                print('resolution',x)
                print('normMP',y)

                iPlot = iPlot + 1

        axes[iVelo].legend(frameon=False, loc=4, fontsize=8, handlelength=4)

        axes[iVelo].set_xlabel("Grid resolution")
        axes[iVelo].set_ylabel(r"$L_2$ error norm")
        axes[iVelo].set_title(velocityLabels[velocity])

        iVelo = iVelo + 1


    plt.tight_layout()
    plt.savefig("velocity_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    velocity_scaling()
