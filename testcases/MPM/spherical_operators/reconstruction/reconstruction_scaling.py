from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

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

def get_norm_vertex(filenameIC, filename, filenameMesh, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    iceAreaVertexAnalytical = fileIC.variables["iceAreaVertex"][:]

    fileIC.close()

    fileMesh = Dataset(filenameMesh, "r")

    nVertices = len(fileMesh.dimensions["nVertices"])

    latVertex = fileMesh.variables["latVertex"][:]

    areaTriangle = fileMesh.variables["areaTriangle"][:]

    fileMesh.close()

    fileMPAS = Dataset(filename, "r")

    iceAreaVertex = fileMPAS.variables["iceAreaVertex"][0,:]

    norm = L2_norm_vertex(iceAreaVertex, iceAreaVertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    fileMPAS.close()

    return norm

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

def get_norm_cell(filenameIC, filename, filenameMesh, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    iceAreaCellAnalytical = fileIC.variables["iceAreaCell"][:]

    fileIC.close()

    fileMesh = Dataset(filenameMesh, "r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]

    areaCell = fileMesh.variables["areaCell"][:]

    fileMesh.close()

    fileMPAS = Dataset(filename, "r")

    iceAreaCell = fileMPAS.variables["iceAreaCell"][0,:]

    norm = L2_norm_cell(iceAreaCell, iceAreaCellAnalytical, nCells, latCell, areaCell, latitudeLimit)

    fileMPAS.close()

    return norm

#--------------------------------------------------------

def get_resolution(filenameMesh, latitudeLimit):

    fileMPAS = Dataset(filenameMesh, "r")

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

def reconstruction_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    tests = ["1", "x", "y", "z", "lat", "lon", "nonlin"]
    #tests = ["x","nonlin"]

    resolutions = [2562,10242,40962,163842]

    lineColours = ["black","black","gray","gray","silver","silver","blue","blue","green","green","cyan","cyan","red","red"]
    lineStyles  = ["-","--","-","--","-","--","-","--","-","--","-","--","-","--"]

    latitudeLimit = 20.0

    xMin = 6e-3
    xMax = 1e-1

    # quadratic scaling
    scale = 1.0e-3 / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    # linear scaling
    scale = 1.0e-3 / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale

    fig, axes = plt.subplots(2, 1, figsize=(4,6))

    axes[0].loglog([xMin, xMax],[scaleMinLin,scaleMaxLin],linestyle=':', color='k')
    axes[0].loglog([xMin, xMax],[scaleMinQuad,scaleMaxQuad],linestyle=':', color='k')
    axes[1].loglog([xMin, xMax],[scaleMinLin,scaleMaxLin],linestyle=':', color='k')
    axes[1].loglog([xMin, xMax],[scaleMinQuad,scaleMaxQuad],linestyle=':', color='k')

    iPlot = 0

    for test in tests:

            x = []
            yVertex = []
            yCell = []

            for resolution in resolutions:

                filenameMesh = "./grid.%i.nc" %(resolution)
                filename = "./output_%s_%i/output.2000.nc" %(test,resolution)
                filenameIC = "./ic_%s_%i.nc" %(test,resolution)

                print(filename, filenameIC)
                if (not os.path.exists(filename)):
                    raise Exception("Missing output file: %s" %(filename))
                if (not os.path.exists(filenameIC)):
                    raise Exception("Missing IC file: %s" %(filenameIC))

                norm = get_norm_vertex(filenameIC, filename, filenameMesh, latitudeLimit)
                norm = get_norm_cell(filenameIC, filename, filenameMesh, latitudeLimit)

                x.append(get_resolution(filenameMesh, latitudeLimit))
                yVertex.append(norm)
                yCell.append(norm)

            axes[0].loglog(x,yVertex, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)
            print('iceAreaVertex', x, yVertex)
            axes[1].loglog(x,yCell, marker='^', color=lineColours[iPlot+1], ls=lineStyles[iPlot+1], markersize=5.0)
            print('iceAreaCell', x, yCell)

            iPlot = iPlot + 2

    legendLabels = ["linear scaling", "quadratic scaling", "1", "x", "y", "z", "lat", "lon", "nonlin"]
    axes[0].legend(legendLabels, frameon=False, loc=3, fontsize=6, handlelength=4)
    axes[1].legend(legendLabels, frameon=False, loc=3, fontsize=6, handlelength=4)
    axes[0].set_title("Vertex Reconstruction")
    axes[1].set_title("Cell Reconstruction")
    axes[0].set_xlabel("Grid resolution")
    axes[1].set_xlabel("Grid resolution")
    axes[0].set_ylabel(r"$L_2$ error norm")
    axes[1].set_ylabel(r"$L_2$ error norm")
    axes[0].set_xlim([xMin, xMax])
    axes[1].set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("reconstruction_scaling.png",dpi=400)
    plt.savefig("reconstruction_scaling.pdf")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_scaling()
