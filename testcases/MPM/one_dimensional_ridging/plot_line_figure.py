from netCDF4 import Dataset
from math import fabs
import matplotlib.pyplot as plt
from netCDF4 import stringtochar
import numpy as np

#-------------------------------------------------------------------------------

def plot_line_figure():

    title = None

    outputFilename = "ridging.png"

    outputDirs = ["./output_orig",
                  "./output_mpm.evp",
                  "./output_mpm.decoh"]

    dirLabel = ["Orig",
                "MPM evp",
                "MPM decohesive"]

    colors = ["black",
              "red",
              "blue",
              "green",
              "orange",
              "pink"]


    times0 = ["2000-01-02_00:00:00",
              "2000-01-06_00:00:00",
              "2000-01-31_00:00:00"]

    timeLabels = ["1 day",
                  "5 days",
                  "30 days"]

    lineStyles = ["dashdot",
                  "dashed",
                  "solid"]


    yVertex0 = 50000
    yCell0   = 55000

    fileMesh = Dataset("grid.nc","r")

    nVertices = len(fileMesh.dimensions["nVertices"])
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]

    nCells = len(fileMesh.dimensions["nCells"])
    xCell = fileMesh.variables["xCell"][:]
    yCell = fileMesh.variables["yCell"][:]

    fileMesh.close()

    vertices = []
    xvs = []
    for iVertex in range(0,nVertices):
        if (fabs(yVertex[iVertex]-yVertex0) < 1e-11):
            vertices.append(iVertex)
            xvs.append(xVertex[iVertex])

    cells = []
    xcs = []
    for iCell in range(0,nCells):
        if (fabs(yCell[iCell]-yCell0) < 1e-11):
            cells.append(iCell)
            xcs.append(xCell[iCell])

    fig, axes = plt.subplots(4,2,figsize=(8,10))

    for iDir, outputDir in enumerate(outputDirs):

        filein = Dataset("%s/output.2000.nc" %(outputDir),"r")
        xtime = filein.variables["xtime"][:]
        times = []
        for iTime in range(0,len(xtime)):
            xtimeStr = b''.join(xtime[iTime]).decode('utf-8').rstrip("\x00")
            times.append(xtimeStr)

        uVelocity = filein.variables["uVelocity"][:]
        icePressure = filein.variables["icePressure"][:]
        iceVolumeCell = filein.variables["iceVolumeCell"][:]
        iceAreaCell = filein.variables["iceAreaCell"][:]
        stressDivergenceU = filein.variables["stressDivergenceU"][:]
        stressDivergenceV = filein.variables["stressDivergenceV"][:]
        divergence = filein.variables["divergence"][:]
        shear      = filein.variables["shear"][:]

        filein.close()


        for iLabel, time0 in enumerate(times0):
            iTime = -1
            for jTime in range(0,len(times)):
                if (time0 == times[jTime]):
                    iTime = jTime

            us = []
            for iVertex in vertices:
                us.append(uVelocity[iTime,iVertex])
            axes[0,0].plot(xvs,us,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            ps = []
            for iCell in cells:
                ps.append(icePressure[iTime,iCell])
            axes[0,1].plot(xcs,ps,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            vs = []
            for iCell in cells:
                vs.append(iceVolumeCell[iTime,iCell])
            axes[1,0].plot(xcs,vs,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            cs = []
            for iCell in cells:
                cs.append(iceAreaCell[iTime,iCell])
            axes[1,1].plot(xcs,cs,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            dus = []
            for iVertex in vertices:
                dus.append(stressDivergenceU[iTime,iVertex])
            axes[2,0].plot(xvs,dus,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            dvs = []
            for iVertex in vertices:
                dvs.append(stressDivergenceV[iTime,iVertex])
            axes[2,1].plot(xvs,dvs,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            ds = []
            for iCell in cells:
                ds.append(divergence[iTime,iCell])
            axes[3,0].plot(xcs, ds,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

            ss = []
            for iCell in cells:
                ss.append(shear[iTime,iCell])
            axes[3,1].plot(xcs,ss,label=timeLabels[iLabel]+" "+dirLabel[iDir],linestyle=lineStyles[iLabel],color=colors[iDir])

    axes[0,0].set_ylabel("uVelocity")
    axes[0,1].set_ylabel("icePressure")
    axes[1,0].set_ylabel("iceVolumeCell")
    axes[1,1].set_ylabel("iceAreaCell")
    axes[2,0].set_ylabel("stressDivergenceU")
    axes[2,1].set_ylabel("stressDivergenceV")
    axes[3,0].set_ylabel("divergence")
    axes[3,1].set_ylabel("shear")

    axes[2,0].set_ylim((-0.3,0.1))

    plt.suptitle(title)

    plt.legend(fontsize=6)
    plt.tight_layout()
    plt.savefig(outputFilename,dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_line_figure()
