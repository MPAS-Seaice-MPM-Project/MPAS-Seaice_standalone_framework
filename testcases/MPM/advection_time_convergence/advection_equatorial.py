from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import glob
from average_particle_ice_area_to_cell import average_particle_ice_area_to_cell

#---------------------------------------------------------

def get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, array):

    # find initial point
    radius = 6371229.0
    xStart = radius
    yStart = 0.0
    zStart = 0.0

    minDistance = 1e30
    iCellStart = -1
    for iCell in range(0,nCells):

        distance = math.sqrt(math.pow((xCell[iCell] - xStart),2) + \
                             math.pow((yCell[iCell] - yStart),2) + \
                             math.pow((zCell[iCell] - zStart),2))

        if (distance < minDistance):
            minDistance = distance
            iCellStart = iCell


    iCell = iCellStart

    lats = []
    lons = []
    y = []
    deltaLon = 0.0

    while True:

        lats.append(latCell[iCell])
        lons.append(lonCell[iCell])
        y.append(array[iCell])

        iCellNext = -1
        nextSmallestLatitude = 1e30

        for iCellOnCell in range(0,nEdgesOnCell[iCell]):

            iCellNeighbour = cellsOnCell[iCell,iCellOnCell] - 1

            lonNeighbour = lonCell[iCellNeighbour]

            if (lonNeighbour < lonCell[iCell] - math.pi):
                lonNeighbour = lonNeighbour + 2.0 * math.pi
            if (lonNeighbour > lonCell[iCell] + math.pi):
                lonNeighbour = lonNeighbour - 2.0 * math.pi

            if (lonNeighbour > lonCell[iCell]):

                if (math.fabs(latCell[iCellNeighbour]) < nextSmallestLatitude):

                    nextSmallestLatitude = math.fabs(latCell[iCellNeighbour])
                    iCellNext = iCellNeighbour
                    lonNeighbourAll = lonNeighbour

        deltaLon = deltaLon + (lonNeighbourAll - lonCell[iCell])

        if (deltaLon > 2.0 * math.pi):
            break

        iCell = iCellNext

    points = zip(lons,y)
    sorted_points = sorted(points)
    lons = [point[0] for point in sorted_points]
    y    = [point[1] for point in sorted_points]

    return lons, y

#---------------------------------------------------------

def advection_equatorial(runtype):

    #res = "2562"
    #res = "10242"
    #res = "40962"
    #res = "163842"
    res = 3600

    experiment1 = "cosine_bell"
    experiment2 = "slotted_cylinder"

    iTime = -1

    # mesh
    filename = "./output_%s_%i_%s/output.2000.nc" %(experiment1,res,runtype)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]
    lonCell = fileMPAS.variables["lonCell"][:]

    xCell = fileMPAS.variables["xCell"][:]
    yCell = fileMPAS.variables["yCell"][:]
    zCell = fileMPAS.variables["zCell"][:]

    cellsOnCell = fileMPAS.variables["cellsOnCell"][:]
    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]

    fileMPAS.close()

    # cosine bell
    filenamesParticleCB = sorted(glob.glob("./output_%s_%i_%s/particles_output*" %(experiment1,res,runtype)))

    iceArea_CB_Initial = average_particle_ice_area_to_cell(filenamesParticleCB[ 0], nCells)
    iceArea_CB         = average_particle_ice_area_to_cell(filenamesParticleCB[-1], nCells)

    # slotted cylinder
    filenamesParticleSC = sorted(glob.glob("./output_%s_%i_%s/particles_output*" %(experiment2,res,runtype)))

    iceArea_SC_Initial = average_particle_ice_area_to_cell(filenamesParticleSC[ 0], nCells)
    iceArea_SC         = average_particle_ice_area_to_cell(filenamesParticleSC[-1], nCells)

    print(np.amin(iceArea_SC_Initial), np.amax(iceArea_SC_Initial))
    print(np.amin(iceArea_SC),         np.amax(iceArea_SC))

    lons, yiceArea_CB_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB_Initial)
    lons, yiceArea_CB = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB)

    lons, yiceArea_SC_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC_Initial)
    lons, yiceArea_SC = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC)

    cm = 1/2.54  # centimeters in inches
    plt.rcParams["font.family"] = "Times New Roman"
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    linewidth = 1.0

    fig, axes = plt.subplots(1, 2, figsize=(15*cm,7*cm))

    # area
    axes[0].plot(lons, yiceArea_CB_Initial, linewidth=linewidth, color="black", ls="-")
    axes[0].plot(lons, yiceArea_CB,      linewidth=linewidth, color="red", ls="--")

    axes[0].legend(["Initial","MPM",""],frameon=False,fontsize=8)
    axes[0].set_xlabel("Longitude (radians)")
    axes[0].set_ylabel("Equatorial ice concentration")

    axes[0].set_xlim([0.5, math.pi-0.5])
    axes[0].set_title("(a) Cosine bell", loc='left')

    axes[1].plot(lons, yiceArea_SC_Initial, linewidth=linewidth, color="black", ls="-")
    axes[1].plot(lons, yiceArea_SC,      linewidth=linewidth, color="red", ls="--")

    axes[1].set_xlabel("Longitude (radians)")
    axes[1].set_xlim([0.5, math.pi-0.5])
    axes[1].set_title("(b) Slotted cylinder", loc='left')

    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("advection_equatorial_%s.png" %(runtype),dpi=300)
    plt.savefig("advection_equatorial_%s.eps" %(runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    advection_equatorial()
