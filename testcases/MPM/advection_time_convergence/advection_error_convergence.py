import sys
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob
from average_particle_ice_area_to_cell import average_particle_ice_area_to_cell
from check_particle_positions_start_end import check_particle_positions_start_end
import argparse

#--------------------------------------------------------

def L2_norm(numerical, analytical, nPoints, areaCell):

    norm  = 0.0
    denom = 0.0

    for iPoint in range(0,nPoints):

        norm  += areaCell[iPoint] * math.pow(numerical[iPoint] - analytical[iPoint],2)

        denom += areaCell[iPoint] * math.pow(analytical[iPoint],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_area(filenameMesh, filenameParticlesTemplate):

    fileMPAS = Dataset(filenameMesh, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    areaCell = fileMPAS.variables["areaCell"][:]

    fileMPAS.close()

    filenamesParticle = sorted(glob.glob(filenameParticlesTemplate))

    iceAreaCellInitial = average_particle_ice_area_to_cell(filenamesParticle[ 0], nCells)
    iceAreaCellFinal   = average_particle_ice_area_to_cell(filenamesParticle[-1], nCells)

    norm = L2_norm(iceAreaCellFinal, iceAreaCellInitial, nCells, areaCell)

    return norm

#--------------------------------------------------------

def get_resolution(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    # mean distance between cells
    dcEdge = fileMPAS.variables["dcEdge"][:]

    resolution = 0.0
    for iEdge in range(0,nEdges):
        resolution = resolution + dcEdge[iEdge]

    resolution = resolution / float(nEdges)

    fileMPAS.close()

    return resolution / 1000.0

#--------------------------------------------------------

def get_grid_size(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    return nCells

#--------------------------------------------------------

def advection_error_convergence(runtype):

    resolutions = [3600, 1800, 900, 450]

    experiments = ["cosine_bell","slotted_cylinder"]

    legendLabels = ["CB", "SC", "", ""]

    markers = ['o','o','^','^']
    linestyles = ["-", "--", "-", "--"]
    dashes=[(1,0),(5,2.5),(1,0),(5,2.5)]

    xMin = 450
    xMax = 3600

    yMin = 3e-2
    yMax = 1

    scaleArea1 = 0.15 / math.pow(xMax,1)
    scaleMinArea1 = math.pow(xMin,1) * scaleArea1
    scaleMaxArea1 = math.pow(xMax,1) * scaleArea1

    scaleArea2 = 0.15 / math.pow(xMax,2)
    scaleMinArea2 = math.pow(xMin,2) * scaleArea2
    scaleMaxArea2 = math.pow(xMax,2) * scaleArea2

    scalePos1 = 0.15e4 / math.pow(xMax,1)
    scaleMinPos1 = math.pow(xMin,1) * scalePos1
    scaleMaxPos1 = math.pow(xMax,1) * scalePos1

    scalePos2 = 0.15e4 / math.pow(xMax,2)
    scaleMinPos2 = math.pow(xMin,2) * scalePos2
    scaleMaxPos2 = math.pow(xMax,2) * scalePos2


    scaleThickness = 10e-1 / math.pow(xMin,1)

    scaleMinThickness = math.pow(xMin,1) * scaleThickness
    scaleMaxThickness = math.pow(xMax,1) * scaleThickness

    #plot
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

    #positions

    fig, axes = plt.subplots(1, 1, figsize=(8*cm,6.5*cm))

    axes.loglog([xMin, xMax], [scaleMinPos1, scaleMaxPos1], linestyle=':', color='k', label="_nolegend_", lw=1)
    axes.loglog([xMin, xMax], [scaleMinPos2, scaleMaxPos2], linestyle=':', color='k', label="_nolegend_", lw=1)

    iPlot = 0
    for experiment in experiments:

         xPos = []
         yPos = []

         r = 0
         for resolution in resolutions:

             filename = "./output_%s_%i_%s/output.2000.nc" %(experiment,resolution,runtype)
             filenameParticlesTemplate = "./output_%s_%i_%s/particles_output*" %(experiment,resolution,runtype)

             #norm = get_norm_area(filename, filenameParticlesTemplate)
             norm = check_particle_positions_start_end(filenameParticlesTemplate)
             xPos.append(resolutions[r])
             yPos.append(norm)
             r = r + 1

         error = float('nan')
         print(xPos, yPos)

         axes.loglog(xPos, yPos, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)

         iPlot = iPlot + 1


    axes.legend(legendLabels, frameon=False, loc=4, fontsize=8, handlelength=4)

    plt.minorticks_off()

    axes.set_xlabel("Time resolution (sec)")
    axes.set_ylabel(r"max error norm")
    axes.set_xticks([450,900,1800,3600])
    axes.set_xticklabels(["450","900","1800","3600"])
    axes.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')

    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
    plt.savefig("advection_error_pos_convergence_%s.png" %(runtype),dpi=300)
    plt.savefig("advection_error_pos_convergence_%s.eps" %(runtype))

    #area
    fig, axes = plt.subplots(1, 1, figsize=(8*cm,6.5*cm))

    axes.loglog([xMin, xMax], [scaleMinArea1, scaleMaxArea1], linestyle=':', color='k', label="_nolegend_", lw=1)
    axes.loglog([xMin, xMax], [scaleMinArea2, scaleMaxArea2], linestyle=':', color='k', label="_nolegend_", lw=1)

    iPlot = 0
    for experiment in experiments:

         xArea = []
         yArea = []

         r = 0
         for resolution in resolutions:

             filename = "./output_%s_%i_%s/output.2000.nc" %(experiment,resolution,runtype)
             filenameParticlesTemplate = "./output_%s_%i_%s/particles_output*" %(experiment,resolution,runtype)

             norm = get_norm_area(filename, filenameParticlesTemplate)
             xArea.append(resolutions[r])
             yArea.append(norm)
             r = r + 1

         error = float('nan')
         print(xArea, yArea)
         error = yArea[-1]

         axes.loglog(xArea, yArea, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)

         iPlot = iPlot + 1


    axes.legend(legendLabels, frameon=False, loc=4, fontsize=8, handlelength=4)

    plt.minorticks_off()

    axes.set_xlabel("Time resolution (sec)")
    axes.set_ylabel(r"$L_2$ error norm")
    axes.set_xticks([450,900,1800,3600])
    axes.set_xticklabels(["450","900","1000","3600"])
    axes.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')

    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
    plt.savefig("advection_error_convergence_%s.png" %(runtype),dpi=300)
    plt.savefig("advection_error_convergence_%s.eps" %(runtype))

    try:
        import colorama
        if (error < 1e-12):
           message = "TEST PASSED: Absolute convergence error: %g" %(error)
           print(colorama.Style.BRIGHT + colorama.Fore.GREEN + message + colorama.Style.RESET_ALL)
        else:
           message = "TEST FAILED: Absolute convergence error: %g" %(error)
           print(colorama.Style.BRIGHT + colorama.Fore.MAGENTA + message + colorama.Style.RESET_ALL)
    except ImportError:
        print(message)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot the advection error as a function of mesh size')

    parser.add_argument('-t', required=True, dest='runtype', help='plot data for polympo or nonpolympo run')

    args = parser.parse_args()

    advection_error_convergence(args.runtype)
