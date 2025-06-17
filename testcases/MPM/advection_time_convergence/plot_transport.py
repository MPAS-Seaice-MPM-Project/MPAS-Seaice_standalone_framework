from netCDF4 import Dataset
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import glob
import argparse

#-------------------------------------------------------------------------------

def plot_file(filenameIn, filenameOut):

    iTime = 0

    filein = Dataset(filenameIn,"r")

    sphere_radius = filein.sphere_radius

    nParticles = len(filein.dimensions["nParticles"])

    posnMP = filein.variables["posnMP"][:]

    procIDParticle = filein.variables["procIDParticle"][:]

    indexToCellIDParticle = filein.variables["indexToCellIDParticle"][:]

    nValid = 0
    valid = []
    colors = []
    for iParticle in range(0,nParticles):
        mag = math.sqrt(math.pow(posnMP[iTime,iParticle,0],2) + \
                        math.pow(posnMP[iTime,iParticle,1],2) + \
                        math.pow(posnMP[iTime,iParticle,2],2))
        if (mag > 1e-8):
            nValid += 1
            valid.append(True)
            colors.append(procIDParticle[iTime,iParticle])
        else:
            valid.append(False)



    filein.close()

    fig, axis = plt.subplots(subplot_kw={'projection':'3d'})

    axis.scatter(posnMP[iTime,valid,0],posnMP[iTime,valid,1],posnMP[iTime,valid,2],c=colors)

    axis.set_xlim((-1.1*sphere_radius,1.1*sphere_radius))
    axis.set_ylim((-1.1*sphere_radius,1.1*sphere_radius))
    axis.set_zlim((-1.1*sphere_radius,1.1*sphere_radius))

    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")

    plt.savefig(filenameOut)
    plt.cla()
    plt.close(fig)


    return nValid

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot the MPM advection test case')

    parser.add_argument('-f', required=True, dest='filenameTemplate', help='filename template for particle files')

    args = parser.parse_args()

    nskip = 1

    filenames = sorted(glob.glob(args.filenameTemplate))

    iFile = 0
    for filenameIn in filenames[::nskip]:
        print(filenameIn)

        filenameOut = "./plots/plot_transport_%4.4i.png" %(iFile)

        nValid = plot_file(filenameIn, filenameOut)

        print(filenameOut, nValid)
        iFile += 1
