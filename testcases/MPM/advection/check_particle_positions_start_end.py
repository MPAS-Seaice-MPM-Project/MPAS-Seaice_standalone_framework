from netCDF4 import Dataset
import sys
import math
import glob
import argparse

#-------------------------------------------------------------------------------

def check_particle_positions_start_end(filenameTemplate):

    filenames = sorted(glob.glob(filenameTemplate))

    filein = Dataset(filenames[0],"r")

    nParticles1 = len(filein.dimensions["nParticles"])

    statusMP1 = filein.variables["statusMP"][0,:]
    cellIDCreationMP1 = filein.variables["cellIDCreationMP"][0,:]
    creationIndexMP1 = filein.variables["creationIndexMP"][0,:]
    posnMP1 = filein.variables["posnMP"][0,:,:]

    filein.close()


    errors = []
    iFile = 0
    for filenameIn in filenames:
        iFile += 1

        filein = Dataset(filenameIn,"r")

        nParticles2 = len(filein.dimensions["nParticles"])

        statusMP2 = filein.variables["statusMP"][0,:]
        cellIDCreationMP2 = filein.variables["cellIDCreationMP"][0,:]
        creationIndexMP2 = filein.variables["creationIndexMP"][0,:]
        posnMP2 = filein.variables["posnMP"][0,:,:]

        filein.close()

        if (nParticles1 != nParticles2):
            raise Exception("nParticles not consistent for " + filenameIn)


        particleIDToIndex = {}
        for iParticle2 in range(0, nParticles2):
            if (statusMP2[iParticle2] == 1):
                particleIDToIndex[(cellIDCreationMP2[iParticle2],creationIndexMP2[iParticle2])] = iParticle2


        maxDistance = -sys.float_info.max
        for iParticle1 in range(0,nParticles1):
            if (statusMP1[iParticle1] == 1):
                iParticle2 = particleIDToIndex[(cellIDCreationMP1[iParticle1],creationIndexMP1[iParticle1])]
                distance = math.sqrt(math.pow(posnMP2[iParticle2,0]-posnMP1[iParticle1,0],2) +
                                     math.pow(posnMP2[iParticle2,1]-posnMP1[iParticle1,1],2) +
                                     math.pow(posnMP2[iParticle2,2]-posnMP1[iParticle1,2],2))
                maxDistance = max(maxDistance,distance)

        print("max distance: ", iFile, maxDistance)

        errors.append(maxDistance)

    print("Final max distance: ", errors[-1])

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Determine error in particle positions after one rotation on sphere')

    parser.add_argument('-f', required=True, dest='filenameTemplate', help='filename template for particle files')

    args = parser.parse_args()

    check_particle_positions_start_end(args.filenameTemplate)
