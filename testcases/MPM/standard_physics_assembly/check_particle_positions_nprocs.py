from netCDF4 import Dataset
import argparse
import glob
import sys

#-------------------------------------------------------------------------------

def check_particle_positions_nprocs(nProcs1, nProcs2):

    filenames1 = sorted(glob.glob("./output_%i/particles*" %(nProcs1)))
    filenames2 = sorted(glob.glob("./output_%i/particles*" %(nProcs2)))

    print("filename for %i: %s" %(nProcs1, filenames1[-1]))
    print("filename for %i: %s" %(nProcs2, filenames2[-1]))

    file1 = Dataset(filenames1[-1], "r")
    file2 = Dataset(filenames2[-1], "r")

    nParticles1 = len(file1.dimensions["nParticles"])
    nParticles2 = len(file2.dimensions["nParticles"])

    posnMP1 = file1.variables["posnMP"][0,:]
    posnMP2 = file2.variables["posnMP"][0,:]

    statusMP1 = file1.variables["statusMP"][0,:]
    statusMP2 = file2.variables["statusMP"][0,:]

    cellIDCreationMP1 = file1.variables["cellIDCreationMP"][0,:]
    cellIDCreationMP2 = file2.variables["cellIDCreationMP"][0,:]

    creationIndexMP1 = file1.variables["creationIndexMP"][0,:]
    creationIndexMP2 = file2.variables["creationIndexMP"][0,:]

    file1.close()
    file2.close()

    error = False
    nActive1 = sum(statusMP1)
    nActive2 = sum(statusMP2)
    if (nActive1 != nActive2):
        error = True

    if (error):
        print("Difference in number of active material points detected for nProcs: %i %i" %(nProcs1, nProcs2))
        print("%i procs has %i acitve mps; %i procs has %i active mps" %(nProcs1, nActive1, nProcs2, nActive2))
        sys.exit()

    cellIDToIndex2 = {}
    for iParticle2 in range(0,nParticles2):
        if (statusMP2[iParticle2] == 1):
            cellIDToIndex2[(cellIDCreationMP2[iParticle2], creationIndexMP2[iParticle2])] = iParticle2


    error = False
    for iParticle1 in range(0,nParticles1):

        if (statusMP1[iParticle1] == 1):

            iParticle2 = cellIDToIndex2[(cellIDCreationMP1[iParticle1], creationIndexMP1[iParticle1])]

            dx = posnMP2[iParticle2,0] - posnMP1[iParticle1,0]
            dy = posnMP2[iParticle2,1] - posnMP1[iParticle1,1]
            dz = posnMP2[iParticle2,2] - posnMP1[iParticle1,2]

            if (dx != 0.0):
                print("dx different: ", dx, iParticle1, iParticle2)
                error = True
            if (dy != 0.0):
                print("dy different: ", dy, iParticle1, iParticle2)
                error = True
            if (dz != 0.0):
                print("dz different: ", dz, iParticle1, iParticle2)
                error = True

    if (error):
        raise Exception("Error: Differences in particle position detected for nProcs: %i %i" %(nProcs1, nProcs2))

    print("No differences detected in particle position for nProcs: %i %i" %(nProcs1, nProcs2))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(help="Check that two different processor counts give the same answer")

    parser.add_argument('--n1', dest='nProcs1', type=int)
    parser.add_argument('--n2', dest='nProcs2', type=int)

    args = parser.parse_args()

    check_particle_positions_nprocs(args.nProcs1, args.nProcs2)
