from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import argparse

#-------------------------------------------------------------------------------

def create_graph_file_basic(nCells, nProcs):

    filein = Dataset("grid.%i.nc" %(nCells),"r")

    nCells = len(filein.dimensions["nCells"])

    lonCell = filein.variables["lonCell"][:]
    latCell = filein.variables["latCell"][:]

    filein.close()

    fileout = open("./graphs/graph.%i.info.part.%i" %(nCells,nProcs),"w")

    dlon = (2.0 * math.pi) / float(nProcs)

    procArray = []
    for iCell in range(0,nCells):
        for iProc in range(0,nProcs):

            lon1 = dlon * float(iProc)
            lon2 = dlon * float(iProc+1)
            if (iProc == 0):
                lon1 -= 1e-4

            if (lonCell[iCell] >  lon1 and
                lonCell[iCell] <= lon2):

                fileout.write("%i\n" %(iProc))
                procArray.append(iProc)


    fileout.close()

    fig, axis = plt.subplots()

    plt.scatter(lonCell,latCell,c=procArray)
    plt.savefig("./graphs/graph.%i.info.part.%i.png" %(nCells, nProcs), dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', dest='nCells', type=int, required=True)
    parser.add_argument('-n', dest='nProcs', type=int, required=True)

    args = parser.parse_args()

    create_graph_file_basic(args.nCells, args.nProcs)
