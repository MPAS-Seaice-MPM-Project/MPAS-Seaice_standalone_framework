from netCDF4 import Dataset
import math, os, sys
import random
import argparse

#--------------------------------------------------------------------

def create_ics(nCells):

    # grid in
    gridFilename = "grid.%i.nc" %(nCells)
    gridFile = Dataset(gridFilename, "r")

    nCells = len(gridFile.dimensions["nCells"])
    maxEdges = len(gridFile.dimensions["maxEdges"])

    gridFile.close()

    # ice out
    icFilename = "ic.%i.nc" %(nCells)
    icFile = Dataset(icFilename, "w", format="NETCDF3_CLASSIC")

    icFile.createDimension("nCells", nCells)
    icFile.createDimension("maxEdges", maxEdges)

    fieldSubAssemblyTest = icFile.createVariable("fieldSubAssemblyTest", 'd', dimensions=("nCells","maxEdges"))

    for iCell in range(0,nCells):
        fieldSubAssemblyTest[iCell] = random.random()

    icFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells")
    args = parser.parse_args()

    create_ics(args.nCells)
