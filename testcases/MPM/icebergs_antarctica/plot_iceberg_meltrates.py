from netCDF4 import Dataset
import glob
import numpy as np
import argparse
from tqdm import tqdm
from iceberg_plot_utils import plot_cell_field

#-------------------------------------------------------------------------------

def iceberg_meltrates(filenameTemplate,
                      location):

    gigaTonneToKg = 1.0e9 * 1.0e3
    yearToSecond = 60.0 * 60.0 * 24.0 * 365.0
    kgToGigaTonne = 1.0 / gigaTonneToKg
    secondToYear = 1.0 / yearToSecond

    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])
    nEdges = len(fileMesh.dimensions["nEdges"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    cellsOnEdge = fileMesh.variables["cellsOnEdge"][:]-1
    verticesOnEdge = fileMesh.variables["verticesOnEdge"][:]-1
    areaCell = fileMesh.variables["areaCell"][:]
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]
    latEdge = fileMesh.variables["latEdge"][:]

    fileMesh.close()

    filenames = sorted(glob.glob(filenameTemplate))

    icebergMeltRateCell = np.zeros(nCells)
    nAvg = 0
    for filename in tqdm(filenames):
        fileIn = Dataset(filename,"r")

        icebergMeltRateCellIn = fileIn.variables["icebergMeltRateCell"][:,:]
        icebergMeltRateCell[:] += np.sum(icebergMeltRateCellIn, axis=0)

        fileIn.close()

        nAvg += icebergMeltRateCellIn.shape[0]

    icebergMeltRateCell[:] /= nAvg

    totalDischarge = 0.0
    for iCell in range(0,nCells):
        totalDischarge += icebergMeltRateCell[iCell] * areaCell[iCell]

    totalDischarge *= kgToGigaTonne * yearToSecond
    print("totalDischarge: ", totalDischarge, " Gt/year")

    plot_cell_field(icebergMeltRateCell,
                    "Iceberg melting",
                    r'Melt flux ($\mathrm{kg}/\mathrm{m}^2/\mathrm{s}$)',
                    "iceberg_melting.png",
                    location,
                    nEdges,
                    nCells,
                    nEdgesOnCell,
                    cellsOnEdge,
                    latEdge,
                    latCell,
                    verticesOnEdge,
                    verticesOnCell,
                    latVertex,
                    lonVertex)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    iceberg_meltrates(args.filenameTemplate,
                      args.location)
