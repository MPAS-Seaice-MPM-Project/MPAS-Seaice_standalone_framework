from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numba import njit

#---------------------------------------------------------

@njit
def average(nCells,
            nParticles,
            statusMP,
            iCellMP,
            iceAreaCellMP):

    iceAreaCell = np.zeros(nCells)
    nParticlesCell = np.zeros(nCells)

    for iParticle in range(0,nParticles):
        if (statusMP[iParticle] == 1):
            iCell = iCellMP[iParticle]
            iceAreaCell[iCell] += iceAreaCellMP[iParticle]
            nParticlesCell[iCell] += 1.0

    for iCell in range(0,nCells):
        if (nParticlesCell[iCell] > 0.0):
            iceAreaCell[iCell] /= nParticlesCell[iCell]

    return iceAreaCell

#---------------------------------------------------------

def average_particle_ice_area_to_cell(filenameParticles, nCells):

    filein = Dataset(filenameParticles,"r")

    nParticles = len(filein.dimensions["nParticles"])

    statusMP      = ma.getdata(filein.variables["statusMP"][0,:])
    iCellMP       = ma.getdata(filein.variables["iCellMP"][0,:])-1
    iceAreaCellMP = ma.getdata(filein.variables["iceAreaCellMP"][0,:])

    filein.close()

    iceAreaCell = average(nCells,
                          nParticles,
                          statusMP,
                          iCellMP,
                          iceAreaCellMP)

    return iceAreaCell

#---------------------------------------------------------
