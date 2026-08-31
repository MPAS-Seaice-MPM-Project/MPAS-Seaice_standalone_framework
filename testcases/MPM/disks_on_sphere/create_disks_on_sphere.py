import sys
import math
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from numba import njit
sys.path.append("../../../utils/MPM/particle_initialization/")
from initial_particle_positions import create_particles_file
from initial_particle_positions import place_particles_in_cell
import argparse

# Check if point is inside disk
#-------------------------------------------------------------------------------
@njit
def in_disk(x, y, z, geomType):
                
    r = math.sqrt(x*x + y*y + z*z)
    lat = math.asin(z/r) * 180.0 / math.pi
    lon =  math.atan2(y, x)  * 180 / math.pi
    iceArea = 0.0
    iceVolume = 0.0
    in_geom = False
    
    if (geomType == "single_particle"):
        if((lat ** 2 + lon ** 2) < 1.0):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0
            geomType = "skip"
    elif (geomType == 'skip'):
        in_geom = False
        iceArea = 0.0
        iceVolume = 0.0
        geomType = "skip"
        
    elif(geomType == 'two_disks'):
        radius = 8
        centerlon = 9
        if(math.sqrt(lat ** 2 + (lon+centerlon) ** 2) < radius):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0
            
        if(math.sqrt(lat ** 2 + (lon-centerlon) ** 2) < radius):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0
            
    elif(geomType == 'many_disks'):
        radius = 2
        latwrap = lat % 10 - 5
        lonwrap = lon % 10 - 5
        
        if(math.sqrt(latwrap ** 2 + lonwrap **2) < radius and abs(lat) < 18 and abs(lon) <18):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    return in_geom, geomType, iceArea, iceVolume

#-------------------------------------------------------------------------------
# Check if particles in this cell are inside disk, create file

def check_cell_particles(xInCell, yInCell, zInCell, geomType, iCell,posnMP, latCellMP,
                         lonCellMP, cellIDCreationMP,creationIndexMP, iceAreaCellMP,
                         iceVolumeCellMP, nParticlesPerCellActual):

    # check if particles in cell are in geometry
    k = 0
    for x, y, z in zip(xInCell, yInCell, zInCell):
        in_geom, geomType, iceArea, iceVolume  = in_disk(x, y, z, geomType)

        if (in_geom):
            k = k + 1
            lat = math.asin(z) 
            lon =  math.atan2(y, x)
            posnMP.append([x, y, z])
            latCellMP.append(lat)
            lonCellMP.append(lon)
            cellIDCreationMP.append(iCell+1)
            creationIndexMP.append(k)
            iceAreaCellMP.append(iceArea)
            iceVolumeCellMP.append(iceVolume)

    nParticlesPerCellActual = k

    return \
        posnMP, \
        latCellMP, \
        lonCellMP, \
        cellIDCreationMP, \
        creationIndexMP, \
        iceAreaCellMP, \
        iceVolumeCellMP, \
        nParticlesPerCellActual, \
        geomType



def create_write_particle_file(filenameMesh,
                               filenameOut,
                               geomType):

    # mesh info
    fileMesh = Dataset(filenameMesh,"r")

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell   = ma.getdata(fileMesh.variables["nEdgesOnCell"][:])
    verticesOnCell = ma.getdata(fileMesh.variables["verticesOnCell"][:])

    areaCell  = ma.getdata(fileMesh.variables["areaCell"][:])
    latVertex = ma.getdata(fileMesh.variables["latVertex"][:])
    lonVertex = ma.getdata(fileMesh.variables["lonVertex"][:])
    xVertex   = ma.getdata(fileMesh.variables["xVertex"][:])
    yVertex   = ma.getdata(fileMesh.variables["yVertex"][:])
    zVertex   = ma.getdata(fileMesh.variables["zVertex"][:])
    xCell     = ma.getdata(fileMesh.variables["xCell"][:])
    yCell     = ma.getdata(fileMesh.variables["yCell"][:])
    zCell     = ma.getdata(fileMesh.variables["zCell"][:])


    fileMesh.close()

    verticesOnCell[:] -= 1

    # average cell size
    averageCellSize = np.mean(areaCell)
    earthRadius = math.sqrt(xVertex[1] ** 2 + yVertex[1] ** 2 + zVertex[1] ** 2 )

    # particle positions
    nParticlesCell = np.zeros(nCells,dtype="i")

    posnMP = []
    latCellMP = []
    lonCellMP = []
    cellIDCreationMP = []
    creationIndexMP = []
    iceAreaCellMP = []
    iceVolumeCellMP = []

    for iCell in range(0, nCells):

        xInCell, \
            yInCell, \
            zInCell = place_particles_in_cell("number","9",areaCell[iCell], averageCellSize,"onePerEdge",True,earthRadius,
                                              nEdgesOnCell[iCell], verticesOnCell[iCell,:], latVertex, lonVertex, 
                                              xVertex, yVertex, zVertex,
                                              xCell[iCell], yCell[iCell], zCell[iCell])

        posnMP, \
            latCellMP, \
            lonCellMP, \
            cellIDCreationMP, \
            creationIndexMP, \
            iceAreaCellMP, \
            iceVolumeCellMP, \
            nParticlesCell[iCell], \
            geomType = check_cell_particles(xInCell, yInCell, zInCell,
                                                     geomType,
                                                     iCell,
                                                     posnMP,
                                                     latCellMP,
                                                     lonCellMP,
                                                     cellIDCreationMP,
                                                     creationIndexMP,
                                                     iceAreaCellMP,
                                                     iceVolumeCellMP,
                                                     nParticlesCell[iCell])
                                                     
    nParticles = len(posnMP)
    posnMP = np.array(posnMP)
    latCellMP = np.array(latCellMP)
    lonCellMP = np.array(lonCellMP)
    cellIDCreationMP = np.array(cellIDCreationMP)
    creationIndexMP = np.array(creationIndexMP)
    nParticlesCell = np.array(nParticlesCell)
    iceAreaCellMP = np.array(iceAreaCellMP)
    iceVolumeCellMP = np.array(iceVolumeCellMP)
    iceAreaCategoryMP = iceAreaCellMP
    iceVolumeCategoryMP = iceVolumeCellMP

    # output
    create_particles_file(filenameOut,
                          nParticles,
                          nCells,
                          1,
                          posnMP,
                          latCellMP,
                          lonCellMP,
                          cellIDCreationMP,
                          creationIndexMP,
                          nParticlesCell,
                          iceAreaCellMP,
                          iceVolumeCellMP,
                          iceAreaCategoryMP,
                          iceVolumeCategoryMP)
                          

    uvVelMP = np.zeros([nParticles, 2])
    for iParticle in range(0,nParticles):
        lat = math.asin(posnMP[iParticle,2]) * 180 / math.pi
        lon = math.atan2(posnMP[iParticle,1],posnMP[iParticle,0])  * 180 / math.pi
        uvVelMP[iParticle, 0] = -math.copysign(0.02,lon)
        if(geomType == 'many_disks' or geomType == 'single_particle'):
            uvVelMP[iParticle, 1] = -math.copysign(0.02,lat)
                  
    # append velocity to particle file
    fileOut = Dataset(filenameOut, "a", format="NETCDF3_CLASSIC")
    var = fileOut.createVariable("uvVelMP","d",dimensions=["nParticles", "TWO"])
    var[:] = uvVelMP[:]
    fileOut.close()
            

#-----------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="filenameMesh", required=True)
    parser.add_argument('-o', dest="filenameOut", required=True)
    parser.add_argument('-i', dest="geomType", required=True, choices=["single_particle","two_disks","many_disks"])

    args = parser.parse_args()

    create_write_particle_file(args.filenameMesh,
                               args.filenameOut,
                               args.geomType)