from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from datetime import datetime, timedelta
import netCDF4

#-------------------------------------------------------------------------------

def matching_function(x, x1, x2, reverse):

    val = 0.0
    if (x < x1):
        val = 0.0
    elif (x > x2):
        val = 1.0
    else:
        val = 0.5 + 0.5 * math.sin(((x-x1)/(x2-x1)) - 0.5)

    if (reverse):
        val = 1.0 - val

    return val

#-------------------------------------------------------------------------------

def create_forcing():

    # mesh data
    filein = Dataset("grid.nc","r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]-1

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    filein.close()

    xMin = np.amin(xCell)
    xMax = np.amax(xCell)
    yMin = np.amin(yCell)
    yMax = np.amax(yCell)
    Lx = xMax - xMin
    Ly = yMax - yMin

    print("Lx: ", Lx, ", Ly: ", Ly)

    dMargin = 0.1 * Ly
    dTransition = 0.05 * Ly
    
    x1 = xMin + dMargin
    x2 = xMin + dMargin + dTransition
    x3 = 0.5*Lx - dMargin - dTransition
    x4 = 0.5*Lx - dMargin

    y1 = yMin + dMargin
    y2 = yMin + dMargin + dTransition
    y3 = yMax - dMargin - dTransition
    y4 = yMax - dMargin

    # spatially varying fields
    m = np.ones(nCells)
    for iCell in range(0,nCells):
        m[iCell] = m[iCell] * \
            matching_function(xCell[iCell], x1, x2, False) * \
            matching_function(xCell[iCell], x3, x4, True) * \
            matching_function(yCell[iCell], y1, y2, False) * \
            matching_function(yCell[iCell], y3, y4, True)

    # sst
    seaSurfaceTemperature1 = 30.0
    seaSurfaceTemperature2 = -1.9045826499242646
    seaSurfaceTemperature = np.ones(nCells)
    for iCell in range(0,nCells):
        seaSurfaceTemperature[iCell] = \
            seaSurfaceTemperature1 + (seaSurfaceTemperature2 - seaSurfaceTemperature1)*m[iCell]

    # sst
    airTemperature1 = 30.0 + 273.15
    airTemperature2 = -30.0 + 273.15
    airTemperature = np.ones(nCells)
    for iCell in range(0,nCells):
        airTemperature[iCell] = \
            airTemperature1 + (airTemperature2 - airTemperature1)*m[iCell]

    patches = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        patches.append(Polygon(vertices,closed=True))

    fig, axis = plt.subplots()
    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(seaSurfaceTemperature)
    axis.add_collection(pc)
    axis.set_xlim((xMin-0.05*Ly,xMax+0.05*Ly))
    axis.set_ylim((yMin-0.05*Ly,yMax+0.05*Ly))
    axis.set_aspect("equal")
    axis.set_title("seaSurfaceTemperature")
    fig.colorbar(pc)
    plt.savefig("seaSurfaceTemperature.png",dpi=600)

    fig, axis = plt.subplots()
    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(airTemperature)
    axis.add_collection(pc)
    axis.set_xlim((xMin-0.05*Ly,xMax+0.05*Ly))
    axis.set_ylim((yMin-0.05*Ly,yMax+0.05*Ly))
    axis.set_aspect("equal")
    axis.set_title("airTemperature")
    fig.colorbar(pc)
    plt.savefig("airTemperature.png",dpi=600)

    # output ic file
    seaSurfaceSalinity = 34.0
    uOceanVelocity = 0.0
    vOceanVelocity = 0.0
    seaSurfaceTiltU = 0.0
    seaSurfaceTiltV = 0.0
    oceanMixedLayerDepth = 50.0
    oceanHeatFluxConvergence = 0.0
    cloudFraction = 0.0
    rainfallRate = 0.0
    airSpecificHumidity = 0.0
    uAirVelocity = 10.0
    vAirVelocity = 0.0

    # atmos six hourly
    #nTimes = 365*4 # don't need a whole year of forcing for run
    nTimes = 20*4
    
    fileout = Dataset("atmosphere_forcing_six_hourly.2000.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("StrLen",64)
    fileout.createDimension("Time",None)

    time = datetime.fromisoformat('2000-01-01T00:00:00')
    times = []
    for iTime in range(0,nTimes):
        time = time + timedelta(hours=6)
        times.append(time)

    varOutTime = fileout.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,nTimes):
        timeStr = times[iTime].strftime('%Y-%m-%d_%H:%M:%S')
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')

    var = fileout.createVariable("airTemperature","d",dimensions=["Time","nCells"])
    for iTime in range(0,nTimes):
        var[iTime,:] = airTemperature[:]

    var = fileout.createVariable("airSpecificHumidity","d",dimensions=["Time","nCells"])
    var[:] = airSpecificHumidity

    var = fileout.createVariable("uAirVelocity","d",dimensions=["Time","nCells"])
    var[:] = uAirVelocity

    var = fileout.createVariable("vAirVelocity","d",dimensions=["Time","nCells"])
    var[:] = vAirVelocity

    fileout.close()

    # atmos monthly climatology
    nTimes = 12
    
    fileout = Dataset("atmosphere_forcing_monthly.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("StrLen",64)
    fileout.createDimension("Time",None)

    times = []
    for iTime in range(0,nTimes):
        times.append(datetime(year=1, month=iTime+1, day=15))

    varOutTime = fileout.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,nTimes):
        timeStr = times[iTime].strftime('%Y-%m-%d_%H:%M:%S')
        timeStr = "0000" + timeStr[4:]
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')

    var = fileout.createVariable("cloudFraction","d",dimensions=["Time","nCells"])
    var[:] = cloudFraction

    var = fileout.createVariable("rainfallRate","d",dimensions=["Time","nCells"])
    var[:] = rainfallRate

    fileout.close()

    # ocean monthly climatology
    nTimes = 12
    
    fileout = Dataset("ocean_forcing_monthly.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("StrLen",64)
    fileout.createDimension("Time",None)

    times = []
    for iTime in range(0,nTimes):
        times.append(datetime(year=1, month=iTime+1, day=15))

    varOutTime = fileout.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,nTimes):
        timeStr = times[iTime].strftime('%Y-%m-%d_%H:%M:%S')
        timeStr = "0000" + timeStr[4:]
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')
            
    var = fileout.createVariable("seaSurfaceTemperature","d",dimensions=["Time","nCells"])
    for iTime in range(0,nTimes):
        var[iTime,:] = seaSurfaceTemperature[:]

    var = fileout.createVariable("seaSurfaceSalinity","d",dimensions=["Time","nCells"])
    var[:] = seaSurfaceSalinity

    var = fileout.createVariable("uOceanVelocity","d",dimensions=["Time","nCells"])
    var[:] = uOceanVelocity

    var = fileout.createVariable("vOceanVelocity","d",dimensions=["Time","nCells"])
    var[:] = vOceanVelocity

    var = fileout.createVariable("seaSurfaceTiltU","d",dimensions=["Time","nCells"])
    var[:] = seaSurfaceTiltU

    var = fileout.createVariable("seaSurfaceTiltV","d",dimensions=["Time","nCells"])
    var[:] = seaSurfaceTiltV

    var = fileout.createVariable("oceanMixedLayerDepth","d",dimensions=["Time","nCells"])
    var[:] = oceanMixedLayerDepth

    var = fileout.createVariable("oceanHeatFluxConvergence","d",dimensions=["Time","nCells"])
    var[:] = oceanHeatFluxConvergence

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_forcing()
