from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from datetime import datetime, timedelta
import netCDF4

#-------------------------------------------------------------------------------

def create_forcing(gridFilename):

    # mesh data
    filein = Dataset(gridFilename,"r")

    nCells = len(filein.dimensions["nCells"])

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    filein.close()

    xMin = np.amin(xCell)
    xMax = np.amax(xCell)
    yMin = np.amin(yCell)
    yMax = np.amax(yCell)
    Lx = xMax - xMin
    Ly = yMax - yMin

    print("Lx: ", Lx, ", Ly: ", Ly)

    # output forcing file
    seaSurfaceTemperature = -1.9045826499242646
    airTemperature = -38.0 + 273.15
    seaSurfaceSalinity = 32.0
    uOceanVelocity = 0.0
    vOceanVelocity = 0.0
    seaSurfaceTiltU = 0.0
    seaSurfaceTiltV = 0.0
    oceanMixedLayerDepth = 50.0
    oceanHeatFluxConvergence = 0.0
    cloudFraction = 0.0
    rainfallRate = 0.0
    airSpecificHumidity = 0.0

    uAirSpeed = 3.0
    vAirSpeed = 0.0

    # atmos six hourly
    nTimes = 365*4

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
        timeStr = times[iTime].isoformat(sep="_", timespec="seconds")
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')

    var = fileout.createVariable("airTemperature","d",dimensions=["Time","nCells"])
    var[:] = airTemperature

    var = fileout.createVariable("airSpecificHumidity","d",dimensions=["Time","nCells"])
    var[:] = airSpecificHumidity

    var = fileout.createVariable("uAirVelocity","d",dimensions=["Time","nCells"])
    for iTime in range(0,nTimes):

       var[iTime:] = -uAirSpeed

       if (iTime <= 20):
          for iCell in range(0,nCells):
            if (xCell[iCell] >= 90000):
               var[iTime,iCell] = uAirSpeed
       elif (iTime > 40 and iTime <= 80):
          var[iTime:] = 0.0

       if (iTime > 80 and iTime <= 100):
          for iCell in range(0,nCells):
            if (xCell[iCell] >= 90000):
               var[iTime,iCell] = uAirSpeed
       elif (iTime > 120 and iTime <= 160):
          var[iTime:] = 0.0

       elif (iTime > 160 and iTime <= 180):
          for iCell in range(0,nCells):
            if (xCell[iCell] >= 90000):
               var[iTime,iCell] = uAirSpeed
       elif (iTime > 200 and iTime <= 240):
          var[iTime:] = 0.0

       if (iTime > 240 and iTime <= 260):
          for iCell in range(0,nCells):
            if (xCell[iCell] >= 90000):
               var[iTime,iCell] = uAirSpeed
       elif (iTime > 280 and iTime <= 320):
          var[iTime:] = 0.0

       if (iTime > 320 and iTime <= 340):
          for iCell in range(0,nCells):
            if (xCell[iCell] >= 90000):
               var[iTime,iCell] = uAirSpeed
       elif (iTime > 360):
          var[iTime:] = 0.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["Time","nCells"])
    var[:] = vAirSpeed

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
        timeStr = times[iTime].isoformat(sep="_", timespec="seconds")
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
        timeStr = times[iTime].isoformat(sep="_", timespec="seconds")
        timeStr = "0000" + timeStr[4:]
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')

    var = fileout.createVariable("seaSurfaceTemperature","d",dimensions=["Time","nCells"])
    for iTime in range(0,nTimes):
        var[iTime,:] = seaSurfaceTemperature

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

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)

    args = parser.parse_args()

    create_forcing(args.gridFilename)
