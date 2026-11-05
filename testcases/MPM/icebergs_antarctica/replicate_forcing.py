import shutil
from netCDF4 import Dataset
import netCDF4
import numpy as np

#-------------------------------------------------------------------------------

def replicate_forcing():

    filenameIn = "forcing_iceberg_monthly.2010.nc"

    years = [x for x in range(2011,2031)]

    for year in years:

        filenameOut = "forcing_iceberg_monthly.%4.4i.nc" %(year)

        shutil.copyfile(filenameIn, filenameOut)

        filein = Dataset(filenameOut,"a")

        nTimes = len(filein.dimensions["Time"])

        times = filein.variables["xtime"][:]

        var = filein.variables["xtime"]

        yearStr = "%4.4i" %(year)
        for iTime in range(0,nTimes):
            timeStr = b''.join(times[iTime]).decode('UTF-8').replace('\x00','')
            timeStr = yearStr + timeStr[4:]
            var[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'), encoding='ascii')

        filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    replicate_forcing()
