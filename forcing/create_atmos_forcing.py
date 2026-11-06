from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import subprocess
import sys
import glob
import configparser
from create_forcing import create_scrip_grid_file, get_mpas_grid_info, create_scrip_file_MPAS, write_scrip_in_file, create_output_times, get_remapping_data
import argparse

#-------------------------------------------------------------------------------

def create_T62_remap_file(filenameT62,
                          filenameScrip,
                          title,
                          latDimname,
                          lonDimname,
                          latVarname,
                          lonVarname):

    fileT62 = Dataset(filenameT62,"r")

    nLat = len(fileT62.dimensions[latDimname])
    nLon = len(fileT62.dimensions[lonDimname])

    LATin = fileT62.variables[latVarname][:]
    LONin = fileT62.variables[lonVarname][:]

    fileT62.close()

    LAT = np.zeros(nLat+2)
    LON = np.zeros(nLon+2)

    LAT[1:-1] = LATin[:]
    LON[1:-1] = LONin[:]

    # center
    latCenter = np.zeros((nLat,nLon))
    lonCenter = np.zeros((nLat,nLon))

    for iLon in range(0,nLon):
       latCenter[:,iLon] = LATin[:]

    for iLat in range(0,nLat):
       lonCenter[iLat,:] = LONin[:]

    latCenter = np.radians(latCenter)
    lonCenter = np.radians(lonCenter)

    # corners
    latCorner = np.zeros((nLat,nLon,4))
    lonCorner = np.zeros((nLat,nLon,4))

    for iLon in range(0,nLon):

        iLon2 = iLon + 1

        lonCorner[:,iLon,0] = 0.5 * (LON[iLon2] + LON[iLon2-1])
        lonCorner[:,iLon,1] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,2] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,3] = 0.5 * (LON[iLon2] + LON[iLon2-1])

    lonCorner = (np.where(lonCorner < 0.0, lonCorner + 360.0, lonCorner))
    lonCorner = np.radians(lonCorner)

    for iLat in range(0,nLat):

        iLat2 = iLat + 1

        latCorner[iLat,:,0] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,1] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,2] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])
        latCorner[iLat,:,3] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])

    latCorner = np.radians(latCorner)

    # create file
    nGridSize = nLat * nLon
    nGridCorners = 4
    gridRank = 2
    gridDims = np.array([nLon,nLat])
    gridImask = np.ones(nGridSize,dtype="i")

    latCornerScrip = np.zeros((nLat*nLon,4))
    lonCornerScrip = np.zeros((nLat*nLon,4))

    for iLat in range(0,nLat):
        for iLon in range(0,nLon):
            for iCorner in range(0,4):
                ij = iLat * nLon + iLon
                latCornerScrip[ij,iCorner] = latCorner[iLat,iLon,iCorner]
                lonCornerScrip[ij,iCorner] = lonCorner[iLat,iLon,iCorner]

    create_scrip_grid_file(filenameScrip,
                           nGridSize,
                           nGridCorners,
                           gridRank,
                           gridDims,
                           latCenter.flatten(),
                           lonCenter.flatten(),
                           gridImask,
                           latCornerScrip,
                           lonCornerScrip,
                           title)

    return nGridSize

#-------------------------------------------------------------------------------

def create_forcing(yearStart,
                   yearStop,
                   timeType,
                   inputDir,
                   filenameTemplatesIn,
                   varnamesIn,
                   outputDir,
                   filenameOutTemplate,
                   varnamesOut,
                   remapMatrix,
                   dstGridSize):

    # loop over years
    for year in range(yearStart,yearStop+1):

        print("  Year: %i of %i to %i" %(year, yearStart, yearStop))

        # create output file
        filenameOut = outputDir+"/"+filenameOutTemplate.replace("$Y",str(year))
        fileForcing = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

        # dimensions
        nCells = fileForcing.createDimension("nCells",dstGridSize)
        StrLen = fileForcing.createDimension("StrLen",64)
        Time   = fileForcing.createDimension("Time",)

        # time
        xtimes = create_output_times(timeType, year)
        nTimes = len(xtimes)
        varXtime = fileForcing.createVariable("xtime","c",dimensions=["Time","StrLen"])
        for iTime in range(0,nTimes):
            varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'), encoding='ascii')
            varXtime[iTime,19:] = " "*45

        # loop over variables
        for iVariable in range(0,len(varnamesIn)):

            print("    Variable: %s to %s" %(varnamesIn[iVariable], varnamesOut[iVariable]))

            # open input file
            filenameTemplate = inputDir+filenameTemplatesIn[iVariable].replace("$Y",str(year))
            filenamesInput = sorted(glob.glob(filenameTemplate))
            if (len(filenamesInput) == 0):
                raise Exception("Empty filenamesInput list: "+filenameTemplate)
            filenameInput = filenamesInput[0]
            fileInput = Dataset(filenameInput,"r")
            arrayIn = fileInput.variables[varnamesIn[iVariable]][:]
            fileInput.close()

            arrayOut = np.zeros((nTimes,dstGridSize))

            # loop over times
            for iTime in range(0,nTimes):

                arrayInTime = arrayIn[iTime,:,:].flatten()
                arrayOut[iTime,:] = remapMatrix.dot(arrayInTime)

            # output variable to netcdf file
            var = fileForcing.createVariable(varnamesOut[iVariable],"d",dimensions=["Time","nCells"])
            var[:] = arrayOut[:]

        # close forcing file
        fileForcing.close()

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remapInputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    gridFile1 = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    gridFile2 = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interpFile1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interpFile2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    mapName1 = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapName2 = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapMethod = 'bilinear'\n")
    scripFile.write("    normalizeOpt = 'fracArea'\n")
    scripFile.write("    outputFormat = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

#-------------------------------------------------------------------------------

def perform_remapping(startYear,
                      endYear,
                      timeType,
                      meshDir,
                      meshFilename,
                      latDimname,
                      lonDimname,
                      latVarname,
                      lonVarname,
                      inputDir,
                      filenameTemplatesIn,
                      varnamesIn,
                      filenameMPASGrid,
                      outputDir,
                      filenameOutTemplate,
                      varnamesOut):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename  = "remap_grid_MPAS_tmp.nc"
    dstGridSize = create_scrip_file_MPAS(filenameMPASGrid,
                                         scripGridFilename)

    # create T62 remapping file
    print("create_T62_remap_file")
    filenameT62 = meshDir+meshFilename
    scripT62Filename = "remap_grid_T62_tmp.nc"
    srcGridSize = create_T62_remap_file(filenameT62,
                                        scripT62Filename,
                                        "T62",
                                        latDimname,
                                        lonDimname,
                                        latVarname,
                                        lonVarname)

    # create input scrip file
    print("write_scrip_in_file")
    write_scrip_in_file("T62")

    # run scrip to generate weights
    print("ESMF_RegridWeightGen")
    process = subprocess.Popen(["ESMF_RegridWeightGen",
                                "--source",     "remap_grid_T62_tmp.nc",
                                "--destination","remap_grid_MPAS_tmp.nc",
                                "--weight",     "remap_T62_to_MPAS_tmp.nc",
                                "--method",     "bilinear",
                                "--weight_only"])
    process.wait()
    if (process.returncode != 0):
        print("ESMF_RegridWeightGen error: ", process.returncode)
        exit(1);

    # get remapping weights
    print("get_remapping_data")
    filenameRemapping = "remap_T62_to_MPAS_tmp.nc"
    remapMatrix = get_remapping_data(filenameRemapping,
                                     srcGridSize,
                                     dstGridSize)

    # combined output file
    print("create_forcing")
    create_forcing(startYear,
                   endYear,
                   timeType,
                   inputDir,
                   filenameTemplatesIn,
                   varnamesIn,
                   outputDir,
                   filenameOutTemplate,
                   varnamesOut,
                   remapMatrix,
                   dstGridSize)

#-------------------------------------------------------------------------------

def create_forcing_from_config(configFilename):

    config = configparser.ConfigParser()
    config.read(configFilename)

    # time
    startYear              = config.getint('time',         'startYear')
    endYear                = config.getint('time',         'endYear')
    timeType               = config.get   ('time',         'timeType')

    # input_mesh
    meshDir                = config.get   ('input_mesh',   'meshDir')
    meshFilename           = config.get   ('input_mesh',   'meshFilename')
    latDimname             = config.get   ('input_mesh',   'latDimname')
    lonDimname             = config.get   ('input_mesh',   'lonDimname')
    latVarname             = config.get   ('input_mesh',   'latVarname')
    lonVarname             = config.get   ('input_mesh',   'lonVarname')

    # input_fields
    inputDir               = config.get   ('input_fields', 'inputDir')
    filenameTemplatesInStr = config.get   ('input_fields', 'filenameTemplatesIn')
    varnamesInStr          = config.get   ('input_fields', 'varnamesIn')

    # output
    filenameMPASGrid       = config.get   ('output',       'filenameMPASGrid')
    outputDir              = config.get   ('output',       'outputDir')
    filenameOutTemplate    = config.get   ('output',       'filenameOutTemplate')
    varnamesOutStr         = config.get   ('output',       'varnamesOut')

    # list conversion
    filenameTemplatesIn = list(filter(None, [x.strip() for x in filenameTemplatesInStr.splitlines()]))
    varnamesIn          = list(filter(None, [x.strip() for x in varnamesInStr.splitlines()]))
    varnamesOut         = list(filter(None, [x.strip() for x in varnamesOutStr.splitlines()]))

    if (len(filenameTemplatesIn) != len(varnamesIn) or
        len(filenameTemplatesIn) != len(varnamesOut)):
        raise Exception("Config lists not same size")

    perform_remapping(startYear,
                      endYear,
                      timeType,
                      meshDir,
                      meshFilename,
                      latDimname,
                      lonDimname,
                      latVarname,
                      lonVarname,
                      inputDir,
                      filenameTemplatesIn,
                      varnamesIn,
                      filenameMPASGrid,
                      outputDir,
                      filenameOutTemplate,
                      varnamesOut)

#-------------------------------------------------------------------------------

'''
create_atmos_forcing.py
=======================

Usage
-----

This script creates atmospheric forcing from input data defined on lat/lon grids.
Examples include CORE-II, AOMIP climatologies, and JRA-3Q.

Usage: python create_atmos_forcing.py -c configFilename

where configFilename is a python config file. The following example gives values
appropriate for CORE-II forcing:

[time]
startYear = 2000
endYear = 2000
timeType = sixhourly_noleap

[input_mesh]
meshDir = data/CORE-II/
meshFilename = /t_10/t_10.2000.nc
latDimname = LAT
lonDimname = LON
latVarname = LAT
lonVarname = LON

[input_fields]
inputDir = data/CORE-II/
filenameTemplatesIn =
                    /t_10/t_10.$Y.*nc
                    /q_10/q_10.$Y.*nc
                    /u_10/u_10.$Y.*nc
                    /v_10/v_10.$Y.*nc
varnamesIn =
           T_10_MOD
           Q_10_MOD
           U_10_MOD
           V_10_MOD

[output]
filenameMPASGrid = domain_QU120km/seaice_QU_120km.nc
outputDir = tmp
filenameOutTemplate = LYq_six_hourly.$Y.nc
varnamesOut =
            airTemperature
            airSpecificHumidity
            uAirVelocity
            vAirVelocity

ESMF_RegridWeightGen
--------------------

This script requires the ESMF_RegridWeightGen utility to be installed.

CORE-II data
------------

Six-hourly air temperature, velocity and specific humidity comes from CORE-II.
Data files can be obtained from
https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html
To generate forcing for a given year YYYY, the following files are required:
${dataDirSixHourly}/t_10/t_10.YYYY.*.nc
${dataDirSixHourly}/q_10/q_10.YYYY.*.nc
${dataDirSixHourly}/u_10/u_10.YYYY.*.nc
${dataDirSixHourly}/v_10/v_10.YYYY.*.nc
where ${dataDirSixHourly} is the local location of the six hourly data.

AOMIP climatologies
-------------------

Monthly climatologies of cloudiness and precipitation comes from AOMIP.
The following data files are required:
${dataDirMonthly}/cldf.omip.nc
${dataDirMonthly}/prec.nmyr.nc
where ${dataDirMonthly} is the local location of the monthly data.
These files can be obtained from:
https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/
MPAS-Seaice_clim_data.tar.gz
'''

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create atmospheric forcing from lat/lon input data')

    parser.add_argument('-c', dest='configFilename', required=True, help='Config filename')

    args = parser.parse_args()

    create_forcing_from_config(args.configFilename)
