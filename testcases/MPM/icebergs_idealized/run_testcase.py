import sys

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from plot_results import plot_testcase
from grounding_tests import groundingTestcases, write_grounding_test_files
from seaice_stress_tests import seaiceStressTestcases, write_seaice_stress_test_files

import numpy as np
from netCDF4 import Dataset
from math import pi
import os
import shutil

#-------------------------------------------------------------------------------

def write_iceberg_file(nIcebergs,
                       latIceberg,
                       lonIceberg,
                       icebergLength,
                       icebergHeight,
                       filename):

    fileout = Dataset(filename,"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",nIcebergs)

    var = fileout.createVariable("cellIDCreationIB", "i", dimensions=["nIcebergs"])
    var[:] = 1

    var = fileout.createVariable("latIceberg", "d", dimensions=["nIcebergs"])
    var[:] = latIceberg[:]

    var = fileout.createVariable("lonIceberg", "d", dimensions=["nIcebergs"])
    var[:] = lonIceberg[:]

    var = fileout.createVariable("icebergLength", "d", dimensions=["nIcebergs"])
    var[:] = icebergLength[:]

    var = fileout.createVariable("icebergHeight", "d", dimensions=["nIcebergs"])
    var[:] = icebergHeight[:]

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_air_drag(nCells):

    fileout = Dataset("ic_air_drag.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_ocean_drag(nCells):

    fileout = Dataset("ic_ocean_drag.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_coriolis():

    fileout = Dataset("ic_coriolis.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",1)

    var = fileout.createVariable("uVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = 0.0

    var = fileout.createVariable("vVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = 1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_surface_tilt(nCells):

    fileout = Dataset("ic_surface_tilt.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_wave_radiation(nCells):

    fileout = Dataset("ic_wave_radiation.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = -1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_seaice_force(nCells,
                               nVertices,
                               concLabel,
                               iceConcentration,
                               icePressure,
                               iceThickness=1.0,
                               seaiceVelocity=1.0,
                               oceanVelocity=None):

    fileout = Dataset("ic_seaice_force_%s.nc" %(concLabel),"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)

    var = fileout.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = 0.0

    var = fileout.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = seaiceVelocity

    var = fileout.createVariable("iceAreaCell","d",dimensions=["nCells"])
    var[:] = iceConcentration

    # cell average sea ice volume (m), sets the sea ice drag contact area
    var = fileout.createVariable("iceVolumeCell","d",dimensions=["nCells"])
    var[:] = iceConcentration * iceThickness

    var = fileout.createVariable("icePressure","d",dimensions=["nCells"])
    var[:] = icePressure

    if (oceanVelocity is not None):

        var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
        var[:] = 0.0

        var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
        var[:] = oceanVelocity

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_grounding(nCells, lat):

    fileout = Dataset("ic_grounding.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 1.0

    # shallow region bathymetry above iceberg draft (~828 m) so iceberg grounds
    lat1 = (3.0*pi)/8.0 - pi/16.0
    lat2 = (3.0*pi)/8.0 + pi/16.0
    bathymetry = np.where((lat > lat1) & (lat < lat2), -700.0, -1100.0)

    var = fileout.createVariable("bathymetry","d",dimensions=["nCells"])
    var[:] = bathymetry

    fileout.close()

#-------------------------------------------------------------------------------

def run_testcase():

    nCells = 2562

    print("Get testcase data")
    print("=================")
    get_testcase_data_spherical()

    print("Grid data")
    print("=========")
    try:
        os.remove("grid.nc")
    except:
        pass
    os.symlink("grid.%i.nc" %(nCells), "grid.nc")

    filenameGrid = "grid.nc"
    fileGrid = Dataset(filenameGrid,"r")
    nVertices = len(fileGrid.dimensions["nVertices"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    latCell = fileGrid.variables["latCell"][:]
    lonCell = fileGrid.variables["lonCell"][:]
    fileGrid.close()

    print("Iceberg particle input")
    print("======================")
    write_iceberg_file(1,
                       np.array([pi/4.0]),
                       np.array([0.0]),
                       np.array([1000.0]),
                       np.array([1000.0]),
                       "icebergs_midlatitude.nc")

    print("Iceberg IC data")
    print("===============")
    write_ic_file_air_drag(nCells)
    write_ic_file_ocean_drag(nCells)
    write_ic_file_surface_tilt(nCells)
    write_ic_file_wave_radiation(nCells)
    write_ic_file_seaice_force(nCells, nVertices, "low_conc",  0.01, 0.0)
    write_ic_file_seaice_force(nCells, nVertices, "mid_conc",  0.5,  0.0)
    # iceberg at rest accelerated by ocean and sea ice drag until the sea ice force
    # needed to move it with the sea ice is below the capture limit (pressure * width)
    write_ic_file_seaice_force(nCells, nVertices, "high_conc", 0.99, 4.5e4,
                               seaiceVelocity=0.5, oceanVelocity=0.4)
    write_ic_file_coriolis()
    write_ic_file_grounding(nCells, latCell)

    print("Grounding test files")
    print("====================")
    write_grounding_test_files(nCells, latCell)

    print("Sea ice stress test files")
    print("=========================")
    write_seaice_stress_test_files(nCells, nVertices, maxEdges)

    testcases = ["air_drag",
                 "ocean_drag",
                 "coriolis",
                 "surface_tilt",
                 "wave_radiation",
                 "seaice_force_low_conc",
                 "seaice_force_mid_conc",
                 "seaice_force_high_conc",
                 "grounding",
                 "grounding_force"] + groundingTestcases + seaiceStressTestcases

    print("Individual test cases")
    print("=====================")
    for testcase in testcases:

        print("   ", testcase)

        try:
            os.remove("namelist.seaice")
        except:
            pass
        try:
            os.remove("streams.seaice")
        except:
            pass

        os.symlink("namelist.seaice.%s" %(testcase),"namelist.seaice")
        os.symlink("streams.seaice.%s" %(testcase),"streams.seaice")

        # remove output from previous runs
        shutil.rmtree("output", ignore_errors=True)
        shutil.rmtree("output_%s" %(testcase), ignore_errors=True)

        cmd = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        os.system(cmd)

        os.rename("output", "output_%s" %(testcase))

        plot_testcase(testcase)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
