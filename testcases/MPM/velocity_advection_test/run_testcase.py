import sys

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical
from create_ics import create_ics
from create_particles import create_particles
from add_uvVelMP_to_particles_file import add_uvVelMP_to_particles_file
from run_model import run_model
from advection_error_convergence import advection_error_convergence

sys.path.append("../advection")
from check_particles_moved import check_particles_moved

sys.path.append("../../../testing")
from testing_utils import create_new_namelist

import math
import f90nml
import argparse

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\n Velocity Advection convergence test case\n")
        logFile.write(  "===============================\n")
        logFile.flush()
    else:
        logFile = None

    reses = ["2562", "10242", "40962", "163842"]
    #reses = ["2562"]
    forceTypes = ["cos_lat"]
    DynamicsTimeStep = [3600.0, 1800.0, 900.0, 450.0]
    rotateCartesianGrid = False
    earthRadius = 6371229.0

    #usePolympos = [False, True]
    usePolympos = [False]

    print("Get testcase data")
    print("=================")
    get_testcase_data_spherical()

    print("\nCreate ICs")
    print("==========")
    create_ics(earthRadius, rotateCartesianGrid)

    print("\nCreate particles")
    print("================")
    create_particles()

    print("\nAdd uvVelMP to particles file")
    print("================")
    add_uvVelMP_to_particles_file(earthRadius, rotateCartesianGrid)

    print("\nCreate namelists")
    print("================")

    for usePolympo in usePolympos:

        print("usePolympo: ", usePolympo)
        if (usePolympo):
            usePolympoStr = "polympo"
        else:
            usePolympoStr = "nonpolympo"

        r = 0
        for res in reses:

            nmlChanges = {"mpm":{"config_use_mpm_polympo":usePolympo},
                          "seaice_model":{"config_dt":DynamicsTimeStep[r]}}

            new_namelist = "namelist.seaice.%s.%s" %(usePolympoStr, res)
            create_new_namelist("namelist.seaice.advection_convergence", new_namelist, nmlChanges)

            r = r + 1

        print("\nRun models")
        print("==========")
        run_model(usePolympo,
                  logFile)

        print("\nCheck particles moved")
        print("=====================")
        particleFilename1 = "particles_output.2000-01-01_00.00.00.nc"
        particleFilename2 = "particles_output.2000-01-01_04.00.00.nc"

        for forceType in forceTypes:
            for res in reses:

                check_particles_moved("./output_"+forceType+"_"+res+"_"+usePolympoStr+"/"+particleFilename1,
                                      "./output_"+forceType+"_"+res+"_"+usePolympoStr+"/"+particleFilename2,
                                      logFile)

        print("\nAdvection error convergence")
        print("===========================")

        advection_error_convergence(usePolympoStr, logFile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.logFilename)
