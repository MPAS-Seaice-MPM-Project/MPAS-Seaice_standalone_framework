import sys

sys.path.append("../../utils/testcases")
from log_messages import log_message

import os
import subprocess
import time
import argparse

from create_final_report import create_final_report

#-------------------------------------------------------------------------------

def run_testcases(optional=False):

    # List of test case directories
    testCaseDirs = [
        "advection",
        "advection_convergence",
        "advection_convergence_mapping",
        "assembly",
        "column_on_particles",
        "polynya",
        "spherical_operators/interpolation",
        "spherical_operators/reconstruction",
        "spherical_operators/interpolation_cell",
        "standard_physics"
    ]
    testCaseDirsOptional = [
        "spherical_operators/strain",
        "spherical_operators/strain_stress_divergence"
    ]

    testCaseDirsUse = testCaseDirs
    if (optional):
        testCaseDirsUse += testCaseDirsOptional

    scriptDir = os.path.dirname(os.path.abspath(__file__))

    logFilename = "%s/log_mpm_testcases.txt" %(scriptDir)
    logFile = open(logFilename, "w")
    logFile.write("MPAS-Seaice-MPM test cases\n")
    logFile.write("##########################\n")

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        message = "\nMPAS_SEAICE_EXECUTABLE not specified. Using standard location.\n"
        logFile.write(message)
        print(message)
    else:
        message = "\nUsing MPAS-Seaice-MPM executable: %s\n" %(MPAS_SEAICE_EXECUTABLE)
        logFile.write(message)
        print(message)

    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND == "mpirun" or
        MPAS_SEAICE_TESTCASES_RUN_COMMAND == "srun"):
        message = "Using MPAS_SEAICE_TESTCASES_RUN_COMMAND: %s\n" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND)
        logFile.write(message)
        print(message)
    else:
        message = "Unsupported MPAS_SEAICE_TESTCASES_RUN_COMMAND type: %s\n" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND)
        logFile.write(message)
        print(message)
        message = "  Must be 'mpirun' or 'srun'\n"
        logFile.write(message)
        print(message)
        sys.exit()

    logFile.flush()
    logFile.close()

    for testDir in testCaseDirs:
        start_time = time.perf_counter()
        result = subprocess.run(["python", "run_testcase.py", "-l", logFilename], cwd=testDir)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        if result.returncode != 0:
            logFile = open(logFilename, "a")
            log_message("   Test case failed", "red", doPrint=False, logFile=logFile)
            logFile.flush()
            logFile.close()

        logFile = open(logFilename, "a")
        logFile.write(f"   Elapsed time: {elapsed_time:.4f} seconds\n")
        logFile.flush()
        logFile.close()

    logFile = open(logFilename, "a")
    logFile.write("\nTesting Completed\n")
    logFile.flush()
    logFile.close()

    create_final_report()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Run MPM test cases")

    parser.add_argument('-o', dest='optional', action='store_true', help="Run the optional test cases")

    args = parser.parse_args()

    run_testcases(args.optional)
