import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical
from log_messages import log_message

from reconstruction_scaling import reconstruction_scaling
from reconstruction_map import reconstruction_map
from create_ic import create_ic

sys.path.append("../reconstruction_cell")
from run_model import run_model
from create_particles import create_particles

import os
import argparse

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nReconstruction test case\n")
        logFile.write(  "========================\n")
        logFile.flush()
    else:
        logFile = None

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    create_particles()

    run_model(logFile)

    reconstruction_map()

    reconstruction_scaling()

    if (logFile is not None):
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        log_message("Check plots", "yellow", logFile=logFile)
        log_message("  "+scriptDir+"/reconstruction_map.png", "magenta", logFile=logFile)
        log_message("  "+scriptDir+"/reconstruction_scaling.png", "magenta", logFile=logFile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='logFilename')
    args = parser.parse_args()

    run_testcase(args.logFilename)
