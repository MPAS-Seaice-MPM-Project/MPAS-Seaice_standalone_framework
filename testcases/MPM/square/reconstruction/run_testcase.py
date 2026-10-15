import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_square import get_testcase_data_square
from log_messages import log_message

sys.path.append("../../../square/operators_strain_stress_divergence")
from create_grids import create_grids

from create_ics import create_ics
from create_particles import create_particles
from run_model import run_model
from reconstruction_scaling import reconstruction_scaling
from reconstruction_map import reconstruction_map

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

    print('create grids')
    create_grids()

    print('create ic')
    create_ics('linear')

    print('create particles')
    create_particles('linear')

    print('run model')
    run_model(logFile)

    print('reconstruction map')
    reconstruction_map()

    print('reconstruction scaling')
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
