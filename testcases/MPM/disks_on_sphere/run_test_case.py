import os, sys

sys.path.append("../../../../utils/testcases")
from create_test_case import create_testcase
from execute_model import execute_model

#-------------------------------------------------------------------------------

def run_model():
    
    if(not (os.path.exists('grid.nc') and os.path.exists('particles.nc'))):
        print("grid.nc or particles.nc not found, running 'create_test_case.py' first")
        create_testcase()


    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    nProcs = 1
    logFile = None
    execute_model(MPAS_SEAICE_EXECUTABLE,
                  nProcs,
                  logFile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
