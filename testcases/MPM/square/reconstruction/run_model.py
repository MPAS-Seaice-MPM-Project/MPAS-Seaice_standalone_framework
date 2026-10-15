import sys

sys.path.append("../../../utils/testcases")
from log_messages import log_message
from execute_model import execute_model

import os

#-------------------------------------------------------------------------------

def run_model(logFile=None):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    gridTypes = ["hex","quad"]
    #gridTypes = ["hex"]

    grids = {"hex" :["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}
    gridTypes = ["hex","quad"]
    #gridTypes = ["hex"]

    for gridType in gridTypes:

       print("Grid type: ", gridType)

       for grid in grids[gridType]:

          print("    Grid: ", grid)

          os.system("rm grid.nc")
          os.system("rm ic.nc")
          os.system("rm particles.nc")
          os.system("ln -s grid_%s_%s.nc grid.nc" %(gridType,grid))
          os.system("ln -s ic_%s_%s.nc ic.nc" %(gridType,grid))
          os.system("ln -s particles_%s_%s.nc particles.nc" %(gridType,grid))

          execute_model(MPAS_SEAICE_EXECUTABLE,
                        1,
                        None)

          os.system("rm -rf output_%s_%s" %(gridType, grid))
          os.system("mv output output_%s_%s" %(gridType, grid))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
