import os
import argparse
import sys
sys.path.append("../../../testing")
from testing_utils import add_pio_namelist_changes, create_new_namelist

#-------------------------------------------------------------------------------

def run_model(nProcs):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    MPAS_SEAICE_METIS_PATH = os.environ.get('MPAS_SEAICE_METIS_PATH')
    if (MPAS_SEAICE_METIS_PATH is None):
        MPAS_SEAICE_METIS_PATH = "gpmetis"

    if (not os.path.isdir("output")):
        os.mkdir("output")

    cmd = "rm -rf output_%i" %(nProcs)
    print(cmd)
    os.system(cmd)

    nmlChanges = {}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)
    create_new_namelist("namelist.seaice.default", "namelist.seaice", nmlChanges)

    cmd = "mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    #cmd = "srun --nodes=1 --cpus-per-task=1 --ntasks-per-node=%i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    print(cmd)
    os.system(cmd)

    cmd = "mv output output_%i" %(nProcs)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest="nProcs", type=int)
    args = parser.parse_args()

    run_model(args.nProcs)
