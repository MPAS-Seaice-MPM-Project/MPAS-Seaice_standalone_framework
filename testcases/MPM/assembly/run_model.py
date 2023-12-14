import os
import argparse

#-------------------------------------------------------------------------------

def run_model(nCells, nProcs):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    MPAS_SEAICE_METIS_PATH = os.environ.get('MPAS_SEAICE_METIS_PATH')
    if (MPAS_SEAICE_METIS_PATH is None):
        MPAS_SEAICE_METIS_PATH = "gpmetis"

    if (not os.path.isdir("output")):
        os.mkdir("output")

    cmd = "rm grid.nc ic.nc graph.info.part.%i" %(nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s grid.%i.nc grid.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s ic.%i.nc ic.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "rm graph.%i.info.part.%i" %(nCells, nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "%s graph.%i.info %i" %(MPAS_SEAICE_METIS_PATH, nCells, nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s graph.%i.info.part.%i graph.info.part.%i" %(nCells, nProcs, nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "rm -rf output_%i" %(nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    print(cmd)
    os.system(cmd)

    cmd = "mv output output_%i" %(nProcs)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells", type=int)
    parser.add_argument('-n', dest="nProcs", type=int)
    args = parser.parse_args()

    run_model(args.nCells, args.nProcs)
