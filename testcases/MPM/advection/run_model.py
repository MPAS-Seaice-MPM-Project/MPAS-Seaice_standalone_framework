import os
import argparse

#-------------------------------------------------------------------------------

def run_model(nCells, nProcs):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    cmd = "rm grid.nc ic.nc log.seaice.*"
    print(cmd)
    os.system(cmd)

    cmd = "ln -s grid.%i.nc grid.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s ic_slotted_cylinder_%i.nc ic.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    if (nProcs > 1):
        os.chdir("./graphs")

        cmd = "ln -s graph.%i.info.part.%i graph.info.part.%i" %(nCells, nProcs, nProcs)
        print(cmd)
        os.system(cmd)

        os.chdir("..")

    cmd = "mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    print(cmd)
    os.system(cmd)

    cmd = "mv output output_%i" %(nProcs)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells", type=int, required=True)
    parser.add_argument('-p', dest="nProcs", type=int, required=True)
    args = parser.parse_args()

    run_model(args.nCells, args.nProcs)
