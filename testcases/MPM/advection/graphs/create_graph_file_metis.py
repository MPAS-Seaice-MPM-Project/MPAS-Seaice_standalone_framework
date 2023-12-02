import argparse
import os

#-------------------------------------------------------------------------------

def create_graph_file_metis(nCells, nProcs):

    MPAS_SEAICE_METIS_PATH = os.environ.get('MPAS_SEAICE_METIS_PATH')
    if (MPAS_SEAICE_METIS_PATH is None):
        MPAS_SEAICE_METIS_PATH = "gpmetis"

    cmd = "rm ./graphs/graph.%i.info.part.%i" %(nCells, nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "%s graph.%i.info %i" %(MPAS_SEAICE_METIS_PATH, nCells, nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "mv graph.%i.info.part.%i ./graphs/" %(nCells, nProcs)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', dest='nCells', type=int, required=True)
    parser.add_argument('-n', dest='nProcs', type=int, required=True)

    args = parser.parse_args()

    create_graph_file_metis(args.nCells, args.nProcs)
