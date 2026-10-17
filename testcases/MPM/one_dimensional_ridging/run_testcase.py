import sys

sys.path.append("../../../utils/testcases/")
from create_square_quad_mesh import create_square_quad_mesh
from create_square_hex_mesh import create_square_hex_mesh
from plot_mesh import plot_mesh

sys.path.append("../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file

from create_ics import create_ics
from run_model import run_model
from plot_heatmap import plot_heatmap

import os
import argparse

lxDefault = 1000000.0
lyDefault = 100000.0
nxDefault = 100
nyDefault = 10
dcDefault = 10000.0

#-------------------------------------------------------------------------------

def run_testcase(runtype,
                 meshtype,
                 lx=lxDefault,
                 ly=lyDefault,
                 nx=nxDefault,
                 ny=nyDefault,
                 dc=dcDefault):

    print("Create grid...")
    x0 = 0.0
    y0 = 0.0

    if (meshtype == "quad"):

        gridFilename = create_square_quad_mesh(nx, ny,
                                               lx, ly,
                                               x0, y0)

    elif (meshtype == "hex"):

        gridFilename = create_square_hex_mesh(dc,
                                              lx, ly,
                                              x0, y0)

    else:
        raise Exception("Unknown mesh type: "+meshtype)

    plot_mesh(gridFilename)

    print("Create ICs...")
    create_ics(gridFilename)

    print("Create particles...")
    particleInitType = "number"
    particleInitNumber = "9"
    particlePositionInitType = "even"
    create_particles_from_cell_file(gridFilename,
                                    "ic.nc",
                                    particleInitType,
                                    particleInitNumber,
                                    particlePositionInitType,
                                    "particles.nc")

    print("Sym link grid...")
    if (os.path.isfile("grid.nc")):
        os.system("rm grid.nc")
    os.symlink(gridFilename,"grid.nc")

    print("Run model...")
    run_model(runtype)

    print("Plot heatmap...")
    plot_heatmap("grid.nc",
                 "output_%s/output.2000.nc" %(runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-r', dest="runtype",  choices=["orig","mpm.decoh","mpm.evp"], default="orig")
    parser.add_argument('-m', dest="meshtype", choices=["quad","hex"], default="quad")
    parser.add_argument('--lx', dest="lx", type=float, default=lxDefault)
    parser.add_argument('--ly', dest="ly", type=float, default=lyDefault)
    parser.add_argument('--nx', dest="nx", type=int,   default=nxDefault)
    parser.add_argument('--ny', dest="ny", type=int,   default=nyDefault)
    parser.add_argument('--dc', dest="dc", type=float, default=dcDefault)

    args = parser.parse_args()

    run_testcase(args.runtype,
                 args.meshtype,
                 args.lx,
                 args.ly,
                 args.nx,
                 args.ny,
                 args.dc)
