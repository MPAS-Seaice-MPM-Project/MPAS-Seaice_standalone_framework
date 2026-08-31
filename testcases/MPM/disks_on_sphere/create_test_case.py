import sys, os
sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical
from create_disks_on_sphere import create_write_particle_file 


def create_testcase():
    resolution = "2562" # {"2562","10242","40962","163842"}
    testcase = "single_particle" # {"single_particle","two_disks","many_disks"}
    
    filenameOut = "particles.nc"
    filenameMesh = "grid.%s.nc" %(resolution)

    print("\nGet Grid files")
    print("=================")
    if (not os.path.exists(filenameMesh)):
        get_testcase_data_spherical()
        os.system("ln -s grid.%s.nc grid.nc" %(resolution))


    print("\nCreate particles file")
    print("================")
    create_write_particle_file(filenameMesh,filenameOut,testcase)

if __name__ == "__main__":
    create_testcase()