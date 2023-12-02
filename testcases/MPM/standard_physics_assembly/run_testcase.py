from run_model import run_model
from check_results import check_results
from check_particle_positions_nprocs import check_particle_positions_nprocs
import sys
sys.path.append("../../../testing")
from compare_mpas_files import compare_files
from testing_utils import get_domain, print_colour, create_test_directory
import os

# domains directory
domainsDir = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
if (domainsDir == None):
        print("Environment variable MPAS_SEAICE_DOMAINS_DIR must be set if no domains directory specified")
        sys.exit()
if (not os.path.exists(domainsDir)):
        print("Requested domains directory does not exist")
        sys.exit()

# get domain
domain="domain_QU120km"
get_domain(domainsDir, domain)

n1 = 1
n16 = 16
n32 = 32

run_model(1)

run_model(n16)

run_model(n32)

file1="output_1/output.2000.nc"

file2="output_16/output.2000.nc"

file3="output_32/output.2000.nc"

# make a test directory
create_test_directory("testDir1_16")
os.chdir("../")

# make log file
logfile = open("log_test.txt", "w")
title = "Test: Parallelism, n1=1, n2=16\n"
print_colour(title, "title")
logfile.write(title)

# run comparison
compare_files(file1,file2,logfile)

if (os.path.isfile("vars_differ.nc")):
     cmd = "mv vars_differ.nc testDir1_16"
     os.system(cmd)

check_particle_positions_nprocs(n1,n16)

# make a test directory
create_test_directory("testDir16_32")
os.chdir("../")

# log file
title = "Test: Parallelism, n1=16, n2=32\n"
print_colour(title, "title")
logfile.write(title)

# run comparison
compare_files(file2,file3,logfile)

if (os.path.isfile("vars_differ.nc")):
     cmd = "mv vars_differ.nc testDir16_32"
     os.system(cmd)

check_particle_positions_nprocs(n16,n32)
