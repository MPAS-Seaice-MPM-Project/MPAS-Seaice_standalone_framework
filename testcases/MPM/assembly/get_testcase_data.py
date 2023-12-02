import subprocess
import os

#-------------------------------------------------------------------------------

def get_testcase_data():

    reses = [2562, 10242, 40962, 163842]

    dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/strain_stress_divergence/"

    for res in reses:

        filename = "grid.%i.nc" %(res)
        if (not os.path.isfile(filename)):

            args = ["wget", dirName+filename]

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            while process.poll() is None:
                line = process.stdout.readline()

        filename = "graph.%i.info" %(res)
        if (not os.path.isfile(filename)):

            args = ["wget", dirName+filename]

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            while process.poll() is None:
                line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data()
