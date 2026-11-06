import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *
from datetime import datetime
import cftime

#-------------------------------------------------------------------------

def parallelism(mpasDevelopmentDir,
                SAFrameDirDev,
                domainsDir,
                domain,
                configuration,
                options,
                check,
                oversubscribe,
                np1,
                np2):

    # start time
    startTimeStr = "2000-01-01_00:00:00"
    if ("startTime" in options.keys()):
        startTimeStr = options["startTime"]
    dt = datetime.strptime(startTimeStr, "%Y-%m-%d_%H:%M:%S")
    startTime = cftime.DatetimeNoLeap(dt.year, dt.month,  dt.day,
                                      dt.hour, dt.minute, dt.second)

    # run duration
    runDurationInterval = "HOURS:24"
    if ("runDuration" in options.keys()):
        runDurationInterval = options["runDuration"]
    runDuration, runDurationStr = run_duration(runDurationInterval)

    # compare time
    compareTime = startTime + runDuration
    compareTimeStr = compareTime.isoformat(sep="_", timespec="seconds").replace(":",".")

    # find available directory name
    iTest = 1
    dirExists = True
    while (dirExists):
        testDir = "parallelism_%i.%s.%s" %(iTest,configuration,domain)
        iTest = iTest + 1
        dirExists = os.path.isdir(testDir)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Parallelism, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    multipleBlocks = False
    if ("multipleBlocks" in options.keys() and options["multipleBlocks"] == "True"):
        multipleBlocks = True

    print("multipleBlocks: ", multipleBlocks)
    logfile.write("multipleBlocks: %s" %(multipleBlocks))

    # first run
    nProcs = np1

    nmlChanges = {"seaice_model": {"config_start_time":startTimeStr,
                                   "config_run_duration":runDurationStr}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":runDurationStr}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("development1",
                  mpasDevelopmentDir,
                  SAFrameDirDev,
                  domainsDir,
                  domain,
                  configuration,
                  nmlChanges,
                  streamChanges,
                  nProcs,
                  logfile,
                  oversubscribe) != 0):
        run_failed("parallelism")
        os.chdir("..")
        return 1

    # second run
    nProcs = np2

    if (not multipleBlocks):
        nmlChanges = {"seaice_model": {"config_start_time":startTimeStr,
                                       "config_run_duration":runDurationStr}}
    else:
        nmlChanges = {"seaice_model": {"config_start_time":startTimeStr,
                                       "config_run_duration":runDurationStr},
                     "decomposition": {"config_block_decomp_file_prefix":'graphs/graph.info.eq.part.',
                                       "config_number_of_blocks": 96,
                                       "config_explicit_proc_decomp": True,
                                       "config_proc_decomp_file_prefix":'graphs/graph.info.eq_block.part.'}}

    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":runDurationStr}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("development2",
                  mpasDevelopmentDir,
                  SAFrameDirDev,
                  domainsDir,
                  domain,
                  configuration,
                  nmlChanges,
                  streamChanges,
                  nProcs,
                  logfile,
                  oversubscribe) != 0):
        run_failed("parallelism")
        os.chdir("..")
        return 1


    # compare
    restart_file = "restart.%s.nc" %(compareTimeStr)

    file1 = "./development1/restarts/%s" %(restart_file)
    file2 = "./development2/restarts/%s" %(restart_file)

    logfile.write("file1: %s\n" %(file1))
    logfile.write("file2: %s\n" %(file2))

    ignoreVarname = ["cellsOnCell","verticesOnCell","edgesOnEdge","edgesOnCell","localCellIDCreationMP"]

    if (check):
        ignoreVarname.append("testArrayReproducibility")
        ignoreVarname.append("testArrayRestartability")

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile,ignoreVarname)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "parallelism")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
