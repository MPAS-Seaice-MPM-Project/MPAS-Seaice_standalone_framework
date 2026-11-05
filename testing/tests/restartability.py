import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *
from datetime import datetime
import cftime

#-------------------------------------------------------------------------

def restartability(mpasDevelopmentDir,
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
    runDuration,  runDurationStr  = run_duration(runDurationInterval)
    runDuration2, runDurationStr2 = run_duration(runDurationInterval, 2)

    compareTime = startTime + runDuration2
    compareTimeStr = compareTime.isoformat(sep="_", timespec="seconds")

    # find available directory name
    iTest = 1
    dirExists = True
    while (dirExists):
        testDir = "restartability_%i.%s.%s" %(iTest,configuration,domain)
        iTest = iTest + 1
        dirExists = os.path.isdir(testDir)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Restartability, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    # base run
    nProcs = np1

    nmlChanges = {"seaice_model": {"config_start_time":startTimeStr,
                                   "config_run_duration":runDurationStr2}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":runDurationStr}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("base",
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
        run_failed("restartability")
        os.chdir("..")
        return 1

    # first restart run
    nProcs = np1

    nmlChanges = {"seaice_model": {"config_start_time":startTimeStr,
                                   "config_run_duration":runDurationStr}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":runDurationStr}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("restart1",
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
        run_failed("restartability")
        os.chdir("..")
        return 1

    # restart
    nProcs = np2

    bgcRestart = False
    if ("bgc" in options.keys() and options["bgc"] == "True"):
        bgcRestart = True

    snowModsRestart = False
    if ("snow_tracer_physics" in options.keys() and options["snow_tracer_physics"] == "True"):
        snowModsRestart = True

    if (not bgcRestart):
        if (not snowModsRestart):
             nmlChanges = {"seaice_model": {"config_start_time":"file"},
                           "restart": {"config_do_restart":True}}
        else:
             nmlChanges = {"seaice_model": {"config_start_time":"file"},
                      "restart": {"config_do_restart":True,
                                  "config_do_restart_snow_density":True,
                                  "config_do_restart_snow_grain_radius":True}}
    else:
        nmlChanges = {"seaice_model": {"config_start_time":"file",
                                       "config_run_duration":runDurationStr},
                      "restart": {"config_do_restart":True,
                                  "config_do_restart_bgc":True,
                                  "config_do_restart_hbrine":True}}

    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = []
    if ("restart2ChangeStreamFilename" in options.keys()):
        streamName  = options["restart2ChangeStreamFilename"].split(":")[0]
        newFileName = options["restart2ChangeStreamFilename"].split(":")[1]
        streamChanges.append({"streamName":streamName,
                              "attributeName":"filename_template",
                              "newValue":newFileName})

    os.system("cp -rf restart1 restart2")
    if (restart_model("restart2",
                      nmlChanges,
                      streamChanges,
                      nProcs,
                      logfile,
                      oversubscribe,
                      ["particles.nc"]) != 0):
        run_failed("restartability")
        os.chdir("..")
        return 1


    # compare
    restart_file = "restart.%s.nc" %(compareTimeStr)
    file1 = "./base/restarts/%s" %(restart_file)
    file2 = "./restart2/restarts/%s" %(restart_file)
    logfile.write("file1: %s\n" %(file1))
    logfile.write("file2: %s\n" %(file2))

    ignoreVarname = ["cellsOnCell","verticesOnCell","edgesOnEdge","edgesOnCell","localCellIDCreationMP"]
    if (check):
        ignoreVarname.append("testArrayParallelism")
        ignoreVarname.append("testArrayReproducibility")

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile,ignoreVarname)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "restartability")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
