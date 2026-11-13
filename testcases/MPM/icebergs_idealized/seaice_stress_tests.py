import glob
import re
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from colorama import Fore, Style

# Idealized tests of the iceberg sea ice stress force
# (config_iceberg_seaice_force_type = "stress"). The variational scheme sea ice
# stresses (stress11var, stress22var, stress12var) are prescribed through the
# initial conditions stream and stay fixed since the velocity solver is off.
# Namelist and streams files are generated from namelist.seaice.seaice_force_mid_conc
# and streams.seaice.seaice_force_mid_conc with only the sea ice force active and
# no sea ice capture, so the iceberg is accelerated only by the sea ice stress force.

#-------------------------------------------------------------------------------

# iceberg size
icebergHeight = 1000.0 # m
icebergLength = 1000.0 # m

# Test definitions
#   stress: (sigma11, sigma22, sigma12) prescribed uniform stress (N/m)
#   pressure: sea ice pressure (N/m), upper limit on the compressive stress
#   seaiceVelocity: (u, v) uniform sea ice velocity (m/s)
seaiceStressTests = {
    # uniaxial compression along the sea ice motion
    "seaice_stress_compression": {"stress": (0.0, -2.0e4, 0.0), "pressure": 3.0e4,
                                  "seaiceVelocity": (0.0, 0.5)},
    # compressive stress limited by the sea ice pressure
    "seaice_stress_pressure_limit": {"stress": (0.0, -2.0e4, 0.0), "pressure": 1.0e4,
                                     "seaiceVelocity": (0.0, 0.5)},
    # tension along the sea ice motion gives no force
    "seaice_stress_tension": {"stress": (0.0, 2.0e4, 0.0), "pressure": 3.0e4,
                              "seaiceVelocity": (0.0, 0.5)},
    # compression perpendicular to the sea ice motion gives no force
    "seaice_stress_perpendicular": {"stress": (-2.0e4, 0.0, 0.0), "pressure": 3.0e4,
                                    "seaiceVelocity": (0.0, 0.5)},
    # isotropic compression with sea ice moving at 45 degrees, tests the stress
    # tensor projection onto a general contact normal
    "seaice_stress_isotropic_diagonal": {"stress": (-2.0e4, -2.0e4, 0.0), "pressure": 3.0e4,
                                         "seaiceVelocity": (0.5/sqrt(2.0), 0.5/sqrt(2.0))},
    # pure shear with sea ice moving at 45 degrees: normal stress -sigma12
    "seaice_stress_shear_diagonal": {"stress": (0.0, 0.0, -2.0e4), "pressure": 3.0e4,
                                     "seaiceVelocity": (0.5/sqrt(2.0), 0.5/sqrt(2.0))},
}

seaiceStressTestcases = list(seaiceStressTests.keys())

#-------------------------------------------------------------------------------

def read_namelist(filename):

    values = {}
    for line in open(filename,"r"):
        match = re.match(r"^\s*(config_\w+)\s*=\s*(.+?)\s*$", line)
        if (match is not None):
            values[match.group(1)] = match.group(2).strip("'\"")

    return values

#-------------------------------------------------------------------------------

def replace_once(text, old, new):

    if (text.count(old) != 1):
        raise Exception("Expected one occurrence of: %s" %(old))
    return text.replace(old, new)

#-------------------------------------------------------------------------------

def write_namelist(testcase):

    settings = {
        "config_run_duration": "'00-00-01_00:00:00'",
        "config_stress_divergence_scheme": "'variational'",
        "config_iceberg_use_ocean_drag_force": "false",
        "config_iceberg_use_air_drag_force": "false",
        "config_iceberg_use_coriolis_force": "false",
        "config_iceberg_use_surface_tilt_force": "false",
        "config_iceberg_use_wave_radiation_force": "false",
        "config_iceberg_use_seaice_force": "true",
        "config_iceberg_seaice_force_type": "\"stress\"",
        "config_iceberg_seaice_capture_type": "\"none\"",
        "config_iceberg_apply_seaice_force_to_seaice": "true",
        "config_iceberg_grounding_type": "\"none\""}

    lines = open("namelist.seaice.seaice_force_mid_conc","r").read().split("\n")

    for key in settings:
        iLines = [iLine for iLine, line in enumerate(lines) if re.match(r"^\s*%s\s*=" %(key), line)]
        if (len(iLines) != 1):
            raise Exception("Expected one occurrence of %s in namelist.seaice.seaice_force_mid_conc" %(key))
        lines[iLines[0]] = "    %s = %s" %(key, settings[key])

    fileout = open("namelist.seaice.%s" %(testcase),"w")
    fileout.write("\n".join(lines))
    fileout.close()

#-------------------------------------------------------------------------------

def write_streams(testcase):

    streams = open("streams.seaice.seaice_force_mid_conc","r").read()

    streams = replace_once(streams,
                           'filename_template="ic_seaice_force_mid_conc.nc"',
                           'filename_template="ic_%s.nc"' %(testcase))
    streams = replace_once(streams,
                           '        <var name="icePressure"/>\n</stream>',
                           '        <var name="icePressure"/>\n'
                           '        <var name="stress11var"/>\n'
                           '        <var name="stress22var"/>\n'
                           '        <var name="stress12var"/>\n'
                           '</stream>')

    # hourly iceberg output so each time step is output, with sea ice stress at the iceberg
    iStart = streams.index('<stream name="icebergsOut"')
    iEnd = streams.index('</stream>', iStart)
    block = replace_once(streams[iStart:iEnd],
                         'output_interval="00-00-01_00:00:00"',
                         'output_interval="00-00-00_01:00:00"')
    for varname in ["seaIceStress11Iceberg", "seaIceStress22Iceberg", "seaIceStress12Iceberg",
                    "seaIcePressureIceberg", "uVelocitySeaiceIceberg", "vVelocitySeaiceIceberg"]:
        block += '        <var name="%s"/>\n' %(varname)
    streams = streams[:iStart] + block + streams[iEnd:]

    fileout = open("streams.seaice.%s" %(testcase),"w")
    fileout.write(streams)
    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file(testcase, test, nCells, nVertices, maxEdges):

    fileout = Dataset("ic_%s.nc" %(testcase),"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)
    fileout.createDimension("maxEdges",maxEdges)

    var = fileout.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = test["seaiceVelocity"][0]

    var = fileout.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = test["seaiceVelocity"][1]

    var = fileout.createVariable("iceAreaCell","d",dimensions=["nCells"])
    var[:] = 0.99

    var = fileout.createVariable("iceVolumeCell","d",dimensions=["nCells"])
    var[:] = 0.99

    var = fileout.createVariable("icePressure","d",dimensions=["nCells"])
    var[:] = test["pressure"]

    # uniform stress at all cell corners
    for varname, stress in zip(["stress11var", "stress22var", "stress12var"], test["stress"]):
        var = fileout.createVariable(varname,"d",dimensions=["nCells","maxEdges"])
        var[:] = stress

    fileout.close()

#-------------------------------------------------------------------------------

def write_seaice_stress_test_files(nCells, nVertices, maxEdges):

    for testcase in seaiceStressTestcases:

        write_namelist(testcase)
        write_streams(testcase)
        write_ic_file(testcase, seaiceStressTests[testcase], nCells, nVertices, maxEdges)

#-------------------------------------------------------------------------------

def expected_force(test, uVelocityIceberg, vVelocityIceberg, icebergMass, icebergWidth, dt):

    # expected sea ice stress force on the iceberg for the start of step iceberg velocity
    uSeaice, vSeaice = test["seaiceVelocity"]
    sigma11, sigma22, sigma12 = test["stress"]

    du = uSeaice - uVelocityIceberg
    dv = vSeaice - vVelocityIceberg
    relativeSpeed = sqrt(du**2 + dv**2)
    if (relativeSpeed <= 1.0e-6):
        return 0.0, 0.0

    nu = du / relativeSpeed
    nv = dv / relativeSpeed
    normalStress = nu*(sigma11*nu + sigma12*nv) + nv*(sigma12*nu + sigma22*nv)
    compressiveStress = min(test["pressure"], max(0.0, -normalStress))
    forceMagnitude = min(icebergWidth * compressiveStress, icebergMass * relativeSpeed / dt)

    return forceMagnitude * nu, forceMagnitude * nv

#-------------------------------------------------------------------------------

def read_output(testcase):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))
    if (len(filenames) == 0):
        return None

    varnames = ["uVelocityIceberg", "vVelocityIceberg", "uSeaiceForceIceberg", "vSeaiceForceIceberg",
                "seaIceStress11Iceberg", "seaIceStress22Iceberg", "seaIceStress12Iceberg",
                "icebergVolume", "icebergHeight", "icebergLength"]

    output = {"times": [], "forceOnSeaiceU": [], "forceOnSeaiceV": []}
    for varname in varnames:
        output[varname] = []

    for filename in filenames:
        filein = Dataset(filename,"r")
        output["times"].append(filein.variables["daysSinceStartOfSim"][0] * 24.0)
        for varname in varnames:
            output[varname].append(filein.variables[varname][0,0])
        areaTriangle = filein.variables["areaTriangle"][:]
        output["forceOnSeaiceU"].append(np.sum(filein.variables["icebergForceU"][0,:] * areaTriangle))
        output["forceOnSeaiceV"].append(np.sum(filein.variables["icebergForceV"][0,:] * areaTriangle))
        filein.close()

    for varname in output:
        output[varname] = np.array(output[varname], dtype=float)

    return output

#-------------------------------------------------------------------------------

def check_seaice_stress_test(testcase, test, physics, output):

    checks = []

    def check(name, passed, detail):
        checks.append((name, bool(passed), detail))

    nSteps = len(output["times"])
    icebergMass = output["icebergVolume"][0] * physics["icebergDensity"]
    icebergWidth = output["icebergVolume"][0] / (output["icebergLength"][0] * output["icebergHeight"][0])

    finite = all(np.all(np.isfinite(output[varname])) for varname in output)
    check("finite", finite, "no NaN or Inf in output")

    # uniform stress averaged to vertices and interpolated to the iceberg is unchanged
    for iComponent, component in enumerate(["11", "22", "12"]):
        stressIceberg = output["seaIceStress%sIceberg" %(component)]
        error = np.max(np.abs(stressIceberg - test["stress"][iComponent]))
        check("stress%s at iceberg" %(component), error <= 1.0e-6 * max(1.0, abs(test["stress"][iComponent])),
              "max |sigma%s - prescribed| = %g N/m" %(component, error))

    # force for each step from the start of step iceberg velocity
    maxError = 0.0
    maxForce = 0.0
    for iStep in range(nSteps):
        if (iStep == 0):
            uInit, vInit = 0.0, 0.0
        else:
            uInit, vInit = output["uVelocityIceberg"][iStep-1], output["vVelocityIceberg"][iStep-1]
        uExpected, vExpected = expected_force(test, uInit, vInit, icebergMass, icebergWidth, physics["dt"])
        error = sqrt((output["uSeaiceForceIceberg"][iStep] - uExpected)**2 +
                     (output["vSeaiceForceIceberg"][iStep] - vExpected)**2)
        maxError = max(maxError, error)
        maxForce = max(maxForce, sqrt(uExpected**2 + vExpected**2))
    check("sea ice stress force", maxError <= 1.0e-6 * max(1.0, maxForce),
          "max |F - expected| = %g N, max expected |F| = %g N" %(maxError, maxForce))

    # iceberg momentum balance: only the sea ice force acts
    uInits = np.concatenate(([0.0], output["uVelocityIceberg"][:-1]))
    vInits = np.concatenate(([0.0], output["vVelocityIceberg"][:-1]))
    uAcceleration = icebergMass * (output["uVelocityIceberg"] - uInits) / physics["dt"]
    vAcceleration = icebergMass * (output["vVelocityIceberg"] - vInits) / physics["dt"]
    error = np.max(np.sqrt((uAcceleration - output["uSeaiceForceIceberg"])**2 +
                           (vAcceleration - output["vSeaiceForceIceberg"])**2))
    check("momentum balance", error <= 1.0e-6 * max(1.0, maxForce),
          "max |M du/dt - F| = %g N" %(error))

    # equal and opposite force on the sea ice
    error = np.max(np.sqrt((output["forceOnSeaiceU"] + output["uSeaiceForceIceberg"])**2 +
                           (output["forceOnSeaiceV"] + output["vSeaiceForceIceberg"])**2))
    check("reaction force on sea ice", error <= 1.0e-6 * max(1.0, maxForce),
          "max |F_seaice + F_iceberg| = %g N" %(error))

    # iceberg never moves faster than the sea ice
    seaiceSpeed = sqrt(test["seaiceVelocity"][0]**2 + test["seaiceVelocity"][1]**2)
    icebergSpeed = np.sqrt(output["uVelocityIceberg"]**2 + output["vVelocityIceberg"]**2)
    check("no overshoot of sea ice velocity", np.max(icebergSpeed) <= seaiceSpeed * (1.0 + 1.0e-10),
          "max iceberg speed %.6f m/s, sea ice speed %.6f m/s" %(np.max(icebergSpeed), seaiceSpeed))

    return checks

#-------------------------------------------------------------------------------

def plot_seaice_stress_test(testcase):

    output = read_output(testcase)
    if (output is None):
        print("    %s: no output" %(testcase))
        return False

    test = seaiceStressTests[testcase]
    namelist = read_namelist("namelist.seaice.%s" %(testcase))
    physics = {"icebergDensity": float(namelist["config_iceberg_density"]),
               "dt": float(namelist["config_dt"])}

    checks = check_seaice_stress_test(testcase, test, physics, output)

    allPassed = all(passed for name, passed, detail in checks)
    for name, passed, detail in checks:
        print("        %s: %s: %s (%s)" %(testcase, Fore.GREEN+"PASS"+Style.RESET_ALL if passed else Fore.RED+"FAIL"+Style.RESET_ALL, name, detail))

    times = output["times"]

    fig, axes = plt.subplots(2, 1, figsize=(7,7), sharex=True)

    axes[0].plot(times, output["uVelocityIceberg"], color="red", marker="x", label="uVelocityIceberg")
    axes[0].plot(times, output["vVelocityIceberg"], color="green", marker="x", label="vVelocityIceberg")
    axes[0].axhline(test["seaiceVelocity"][0], color="red", linestyle="dotted", label="uVelocity sea ice")
    axes[0].axhline(test["seaiceVelocity"][1], color="green", linestyle="dotted", label="vVelocity sea ice")
    axes[0].set_ylabel("Velocity (m/s)")
    axes[0].legend()

    axes[1].plot(times, output["uSeaiceForceIceberg"], color="red", marker="x", label="u sea ice force on iceberg")
    axes[1].plot(times, output["vSeaiceForceIceberg"], color="green", marker="x", label="v sea ice force on iceberg")
    axes[1].plot(times, -output["forceOnSeaiceU"], color="red", linestyle="dashed", label="-u force on sea ice")
    axes[1].plot(times, -output["forceOnSeaiceV"], color="green", linestyle="dashed", label="-v force on sea ice")
    axes[1].set_ylabel("Force (N)")
    axes[1].set_xlabel("Time (hours)")
    axes[1].legend()

    axes[0].set_title("%s: %s" %(testcase, "PASS" if allPassed else "FAIL"))

    plt.tight_layout()
    plt.savefig("results_%s.png" %(testcase),dpi=300)
    plt.close(fig)

    return allPassed
