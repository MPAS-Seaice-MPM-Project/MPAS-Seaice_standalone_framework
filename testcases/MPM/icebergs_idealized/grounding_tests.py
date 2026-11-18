import glob
import re
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt, atan, sin, cos
from colorama import Fore, Style

# Idealized tests of the individual iceberg grounding forces
# (config_iceberg_grounding_type = "force"). Namelist and streams files are
# generated from namelist.seaice.grounding_force and streams.seaice.grounding_force
# with only ocean drag and the selected grounding forces active.

#-------------------------------------------------------------------------------

# constants matching the standalone model (mpas_seaice_constants.F)
densitySeaWater = 1026.0 # kg/m^3
gravity = 9.80616 # m/s^2

# iceberg size
icebergHeight = 1000.0 # m
icebergLength = 1000.0 # m

# iceberg start latitude (icebergs_midlatitude.nc)
latIcebergStart = pi/4.0

# Test definitions
#   forces: (sediment, gravity, solid) force switches
#   bathymetry: "deep" (no grounding), "flat" or "ramp" (depth increasing northwards)
#   penetration: iceberg draft minus ramp/flat ocean depth at the iceberg start latitude (m)
#   slope: ramp ocean depth gradient (m/m)
#   oceanVelocity: uniform northward ocean velocity (m/s)
#   icebergVelocity: initial northward iceberg velocity (m/s)
groundingTests = {
    "grnd_none_ref":         {"groundingType": "none",  "forces": (True,  True,  True),  "bathymetry": "deep", "oceanVelocity": 0.5},
    "grnd_none":             {"groundingType": "force", "forces": (True,  True,  True),  "bathymetry": "deep", "oceanVelocity": 0.5},
    "grnd_rest":             {"groundingType": "force", "forces": (True,  True,  True),  "bathymetry": "flat", "penetration": 20.0, "oceanVelocity": 0.0},
    "grnd_sediment_stuck":   {"groundingType": "force", "forces": (True,  False, False), "bathymetry": "flat", "penetration": 4.0,  "oceanVelocity": 0.1},
    "grnd_sediment_stop":    {"groundingType": "force", "forces": (True,  False, False), "bathymetry": "flat", "penetration": 20.0, "oceanVelocity": 0.0, "icebergVelocity": 0.1},
    "grnd_sediment_partial": {"groundingType": "force", "forces": (True,  False, False), "bathymetry": "flat", "penetration": 4.0,  "oceanVelocity": 0.5},
    "grnd_sediment_full_10": {"groundingType": "force", "forces": (True,  False, False), "bathymetry": "flat", "penetration": 10.0, "oceanVelocity": 0.5},
    "grnd_sediment_full_20": {"groundingType": "force", "forces": (True,  False, False), "bathymetry": "flat", "penetration": 20.0, "oceanVelocity": 0.5},
    "grnd_solid":            {"groundingType": "force", "forces": (False, False, True),  "bathymetry": "flat", "penetration": 8.01, "oceanVelocity": 0.5},
    "grnd_gravity":          {"groundingType": "force", "forces": (False, True,  False), "bathymetry": "ramp", "penetration": 13.0, "slope": 2.0e-4, "oceanVelocity": 0.0},
    "grnd_solid_upslope":    {"groundingType": "force", "forces": (False, False, True),  "bathymetry": "ramp", "penetration": 8.15, "slope": 2.0e-6, "oceanVelocity": -2.0},
    "grnd_solid_downslope":  {"groundingType": "force", "forces": (False, False, True),  "bathymetry": "ramp", "penetration": 8.15, "slope": 2.0e-6, "oceanVelocity": 2.0},
    # interpolation of the ramp to the iceberg increases the start penetration by ~0.47 m
    # for a 2e-4 slope, so the nominal penetration gives ~8.05 m at the iceberg
    "grnd_all_ramp":         {"groundingType": "force", "forces": (True,  True,  True),  "bathymetry": "ramp", "penetration": 7.58, "slope": 2.0e-4, "oceanVelocity": 1.0},
}

groundingTestcases = list(groundingTests.keys())

#-------------------------------------------------------------------------------

def read_namelist(filename):

    values = {}
    for line in open(filename,"r"):
        match = re.match(r"^\s*(config_\w+)\s*=\s*(.+?)\s*$", line)
        if (match is not None):
            values[match.group(1)] = match.group(2).strip("'\"")

    return values

#-------------------------------------------------------------------------------

class GroundingPhysics:

    # Iceberg properties and forces as calculated by the model

    def __init__(self, namelist):

        self.icebergDensity = float(namelist["config_iceberg_density"])
        self.earthRadius = float(namelist["config_earth_radius"])
        self.dt = float(namelist["config_dt"])

        self.sedimentDepth = float(namelist["config_iceberg_sediment_depth"])
        self.scourWidth = float(namelist["config_iceberg_scour_width"])
        self.sedimentShearStrength = float(namelist["config_iceberg_sediment_shear_strength"])
        self.sedimentDensity = float(namelist["config_iceberg_sediment_density"])
        self.staticSpeed = float(namelist["config_iceberg_static_speed"])

        icebergWidth = icebergLength / float(namelist["config_iceberg_initial_length_to_width_ratio"])
        icebergVolume = icebergLength * icebergHeight * icebergWidth
        self.icebergMass = icebergVolume * self.icebergDensity

        if (namelist["config_iceberg_size_type"] == "angle_integral"):
            icebergSize = (2.0 / pi) * (icebergLength + icebergWidth)
        else:
            icebergSize = 0.5 * (icebergLength + icebergWidth)

        self.draft = icebergHeight * (self.icebergDensity / densitySeaWater)

        formArea = self.draft * icebergSize
        skinArea = icebergLength * icebergWidth
        self.oceanDragCoefficient = densitySeaWater * (
            float(namelist["config_iceberg_ocean_form_drag_coefficient"]) * formArea * 0.5 +
            float(namelist["config_iceberg_ocean_skin_drag_coefficient"]) * skinArea)

    def scour_depth(self, penetration):
        return np.clip(penetration, 0.0, self.sedimentDepth)

    def freeboard_anomaly(self, penetration):
        return np.maximum(0.0, penetration - self.sedimentDepth)

    def sediment_force(self, scourDepth):
        sedimentUnitWeight = gravity * (self.sedimentDensity - densitySeaWater)
        return 0.5 * sedimentUnitWeight * scourDepth**2 * self.scourWidth + \
               2.0 * self.sedimentShearStrength * scourDepth * self.scourWidth + \
               0.5 * sqrt(2.0) * self.sedimentShearStrength * scourDepth**2

    def buoyancy_deficit_force(self, freeboardAnomaly):
        return self.icebergMass * (densitySeaWater / self.icebergDensity) * \
            (freeboardAnomaly / icebergHeight) * gravity

    def gravity_force(self, freeboardAnomaly, slope):
        angleGradient = atan(slope)
        return self.buoyancy_deficit_force(freeboardAnomaly) * sin(angleGradient) * cos(angleGradient)

    def ocean_drag(self, relativeSpeed):
        return self.oceanDragCoefficient * relativeSpeed**2

    def steady_speed(self, oceanVelocity, resistingForce):
        # iceberg speed where ocean drag balances a constant resisting force
        return abs(oceanVelocity) - sqrt(resistingForce / self.oceanDragCoefficient)

#-------------------------------------------------------------------------------

def replace_once(text, old, new):

    if (text.count(old) != 1):
        raise Exception("Expected one occurrence of: %s" %(old))
    return text.replace(old, new)

#-------------------------------------------------------------------------------

def write_namelist(testcase, test):

    settings = {
        "config_dt": "3600.0",
        "config_run_duration": "'00-00-05_00:00:00'",
        "config_iceberg_use_surface_wave_ablation": "false",
        "config_iceberg_use_bouyant_convection_ablation": "false",
        "config_iceberg_use_basal_melt_ablation": "false",
        "config_iceberg_use_rollover": "false",
        "config_iceberg_use_coriolis_force": "false",
        "config_iceberg_use_air_drag_force": "false",
        "config_iceberg_use_ocean_drag_force": "true",
        "config_iceberg_use_surface_tilt_force": "false",
        "config_iceberg_use_wave_radiation_force": "false",
        "config_iceberg_use_seaice_force": "false",
        "config_iceberg_grounding_type": "\"%s\"" %(test["groundingType"]),
        "config_iceberg_use_grounding_sediment_force": str(test["forces"][0]).lower(),
        "config_iceberg_use_grounding_gravity_force": str(test["forces"][1]).lower(),
        "config_iceberg_use_grounding_solid_force": str(test["forces"][2]).lower()}

    lines = open("namelist.seaice.grounding_force","r").read().split("\n")

    for key in settings:
        iLines = [iLine for iLine, line in enumerate(lines) if re.match(r"^\s*%s\s*=" %(key), line)]
        if (len(iLines) != 1):
            raise Exception("Expected one occurrence of %s in namelist.seaice.grounding_force" %(key))
        lines[iLines[0]] = "    %s = %s" %(key, settings[key])

    fileout = open("namelist.seaice.%s" %(testcase),"w")
    fileout.write("\n".join(lines))
    fileout.close()

#-------------------------------------------------------------------------------

def write_streams(testcase):

    streams = open("streams.seaice.grounding_force","r").read()

    streams = replace_once(streams,
                           'filename_template="ic_grounding.nc"',
                           'filename_template="ic_%s.nc"' %(testcase))
    streams = replace_once(streams,
                           '        <var name="bathymetry"/>\n</stream>',
                           '        <var name="bathymetry"/>\n'
                           '        <var name="uVelocityIceberg"/>\n'
                           '        <var name="vVelocityIceberg"/>\n'
                           '</stream>')

    # hourly iceberg output so each time step is output
    iStart = streams.index('<stream name="icebergsOut"')
    iEnd = streams.index('</stream>', iStart)
    block = replace_once(streams[iStart:iEnd],
                         'output_interval="00-00-01_00:00:00"',
                         'output_interval="00-00-00_01:00:00"')
    streams = streams[:iStart] + block + streams[iEnd:]

    fileout = open("streams.seaice.%s" %(testcase),"w")
    fileout.write(streams)
    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file(testcase, test, physics, nCells, latCell):

    if (test["bathymetry"] == "deep"):
        depth = np.full(nCells, 1100.0)
    elif (test["bathymetry"] == "flat"):
        depth = np.full(nCells, physics.draft - test["penetration"])
    elif (test["bathymetry"] == "ramp"):
        depth = physics.draft - test["penetration"] + \
            test["slope"] * physics.earthRadius * (latCell - latIcebergStart)
        depth = np.maximum(depth, 1.0)

    fileout = Dataset("ic_%s.nc" %(testcase),"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nIcebergs",1)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = test["oceanVelocity"]

    var = fileout.createVariable("bathymetry","d",dimensions=["nCells"])
    var[:] = -depth

    var = fileout.createVariable("uVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = 0.0

    var = fileout.createVariable("vVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = test.get("icebergVelocity", 0.0)

    fileout.close()

#-------------------------------------------------------------------------------

def write_grounding_test_files(nCells, latCell):

    for testcase in groundingTestcases:

        test = groundingTests[testcase]

        write_namelist(testcase, test)
        write_streams(testcase)

        physics = GroundingPhysics(read_namelist("namelist.seaice.%s" %(testcase)))
        write_ic_file(testcase, test, physics, nCells, latCell)

#-------------------------------------------------------------------------------

def read_iceberg_output(testcase):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))
    if (len(filenames) == 0):
        return None

    varnames = ["latIceberg", "lonIceberg", "uVelocityIceberg", "vVelocityIceberg",
                "icebergHeight", "icebergBathymetry"]
    forceNames = ["Grounding", "GrndSediment", "GrndGravity", "GrndSolid", "Ocean"]
    for forceName in forceNames:
        varnames += ["u%sForceIceberg" %(forceName), "v%sForceIceberg" %(forceName)]

    output = {"times": [], "icebergBathymetryGradient": []}
    for varname in varnames:
        output[varname] = []

    for filename in filenames:
        filein = Dataset(filename,"r")
        output["times"].append(filein.variables["daysSinceStartOfSim"][0])
        for varname in varnames:
            output[varname].append(filein.variables[varname][0,0])
        output["icebergBathymetryGradient"].append(filein.variables["icebergBathymetryGradient"][0,0,:])
        filein.close()

    for varname in output:
        output[varname] = np.array(output[varname], dtype=float)

    for forceName in forceNames:
        output["mag%s" %(forceName)] = np.sqrt(output["u%sForceIceberg" %(forceName)]**2 +
                                               output["v%sForceIceberg" %(forceName)]**2)

    output["speed"] = np.sqrt(output["uVelocityIceberg"]**2 + output["vVelocityIceberg"]**2)

    return output

#-------------------------------------------------------------------------------

def angle_between(u1, v1, u2, v2):

    return np.degrees(np.arccos(np.clip((u1*u2 + v1*v2) /
                                        (np.sqrt(u1**2 + v1**2) * np.sqrt(u2**2 + v2**2)), -1.0, 1.0)))

#-------------------------------------------------------------------------------

def relative_error(value, expected):

    return abs(value - expected) / abs(expected)

#-------------------------------------------------------------------------------

def check_grounding_test(testcase, test, physics, output):

    checks = []

    def check(name, passed, detail):
        checks.append((name, bool(passed), detail))

    v = output["vVelocityIceberg"]
    speed = output["speed"]
    vSteady = np.mean(v[-24:])

    # penetration (draft minus depth) from the model bathymetry at the iceberg
    draft = output["icebergHeight"] * (physics.icebergDensity / densitySeaWater)
    penetration = draft + output["icebergBathymetry"]
    freeboardAnomaly = physics.freeboard_anomaly(penetration)

    # tests common to all
    finite = all(np.all(np.isfinite(output[varname])) for varname in output)
    check("finite", finite, "no NaN or Inf in output")

    if (test["groundingType"] == "none" or test["bathymetry"] == "deep"):
        check("zero grounding force", np.max(output["magGrounding"]) == 0.0,
              "max |F_grounding| = %g N" %(np.max(output["magGrounding"])))

    # dissipative forces oppose the velocity at the start of each time step
    uDissipative = output["uGrndSedimentForceIceberg"] + output["uGrndSolidForceIceberg"]
    vDissipative = output["vGrndSedimentForceIceberg"] + output["vGrndSolidForceIceberg"]
    moving = (speed[:-1] > physics.staticSpeed) & (np.sqrt(uDissipative[1:]**2 + vDissipative[1:]**2) > 0.0)
    if (np.any(moving)):
        angles = angle_between(uDissipative[1:][moving], vDissipative[1:][moving],
                               output["uVelocityIceberg"][:-1][moving], output["vVelocityIceberg"][:-1][moving])
        check("dissipative opposes motion", np.min(angles) > 179.0,
              "min angle to velocity = %.3f deg" %(np.min(angles)))

    # no reversal of the iceberg velocity
    if (np.max(np.abs(v)) > 0.0):
        sign = np.sign(v[np.abs(v) > 1.0e-12])
        check("no velocity reversal", np.all(sign == sign[0]) if len(sign) > 0 else True,
              "northward velocity keeps the same sign")

    # test specific checks
    if (testcase == "grnd_none"):

        reference = read_iceberg_output("grnd_none_ref")
        if (reference is None):
            check("matches grounding none", False, "output_grnd_none_ref not found")
        else:
            difference = max(np.max(np.abs(output["vVelocityIceberg"] - reference["vVelocityIceberg"])),
                             np.max(np.abs(output["latIceberg"] - reference["latIceberg"])))
            check("matches grounding none", difference == 0.0, "max difference = %g" %(difference))

    elif (testcase in ["grnd_rest", "grnd_sediment_stuck"]):

        check("iceberg stays at rest", np.max(speed) <= 1.0e-12, "max speed = %g m/s" %(np.max(speed)))

        if (testcase == "grnd_sediment_stuck"):
            # static sediment force limited to the ocean drag on the static iceberg
            expected = physics.ocean_drag(test["oceanVelocity"])
            error = np.max([relative_error(force, expected) for force in output["magGrndSediment"][1:]])
            check("static sediment force balances drivers", error < 1.0e-6,
                  "|F_sediment| = %g N, expected %g N" %(output["magGrndSediment"][1], expected))
            check("unlimited sediment force exceeds drivers",
                  physics.sediment_force(physics.scour_depth(test["penetration"])) > expected,
                  "F_sediment = %g N, drivers %g N" %(physics.sediment_force(physics.scour_depth(test["penetration"])), expected))

    elif (testcase == "grnd_sediment_stop"):

        # grounded iceberg velocity set to zero in a single time step
        check("iceberg stops in one time step", v[0] > 0.0 and np.max(speed[1:]) == 0.0,
              "initial v = %g m/s, max subsequent speed = %g m/s" %(v[0], np.max(speed[1:])))
        expectedForce = physics.sediment_force(physics.scour_depth(test["penetration"]))
        drivers = physics.icebergMass * test["icebergVelocity"] / physics.dt
        check("sediment force exceeds stopping force", expectedForce > drivers,
              "F_sediment = %g N, M v/dt = %g N" %(expectedForce, drivers))

    elif (testcase in ["grnd_sediment_partial", "grnd_sediment_full_10", "grnd_sediment_full_20"]):

        expectedForce = physics.sediment_force(physics.scour_depth(test["penetration"]))
        error = relative_error(output["magGrndSediment"][-1], expectedForce)
        check("sediment force", error < 1.0e-6,
              "|F_sediment| = %g N, expected %g N" %(output["magGrndSediment"][-1], expectedForce))

        expectedSpeed = physics.steady_speed(test["oceanVelocity"], expectedForce)
        error = relative_error(vSteady, expectedSpeed)
        check("steady drift speed", error < 1.0e-3,
              "v = %.6f m/s, expected %.6f m/s" %(vSteady, expectedSpeed))

        if (testcase == "grnd_sediment_full_20"):
            other = read_iceberg_output("grnd_sediment_full_10")
            if (other is None):
                check("sediment force capped", False, "output_grnd_sediment_full_10 not found")
            else:
                error = relative_error(output["magGrndSediment"][-1], other["magGrndSediment"][-1])
                check("sediment force capped", error < 1.0e-10,
                      "penetration 20 m vs 10 m relative difference = %g" %(error))

    elif (testcase == "grnd_solid"):

        expectedForce = physics.buoyancy_deficit_force(physics.freeboard_anomaly(test["penetration"]))
        error = relative_error(output["magGrndSolid"][-1], expectedForce)
        check("solid force", error < 1.0e-6,
              "|F_solid| = %g N, expected %g N" %(output["magGrndSolid"][-1], expectedForce))

        expectedSpeed = physics.steady_speed(test["oceanVelocity"], expectedForce)
        error = relative_error(vSteady, expectedSpeed)
        check("steady drift speed", error < 1.0e-3,
              "v = %.6f m/s, expected %.6f m/s" %(vSteady, expectedSpeed))

    elif (testcase == "grnd_gravity"):

        # bathymetry gradient points southwards (ocean deepens northwards)
        gradient = output["icebergBathymetryGradient"][1]
        magGradient = sqrt(gradient[0]**2 + gradient[1]**2)
        error = relative_error(magGradient, test["slope"])
        check("bathymetry gradient magnitude", error < 0.05,
              "|grad b| = %g, expected %g" %(magGradient, test["slope"]))
        angle = angle_between(gradient[0], gradient[1], 0.0, -1.0)
        check("bathymetry gradient direction", angle < 5.0,
              "angle to south = %.3f deg" %(angle))

        angle = angle_between(output["uGrndGravityForceIceberg"][1], output["vGrndGravityForceIceberg"][1], 0.0, 1.0)
        check("gravity force downslope", angle < 5.0, "angle to north = %.3f deg" %(angle))

        expectedForce = physics.gravity_force(freeboardAnomaly[1], test["slope"])
        error = relative_error(output["magGrndGravity"][1], expectedForce)
        check("gravity force magnitude", error < 0.05,
              "|F_gravity| = %g N, expected %g N" %(output["magGrndGravity"][1], expectedForce))

        check("iceberg slides downslope", output["latIceberg"][-1] > output["latIceberg"][0],
              "latitude change = %g rad" %(output["latIceberg"][-1] - output["latIceberg"][0]))
        check("freeboard anomaly decreases", np.all(np.diff(freeboardAnomaly[1:]) <= 1.0e-12),
              "freeboard anomaly %.3f m -> %.3f m" %(freeboardAnomaly[1], freeboardAnomaly[-1]))

    elif (testcase == "grnd_solid_upslope"):

        criticalAnomaly = physics.ocean_drag(test["oceanVelocity"]) / physics.buoyancy_deficit_force(1.0)
        check("iceberg moves upslope", v[1] < 0.0, "initial v = %.4f m/s" %(v[1]))
        check("freeboard anomaly increases", freeboardAnomaly[-1] > freeboardAnomaly[1],
              "freeboard anomaly %.4f m -> %.4f m" %(freeboardAnomaly[1], freeboardAnomaly[-1]))
        check("freeboard anomaly below critical", freeboardAnomaly[-1] <= criticalAnomaly * (1.0 + 1.0e-3),
              "final %.4f m, critical %.4f m" %(freeboardAnomaly[-1], criticalAnomaly))
        check("iceberg slows", speed[-1] < 0.25 * np.max(speed),
              "final speed %.4f m/s, max %.4f m/s" %(speed[-1], np.max(speed)))

    elif (testcase == "grnd_solid_downslope"):

        check("iceberg moves downslope", v[1] > 0.0, "initial v = %.4f m/s" %(v[1]))
        check("freeboard anomaly reaches zero", freeboardAnomaly[-1] == 0.0,
              "final freeboard anomaly %.4f m" %(freeboardAnomaly[-1]))
        error = relative_error(v[-1], test["oceanVelocity"])
        check("iceberg free drifts", error < 0.01,
              "final v = %.4f m/s, ocean %.4f m/s" %(v[-1], test["oceanVelocity"]))

    elif (testcase == "grnd_all_ramp"):

        iSolidZero = np.argmax((output["magGrndSolid"] == 0.0) & (np.arange(len(v)) > 0))
        iSedimentZero = np.argmax((output["magGrndSediment"] == 0.0) & (np.arange(len(v)) > 0))
        check("solid force ends before sediment force", 0 < iSolidZero < iSedimentZero,
              "solid zero at step %i, sediment zero at step %i" %(iSolidZero, iSedimentZero))
        error = relative_error(v[-1], test["oceanVelocity"])
        check("iceberg free drifts", error < 0.01,
              "final v = %.4f m/s, ocean %.4f m/s" %(v[-1], test["oceanVelocity"]))

    return checks

#-------------------------------------------------------------------------------

def plot_grounding_test(testcase):

    output = read_iceberg_output(testcase)
    if (output is None):
        print("    %s: no output" %(testcase))
        return

    test = groundingTests[testcase]
    physics = GroundingPhysics(read_namelist("namelist.seaice.%s" %(testcase)))

    checks = check_grounding_test(testcase, test, physics, output)

    allPassed = all(passed for name, passed, detail in checks)
    for name, passed, detail in checks:
        print("        %s: %s: %s (%s)" %(testcase, Fore.GREEN+"PASS"+Style.RESET_ALL if passed else Fore.RED+"FAIL"+Style.RESET_ALL, name, detail))

    times = output["times"]
    draft = output["icebergHeight"] * (physics.icebergDensity / densitySeaWater)
    penetration = draft + output["icebergBathymetry"]
    #print(draft,output["icebergBathymetry"],penetration,physics.freeboard_anomaly(penetration))

    fig, axes = plt.subplots(3, 1, figsize=(7,9), sharex=True)

    axes[0].plot(times, output["uVelocityIceberg"], label="uVelocityIceberg")
    axes[0].plot(times, output["vVelocityIceberg"], label="vVelocityIceberg")
    axes[0].axhline(test["oceanVelocity"], color="grey", linestyle="dotted", label="vOceanVelocity")
    axes[0].set_ylabel("Velocity (m/s)")
    axes[0].legend()

    axes[1].plot(times, physics.scour_depth(penetration), label="scour depth")
    axes[1].plot(times, physics.freeboard_anomaly(penetration), label="freeboard anomaly")
    axes[1].set_ylabel("Depth (m)")
    axes[1].legend()
    ax2 = axes[1].twinx()
    ax2.plot(times, output["icebergBathymetry"], label="bathymetry", color="red")
    ax2.set_ylabel("Bathymetry (m)")
    ax2.legend()

    for forceName, label in [("GrndSediment", "sediment"), ("GrndGravity", "gravity"),
                             ("GrndSolid", "solid"), ("Ocean", "ocean drag")]:
        axes[2].plot(times, output["mag%s" %(forceName)], label=label)
    axes[2].set_yscale("symlog", linthresh=1.0)
    axes[2].set_ylabel("Force magnitude (N)")
    axes[2].set_xlabel("Time (days)")
    axes[2].legend()

    axes[0].set_title("%s: %s" %(testcase, "PASS" if allPassed else "FAIL"))

    plt.tight_layout()
    plt.savefig("results_%s.png" %(testcase),dpi=300)
    plt.close(fig)

    return allPassed
