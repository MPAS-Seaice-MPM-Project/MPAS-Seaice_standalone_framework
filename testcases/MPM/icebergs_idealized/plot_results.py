import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import argparse
from grounding_tests import groundingTestcases, plot_grounding_test
from seaice_stress_tests import seaiceStressTestcases, plot_seaice_stress_test

#-------------------------------------------------------------------------------

def plot_figure(testType, uColor, vColor):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testType)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    latIcebergs = []
    lonIcebergs = []

    uVelocityIcebergs = []
    vVelocityIcebergs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        latIceberg = filein.variables["latIceberg"][0,0]
        lonIceberg = filein.variables["lonIceberg"][0,0]

        uVelocityIceberg = filein.variables["uVelocityIceberg"][0,0]
        vVelocityIceberg = filein.variables["vVelocityIceberg"][0,0]

        filein.close()

        latIcebergs.append(latIceberg)
        lonIcebergs.append(lonIceberg)

        uVelocityIcebergs.append(uVelocityIceberg)
        vVelocityIcebergs.append(vVelocityIceberg)

    fig, axis = plt.subplots()
    axis2 = axis.twinx()

    axis.plot(times, lonIcebergs, color=uColor, linestyle="dashed", label="lonIceberg")
    axis.plot(times, latIcebergs, color=vColor, linestyle="dashed", label="latIceberg")

    axis2.plot(times, uVelocityIcebergs, color=uColor, linestyle="solid", label="uVelocityIceberg")
    axis2.plot(times, vVelocityIcebergs, color=vColor, linestyle="solid", label="vVelocityIceberg")

    axis.legend(loc="upper left")
    axis2.legend(loc="lower right")

    axis.set_title("Position/velocity %s" %(testType))

    axis.set_xlabel("Time (days)")

    axis.set_ylabel("Position (rads)")
    axis2.set_ylabel("Velocity (m/s)")

    plt.tight_layout()

    plt.savefig("results_%s.png" %(testType),dpi=300)

#-------------------------------------------------------------------------------

def plot_coriolis():

    filenames = sorted(glob.glob("./output_coriolis/icebergs_output.*"))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    xs = []
    ys = []
    zs = []
    vs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        posnIBGeo = filein.variables["posnIBGeo"][0,0,:]

        uVelocityIceberg = filein.variables["uVelocityIceberg"][0,0]
        vVelocityIceberg = filein.variables["vVelocityIceberg"][0,0]
        v = sqrt(pow(uVelocityIceberg,2)+
                 pow(vVelocityIceberg,2))
        vs.append(v)

        filein.close()

        xs.append(posnIBGeo[0])
        ys.append(posnIBGeo[1])
        zs.append(posnIBGeo[2])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    vs = np.array(vs)

    # Create sphere parameterization
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Sphere of radius R centered at (cx, cy, cz)
    R = sqrt(xs[0]*xs[0]+ys[0]*ys[0]+zs[0]*zs[0])
    cx, cy, cz = 0.0, 0.0, 0.0

    x = cx + R * np.outer(np.cos(u), np.sin(v))
    y = cy + R * np.outer(np.sin(u), np.sin(v))
    z = cz + R * np.outer(np.ones_like(u), np.cos(v))

    fig, axis = plt.subplots(subplot_kw={"projection": "3d"})

    # Plot transparent surface
    axis.plot_surface(
        x, y, z,
        alpha=0.3,          # transparency (0=fully transparent, 1=opaque)
        linewidth=0,
        antialiased=True
    )

    lats = np.deg2rad(np.arange(-60, 61, 30))
    phi = np.linspace(0, 2*np.pi, 200)

    for lat in lats:
        z = cz + R * np.sin(lat)
        r_xy = R * np.cos(lat)
        x = cx + r_xy * np.cos(phi)
        y = cy + r_xy * np.sin(phi)
        axis.plot(x, y, z, color="grey", linewidth=0.4)

    lons = np.deg2rad(np.arange(0, 360, 30))
    theta = np.linspace(0, np.pi, 200)

    for lon in lons:
        x = cx + R * np.sin(theta) * np.cos(lon)
        y = cy + R * np.sin(theta) * np.sin(lon)
        z = cz + R * np.cos(theta)
        axis.plot(x, y, z, color="grey", linewidth=0.4)

    axis.plot(xs, ys, zs, linewidth=0.5)
    axis.view_init(elev=20, azim=45)

    axis.set_title("Coriolis trajectory")

    plt.tight_layout()

    plt.savefig("results_coriolis.png",dpi=300)

    fig, axis = plt.subplots()

    axis.plot(times, vs)

    axis.set_xlabel("Time (days)")

    axis.set_title("Coriolis speed")

    plt.tight_layout()

    plt.savefig("results_coriolis_v.png",dpi=300)

#-------------------------------------------------------------------------------

def plot_seaice_forces(testcase):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    fius = []
    fivs = []
    fsus = []
    fsvs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        nVertices = len(filein.dimensions["nVertices"])

        uSeaiceForceIceberg = filein.variables["uSeaiceForceIceberg"][0,0]
        vSeaiceForceIceberg = filein.variables["vSeaiceForceIceberg"][0,0]

        icebergForceU = filein.variables["icebergForceU"][0,:]
        icebergForceV = filein.variables["icebergForceV"][0,:]

        areaTriangle = filein.variables["areaTriangle"][:]

        filein.close()

        icebergForceUTotal = 0.0
        icebergForceVTotal = 0.0
        for iVertex in range(0,nVertices):
            icebergForceUTotal += icebergForceU[iVertex] * areaTriangle[iVertex]
            icebergForceVTotal += icebergForceV[iVertex] * areaTriangle[iVertex]

        fius.append(uSeaiceForceIceberg)
        fivs.append(vSeaiceForceIceberg)
        fsus.append(icebergForceUTotal)
        fsvs.append(icebergForceVTotal)

    fig, axis = plt.subplots()

    axis.plot(times, fius, label="u iceberg", marker="x")
    axis.plot(times, fivs, label="v iceberg", marker="+")
    axis.plot(times, fsus, label="u seaice", marker="x")
    axis.plot(times, fsvs, label="v seaice", marker="+")

    axis.set_title("Sea ice force")

    axis.legend()
    axis.set_xlabel("Time (days)")
    axis.set_ylabel("Force (N)")

    plt.tight_layout()

    plt.savefig("results_seaice_force_%s.png" %(testcase),dpi=300)

#-------------------------------------------------------------------------------

def plot_seaice_capture(testcase):

    # Transition from uncaptured to captured iceberg: capture when the sea ice force
    # needed to move the iceberg with the sea ice is below pressure * width.
    # The implied capture force is recalculated here assuming ocean drag is the
    # only other active force and hourly output at the model time step.

    from grounding_tests import GroundingPhysics, read_namelist

    physics = GroundingPhysics(read_namelist("namelist.seaice.%s" %(testcase)))

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))

    times = []
    vIcebergs = []
    vSeaices = []
    vOceans = []
    captured = []
    seaiceForces = []
    captureLimits = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        times.append(filein.variables["daysSinceStartOfSim"][0] * 24.0)
        vIcebergs.append(filein.variables["vVelocityIceberg"][0,0])
        vSeaices.append(filein.variables["vVelocitySeaiceIceberg"][0,0])
        vOceans.append(filein.variables["vOceanVelocityIceberg"][0,0])
        captured.append(filein.variables["icebergCapturedBySeaice"][0,0])
        seaiceForces.append(sqrt(filein.variables["uSeaiceForceIceberg"][0,0]**2 +
                                 filein.variables["vSeaiceForceIceberg"][0,0]**2))
        icebergWidth = filein.variables["icebergVolume"][0,0] / \
            (filein.variables["icebergLength"][0,0] * filein.variables["icebergHeight"][0,0])
        captureLimits.append(filein.variables["seaIcePressureIceberg"][0,0] * icebergWidth)

        filein.close()

    captured = np.array(captured)

    # sea ice force needed for the iceberg to move with the sea ice over each step
    vIcebergsInit = np.concatenate(([0.0], vIcebergs[:-1]))
    vSeaices = np.array(vSeaices)
    vOceans = np.array(vOceans)
    impliedForces = np.abs(physics.icebergMass * (vSeaices - vIcebergsInit) / physics.dt -
                           physics.oceanDragCoefficient * np.abs(vOceans - vSeaices) * (vOceans - vSeaices))

    fig, axes = plt.subplots(2, 1, figsize=(7,7), sharex=True)

    axes[0].plot(times, vIcebergs, color="red", marker="x", label="vVelocityIceberg")
    axes[0].plot(times, vSeaices, color="blue", linestyle="dashed", label="vVelocitySeaiceIceberg")
    axes[0].plot(times, vOceans, color="green", linestyle="dotted", label="vOceanVelocityIceberg")
    axes[0].set_ylabel("Velocity (m/s)")
    axes[0].legend(loc="lower right")

    axes[1].plot(times, impliedForces, color="blue", marker="+", label="implied capture force")
    axes[1].plot(times, seaiceForces, color="red", marker="x", label="|sea ice force| applied")
    axes[1].plot(times, captureLimits, color="black", linestyle="dashed", label="capture limit")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Force (N)")
    axes[1].set_xlabel("Time (hours)")
    axes[1].legend(loc="center right")

    for axis in axes:
        if (np.any(captured == 1)):
            axis.axvspan(times[np.argmax(captured == 1)], times[-1], color="grey", alpha=0.2)

    axes[0].set_title("Sea ice capture (shaded: captured)")

    plt.tight_layout()

    plt.savefig("results_%s_capture.png" %(testcase),dpi=300)

#-------------------------------------------------------------------------------

def plot_grounding(testcase, uColor, vColor):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    latIcebergs = []
    lonIcebergs = []

    uVelocityIcebergs = []
    vVelocityIcebergs = []

    icebergHeights = []
    icebergBathymetrys = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        latIceberg = filein.variables["latIceberg"][0,0]
        lonIceberg = filein.variables["lonIceberg"][0,0]

        uVelocityIceberg = filein.variables["uVelocityIceberg"][0,0]
        vVelocityIceberg = filein.variables["vVelocityIceberg"][0,0]

        icebergHeight = filein.variables["icebergHeight"][0,0]
        icebergBathymetry = filein.variables["icebergBathymetry"][0,0]
        icebergBathymetry *= -1


        filein.close()

        latIcebergs.append(latIceberg)
        lonIcebergs.append(lonIceberg)

        uVelocityIcebergs.append(uVelocityIceberg)
        vVelocityIcebergs.append(vVelocityIceberg)

        icebergHeights.append(icebergHeight)
        icebergBathymetrys.append(icebergBathymetry)

    fig, axis = plt.subplots()
    axis2 = axis.twinx()
    axis3 = axis.twinx()
    axis3.spines['right'].set_position(('outward', 60))  # offset in points

    axis.plot(times, lonIcebergs, color=uColor, linestyle="dashed", label="lonIceberg")
    axis.plot(times, latIcebergs, color=vColor, linestyle="dashed", label="latIceberg")

    axis2.plot(times, uVelocityIcebergs, color=uColor, linestyle="solid", label="uVelocityIceberg")
    axis2.plot(times, vVelocityIcebergs, color=vColor, linestyle="solid", label="vVelocityIceberg")

    axis3.plot(times, icebergHeights,     color="blue", linestyle="solid", label="icebergHeight")
    axis3.plot(times, icebergBathymetrys, color="teal", linestyle="solid", label="icebergBathymetry")

    axis.legend(loc="upper left")
    axis2.legend(loc="lower right")
    axis3.legend(loc="lower left")

    axis.set_title("Position/velocity %s" %(testcase))

    axis.set_xlabel("Time (days)")

    axis.set_ylabel("Position (rads)")
    axis2.set_ylabel("Velocity (m/s)")
    axis3.set_ylabel("Height (m)")

    plt.tight_layout()

    plt.savefig("results_%s.png" %(testcase),dpi=300)

#-------------------------------------------------------------------------------

def plot_grounding_forces(testcase):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    forceNames = ["Grounding", "GrndSediment", "GrndGravity", "GrndSolid"]

    forces = {}
    for forceName in forceNames:
        forces[forceName] = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        for forceName in forceNames:
            uForce = filein.variables["u%sForceIceberg" %(forceName)][0,0]
            vForce = filein.variables["v%sForceIceberg" %(forceName)][0,0]
            forces[forceName].append(np.sqrt(uForce**2 + vForce**2))

        filein.close()

    fig, axis = plt.subplots()

    for forceName in forceNames:
        axis.plot(times, forces[forceName], label=forceName)

    axis.set_title("Grounding force magnitudes")

    axis.legend()
    axis.set_xlabel("Time (days)")
    axis.set_ylabel("Force (N)")

    plt.tight_layout()

    plt.savefig("results_%s_forces.png" %(testcase),dpi=300)

#-------------------------------------------------------------------------------

def plot_testcase(testcase):

    if (testcase == "air_drag"):
        plot_figure("air_drag","red","green")
    elif (testcase == "ocean_drag"):
        plot_figure("ocean_drag","red","green")
    elif (testcase == "surface_tilt"):
        plot_figure("surface_tilt","red","green")
    elif (testcase == "wave_radiation"):
        plot_figure("wave_radiation","red","green")
    elif (testcase == "coriolis"):
        plot_coriolis()
    elif (testcase == "seaice_force_low_conc"):
        plot_figure("seaice_force_low_conc","red","green")
        plot_seaice_forces("seaice_force_low_conc")
    elif (testcase == "seaice_force_mid_conc"):
        plot_figure("seaice_force_mid_conc","red","green")
        plot_seaice_forces("seaice_force_mid_conc")
    elif (testcase == "seaice_force_high_conc"):
        plot_figure("seaice_force_high_conc","red","green")
        plot_seaice_forces("seaice_force_high_conc")
        plot_seaice_capture("seaice_force_high_conc")
    elif (testcase == "grounding"):
        plot_grounding("grounding","red","green")
    elif (testcase == "grounding_force"):
        plot_grounding("grounding_force","red","green")
        plot_grounding_forces("grounding_force")
    elif (testcase in groundingTestcases):
        plot_grounding_test(testcase)
    elif (testcase in seaiceStressTestcases):
        plot_seaice_stress_test(testcase)

#-------------------------------------------------------------------------------

def plot_results(testcase):

    if (testcase is not None):

        plot_testcase(testcase)

    else:

        testcases = ["air_drag",
                     "ocean_drag",
                     "surface_tilt",
                     "wave_radiation",
                     "coriolis",
                     "seaice_force_low_conc",
                     "seaice_force_mid_conc",
                     "seaice_force_high_conc",
                     "grounding",
                     "grounding_force"] + groundingTestcases + seaiceStressTestcases

        for testcase in testcases:
            plot_testcase(testcase)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='testcase', default=None)

    args = parser.parse_args()

    plot_results(args.testcase)
