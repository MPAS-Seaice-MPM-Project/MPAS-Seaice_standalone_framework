import sys
from math import fabs, sqrt, pow, pi, sin, cos, asin, acos, atan2, tan, atan
import numpy as np
import numpy.ma as ma
import argparse
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from numba import njit
from shapely.geometry import Point, Polygon

#-------------------------------------------------------------------------------

@njit
def mpas_arc_length(ax, ay, az, bx, by, bz):

    cx = bx - ax
    cy = by - ay
    cz = bz - az

    #      r = ax*ax + ay*ay + az*az
    #      c = cx*cx + cy*cy + cz*cz
    #
    #      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

    r = sqrt(ax*ax + ay*ay + az*az)
    c = sqrt(cx*cx + cy*cy + cz*cz)
    #      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
    return r * 2.0 * asin(c/(2.0*r))

#-------------------------------------------------------------------------------

@njit
def mpas_sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz):

    a = acos(max(min(bx*cx + by*cy + bz*cz,1.0),-1.0))      # Eqn. (3)
    b = acos(max(min(ax*cx + ay*cy + az*cz,1.0),-1.0))      # Eqn. (2)
    c = acos(max(min(ax*bx + ay*by + az*bz,1.0),-1.0))      # Eqn. (1)

    ABx = bx - ax
    ABy = by - ay
    ABz = bz - az

    ACx = cx - ax
    ACy = cy - ay
    ACz = cz - az

    Dx =   (ABy * ACz) - (ABz * ACy)
    Dy = -((ABx * ACz) - (ABz * ACx))
    Dz =   (ABx * ACy) - (ABy * ACx)

    s = 0.5*(a + b + c)
    #      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   # Eqn. (28)
    sin_angle = sqrt(min(1.0,max(0.0,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   # Eqn. (28)

    if ((Dx*ax + Dy*ay + Dz*az) >= 0.0):
        return  2.0 * asin(max(min(sin_angle,1.0),-1.0))
    else:
        return -2.0 * asin(max(min(sin_angle,1.0),-1.0))

#-------------------------------------------------------------------------------

@njit
def mpas_rotate_about_vector(x, y, z, theta, a, b, c, u, v, w):

    vw2 = v*v + w*w
    uw2 = u*u + w*w
    uv2 = u*u + v*v
    m = sqrt(u*u + v*v + w*w)

    xp = (a*vw2 + u*(-b*v-c*w+u*x+v*y+w*z) + ((x-a)*vw2+u*(b*v+c*w-v*y-w*z))*cos(theta) + m*(-c*v+b*w-w*y+v*z)*sin(theta))/(m*m)
    yp = (b*uw2 + v*(-a*u-c*w+u*x+v*y+w*z) + ((y-b)*uw2+v*(a*u+c*w-u*x-w*z))*cos(theta) + m*( c*u-a*w+w*x-u*z)*sin(theta))/(m*m)
    zp = (c*uv2 + w*(-a*u-b*v+u*x+v*y+w*z) + ((z-c)*uv2+w*(a*u+b*v-u*x-v*y))*cos(theta) + m*(-b*u+a*v-v*x+u*y)*sin(theta))/(m*m)

    return xp, yp, zp

#-------------------------------------------------------------------------------

@njit
def mpas_mirror_point(xPoint,
                      yPoint,
                      zPoint,
                      xA,
                      yA,
                      zA,
                      xB,
                      yB,
                      zB):

    #
    # Find the spherical angle between arcs AP and AB (where P is the point to be reflected)
    #
    alpha = mpas_sphere_angle(xA, yA, zA, xPoint, yPoint, zPoint, xB, yB, zB)

    #
    # Rotate the point to be reflected by twice alpha about the vector from the origin to A
    #
    xMirror, yMirror, zMirror = mpas_rotate_about_vector(xPoint,
                                                         yPoint,
                                                         zPoint,
                                                         2.0 * alpha,
                                                         0.0,
                                                         0.0,
                                                         0.0,
                                                         xA,
                                                         yA,
                                                         zA)

    return xMirror, yMirror, zMirror

#-------------------------------------------------------------------------------

@njit
def mpas_signed_triangle_area_sphere(ax, ay, az, bx, by, bz, cx, cy, cz, radius):

     ab = mpas_arc_length(ax, ay, az, bx, by, bz)/radius
     bc = mpas_arc_length(bx, by, bz, cx, cy, cz)/radius
     ca = mpas_arc_length(cx, cy, cz, ax, ay, az)/radius
     semiperim = 0.5 * (ab + bc + ca)

     tanqe = sqrt(max(0.0,tan(0.5 * semiperim) * tan(0.5 * (semiperim - ab)) \
                  * tan(0.5* (semiperim - bc)) * tan(0.5 * (semiperim - ca))))

     area = 4.0 * radius * radius * atan(tanqe)

     ablenx = bx - ax
     ableny = by - ay
     ablenz = bz - az

     aclenx = cx - ax
     acleny = cy - ay
     aclenz = cz - az

     dlenx =   (ableny * aclenz) - (ablenz * acleny)
     dleny = -((ablenx * aclenz) - (ablenz * aclenx))
     dlenz =   (ablenx * acleny) - (ableny * aclenx)

     if ((dlenx*ax + dleny*ay + dlenz*az) < 0.0):
          return -area
     else:
          return area

#-------------------------------------------------------------------------------

@njit
def mpas_signed_triangle_area_plane(ax, ay, az, bx, by, bz, cx, cy, cz):

     area = 0.5*(ax*(by-cy) - bx*(ay-cy) + cx*(ay-by))

     return area

#-------------------------------------------------------------------------------

@njit
def mpas_in_cell(xPoint,
                 yPoint,
                 zPoint,
                 xCell,
                 yCell,
                 zCell,
                 nEdgesOnCell,
                 verticesOnCell,
                 xVertex,
                 yVertex,
                 zVertex):

    radius = sqrt(xCell * xCell + yCell * yCell + zCell * zCell)
    radius_inv = 1.0 / radius

    inDist = mpas_arc_length(xPoint, yPoint, zPoint, xCell, yCell, zCell)

    mpas_in_cell = True

    for i in range(0,nEdgesOnCell):
        vtx1 = verticesOnCell[i]
        vtx2 = verticesOnCell[(i+1) % nEdgesOnCell]

        xNeighbor, yNeighbor, zNeighbor = mpas_mirror_point(xCell*radius_inv,
                                                            yCell*radius_inv,
                                                            zCell*radius_inv,
                                                            xVertex[vtx1]*radius_inv,
                                                            yVertex[vtx1]*radius_inv,
                                                            zVertex[vtx1]*radius_inv,
                                                            xVertex[vtx2]*radius_inv,
                                                            yVertex[vtx2]*radius_inv,
                                                            zVertex[vtx2]*radius_inv)

        xNeighbor = xNeighbor * radius
        yNeighbor = yNeighbor * radius
        zNeighbor = zNeighbor * radius

        outDist = mpas_arc_length(xPoint, yPoint, zPoint, xNeighbor, yNeighbor, zNeighbor)

        if (outDist < inDist):
            mpas_in_cell = False
            return mpas_in_cell

    return mpas_in_cell

#-------------------------------------------------------------------------------

@njit
def in_geom(x,
            y,
            z,
            geomType):

    in_geom = False

    if (geomType == "slotted_cylinder"):

        iceArea = 0.0
        iceVolume = 0.0

        circleRadius = 0.5
        r = sqrt(pow(z,2) + pow(x,2))

        if (r < circleRadius and y > 0.0):

            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

        if (fabs(x) < 1.0/12.0 and z > -2.0/6.0):

            in_geom = False

    elif (geomType == 'cosine_bell'):

        iceArea = 0.0
        iceVolume = 0.0

        circleRadius = 1.0/3.0
        r = sqrt(pow(z,2) + pow(x,2))

        if (r < circleRadius and y > 0.0):

            in_geom = True
            iceArea = 0.5 * (1.0 + cos((pi * r) / circleRadius))
            iceVolume = 1.0

    elif (geomType == 'uniform'):

        in_geom = True
        iceArea = 1.0
        iceVolume = 1.0

    elif (geomType == 'cap'):

        lat = asin(z)
        if (lat > 1.22):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    elif (geomType == 'ring'):

        lat = asin(z)
        if (lat > 1.0 and lat < 1.25):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    elif (geomType == 'bar'):
        if (x < 160000 and x > 0 and y < 160000 and y > 80000):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    elif (geomType == 'square'):
        if (x < 100000 and x > 0 and y < 120000 and y > 20000):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    elif (geomType == 'BBMsquare'):
        if (x < 1000000 and x > 0 and y < 1020000 and y > 20000):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    elif (geomType == 'disks'):
        if ((x-0.25)*(x-0.25)+(y-0.25)*(y-0.25) < 0.04 or (x-0.75)*(x-0.75)+(y-0.75)*(y-0.75) < 0.04):
            in_geom = True
            iceArea = 1.0
            iceVolume = 1.0

    return in_geom, iceArea, iceVolume

#-------------------------------------------------------------------------------

def place_particles_in_cell(particleInitType,
                            particleInitNumber,
                            areaCell,
                            averageCellSize,
                            particlePositionInitType,
                            on_a_sphere,
                            earthRadius,
                            nEdgesOnCell,
                            verticesOnCell,
                            latVertex,
                            lonVertex,
                            xVertex,
                            yVertex,
                            zVertex,
                            xCell,
                            yCell,
                            zCell):

    # TODO add some checks for minimum sizes, possibly reduce
    # number of added material points to consolidate
    if (particleInitType.strip() == 'number'):
        nParticlesPerCellDesired = particleInitNumber
    elif (particleInitType.strip() == 'area'):
        nParticlesPerCellDesired = particleInitNumber * int(areaCell / averageCellSize)
        nParticlesPerCellDesired = max(1, nParticlesPerCellDesired)
    else:
        # error not a valid option
        raise Exception("Invalid particle_init_type")

    xTrial = []
    yTrial = []
    zTrial = []
    areaTrial = []

    if (particlePositionInitType.strip() == 'even'):

        # place points evenly in circumscribing rectangle
        # WARNING: for 'even' to work, numberToPlace must be a perfect square
        npart = int(sqrt(float(nParticlesPerCellDesired)))
        #if ( (npart*npart) \= nParticlesCell(iCell)) then
        #   print*,'error in placing material points'
        #endif

        if (on_a_sphere):

            # determine rectangle
            latmin =  sys.float_info.max
            latmax = -sys.float_info.max
            lonmin =  sys.float_info.max
            lonmax = -sys.float_info.max
            for iVertexOnCell in range(0, nEdgesOnCell):
                iVertex = verticesOnCell[iVertexOnCell]
                latmin = min(latmin,latVertex[iVertex])
                latmax = max(latmax,latVertex[iVertex])
                lonmin = min(lonmin,lonVertex[iVertex])
                lonmax = max(lonmax,lonVertex[iVertex])

            dx = (latmax - latmin) / float(npart) / 2.0
            dy = (lonmax - lonmin) / float(npart) / 2.0
            # special case for lon near -pi and pi
            if (lonmin * lonmax < 0.0 and lonmax > pi / 2.0):
                dy = (lonmax - lonmin - 2.0 * pi) / float(npart) / 2.0

            k = 0
            for i in range(0, npart):
                for j in range(0, npart):

                    lat = latmin + dx * float(2 * i + 1)
                    lon = lonmin + dy * float(2 * j + 1)
                    if (lon < -pi): lon = 2.0 * pi + lon

                    x = earthRadius * cos(lon) * cos(lat)
                    y = earthRadius * sin(lon) * cos(lat)
                    z = earthRadius * sin(lat)

                    area = earthRadius*earthRadius*((sin(lat) - sin(lat-dx))*(dy))

                    xTrial.append(x)
                    yTrial.append(y)
                    zTrial.append(z)
                    areaTrial.append(area)

        else: # not on a sphere

            # determine rectangle
            xmin =  sys.float_info.max
            xmax = -sys.float_info.max
            ymin =  sys.float_info.max
            ymax = -sys.float_info.max
            for iVertexOnCell in range(0, nEdgesOnCell):
                iVertex = verticesOnCell[iVertexOnCell]
                xmin = min(xmin, xVertex[iVertex])
                xmax = max(xmax, xVertex[iVertex])
                ymin = min(ymin, yVertex[iVertex])
                ymax = max(ymax, yVertex[iVertex])

            # placement limits decrease in size by 5% so particles not on boundary
            # of regular planar polygons
            xminUse = xmin - 0.05 * (xmax-xmin)
            xmaxUse = xmax + 0.05 * (xmax-xmin)
            yminUse = ymin - 0.05 * (ymax-ymin)
            ymaxUse = ymax + 0.05 * (ymax-ymin)

            dx = (xmaxUse-xminUse) / float(npart) / 2.0
            dy = (ymaxUse-yminUse) / float(npart) / 2.0

            k = 0
            for i in range(0, npart):
                for j in range(0, npart):

                    x = xminUse + dx * float(2 * i + 1)
                    y = yminUse + dy * float(2 * j + 1)
                    z = 0.0
                    area = dx*dy

                    xTrial.append(x)
                    yTrial.append(y)
                    zTrial.append(z)
                    areaTrial.append(area)

    elif (particlePositionInitType.strip() == 'onePerEdge'):
        # place one particle per edge in the polygon

        v1 = [0.,0.,0.]
        v2 = [0.,0.,0.]
        v3 = [0.,0.,0.]

        xvert = xVertex[verticesOnCell[:]]
        yvert = yVertex[verticesOnCell[:]]
        zvert = zVertex[verticesOnCell[:]]
        k = 0

        for iVertexOnCell in range (0, nEdgesOnCell):

            v1[0] = xvert[iVertexOnCell]
            if (iVertexOnCell < nEdgesOnCell - 1):
                v2[0] = xvert[iVertexOnCell + 1]
            else:
                v2[0] = xvert[0]
            v3[0] = xCell
            x = (v1[0] + v2[0] + v3[0])/3.0

            v1[1] = yvert[iVertexOnCell]
            if (iVertexOnCell < nEdgesOnCell - 1):
                v2[1] = yvert[iVertexOnCell + 1]
            else:
                v2[1] = yvert[0]
            v3[1] = yCell
            y = (v1[1] + v2[1] + v3[1])/3.0

            v1[2] = zvert[iVertexOnCell]
            if (iVertexOnCell < nEdgesOnCell - 1):
                v2[2] = zvert[iVertexOnCell + 1]
            else:
                v2[2] = zvert[0]
            v3[2] = zCell
            z = (v1[2] + v2[2] + v3[2])/3.0

            if (on_a_sphere):

                R = sqrt(x*x + y*y + z*z)
                x = x*earthRadius/R
                y = y*earthRadius/R
                z = z*earthRadius/R
                v1[0]=v1[0]*earthRadius
                v1[1]=v1[1]*earthRadius
                v1[2]=v1[2]*earthRadius
                v2[0]=v2[0]*earthRadius
                v2[1]=v2[1]*earthRadius
                v2[2]=v2[2]*earthRadius
                v3[0]=v3[0]*earthRadius
                v3[1]=v3[1]*earthRadius
                v3[2]=v3[2]*earthRadius
                area = mpas_signed_triangle_area_sphere(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2], earthRadius)

            else:
                area = mpas_signed_triangle_area_plane(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2])

            xTrial.append(x)
            yTrial.append(y)
            zTrial.append(z)
            areaTrial.append(area)

    elif (particlePositionInitType.strip() == 'poisson'):
        # place points using a poisson distribution
        raise Exception("Poisson position init type not implemented")
    elif (particlePositionInitType.strip() == 'random'):
        # place points using a random distribution
        raise Exception("Random position init type not implemented")
    else:
        raise Exception("Invalid particle initType")

    # check if trial particles are in cell
    xInCell = []
    yInCell = []
    zInCell = []
    areaInCell = []

    for x, y, z, area in zip(xTrial, yTrial, zTrial, areaTrial):

        if (on_a_sphere):
            inCell = mpas_in_cell(x,
                                  y,
                                  z,
                                  xCell,
                                  yCell,
                                  zCell,
                                  nEdgesOnCell,
                                  verticesOnCell[:],
                                  xVertex,
                                  yVertex,
                                  zVertex)
        else:
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell):
                iVertex = verticesOnCell[iVertexOnCell]
                vertices.append((xVertex[iVertex],yVertex[iVertex]))
            polygon = Polygon(vertices)
            inCell = polygon.contains(Point(x,y))

        if (inCell):
            xInCell.append(x)
            yInCell.append(y)
            zInCell.append(z)
            areaInCell.append(area)

    return \
        xInCell, \
        yInCell, \
        zInCell, \
        areaInCell

#-------------------------------------------------------------------------------

def check_cell_particles_in_geom(xInCell,
                                 yInCell,
                                 zInCell,
                                 areaInCell,
                                 on_a_sphere,
                                 earthRadius,
                                 geomType,
                                 iCell,
                                 posnMP,
                                 latCellMP,
                                 lonCellMP,
                                 areaMP,
                                 cellIDCreationMP,
                                 creationIndexMP,
                                 iceAreaCellMP,
                                 iceVolumeCellMP,
                                 nParticlesPerCellActual):

    # check if cell particles are in geometry
    k = 0
    for x, y, z, area in zip(xInCell, yInCell, zInCell, areaInCell):

        if (on_a_sphere):

            inGeom, iceArea, iceVolume  = in_geom(x/earthRadius,
                                                  y/earthRadius,
                                                  z/earthRadius,
                                                  geomType)

            lon = atan2(y, x)
            lat = asin(z/earthRadius)

        else:

            inGeom, iceArea, iceVolume  = in_geom(x,
                                                  y,
                                                  z,
                                                  geomType)

            lon = 0.0
            lat = 0.0

        if (inGeom):
            k = k + 1
            posnMP.append([x, y, z])
            latCellMP.append(lat)
            lonCellMP.append(lon)
            areaMP.append(area)
            cellIDCreationMP.append(iCell+1)
            creationIndexMP.append(k)
            iceAreaCellMP.append(iceArea)
            iceVolumeCellMP.append(iceVolume)

    nParticlesPerCellActual = k

    return \
        posnMP, \
        latCellMP, \
        lonCellMP, \
        areaMP, \
        cellIDCreationMP, \
        creationIndexMP, \
        iceAreaCellMP, \
        iceVolumeCellMP, \
        nParticlesPerCellActual

#-------------------------------------------------------------------------------

def create_particles_file(filenameOut,
                          nParticles,
                          nCells,
                          nCategories,
                          posnMP,
                          latCellMP,
                          lonCellMP,
                          areaMP,
                          cellIDCreationMP,
                          creationIndexMP,
                          nParticlesCell,
                          iceAreaCellMP,
                          iceVolumeCellMP,
                          iceAreaCategoryMP,
                          iceVolumeCategoryMP):

    # output
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("nCategories", nCategories)
    fileOut.createDimension("THREE", 3)
    fileOut.createDimension("TWO", 2)
    fileOut.createDimension("ONE", 1)

    var = fileOut.createVariable("posnMP","d",dimensions=["nParticles","THREE"])
    var[:] = posnMP[:]

    var = fileOut.createVariable("latCellMP","d",dimensions=["nParticles"])
    var[:] = latCellMP[:]

    var = fileOut.createVariable("lonCellMP","d",dimensions=["nParticles"])
    var[:] = lonCellMP[:]

    var = fileOut.createVariable("areaMP","d",dimensions=["nParticles"])
    var[:] = areaMP[:]

    var = fileOut.createVariable("cellIDCreationMP","i",dimensions=["nParticles"])
    var[:] = cellIDCreationMP[:]

    var = fileOut.createVariable("creationIndexMP","i",dimensions=["nParticles"])
    var[:] = creationIndexMP[:]

    var = fileOut.createVariable("nParticlesCell","i",dimensions=["nCells"])
    var[:] = nParticlesCell[:]

    var = fileOut.createVariable("iceAreaCellMP","d",dimensions=["nParticles"])
    var[:] = iceAreaCellMP[:]

    var = fileOut.createVariable("iceVolumeCellMP","d",dimensions=["nParticles"])
    var[:] = iceVolumeCellMP[:]

    var = fileOut.createVariable("iceAreaCategoryMP","d",dimensions=["nParticles","nCategories","ONE"])
    var[:] = iceAreaCategoryMP[:]

    var = fileOut.createVariable("iceVolumeCategoryMP","d",dimensions=["nParticles","nCategories","ONE"])
    var[:] = iceVolumeCategoryMP[:]

    fileOut.close()

#-------------------------------------------------------------------------------

def plot_particles(on_a_sphere,
                   posnMP,
                   x,
                   y,
                   z,
                   earthRadius,
                   filenameOut):

    if (on_a_sphere):

        fig = plt.figure(figsize=(10,10))
        axis = fig.add_subplot(projection='3d')
        axis.scatter(posnMP[:,0],posnMP[:,1],posnMP[:,2])
        axis.set_xlim(np.min([-earthRadius,np.amin(x)])*1.1,np.max([np.amax(x),earthRadius])*1.1)
        axis.set_ylim(np.min([-earthRadius,np.amin(y)])*1.1,np.max([np.amax(y),earthRadius])*1.1)
        axis.set_zlim(np.min([-earthRadius,np.amin(z)])*1.1,np.max([np.amax(z),earthRadius])*1.1)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_zlabel("z")

    else:

        fig, axis = plt.subplots(figsize=(10,10))
        axis.scatter(posnMP[:,0],posnMP[:,1])
        axis.set_xlim(np.amin(x),np.amax(x))
        axis.set_ylim(np.amin(y),np.amax(y))
        axis.set_xlabel("x")
        axis.set_ylabel("y")

    plt.tight_layout()
    plt.savefig(filenameOut)
    plt.close()

#-------------------------------------------------------------------------------

def initial_particle_positions(filenameMesh,
                               filenameOut,
                               particleInitType,
                               particleInitNumber,
                               particlePositionInitType,
                               geomType,
                               earthRadius):

    # mesh info
    fileMesh = Dataset(filenameMesh,"r")

    on_a_sphere = fileMesh.on_a_sphere
    sphere = fileMesh.on_a_sphere
    if (sphere == "NO"):
        on_a_sphere = False
    elif (sphere == "YES"):
        on_a_sphere = True

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell   = ma.getdata(fileMesh.variables["nEdgesOnCell"][:])
    verticesOnCell = ma.getdata(fileMesh.variables["verticesOnCell"][:])

    areaCell  = ma.getdata(fileMesh.variables["areaCell"][:])
    latVertex = ma.getdata(fileMesh.variables["latVertex"][:])
    lonVertex = ma.getdata(fileMesh.variables["lonVertex"][:])
    xVertex   = ma.getdata(fileMesh.variables["xVertex"][:])
    yVertex   = ma.getdata(fileMesh.variables["yVertex"][:])
    zVertex   = ma.getdata(fileMesh.variables["zVertex"][:])
    xCell     = ma.getdata(fileMesh.variables["xCell"][:])
    yCell     = ma.getdata(fileMesh.variables["yCell"][:])
    zCell     = ma.getdata(fileMesh.variables["zCell"][:])

    fileMesh.close()

    verticesOnCell[:] -= 1

    # average cell size
    averageCellSize = np.mean(areaCell)

    # particle positions
    nParticlesCell = np.zeros(nCells,dtype="i")

    posnMP = []
    latCellMP = []
    lonCellMP = []
    areaMP = []
    cellIDCreationMP = []
    creationIndexMP = []
    iceAreaCellMP = []
    iceVolumeCellMP = []

    for iCell in range(0, nCells):

        xInCell, \
            yInCell, \
            zInCell, \
            areaInCell = place_particles_in_cell(particleInitType,
                                              particleInitNumber,
                                              areaCell[iCell],
                                              averageCellSize,
                                              particlePositionInitType,
                                              on_a_sphere,
                                              earthRadius,
                                              nEdgesOnCell[iCell],
                                              verticesOnCell[iCell,:],
                                              latVertex,
                                              lonVertex,
                                              xVertex,
                                              yVertex,
                                              zVertex,
                                              xCell[iCell],
                                              yCell[iCell],
                                              zCell[iCell])

        posnMP, \
            latCellMP, \
            lonCellMP, \
            areaMP, \
            cellIDCreationMP, \
            creationIndexMP, \
            iceAreaCellMP, \
            iceVolumeCellMP, \
            nParticlesCell[iCell] = check_cell_particles_in_geom(xInCell,
                                                                 yInCell,
                                                                 zInCell,
                                                                 areaInCell,
                                                                 on_a_sphere,
                                                                 earthRadius,
                                                                 geomType,
                                                                 iCell,
                                                                 posnMP,
                                                                 latCellMP,
                                                                 lonCellMP,
                                                                 areaMP,
                                                                 cellIDCreationMP,
                                                                 creationIndexMP,
                                                                 iceAreaCellMP,
                                                                 iceVolumeCellMP,
                                                                 nParticlesCell[iCell])

    nParticles = len(posnMP)
    posnMP = np.array(posnMP)
    latCellMP = np.array(latCellMP)
    lonCellMP = np.array(lonCellMP)
    areaMP = np.array(areaMP)
    cellIDCreationMP = np.array(cellIDCreationMP)
    creationIndexMP = np.array(creationIndexMP)
    nParticlesCell = np.array(nParticlesCell)
    iceAreaCellMP = np.array(iceAreaCellMP)
    iceVolumeCellMP = np.array(iceVolumeCellMP)
    iceAreaCategoryMP = iceAreaCellMP
    iceVolumeCategoryMP = iceVolumeCellMP

    # output
    create_particles_file(filenameOut,
                          nParticles,
                          nCells,
                          1,
                          posnMP,
                          latCellMP,
                          lonCellMP,
                          areaMP, 
                          cellIDCreationMP,
                          creationIndexMP,
                          nParticlesCell,
                          iceAreaCellMP,
                          iceVolumeCellMP,
                          iceAreaCategoryMP,
                          iceVolumeCategoryMP)

    # plot
    plot_particles(on_a_sphere,
                   posnMP,
                   xCell,
                   yCell,
                   zCell,
                   earthRadius,
                   "particles.png")

#-------------------------------------------------------------------------------

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="filenameMesh", required=True)
    parser.add_argument('-o', dest="filenameOut", required=True)
    parser.add_argument('-t', dest="particleInitType", required=True, choices=["number","area"])
    parser.add_argument('-n', dest="particleInitNumber", required=True, type=int)
    parser.add_argument('-p', dest="particlePositionInitType", required=True, choices=["even","onePerEdge","poisson","random"])
    parser.add_argument('-i', dest="geomType", required=True, choices=["cosine_bell","slotted_cylinder","uniform"])
    parser.add_argument('-r', dest="earthRadius", type=float, default=6371229.0)
    parser.add_argument('--concIndex', dest="iceConcTimeIndex", default=-1)

    args = parser.parse_args()

    initial_particle_positions(args.filenameMesh,
                               args.filenameOut,
                               args.particleInitType,
                               args.particleInitNumber,
                               args.particlePositionInitType,
                               args.geomType,
                               args.earthRadius)
