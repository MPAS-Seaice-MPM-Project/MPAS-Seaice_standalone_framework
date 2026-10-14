"""
Extract a field of particles from a time series of NetCDF files as VTK files
for plotting in ParaView.

The times in pvd file are NOT the time extracted from netCDF output files,
they are simply the order of the output files.

To run the script:

$ python paraview_vtk_particle_extractor.py -f "[path]/particle*" -v "iceAreaMP, iceVolMP"

where:
[path] is path of the directory holding netCDF output particle data files,
whose pattern are particles*;
"iceAreaMP, iceVolMP, ..." is the list of extracted variables separated
by comma (can contain space).

Notes:
(1) If any extracted variable is not in the netCDF files, the script prints
error message and exit.
(2) The .vtp files are in ./vtk_files/time_series
(3) The .pvd file (contains the filenames of .vtp files) is in ./vtk_files
(4) If ./vtk_files does not exist, then the script will create it. If it already
exists, then the script adds files to it.
"""

from netCDF4 import Dataset
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import chartostring
import glob
import math
import argparse
import os
from datetime import datetime

#-------------------------------------------------------------------------------

def fractional_year_noleap(time_str):

    dt = datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")

    day_of_year = dt.timetuple().tm_yday

    return (
        dt.year
        + (
            (day_of_year - 1)
            + dt.hour / 24
            + dt.minute / 1440
            + dt.second / 86400
        ) / 365.0
    )

#-------------------------------------------------------------------------------

def print_vtp_file(filenameIn,
                   filenameOut,
                   variable_names,
                   nParticlesDimName,
                   posnMPVarName):

    iTime = 0
    filein = Dataset(filenameIn,"r")
    sphere_radius = filein.sphere_radius
    nParticles = len(filein.dimensions[nParticlesDimName])
    posnMP = filein.variables[posnMPVarName][:]

    try:
        xtime = str(chartostring(filein.variables["xtime"][iTime,:]))
        xtimeYear = fractional_year_noleap(xtime)
    except:
        xtime = None
        xtimeYear = None

    # number of valid MPs
    nValid = 0

    # mask of valid (true) and invalid (false) MPs
    valid = []

    # loop over all MPs, check if it is a "valid" MP (i.e. it is "in use") then count it
    for iParticle in range(0,nParticles):
        # distance from globe center to MP
        mag = math.sqrt(math.pow(posnMP[iTime,iParticle,0],2) + \
                        math.pow(posnMP[iTime,iParticle,1],2) + \
                        math.pow(posnMP[iTime,iParticle,2],2))
        # if distance > 0, then it is a valid MP
        if (mag > 1e-8):
            nValid += 1
            valid.append(True)
        else:
            valid.append(False)

    fileout = open(filenameOut, "w")
    fileout.write("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n")
    if (xtime is not None):
        fileout.write("   <metadata>\n")
        fileout.write("      <xtime>\n")
        fileout.write("         %s\n" %(xtime))
        fileout.write("      </xtime>\n")
        fileout.write("   </metadata>\n")
    fileout.write("   <PolyData>\n")
    fileout.write("      <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n" %(nValid))

    # write MP position
    fileout.write("         <Points>\n")
    fileout.write("            <DataArray type=\"Float64\" Name=\"MaterialPoints\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    for iParticle in range(0,nParticles):
        if (valid[iParticle]):
            fileout.write("               %.13E %.13E %.13E\n" %(posnMP[iTime, iParticle, 0], posnMP[iTime, iParticle, 1], posnMP[iTime, iParticle, 2]))
    fileout.write("            </DataArray>\n")
    fileout.write("         </Points>\n")

    # write MP data
    fileout.write("         <PointData Scalars=\"MPsData\">\n")
    for var in variable_names:
        variableData = filein.variables[var][:]
        dataType = filein[var].dtype
        if (dataType == 'float64'):
            dataTypeVtk = 'Float64'
        elif (dataType == 'int32'):
            dataTypeVtk = 'Int32'
        else:
            print("Variable %s: dtype= %s is not supported" % var, filein[var].dtype)
            exit(1)
        nDims = filein[var].ndim
        if (nDims > 2):
            nComponents = filein[var].shape[nDims - 1]
        else:
            nComponents = 1

        fileout.write("            <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n" % (dataTypeVtk, var, nComponents))
        for iParticle in range(0,nParticles):
            if valid[iParticle]:
                if (nComponents == 1):
                    if (dataTypeVtk == 'Float64'):
                        fileout.write("               %.13E\n" % variableData[iTime,iParticle])
                    elif (dataTypeVtk == 'Int32'):
                        fileout.write("               %d\n" % variableData[iTime,iParticle])
                elif (nComponents == 2):
                    if (dataTypeVtk == 'Float64'):
                        fileout.write("               %.13E %.13E\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1]))
                    elif (dataTypeVtk == 'Int32'):
                        fileout.write("               %d %d\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1]))
                elif (nComponents == 3):
                    if (dataTypeVtk == 'Float64'):
                        fileout.write("               %.13E %.13E %.13E\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2]))
                    elif (dataTypeVtk == 'Int32'):
                        fileout.write("               %d %d %d\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2]))
                elif (nComponents == 4):
                    if (dataTypeVtk == 'Float64'):
                        fileout.write("               %.13E %.13E %.13E %.13E\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2], variableData[iTime,iParticle,3]))
                    elif (dataTypeVtk == 'Int32'):
                        fileout.write("               %d %d %d %d\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2], variableData[iTime,iParticle,3]))
                elif (nComponents == 6):
                    if (dataTypeVtk == 'Float64'):
                        fileout.write("               %.13E %.13E %.13E %.13E %.13E %.13E\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2], variableData[iTime,iParticle,3], variableData[iTime,iParticle,4], variableData[iTime,iParticle,5]))
                    elif (dataTypeVtk == 'Int32'):
                        fileout.write("               %d %d %d %d %d %d\n" % (variableData[iTime,iParticle,0], variableData[iTime,iParticle,1], variableData[iTime,iParticle,2], variableData[iTime,iParticle,3], variableData[iTime,iParticle,4], variableData[iTime,iParticle,5]))
                else:
                    print("nComponents = %d, is not supported yet" % nComponents)
        fileout.write("            </DataArray>\n")

    fileout.write("         </PointData>\n")
    fileout.write("      </Piece>\n")
    fileout.write("   </PolyData>\n")
    fileout.write("</VTKFile>\n")
    fileout.close()
    filein.close()

    return nValid, xtimeYear

#-------------------------------------------------------------------------------

def makeFolder(folderName):
    try:
        os.makedirs(folderName)
    except OSError:
        pass

#-------------------------------------------------------------------------------

def write_pvd_header(path, prefix):  # {{{
    pvd_file = open('{}/{}.pvd'.format(path, prefix), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    return pvd_file  # }}}

#-------------------------------------------------------------------------------

def _expand_variable_list(variable_list):
    if isinstance(variable_list, str):
        variable_names = variable_list.replace(' ','').split(',')
    else:
        variable_names = variable_list

    variable_names.sort()
    return variable_names

#-------------------------------------------------------------------------------

def open_netcdf(file_name):
    nc_file = NetCDFFile(file_name, 'r')
    # turn off auto mask (if applicable)
    try:
        nc_file.set_auto_mask(False)
    except AttributeError:
        pass
    return nc_file

#-------------------------------------------------------------------------------

def _add_var(variables, var_name, inc_dims, variable_names, exc_dims=None):
    if var_name in variable_names:
        return

    dims = variables[var_name].dimensions
    supported = False
    for d in inc_dims:
        if d in dims:
            supported = True
    if exc_dims is not None:
        for d in exc_dims:
            if d in dims:
                supported = False
    if supported:
        variable_names.append(var_name)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Print vtk files for particle visualization')
    parser.add_argument('-f', "--file_pattern", required=True, dest='filenameTemplate',
                        help='filename template for particle files')
    parser.add_argument("-v", "--variable_list", required=True, dest="variable_list",
                        help="List of variables to extract"),
    parser.add_argument("-d", "--dim_name", default="nParticles", dest="nParticlesDimName",
                        help="Name of the particle dimension")
    parser.add_argument("-p", "--pos_name", default="posnMP", dest="posnMPVarName",
                        help="Name of the particle position variable name")
    args = parser.parse_args()

    # all time-series files received from command line after -f
    all_time_series_filenames = sorted(glob.glob(args.filenameTemplate))

    # filter to exclude files containing '.lock'
    time_series_filenames = []
    for f in all_time_series_filenames:
        if '.lock' not in f:
            time_series_filenames.append(f)

    # take the first file from the time-series files to get the array of variables listed in the time-series files
    time_series_file = open_netcdf(time_series_filenames[0])
    time_series_variables = time_series_file.variables

    # print variables existing in the output results
    print("==================================================")
    print("Variables existing in the time-series files:")
    count = 0
    for var1 in time_series_variables:
        count = count + 1
        print("%d." % count, var1)
    print("==================================================")

    # list of variables received from command line after -v
    variable_list = args.variable_list

    # split variable names in variable_list and put them in sorted array
    variable_names = _expand_variable_list(variable_list)

    # check if variables received from command line after -v are existing in the time-series files
    for var in variable_names:
        if var not in time_series_variables :
            print("\nError: variable is not in the time-series files: %s\n" % var)
            exit(1)

    nskip = 1

    # testing, do not delete as we may use it for later development
    # for var in variable_names:
    #     print("---------------")
    #     print("var name = %s" % var)
    #     print("dimension (tuple): ")
    #     print(time_series_file[var].dimensions)
    #     print("shape (tuple):")
    #     print(time_series_file[var].shape)
    #     nDims = time_series_file[var].ndim
    #     print("nDims = %d" % nDims)
    #     print("nComps = %d" % time_series_file[var].shape[nDims-1])
    #     print("type= %s" % time_series_file[var].dtype)
    #     dataType = time_series_file.variables[var].dtype
    #     print(type(dataType))
    #     if (dataType == 'float64'):
    #         print("Float64")
    #     elif (dataType == 'int32'):
    #         print("Int32")
    # exit(0)

    # create 2 directories
    makeFolder("./vtk_files")
    makeFolder("./vtk_files/time_series")

    # start the pvd file in ./vtk_files
    out_prefix = "fieldsOnParticles"
    out_dir = "./vtk_files"
    pvd_file = write_pvd_header(out_dir, out_prefix)
    pvd_file.write('<Collection>\n')

    iFile = 0
    for filenameIn in time_series_filenames[::nskip]:
        print(filenameIn)
        filenameOut = "./vtk_files/time_series/fieldsOnParticles_%d.vtp" %(iFile)

        # write vtp files in vtk_files/time_series/
        nValid, xtimeYear = print_vtp_file(filenameIn,
                                           filenameOut,
                                           variable_names,
                                           args.nParticlesDimName,
                                           args.posnMPVarName)
        print(filenameOut, nValid, xtimeYear)

        # write the header for the vtp files
        time_index = iFile
        vtp_file_prefix = "time_series/{}_{:d}".format(out_prefix,
                                                       time_index)
        file_name = '{}/{}.vtp'.format(out_dir, vtp_file_prefix)
        if (xtimeYear is not None):
            real_time = xtimeYear
        else:
            real_time = iFile
        pvd_file.write('<DataSet timestep="{:.16f}" group="" '
                        'part="0"\n'.format(real_time))
        pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))

        iFile += 1

    # finish the pvd file
    pvd_file.write('</Collection>\n')
    pvd_file.write('</VTKFile>\n')
    pvd_file.close()
