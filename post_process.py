#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import random
import sys
import os
from pathos.multiprocessing import ThreadingPool
from io import StringIO
import argparse
import datetime


def read_arguments():
    parser = argparse.ArgumentParser(description="Input data", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-P",
        "--plot",
        default="False",
        help="True: Produce plots of the solutions. False: Do not produce plots",
    )
    parser.add_argument(
        "-PE",
        "--plot_ex_prob",
        default="False",
        help="True: Produce plots of the specified exceedance probabilities. False: Do not produce plots",
    )
    parser.add_argument(
        "-EX",
        "--ex_prob",
        default='',
        help="List of exceedence probabilities to be used for graphical output",
    )
    parser.add_argument(
        "-T",
        "--time_steps",
        default='',
        help="List of time steps to plot (integer >= 0). Type all to plot all the time steps",
    )
    parser.add_argument(
        "-L",
        "--levels",
        default='',
        help="List of vertical levels (integer >= 1) to plot. Type all to plot all the levels",
    )
    parser.add_argument(
        "-D",
        "--days_plot",
        default='',
        help="List of days to plot (YYYYMMDD). Type all to plot all the days",
    )
    parser.add_argument(
        "-C",
        "--convert",
        default="False",
        help="If True, convert output concentration into other species listed with the command -S (--species)",
    )
    parser.add_argument(
        "-S", "--species", default='', help="List of gas species (e.g. CO2)"
    )
    parser.add_argument(
        "-TS", "--tracking_specie", default=None, help="The original emitted specie that is tracked in the simulation"
    )
    parser.add_argument(
        "-PER",
        "--persistence",
        default="False",
        help="If True, calculate the persistence of the gas specie, i.e. the probability to be exposed to a gas species"
             " above specified concentration thresholds for times longer than the specified exposure times for those"
             " thresholds.\n" + "Concentration thresholds and exposure times should be provided in gas_properties.csv",
    )
    parser.add_argument(
        "-N",
        "--nproc",
        default=1,
        help="Maximum number of allowed simultaneous processes",
    )
    parser.add_argument(
        "-M",
        "--models",
        default="all",
        help="Model outputs to post-process. Options: disgas, twodee, all",
    )
    parser.add_argument(
        "-U",
        "--units",
        default=None,
        help="Gas concentration units. Possible options are: ppm, kg/m3",
    )
    parser.add_argument(
        "-PL",
        "--plot_limits",
        default='',
        help="Minimum and maximum value of concentration to display. If unspecified, they are obtained from all "
             "the outputs",
    )
    parser.add_argument("-PI",
        "--plot_isolines",
        default='',
        help="List of gas concentrations values to be used to draw isolines. Optional"
    )
    parser.add_argument(
        "-TA",
        "--time_av",
        default=None,
        help="Generate time-averaged outputs. Specify the time-averaging interval (in hours), or 0 for averaging "
             "over the whole duration",
    )
    parser.add_argument(
        "-OF",
        "--output_format",
        default="GRD",
        help="Select format of the processed output files. Valid options are: GRD",
    )
    parser.add_argument(
        "-PT",
        "--plot_topography",
        default="False",
        help="Plot topography layer (True or False). Warning, it can be time-consuming!",
    )
    parser.add_argument(
        "-TI",
        "--topography_isolines",
        default=100,
        help="Topography height contour lines spatial resolution (in m a.s.l.). Used only if -PT True",
    )
    parser.add_argument(
        "-PR", "--plot_resolution", default=600, help="Specify plot resolution in dpi"
    )
    parser.add_argument(
        "-TP", "--tracking_points", default="False", help="Extrapolate gas concentration at locations specified in the "
                                                          "file tracking_points.txt"
    )
    args = parser.parse_args()
    plot = args.plot
    plot_ex_prob = args.plot_ex_prob
    ex_prob_in = args.ex_prob
    time_steps_in = args.time_steps
    levels_in = args.levels
    days_plot_in = args.days_plot
    species_in = args.species
    original_specie = args.tracking_specie
    nproc = args.nproc
    convert = args.convert
    persistence = args.persistence
    models = args.models
    units = args.units
    plot_limits_in = args.plot_limits
    plot_isolines_in = args.plot_isolines
    time_av = args.time_av
    output_format = args.output_format
    plot_topography = args.plot_topography
    dz_lines_res = args.topography_isolines
    plot_resolution = args.plot_resolution
    tracking_points = args.tracking_points
    ex_prob = ex_prob_in.split(',')
    time_steps = time_steps_in.split(',')
    levels = levels_in.split(',')
    days_plot = days_plot_in.split(',')
    plot_limits = plot_limits_in.split(',')
    plot_isolines_s = plot_isolines_in.split(',')
    plot_isolines = []
    species = species_in.split(',')
    days_to_plot_in = []
    try:
        max_number_processes = int(os.environ["SLURM_NTASKS"])
    except (ValueError, KeyError) as e:
        try:
            max_number_processes = int(nproc)
        except ValueError:
            print("Please provide a valid number for the maximum number of process")
            sys.exit()
    if plot.lower() == "true":
        plot = True
        if days_plot_in == '':
            print("ERROR. Please specify at least one day to plot when --plot==True")
            sys.exit()
        else:
            for day in days_plot:
                if day == 'all':
                    days_to_plot_in.append(day)
                else:
                    try:
                        day_datetime = datetime.datetime.strptime(day, '%Y%m%d')
                        days_to_plot_in.append(day_datetime.strftime('%Y%m%d'))
                    except ValueError:
                        print('ERROR. Wrong format for -D -days_plot')
                        sys.exit()
    elif plot.lower() == "false":
        plot = False
    else:
        print("ERROR. Wrong value for variable -P --plot")
        sys.exit()
    if plot_ex_prob.lower() == "true":
        plot_ex_prob = True
        if ex_prob_in == '':
            print(
                "ERROR. Please specify at least one exceedance probability to plot when --plot_ex_prob==True"
            )
            sys.exit()
    elif plot_ex_prob.lower() == "false":
        plot_ex_prob = False
    else:
        print("ERROR. Wrong value for variable -PE --plot_ex_prob")
        sys.exit()
    if plot or plot_ex_prob:
        if time_steps_in == '':
            print("ERROR. Please specify at least one time step to plot")
            sys.exit()
        if levels_in == '':
            print("ERROR. Please specify at least one level to plot")
            sys.exit()
    if original_specie == None:
        print('ERROR. Please specify the name of the tracked specie')
        sys.exit()
    if species_in == '':
        print("ERROR. Please specify at least one gas specie name")
        sys.exit()
    if convert.lower() == "true":
        convert = True
    elif convert.lower() == "false":
        convert = False
    else:
        print("ERROR. Wrong value for variable -C --convert")
        sys.exit()
    if persistence.lower() == "true":
        persistence= True
    elif persistence.lower() == "false":
        persistence = False
    else:
        print("ERROR. Wrong value for variable -PER --persistence")
        sys.exit()
    exceedance_probabilities = []
    if ex_prob_in != '':
        for prob in ex_prob:
            exceedance_probabilities.append(float(prob))
    if (
        models.lower() != "disgas"
        and models.lower() != "twodee"
        and models.lower() != "all"
    ):
        print("ERROR. Wrong value for variable -M --models")
        sys.exit()
    try:
        units = units.lower()
    except AttributeError:
        print("Please provide an option for -U --units")
        sys.exit()
    if units != "ppm" and units != "kg/m3":
        print("ERROR. Wrong value for variable -U --units")
        sys.exit()
    try:
        time_av = int(time_av)
    except TypeError:
        time_av = None
    except ValueError:
        try:
            time_av = int(float(time_av))
        except ValueError:
            print("ERROR. Please specify a valid time-averaging interval")
            sys.exit()
    min_con = max_con = -1.0
    if len(plot_limits) > 1:
        try:
            min_con = float(plot_limits[0])
            max_con = float(plot_limits[1])
        except ValueError:
            print(
                "ERROR. Please specify valid minimum and maximum concentration -PL --plot_limits"
            )
            sys.exit()
    if len(plot_isolines_s) >= 1:
        for isoline in plot_isolines_s:
            if isoline != '':
                try:
                    plot_isolines.append(float(isoline))
                except ValueError:
                    print("WARNING. Wrong entry for -PI --plot_isolines. Continuing discarding concentration contour lines")
                    plot_isolines = []
                    break
    if output_format.lower() != "grd":
        print(
            "ERROR. Please specify a valid output format. Current valid options are: GRD"
        )
        sys.exit()
    else:
        output_format = "grd"
    if plot_topography.lower() == "true":
        plot_topography_layer = True
    elif plot_topography.lower() == "false":
        plot_topography_layer = False
    else:
        print("ERROR. Wrong value for variable -PT --plot_topography")
        sys.exit()
    if plot_topography_layer:
        try:
            dz_lines_res = float(dz_lines_res)
        except ValueError:
            dz_lines_res = 100
    try:
        plot_resolution = int(plot_resolution)
    except ValueError:
        print("ERROR. Please provide a valid number for -PR --plot_resolution")
        sys.exit()
    if tracking_points.lower() == 'true':
        tracking_points = True
        if not os.path.isfile('tracking_points.txt'):
            print("WARNING. Tracking points option activated but file tracking_points.txt not found. Continuing "
                  "without this option")
            tracking_points = False
    elif tracking_points.lower() == 'false':
        tracking_points = False
    else:
        print("ERROR. Wrong value for variable -TP --tracking_points")
        sys.exit()
    return (
        plot,
        plot_ex_prob,
        time_steps,
        levels,
        days_to_plot_in,
        species,
        original_specie,
        exceedance_probabilities,
        max_number_processes,
        convert,
        persistence,
        models,
        units,
        time_av,
        min_con,
        max_con,
        plot_isolines,
        output_format,
        plot_topography_layer,
        dz_lines_res,
        plot_resolution,
        tracking_points
    )


def folder_structure():
    original_output_folder_name = "simulations"
    post_processing = "post_processing"
    processed_output_folder_name = original_output_folder_name + "_processed"
    ecdf_folder_name = "output_ecdf"
    persistence_folder_name = "output_persistence"
    disgas_outputs = os.path.join(root, post_processing, "disgas")
    twodee_outputs = os.path.join(root, post_processing, "twodee")
    disgas_original_output_folder = os.path.join(
        root, original_output_folder_name, "disgas"
    )
    twodee_original_output_folder = os.path.join(
        root, original_output_folder_name, "twodee"
    )
    disgas_processed_output_folder = os.path.join(
        disgas_outputs, processed_output_folder_name
    )
    twodee_processed_output_folder = os.path.join(
        twodee_outputs, processed_output_folder_name
    )
    disgas_ecdf = os.path.join(disgas_outputs, ecdf_folder_name)
    twodee_ecdf = os.path.join(twodee_outputs, ecdf_folder_name)
    disgas_ecdf_tracking_points = os.path.join(disgas_ecdf, 'tracking_points')
    twodee_ecdf_tracking_points = os.path.join(twodee_ecdf, 'tracking_points')
    disgas_persistence = os.path.join(disgas_outputs, persistence_folder_name)
    twodee_persistence = os.path.join(twodee_outputs, persistence_folder_name)
    try:
        os.mkdir(post_processing)
    except FileExistsError:
        print("Folder post_processing already exists")
    if models == "disgas" or models == "all":
        try:
            os.mkdir(disgas_outputs)
        except FileExistsError:
            print("Folder " + disgas_outputs + " already exists")
        try:
            os.mkdir(disgas_processed_output_folder)
        except FileExistsError:
            print("Folder " + disgas_processed_output_folder + " already exists")
        try:
            os.mkdir(disgas_ecdf)
        except FileExistsError:
            print("Folder " + disgas_ecdf + " already exists")
        try:
            os.mkdir(disgas_ecdf_tracking_points)
        except FileExistsError:
            print("Folder " + disgas_ecdf_tracking_points + " already exists")
        try:
            os.mkdir(disgas_persistence)
        except FileExistsError:
            print("Folder " + disgas_persistence + " already exists")
    if models == "twodee" or models == "all":
        try:
            os.mkdir(twodee_outputs)
        except FileExistsError:
            print("Folder " + twodee_outputs + " already exists")
        try:
            os.mkdir(twodee_processed_output_folder)
        except FileExistsError:
            print("Folder " + twodee_processed_output_folder + " already exists")
        try:
            os.mkdir(twodee_ecdf)
        except FileExistsError:
            print("Folder " + twodee_ecdf + " already exists")
        try:
            os.mkdir(twodee_ecdf_tracking_points)
        except FileExistsError:
            print("Folder " + twodee_ecdf_tracking_points + " already exists")
        try:
            os.mkdir(twodee_persistence)
        except FileExistsError:
            print("Folder " + twodee_persistence + " already exists")
    twodee_input_file = os.path.join(root, "twodee.inp")
    twodee_output_time_step = 0
    if models == "all":
        models_to_elaborate = ["disgas", "twodee"]
    elif models == "disgas":
        models_to_elaborate = ["disgas"]
    else:
        models_to_elaborate = ["twodee"]
        # Read the output time interval from the twodee input file
        with open(twodee_input_file, "r") as twodee_file:
            for line in twodee_file:
                if "OUTPUT_INTERVAL_(SEC)" in line:
                    twodee_output_time_step = float(line.split("=")[1])
        if twodee_output_time_step == 0:
            print("Unable to read the Twodee output time step")
            sys.exit()
    for model in models_to_elaborate:
        if model == "disgas":
            model_outputs = disgas_outputs
            model_processed_output_folder = disgas_processed_output_folder
            ecdf_outputs = disgas_ecdf
            persistence_outputs = disgas_persistence
        else:
            model_outputs = twodee_outputs
            model_processed_output_folder = twodee_processed_output_folder
            ecdf_outputs = twodee_ecdf
            persistence_outputs = twodee_persistence
        graphical_outputs = os.path.join(model_outputs, "graphical_outputs")
        graphical_outputs_simulations = os.path.join(graphical_outputs, "simulations")
        graphical_outputs_ecdf = os.path.join(graphical_outputs, "ecdf")
        graphical_outputs_ecdf_tracking_points = os.path.join(graphical_outputs_ecdf, "tracking_points")
        try:
            os.mkdir(graphical_outputs)
        except FileExistsError:
            print("Folder " + graphical_outputs + " already exists")
        try:
            os.mkdir(graphical_outputs_simulations)
        except FileExistsError:
            print("Folder " + graphical_outputs_simulations + " already exists")
        try:
            os.mkdir(graphical_outputs_ecdf)
        except FileExistsError:
            print("Folder " + graphical_outputs_ecdf + " already exists")
        try:
            os.mkdir(graphical_outputs_ecdf_tracking_points)
        except FileExistsError:
            print("Folder " + graphical_outputs_ecdf_tracking_points + " already exists")

    return (
        disgas_outputs,
        disgas_original_output_folder,
        disgas_processed_output_folder,
        ecdf_folder_name,
        disgas_ecdf,
        disgas_ecdf_tracking_points,
        disgas_persistence,
        twodee_outputs,
        twodee_original_output_folder,
        twodee_processed_output_folder,
        twodee_ecdf,
        twodee_ecdf_tracking_points,
        twodee_persistence,
        models_to_elaborate,
        twodee_output_time_step,
        model_processed_output_folder,
        ecdf_outputs,
        persistence_outputs,
        graphical_outputs,
        graphical_outputs_simulations,
        graphical_outputs_ecdf,
        graphical_outputs_ecdf_tracking_points
    )


def gas_properties():
    def extract_gas_properties(specie):
        data = pd.read_csv(gas_properties_file, on_bad_lines='skip')
        molar_ratio = None
        conc_thresholds = []
        exp_times = []
        if convert:
            if specie == original_specie:
                molar_ratio = 1
            else:
                try:
                    x = np.sort(data[specie + '/' + original_specie])
                    list_x = list(x)
                    samples = random.sample(list_x, 1)
                    molar_ratio = samples[0]
                except KeyError:
                    print('ERROR. Molar ratio ' + specie + '/' + original_specie + ' not found in gas_properties.csv')
                    exit()
        try:
            y = np.sort(data['M_' + specie])
            molar_weight = list(y)[0]
            if molar_weight != molar_weight:
                print('ERROR. Molar weight of ' + specie + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Molar weight of ' + specie + ' not found in gas_properties.csv')
            sys.exit()
        try:
            y = data['CT_' + specie]
            conc_thresholds_temp = list(y)
            conc_thresholds = [x for x in conc_thresholds_temp if x == x]
            if len(conc_thresholds) == 0:
                print('WARNING. Concentration thresholds of ' + specie + ' not found in gas_properties.csv. Gas '
                      'persistence calculation not possible')
        except KeyError:
            print('WARNING. Concentration thresholds of ' + specie + ' not found in gas_properties.csv. Gas '
                  'persistence calculation not possible')
        try:
            y = data['ET_' + specie]
            exp_times_temp = list(y)
            exp_times = [x for x in exp_times_temp if x == x]
            if len(exp_times) == 0:
                print('WARNING. Exposure times of ' + specie + ' not found in gas_properties.csv. Gas '
                      'persistence calculation not possible for gas specie ' + specie)
        except KeyError:
            print('WARNING. Exposure times of ' + specie + ' not found in gas_properties.csv. Gas '
                  'persistence calculation not possible for gas specie' + specie)
        return molar_ratio, molar_weight, conc_thresholds, exp_times

    gas_properties_file = os.path.join(root, "gas_properties.csv")
    try:
        open(gas_properties_file, "r")
    except FileNotFoundError:
        print("File " + gas_properties_file + " not present")
        sys.exit()
    molar_ratios = []
    molar_weights = []
    concentration_thresholds = []
    exposure_times = []
    for specie in species:
        molar_ratio, molar_weight, concentration_thresholds_specie, exposure_times_specie = extract_gas_properties(
            specie)
        molar_ratios.append(molar_ratio)
        molar_weights.append(molar_weight)
        concentration_thresholds.append(concentration_thresholds_specie)
        exposure_times.append(exposure_times_specie)
    molar_ratio, molar_weight, concentration_thresholds_specie, exposure_times_specie = extract_gas_properties(
        original_specie)
    molar_ratios_tracking_specie = molar_ratio
    molar_weights_tracking_specie = molar_weight
    concentration_thresholds_tracking_specie = concentration_thresholds_specie
    exposure_times_tracking_specie = exposure_times_specie
    species_properties = []
    for i in range(0, len(species)):
        gas_specie = {}
        gas_specie["specie_name"] = species[i]
        gas_specie["molar_ratio"] = molar_ratios[i]
        gas_specie["molar_weight"] = molar_weights[i]
        gas_specie["concentration_thresholds"] = concentration_thresholds[i]
        gas_specie["exposure_times"] = exposure_times[i]
        species_properties.append(gas_specie)
    gas_specie = {}
    gas_specie["specie_name"] = original_specie
    gas_specie["molar_ratio"] = molar_ratios_tracking_specie
    gas_specie["molar_weight"] = molar_weights_tracking_specie
    gas_specie["concentration_thresholds"] = concentration_thresholds_tracking_specie
    gas_specie["exposure_times"] = exposure_times_tracking_specie
    species_properties.append(gas_specie)
    return species_properties


def domain(model):
    output_levels = []
    if model == "disgas":
        with open(file="disgas.inp") as input_file:
            for record in input_file:
                try:
                    record_splitted = record.split("=")
                    temp = record_splitted[1].split("(")
                    if "SIMULATION_INTERVAL_(SEC)" in record_splitted[0]:
                        tot_time = int(temp[0])
                    elif "NX" in record_splitted[0]:
                        nx = int(temp[0])
                    elif "NY" in record_splitted[0]:
                        ny = int(temp[0])
                    elif "NZ" in record_splitted[0]:
                        nz = int(temp[0])
                    elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                        dt = int(temp[0])
                    elif "DX_(M)" in record_splitted[0]:
                        dx = float(temp[0])
                    elif "DY_(M)" in record_splitted[0]:
                        dy = float(temp[0])
                    elif "X_ORIGIN_(UTM_M)" in record_splitted[0]:
                        x0 = float(temp[0])
                    elif "Y_ORIGIN_(UTM_M)" in record_splitted[0]:
                        y0 = float(temp[0])
                    elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                        dt = float(temp[0])
                    #New
                    elif "HOUR" in record_splitted[0]:
                        hour_start = int(temp[0])
                    elif "MINUTE" in record_splitted[0]:
                        minute_start = int(temp[0])
                    elif "Z_LAYERS_(M)" in record_splitted[0]:
                        heights = temp[0]
                        heights_list = heights.split(' ')
                        for height in heights_list:
                            try:
                                output_levels.append(float(height))
                            except ValueError:
                                continue
                        output_levels = sorted(output_levels)
                except (IndexError, ValueError):
                    continue
    else:
        with open(file="twodee.inp") as input_file:
            for record in input_file:
                try:
                    record_splitted = record.split("=")
                    temp = record_splitted[1].split("(")
                    if "SIMULATION_INTERVAL_(SEC)" in record_splitted[0]:
                        tot_time = int(temp[0])
                    elif "NX" in record_splitted[0]:
                        nx = int(temp[0])
                    elif "NY" in record_splitted[0]:
                        ny = int(temp[0])
                    elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                        dt = int(temp[0])
                    elif "DX_(M)" in record_splitted[0]:
                        dx = float(temp[0])
                    elif "DY_(M)" in record_splitted[0]:
                        dy = float(temp[0])
                    elif "X_ORIGIN_(UTM_M)" in record_splitted[0]:
                        x0 = float(temp[0])
                    elif "Y_ORIGIN_(UTM_M)" in record_splitted[0]:
                        y0 = float(temp[0])
                    elif "HEIGHTS_(M)" in record_splitted[0]:
                        heights = temp[0]
                        heights_list = heights.split(' ')
                        for height in heights_list:
                            try:
                                output_levels.append(float(height))
                            except ValueError:
                                continue
                        output_levels = sorted(output_levels)
                    elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                        dt = float(temp[0])
                    # New
                    elif "HOUR" in record_splitted[0]:
                        hour_start = int(temp[0])
                    elif "MINUTE" in record_splitted[0]:
                        minute_start = int(temp[0])
                except (IndexError, ValueError):
                    continue
    yf = y0 + (ny - 1) * dy
    xf = x0 + (nx - 1) * dx
    n_time_steps = int(tot_time / dt)
    nz = len(output_levels)
    return x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, output_levels, hour_start, minute_start


def elaborate_tracking_points():
    import utm
    stations = []
    tracking_points_file = os.path.join(root, 'tracking_points.txt')
    station_id = 0
    with open(tracking_points_file, 'r') as tracking_points_file_read:
        for line in tracking_points_file_read:
            x = line.split('\t')[0]
            y = line.split('\t')[1]
            z = line.split('\t')[2]
            try:
                station_x = float(x)
                station_y = float(y)
                station_z = float(z)
            except ValueError:
                continue
            if -90 <= station_y <= 90 and -180 <= station_x <= 180:  # Input is in lat-lon
                try:
                    out_utm = utm.from_latlon(station_x, station_y)
                    station_easting = float(out_utm[0])
                    station_northing = float(out_utm[1])
                except ValueError:
                    print(
                        "WARNING. Invalide coordinate of the tracking point"
                    )
            else:
                station_northing = station_y
                station_easting = station_x
            if y0 <= station_northing <= yf and x0 <= station_easting <= xf and \
                    min(output_levels) <= station_z <= max(output_levels):
                station_id += 1
                stations.append({'station_id': station_id, 'easting': station_easting,
                                        'northing': station_northing, 'elevation': station_z})
            else:
                continue
    return stations


def extract_tracking_points(files_to_interpolate, j):
    def interpolate(x, y, z, levels_interpolation, files):
        from scipy import interpolate
        x_array = np.linspace(x0, xf, nx, endpoint=True)
        y_array = np.linspace(y0, yf, ny, endpoint=True)
        Z1 = np.loadtxt(files[0], skiprows=5)
        if len(levels_interpolation) == 2:
            z_array = np.linspace(levels_interpolation[0], levels_interpolation[1], 2)
            Z2 = np.loadtxt(files[1], skiprows=5)
            Z = np.array([Z1, Z2])
            my_interpolating_function = interpolate.RegularGridInterpolator((x_array, y_array, z_array), Z.T)
            pt = np.array([x, y, z])
        else:
            Z = Z1
            my_interpolating_function = interpolate.RegularGridInterpolator((x_array, y_array), Z.T)
            pt = np.array([x, y])
        return my_interpolating_function(pt)

    files_time_steps = []
    files_time_averaging_steps = []
    for file in files_to_interpolate:
        file_name = file.split(os.sep)[-1]
        file_time_step = file_name.split('_')[2]
        try:
            file_time_step = int(file_time_step.split('.grd')[0])
            if file_time_step not in files_time_steps:
                files_time_steps.append(file_time_step)
        except ValueError:
            file_time_step_avg = file_time_step.split('.grd')[0]
            if file_time_step_avg not in files_time_averaging_steps:
                files_time_averaging_steps.append(file_time_step_avg)
    files_time_steps = sorted(files_time_steps)
    files_time_averaging_steps = sorted(files_time_averaging_steps)
    levels_for_interpolation = []
    for l in range(0, len(species)):
        k = 0
        for station in stations:
            if min(output_levels) <= station['elevation'] <= max(output_levels):
                for i in range(1, len(output_levels) + 1):
                    levels_for_interpolation = []
                    if output_levels[i - 1] == station['elevation']:
                        levels_for_interpolation.append(output_levels[i - 1])
                        break
                    elif output_levels[i - 1] < station['elevation'] < output_levels[i]:
                        levels_for_interpolation.append(output_levels[i - 1])
                        levels_for_interpolation.append(output_levels[i])
                        break
                    else:
                        continue
            i = 0
            for time_step in files_time_steps:
                files_to_use = []
                files_levels = []
                file_directory = ''
                for file in files_to_interpolate:
                    file_name = file.split(os.sep)[-1]
                    file_directory = file.split(file_name)[0]
                    file_level = file_name.split('_')[1]
                    file_level = float(file_level.split('mabg')[0])
                    file_time_step = file_name.split('_')[2]
                    try:
                        file_time_step = int(file_time_step.split('.grd')[0])
                    except ValueError:
                        continue
                    if file_level in levels_for_interpolation and file_time_step == time_step:
                        files_to_use.append(file)
                        files_levels.append(file_level)
                files_levels = sorted(files_levels)
                c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'], files_levels,
                                             files_to_use)
                c[l][k][j][i] = c_interpolated[0]
                tracking_point_file = os.path.join(file_directory, 'TP_' + str(station['station_id']) + '.txt')
                tracking_points_files.append(tracking_point_file)
                with open(tracking_point_file, 'a') as tracking_point_file:
                    tracking_point_file.write(str(time_step) + '\t' + "{0:.2e}".format(c_interpolated[0]) + '\n')
                i += 1
            for time_step in files_time_averaging_steps:
                files_to_use = []
                files_levels = []
                file_directory = ''
                for file in files_to_interpolate:
                    file_name = file.split(os.sep)[-1]
                    file_directory = file.split(file_name)[0]
                    file_level = file_name.split('_')[1]
                    file_level = float(file_level.split('mabg')[0])
                    file_time_step = file_name.split('_')[2]
                    try:
                        file_time_step = int(file_time_step.split('.grd')[0])
                    except ValueError:
                        file_time_step = file_time_step.split('.grd')[0]
                        if file_level in levels_for_interpolation and file_time_step == time_step:
                            files_to_use.append(file)
                            files_levels.append(file_level)
                files_levels = sorted(files_levels)
                c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'], files_levels,
                                             files_to_use)
                c[l][k][j][i] = c_interpolated[0]
                tracking_point_file = os.path.join(file_directory, 'TP_' + str(station['station_id']) + '.txt')
                tracking_points_files.append(tracking_point_file)
                with open(tracking_point_file, 'a') as tracking_point_file:
                    tracking_point_file.write(str(time_step) + '\t' + "{0:.2e}".format(c_interpolated[0]) + '\n')
                i += 1
            k += 1
    return files_time_steps + files_time_averaging_steps


def probabilistic_tracking_points():
    def plot_hazard_curves(file_input, folder):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        c_ecdf = []
        with open(file_input, 'r') as tp_ecdf_file:
            lines = tp_ecdf_file.readlines()[1:]
            time_steps = []
            for line in lines:
                entries = line.split('\t')[1:]
                time_steps.append(line.split('\t')[0])
                for i in range(0, len(entries)):
                    entries[i] = float(entries[i])
                c_ecdf.append(entries)
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        tp_file = file_input.split(os.sep)[-1]
        tp_file = tp_file.split('.txt')[0]
        specie_name = file_input.split(os.sep)[-2]
        specie_name = specie_name.translate(SUB)
        for i in range(0, len(c_ecdf)):
            ii = 1
            if 'tavg' in str(time_steps[i]):
                output_plot_file = os.path.join(folder, tp_file + '_tavg-time_interval' + str(ii) + '.png')
                ii += 1
            else:
                output_plot_file = os.path.join(folder, tp_file + '_time_step' + str(i + 1) + '.png')
            plt.plot(c_ecdf[i], ex_probabilities)
            if units == "ppm":
                plt.title(specie_name + " concentration [ppm]")
                plt.xlabel("C [ppm]")
            else:
                plt.title(specie_name + " concentration [kg m$\mathregular{^{-3}}$]")
                plt.xlabel("C [kg m$\mathregular{^{-3}}$]")
            plt.xlim(min_con_tp, max_con_tp)
            plt.ylabel("Exceedance probability")
            image_buffer = StringIO()
            plt.tight_layout()
            plt.savefig(output_plot_file, dpi=plot_resolution)
            image_buffer.close()
            plt.close()


    delta_quantile = 0.01
    quantiles = []
    ex_probabilities = []
    quantile_string = ''
    for m in range(0, int(1 / delta_quantile + 1)):
        quantiles.append(m * delta_quantile)
        ex_probabilities.append(1 - quantiles[-1])
        quantile_string += '\t' + "{0:.3f}".format(quantiles[-1])
    output_quantile = [[[[0 for i in range(0, len(all_time_steps))] for m in range(0, len(quantiles))]
                       for k in range(0, len(stations))] for l in range (0, len(stations))]
    c_list = []
    for specie in species:
        graphical_outputs_ecdf_tracking_points_specie = os.path.join(graphical_outputs_ecdf_tracking_points, specie)
        try:
            os.mkdir(graphical_outputs_ecdf_tracking_points_specie)
        except FileExistsError:
            print('Folder ' + graphical_outputs_ecdf_tracking_points_specie + ' already exists')
        if model == 'disgas':
            disgas_ecdf_tracking_points_specie = os.path.join(disgas_ecdf_tracking_points, specie)
            try:
                os.mkdir(disgas_ecdf_tracking_points_specie)
            except FileExistsError:
                print('Folder ' + disgas_ecdf_tracking_points_specie + ' already exists')
        else:
            twodee_ecdf_tracking_points_specie = os.path.join(twodee_ecdf_tracking_points, specie)
            try:
                os.mkdir(twodee_ecdf_tracking_points_specie)
            except FileExistsError:
                print('Folder ' + twodee_ecdf_tracking_points_specie + ' already exists')
    min_con_tp = 1000000000000000
    max_con_tp = 0
    for l in range(0, len(species)):
        for k in range(0, len(stations)):
            for i in range(0, len(all_time_steps)):
                for j in range(0, len(days)):
                    c_list.append(c[l][k][j][i])
                for m in range(0, len(quantiles)):
                    output_quantile[l][k][m][i] = np.quantile(c_list, q=quantiles[m])
                    if output_quantile[l][k][m][i] <= min_con_tp:
                        min_con_tp = output_quantile[l][k][m][i]
                    if output_quantile[l][k][m][i] >= max_con_tp:
                        max_con_tp = output_quantile[l][k][m][i]
                c_list = []
    for l in range(0, len(species)):
        for k in range(0, len(stations)):
            station = stations[k]
            if model == 'disgas':
                ecdf_tracking_point_file = os.path.join(disgas_ecdf_tracking_points, species[l],
                                                        'TP_' + str(station['station_id']) + '_ecdf.txt')
            else:
                ecdf_tracking_point_file = os.path.join(twodee_ecdf_tracking_points, species[l],
                                                        'TP_' + str(station['station_id']) + '_ecdf.txt')
            with open(ecdf_tracking_point_file, 'w') as output_file:
                output_file.write(quantile_string + '\n')
                for i in range(0, len(all_time_steps)):
                    ii = 1
                    if 'tavg' in str(all_time_steps[i]):
                        output_quantile_string = 'tavg-time_interval_' + str(ii)
                        ii += 1
                    else:
                        output_quantile_string = 'time_step_' + str(i + 1)
                    for m in range(0, len(quantiles)): \
                            output_quantile_string += '\t' + "{0:.2e}".format(output_quantile[l][k][m][i])
                    output_file.write(output_quantile_string + '\n' )
            plot_file_folder = os.path.join(graphical_outputs_ecdf_tracking_points, specie)
            plot_hazard_curves(ecdf_tracking_point_file, plot_file_folder)


def extract_days():
    days = []
    days_list_path = os.path.join(root, "days_list.txt")
    days_to_plot = []
    with open(days_list_path, "r") as days_list:
        for line in days_list:
            day_temp = line.split(" ")[0]
            day_temp = day_temp.split("-")
            day = day_temp[0] + day_temp[1] + day_temp[2]
            days.append(day)
            try:
                if days_to_plot_in[0] == "all":
                    days_to_plot.append(day)
                else:
                    for day_to_plot in days_to_plot_in:
                        if day_to_plot == day:
                            days_to_plot.append(day_to_plot)
            except IndexError:
                print('No days to plot')
    return days, days_to_plot


def converter(input_file, processed_file, specie_input, model):
    Z = np.loadtxt(input_file, skiprows=5)
    Z[Z < 0] = 0
    if units == "ppm":
        if model == "disgas":
            file_time_step = os.path.split(processed_file)[1]
            file_time_step = file_time_step.split("_")[2]
            file_time_step = file_time_step.split(".grd")[0]
            file_time_h = file_time_step[-4:]
            file_name = input_file.split(os.sep)[-1]
            file_folder = input_file.split(file_name)[0]
            file_folder_daily = file_folder.split("outfiles")[0]
            surface_data = os.path.join(file_folder_daily, "surface_data.txt")
            with open(surface_data) as surface_data_file:
                for line in surface_data_file:
                    try:
                        records = line.split("\t")
                    except BaseException:
                        continue
                    if file_time_h == records[0]:
                        t2m = float(records[2])
                        p2m = float(records[3]) / 100  # in hPa for this conversion
                        break
            for specie in species_properties:
                if specie["specie_name"] == specie_input:
                    molar_weight = specie["molar_weight"]
            conversion_factor = (
                (22.4 / molar_weight) * (t2m / 273) * (1013 / p2m)
            ) * 1000000
            Z_converted = np.multiply(Z, conversion_factor)  # convert kg/m3 to ppm
        else:
            Z_converted = Z
    else:
        if model == "twodee":
            file_time_step = os.path.split(processed_file)[1]
            file_time_step = file_time_step.split("_")[2]
            file_time_step = int(file_time_step.split(".grd")[0])
            file_time_h = file_time_step[-4:]
            file_name = input_file.split(os.sep)[-1]
            file_folder = input_file.split(file_name)[0]
            file_folder_daily = file_folder.split("outfiles")[
                0
            ]  # To be changed when folder structures will change
            surface_data = os.path.join(file_folder_daily, "surface_data.txt")
            with open(surface_data) as surface_data_file:
                for line in surface_data_file:
                    try:
                        records = line.split("\t")
                    except BaseException:
                        continue
                    if file_time_h == records[0]:
                        t2m = float(records[2])
                        p2m = float(records[3]) / 100  # in hPa for this conversion
                        break
            for specie in species_properties:
                if specie["specie_name"] == specie_input:
                    molar_weight = specie["molar_weight"]
            conversion_factor = (
                (molar_weight / 22.4) * (273 / t2m) * (p2m / 1013)
            ) / 1000000
            Z_converted = np.multiply(Z, conversion_factor)  # convert ppm to kg/m3
        else:
            Z_converted = Z
    try:
        np.loadtxt(processed_file, skiprows=5)
    except OSError:
        with open(processed_file, "a") as processed_file:
            if output_format == "grd":
                processed_file.write("DSAA\n")
                processed_file.write(str(nx) + "  " + str(ny) + "\n")
                processed_file.write(str(x0) + "  " + str(xf) + "\n")
                processed_file.write(str(y0) + "  " + str(yf) + "\n")
            if not convert:
                processed_file.write(
                    str(np.amin(Z_converted)) + "  " + str(np.amax(Z_converted)) + "\n"
                )
                np.savetxt(processed_file, Z_converted, fmt="%.2e")
            else:
                for specie in species_properties:
                    if specie["specie_name"] == specie_input:
                        molar_ratio = specie["molar_ratio"]
                        if units == 'ppm':
                            species_conversion_factor = molar_ratio
                        elif units == 'kg/m3':
                            for specie in species_properties:
                                if specie["specie_name"] == specie_input:
                                    molar_weight = specie["molar_weight"]
                                if specie["specie_name"] == original_specie:
                                    molar_weight_tracking_specie = specie["molar_weight"]
                            species_conversion_factor = molar_ratio * (molar_weight / molar_weight_tracking_specie)
                Z_converted = np.multiply(Z_converted, species_conversion_factor)
                processed_file.write(
                    str(np.amin(Z_converted)) + "  " + str(np.amax(Z_converted)) + "\n"
                )
                np.savetxt(processed_file, Z_converted, fmt="%.2e")


def time_average(files_to_average, outfile):
    Z_sum = 0
    for file in files_to_average:
        if output_format == "grd":
            Z = np.loadtxt(file, skiprows=5)
        Z_sum += Z
    Z_avg = np.divide(Z_sum, len(files_to_average))
    # Create header of the processed file
    with open(outfile, "a") as processed_file:
        if output_format == "grd":
            processed_file.write("DSAA\n")
            processed_file.write(str(nx) + "  " + str(ny) + "\n")
            processed_file.write(str(x0) + "  " + str(xf) + "\n")
            processed_file.write(str(y0) + "  " + str(yf) + "\n")
            processed_file.write(
                str(np.amin(Z_avg)) + "  " + str(np.amax(Z_avg)) + "\n"
            )
        np.savetxt(processed_file, Z_avg, fmt="%.2e")


def elaborate_day(day_input, model):
    if model == "disgas":
        model_output_folder = os.path.join(
            disgas_original_output_folder, day_input, "outfiles"
        )
        model_processed_output_folder_daily = os.path.join(
            disgas_processed_output_folder, day_input
        )
    else:
        model_output_folder = os.path.join(
            twodee_original_output_folder, day_input, "outfiles"
        )
        model_processed_output_folder_daily = os.path.join(
            twodee_processed_output_folder, day_input
        )
    try:
        os.mkdir(model_processed_output_folder_daily)
    except FileExistsError:
        print("Folder " + model_processed_output_folder_daily + " already exists")
    for specie in species:
        model_processed_output_folder_specie = os.path.join(
            model_processed_output_folder_daily, specie
        )
        try:
            os.mkdir(model_processed_output_folder_specie)
        except FileExistsError:
            print("Folder " + model_processed_output_folder_specie + " already exists")
        except PermissionError:  # retry
            try:
                os.mkdir(model_processed_output_folder_specie)
            except FileExistsError:
                print(
                    "Folder " + model_processed_output_folder_specie + " already exists"
                )
    files_list_temp = os.listdir(model_output_folder)
    files_list_path = []
    files_list = []
    models = []
    for file in files_list_temp:
        if file[0:2] == "c_":
            files_list.append(file)
            for specie in species:
                files_list_path.append(os.path.join(model_output_folder, file))  #
                models.append(model)
    time_steps_disgas = []
    if model == 'disgas':
        for file in files_list:
            file_name_splitted = file.split("_")
            file_time_step = file_name_splitted[2]
            file_time_step = int(file_time_step.split(".")[0])
            if file_time_step not in time_steps_disgas:
                time_steps_disgas.append(file_time_step)
    time_steps_disgas = sorted(time_steps_disgas)
    converted_files = []
    processed_files = []
    species_list = []
    processed_files_species = []
    time_steps = []
    time_start = datetime.datetime.strptime(day_input + str(hour_start).zfill(2) + str(minute_start).zfill(2),
                                            '%Y%m%d%H%M')
    for specie in species:
        processed_files_specie = []
        for file in files_list:
            species_list.append(specie)
            if model == "twodee":
                file_name_splitted = file.split("_")
                file_level = file_name_splitted[1]
                file_time_step = int(file_name_splitted[2].split(".")[0])
                file_level = float(file_level.split("cm")[0]) / 100
                file_level_s = "{0:.3f}".format(file_level) + 'mabg'
            else:
                file_name_splitted = file.split("_")
                file_level = int(file_name_splitted[1])
                file_level_s = "{0:.3f}".format(output_levels[file_level - 1]) + 'mabg'
                file_time_step = file_name_splitted[2]
                file_time_step = int(file_time_step.split(".")[0])
                file_time_step = time_steps_disgas.index(file_time_step) * dt
            if file_level_s not in processed_files_levels:
                processed_files_levels.append(file_level_s)
            time_validity = time_start + datetime.timedelta(seconds=file_time_step)
            file_validity = datetime.datetime.strftime(time_validity, '%Y%m%d%H%M')
            if time_validity not in time_steps:
                time_steps.append(time_validity)
            if file_validity not in processed_files_steps:
                processed_files_steps.append(file_validity)
            file = "c_" + file_level_s + "_" + file_validity + ".grd"
            converted_file = file
            converted_files.append(converted_file)
            processed_files.append(
                os.path.join(
                    os.path.join(model_processed_output_folder_daily, specie),
                    converted_file,
                )
            )
            processed_files_specie.append(
                os.path.join(
                    os.path.join(model_processed_output_folder_daily, specie),
                    converted_file,
                )
            )
        processed_files_species.append(processed_files_specie)
    n_elaborated_files = 0
    while n_elaborated_files < len(files_list_path):
        start = n_elaborated_files
        end = n_elaborated_files + max_number_processes
        if end > len(files_list_path):
            end = len(files_list_path)
        try:
            pool_files = ThreadingPool(max_number_processes)
            pool_files.map(
                converter,
                files_list_path[start:end],
                processed_files[start:end],
                species_list[start:end],
                models[start:end],
            )
        except BaseException:
            print("Unable to convert files")
            sys.exit()
        n_elaborated_files = end
        if n_elaborated_files == len(files_list_path):
            break
    if time_av is not None:
        averaged_files = []
        time_min = min(time_steps)
        if time_av == 0:
            time_max = max(time_steps)
        else:
            time_max = time_min + datetime.timedelta(seconds=time_av * 3600)
            if time_max > max(time_steps):
                time_max = max(time_steps)
        while time_max <= max(time_steps):
            time_max_s = "{:02d}".format(int(datetime.datetime.strftime(time_min, '%H')) +
                                         int((time_max - time_min).seconds / 3600))
            if datetime.datetime.strftime(time_min, '%H') + "-" + time_max_s + "-tavg" not in tavg_intervals:
                tavg_intervals.append(datetime.datetime.strftime(time_min, '%H') + "-" + time_max_s + "-tavg")
            for i in range(0, len(species)):
                files_to_average = []
                for level in processed_files_levels:
                    files_in_level = []
                    for file in processed_files_species[i]:
                        file_level = os.path.split(file)[1]
                        file_level = file_level.split("_")[1]
                        if file_level == level:
                            files_in_level.append(file)
                    time_averaged_file = os.path.join(
                        os.path.join(model_processed_output_folder_daily, species[i]),
                        "c_"
                        + level
                        + "_"
                        + datetime.datetime.strftime(time_min, '%Y%m%d%H%M')
                        + "-"
                        + datetime.datetime.strftime(time_max, '%Y%m%d%H%M')
                        + "-tavg.grd",
                    )
                    for file in processed_files_species[i]:
                        file_name = os.path.split(file)[1]
                        file_level = file_name.split("_")[1]
                        file_validity_s = file_name.split("_")[2]
                        file_validity_s = file_validity_s.split(".")[0]
                        file_validity = datetime.datetime.strptime(file_validity_s, '%Y%m%d%H%M')
                        if file_level == level:
                            if time_min <= file_validity <= time_max:
                                files_to_average.append(file)
                                averaged_files.append(file)
                    processed_files.append(time_averaged_file)
                    time_average(files_to_average, time_averaged_file)
                    files_to_average = []
            if time_av == 0:
                break
            else:
                time_min = time_max
                if time_min == max(time_steps):
                    break
                time_max = time_min + datetime.timedelta(seconds=time_av * 3600)
                if time_max > max(time_steps):
                    time_max = max(time_steps)
                continue
    return processed_files


def sort_levels(input_array):
    output_array = []
    for level in input_array:
        level_float = float(level.split('mabg')[0])
        output_array.append(level_float)
    output_array = sorted(output_array)
    for i in range (0, len(output_array)):
        output_array[i] = "{0:.3f}".format(output_array[i]) + 'mabg'
    return output_array


def probabilistic_output(model):
    def calculate_quantiles():
        def ecdf(index):
            specie = index[1]
            file_level_s = index[2]
            time_step = index[3]
            ex_prob = index[0]
            quantile = 1 - ex_prob
            output_files = []
            for day in days:
                try:
                    file_time_step = int(time_step)
                    file_time_step = file_time_step * dt
                    time_start = datetime.datetime.strptime(
                        day + str(hour_start).zfill(2) + str(minute_start).zfill(2),
                        '%Y%m%d%H%M')
                    time_validity = time_start + datetime.timedelta(seconds=file_time_step)
                    file_validity = datetime.datetime.strftime(time_validity, '%Y%m%d%H%M')
                    time_step_s = "{:06d}".format(int(time_step))
                except ValueError:
                    tavg_interval_start_s = day + time_step.split('-')[0] + '00'
                    tavg_interval_end_s = day + time_step.split('-')[1] + '00'
                    tavg_interval_start = datetime.datetime.strptime(tavg_interval_start_s, '%Y%m%d%H%M')
                    try:
                        datetime.datetime.strptime(tavg_interval_end_s, '%Y%m%d%H%M')
                    except ValueError:
                        tavg_interval_end = tavg_interval_start + datetime.timedelta(hours=(int(time_step.split('-')[1]) -
                                                                                            int(time_step.split('-')[0])))
                        tavg_interval_end_s = datetime.datetime.strftime(tavg_interval_end, '%Y%m%d%H%M')
                    day_interval = datetime.datetime.strftime(tavg_interval_start, '%Y%m%d')
                    if day != day_interval:
                        continue
                    file_validity = tavg_interval_start_s + '-' + tavg_interval_end_s + '-tavg'
                    time_step_s = time_step
                file_name = "c_" + file_level_s + "_" + file_validity + ".grd"
                output_folder = os.path.join(model_processed_output_folder, day, specie)
                output_files.append(os.path.join(output_folder, file_name))
            ecdf_output_file = os.path.join(
                ecdf_folder,
                str(ex_prob),
                specie,
                "c_"
                + file_level_s
                + "_"
                + time_step_s
                + ".grd",
            )
            output_quantile = np.zeros((ny, nx))
            c_arrays = []
            files_not_available = []
            for file in output_files:
                try:
                    input_file = open(file)
                except FileNotFoundError:
                    print("File " + file + " not found")
                    files_not_available.append(file)
                    continue
                records = []
                nline = 1
                for line in input_file:
                    if nline > 5:
                        records.append(line.split(" "))
                    nline += 1
                c_arrays.append(records)
            for file in files_not_available:
                output_files.remove(file)
            for j in range(0, ny):
                for i in range(0, nx):
                    c_list = []
                    for k in range(0, len(output_files)):
                        c_list.append(float(c_arrays[k][j][i]))
                    output_quantile[j, i] = np.quantile(c_list, q=quantile)
            # Create header of the processed file
            with open(ecdf_output_file, "a") as processed_file:
                if output_format == "grd":
                    processed_file.write("DSAA\n")
                    processed_file.write(str(nx) + "  " + str(ny) + "\n")
                    processed_file.write(str(x0) + "  " + str(xf) + "\n")
                    processed_file.write(str(y0) + "  " + str(yf) + "\n")
                    processed_file.write(
                        str(np.amin(output_quantile))
                        + "  "
                        + str(np.amax(output_quantile))
                        + "\n"
                    )
                np.savetxt(processed_file, output_quantile, fmt="%.2e")

        if model == "disgas":
            ecdf_folder = disgas_ecdf
            model_processed_output_folder = disgas_processed_output_folder
        else:
            ecdf_folder = twodee_ecdf
            model_processed_output_folder = twodee_processed_output_folder
        for probability in exceedance_probabilities:
            prob_folder = os.path.join(ecdf_folder, str(probability))
            try:
                os.mkdir(prob_folder)
            except FileExistsError:
                print("Folder " + prob_folder + " already exists")
            for specie in species:
                specie_folder = os.path.join(ecdf_folder, prob_folder, specie)
                try:
                    os.mkdir(specie_folder)
                except FileExistsError:
                    print("Folder " + specie_folder + " already exists")
        indexes = []
        pools_ecdfs = []
        indexes_tavg = []
        pools_ecdfs_tavg = []
        n_pool = 0
        n_pool_tavg = 0
        for probability in exceedance_probabilities:
            for specie in species:
                pools_ecdfs.append(n_pool)
                if len(tavg_intervals) > 0:
                    pools_ecdfs_tavg.append(n_pool_tavg)
                    n_pool_tavg += 1
                if levels[0] == "all":
                    for i in range(0, nz):
                        if time_steps[0] == "all":
                            for j in range(0, n_time_steps + 1):
                                indexes.append([probability, specie, processed_files_levels[i], j])
                        else:
                            for time_step in time_steps:
                                indexes.append([probability, specie, processed_files_levels[i], time_step])
                        if len(tavg_intervals) > 0:
                            for k in range(0, len(tavg_intervals)):
                                indexes_tavg.append(
                                    [probability, specie, processed_files_levels[i], tavg_intervals[k]]
                                )
                else:
                    for level in processed_files_levels:
                        if time_steps[0] == "all":
                            for j in range(0, n_time_steps + 1):
                                indexes.append([probability, specie, level, j])
                        else:
                            for time_step in time_steps:
                                indexes.append([probability, specie, level, time_step])
                        if len(tavg_intervals) > 0:
                            for k in range(0, len(tavg_intervals)):
                                indexes_tavg.append(
                                    [probability, specie, level, tavg_intervals[k]]
                                )
                n_pool += 1
        n_pool = 0
        n_completed_processes = 0
        while n_completed_processes <= len(indexes):
            start = n_completed_processes
            end = n_completed_processes + max_number_processes
            if end > len(indexes):
                end = len(indexes)
            pools_ecdfs[n_pool] = ThreadingPool(max_number_processes)
            pools_ecdfs[n_pool].map(ecdf, indexes[start:end])
            n_completed_processes = end
            if n_completed_processes == len(indexes):
                break
        n_pool += 1
        if len(tavg_intervals) > 0:
            n_pool_tavg = 0
            n_completed_processes = 0
            while n_completed_processes <= len(indexes_tavg):
                start = n_completed_processes
                end = n_completed_processes + max_number_processes
                if end > len(indexes_tavg):
                    end = len(indexes_tavg)
                try:
                    pools_ecdfs_tavg[n_pool_tavg] = ThreadingPool(max_number_processes)
                    pools_ecdfs_tavg[n_pool_tavg].map(ecdf, indexes_tavg[start:end])
                except BaseException:
                    print("Unable to elaborate days")
                    sys.exit()
                n_completed_processes = end
                if n_completed_processes == len(indexes_tavg):
                    break
            n_pool_tavg += 1

    def calculate_persistence():
        if model == "disgas":
            persistence_folder = disgas_persistence
            model_processed_output_folder = disgas_processed_output_folder
        else:
            persistence_folder = twodee_persistence
            model_processed_output_folder = twodee_processed_output_folder
        for specie in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie:
                    specie_folder = os.path.join(persistence_folder, specie_dict["specie_name"])
                    try:
                        os.mkdir(specie_folder)
                    except FileExistsError:
                        print("Folder " + specie_folder + " already exists")
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    for threshold in concentration_thresholds:
                        threshold_folder = os.path.join(specie_folder, str(threshold))
                        try:
                            os.mkdir(threshold_folder)
                        except FileExistsError:
                            print("Folder " + threshold_folder + " already exists")

        weight = 1 / len(days)
        #for specie in species:
        #    if levels[0] == "all":
        #        for i in range(0, nz):
        #            for j in range(0, len())
        #    else:
        #        for level in processed_files_levels:



    if plot_ex_prob:
        calculate_quantiles()
    if persistence:
        calculate_persistence()


def save_plots(model, min_con, max_con):
    import re

    def plot_file(input, output, dz_lines_res):
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

        def myround(x, prec=2, base=100):
            return round(base * round(float(x) / base), prec)

        def resize_topography(bottom_left_easting, top_right_easting, bottom_left_northing, top_right_northing,
                              topography):
            import numpy as np
            info_records = []
            with open(topography) as original_topography:
                lines_to_save = []
                for line in original_topography:
                    if len(line.split(' ')) <= 2:
                        info_records.append(line.split(' '))
                        lines_to_save.append(line)
            nx_or_top = float(info_records[1][0])
            ny_or_top = float(info_records[1][1])
            x0_or_top = float(info_records[2][0])
            xf_or_top = float(info_records[2][1])
            y0_or_top = float(info_records[3][0])
            yf_or_top = float(info_records[3][1])
            dx = (xf_or_top - x0_or_top) / (nx_or_top - 1)
            dy = (yf_or_top - y0_or_top) / (ny_or_top - 1)
            i_bottom_left = 0
            j_bottom_left = 0
            i_top_right = 0
            j_top_right = 0
            x = x0_or_top
            y = y0_or_top
            if bottom_left_easting < x0_or_top or bottom_left_northing < y0_or_top or top_right_easting > xf_or_top or \
                    top_right_northing > yf_or_top:
                print('ERROR. Specified domain is not consistent with the topography.grd file')
                print('Topography domain limits')
                print('(X0,XF),(Y0,YF)')
                print('(' + str(x0_or_top) + ',' + str(xf_or_top) + '),(' + str(y0_or_top) + ',' + str(yf_or_top) + ')')
                print('Specified domain limits')
                print('(X0,XF),(Y0,YF)')
                print('(' + str(bottom_left_easting) + ',' + str(top_right_easting) + '),(' + str(bottom_left_northing)
                      + ',' + str(top_right_northing) + ')')
                sys.exit()
            while x <= bottom_left_easting:
                i_bottom_left += 1
                i_top_right += 1
                x += dx
            while x <= top_right_easting:
                i_top_right += 1
                x += dx
            while y <= bottom_left_northing:
                j_bottom_left += 1
                j_top_right += 1
                y += dy
            while y <= top_right_northing:
                j_top_right += 1
                y += dy
            Z_top = np.loadtxt("topography.grd", skiprows=5)
            Z_resized = Z_top[j_bottom_left:j_top_right, i_bottom_left:i_top_right]
            nx_resized = i_top_right - i_bottom_left
            ny_resized = j_top_right - j_bottom_left
            z_min = np.amin(Z_resized)
            z_max = np.amax(Z_resized)
            return nx_resized, ny_resized, z_min, z_max, Z_resized

        if plot_topography_layer:
            nx_top, ny_top, min_z, max_z, Z_top = resize_topography(x0, xf, y0, yf, "topography.grd")
            X_top = np.linspace(x0, xf, num=nx_top)
            Y_top = np.linspace(y0, yf, num=ny_top)
            n_levels = 100
            dz = (max_z - min_z) / n_levels
            if dz_lines_res >= max_z:
                dz_lines_res = max_z
            n_levels_lines = int((max_z - min_z) / dz_lines_res)
            dz_lines = myround((max_z - min_z) / (n_levels_lines))
            levels_top = np.arange(min_z + 0.0000001, max_z, dz)
            levels_top_lines = np.arange(min_z, max_z, dz_lines)
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        specie_name = input.split(os.sep)[-2]
        specie_name = specie_name.translate(SUB)
        with open(input) as input_file:
            if output_format == "grd":
                Z = np.loadtxt(input, skiprows=5)
            X = np.linspace(x0, xf, num=nx)
            Y = np.linspace(y0, yf, num=ny)
            n_levels = 10
            dc = (max_con - min_con) / n_levels
            levels = np.arange(min_con + 0.0000001, max_con, dc)
        fig, ax = plt.subplots(figsize=(6, 5), dpi=plot_resolution)
        if plot_topography_layer:
            top = ax.contourf(
                X_top, Y_top, Z_top, levels_top, cmap="Greys", extend="max"
            )
            top_lines = ax.contour(top, levels=levels_top_lines, colors='black', linewidths=0.05)
            ax.clabel(top_lines, inline=True, fontsize=2, fmt='%1.0f')
            top_cbar = fig.colorbar(
                top, orientation="horizontal", format="%.1f", shrink=0.75
            )
            top_cbar.ax.tick_params(labelsize=6)
            top_cbar.set_label("m a.s.l.")
        c_field = plt.contourf(X, Y, Z, levels, cmap="Reds", alpha=0.9, extend="max")
        if len(plot_isolines) != 0:
            specie_isolines = ax.contour(X, Y, Z, levels=plot_isolines, colors='black', linewidths=0.2)
            ax.clabel(specie_isolines, inline=True, fontsize=4, fmt='%1.1f')
        aspect = 20
        pad_fraction = 0.5
        divider = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1.0 / aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax_c = divider.append_axes("right", size=width, pad=pad)
        cbar = fig.colorbar(c_field, cax=cax_c, orientation="vertical", format="%.1e")
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label("ppm")
        if units == "ppm":
            ax.set_title(specie_name + " concentration [ppm]")
        else:
            ax.set_title(specie_name + " concentration [kg m$\mathregular{^{-3}}$]")
        ax.set_aspect("equal")
        ax.set_xlim(x0, xf)
        ax.set_ylim(y0, yf)
        ax.ticklabel_format(style="plain")
        ax.set_xlabel("X_UTM [m]")
        ax.tick_params(labelsize=8)
        ax.set_ylabel("Y_UTM [m]")
        image_buffer = StringIO()
        plt.tight_layout()
        plt.savefig(output)
        image_buffer.close()
        plt.close(fig)
        input_file.close()

    files_to_plot = []
    output_files = []
    if plot:
        if model == "disgas":
            model_outputs = disgas_outputs
            model_processed_output_folder = disgas_processed_output_folder
            ecdf_outputs = disgas_ecdf
        else:
            model_outputs = twodee_outputs
            model_processed_output_folder = twodee_processed_output_folder
            ecdf_outputs = twodee_ecdf
        graphical_outputs = os.path.join(model_outputs, "graphical_outputs")
        graphical_outputs_simulations = os.path.join(graphical_outputs, "simulations")
        graphical_outputs_ecdf = os.path.join(graphical_outputs, "ecdf")
        for day in days_to_plot:
            graphical_outputs_daily = os.path.join(graphical_outputs_simulations, day)
            try:
                os.mkdir(graphical_outputs_daily)
            except FileExistsError:
                print("Folder " + graphical_outputs_daily + " already exists")
            model_processed_output_folder_daily = os.path.join(
                model_processed_output_folder, day
            )
            model_processed_output_folder_species = []
            for specie in species:
                model_processed_output_folder_species.append(
                    os.path.join(model_processed_output_folder_daily, specie)
                )
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_daily, specie))
                except FileExistsError:
                    print(
                        "Folder "
                        + os.path.join(graphical_outputs_daily, specie)
                        + " already exists"
                    )
            files_list_path = []
            files_list = []
            for folder in model_processed_output_folder_species:
                files_list_temp = os.listdir(folder)
                for file in files_list_temp:
                    if 'TP' in file:
                        continue
                    files_list.append(file)
                    files_list_path.append(os.path.join(folder, file))
            i = 0
            for file in files_list_path:
                file_specie = file.split(model_processed_output_folder_daily)
                file_specie = file_specie[1].split(files_list[i])
                file_specie = re.sub("\W+", "", file_specie[0])
                file_name_splitted = files_list[i].split("_")
                file_level = file_name_splitted[1]
                file_time_step = file_name_splitted[2].split(".")[0]
                try:
                    file_time_step_datetime = datetime.datetime.strptime(file_time_step, '%Y%m%d%H%M')
                except ValueError:
                    file_time_step_datetime = datetime.datetime.strptime('999912310000', '%Y%m%d%H%M')
                simulation_start = datetime.datetime.strptime(day + "{:02d}".format(hour_start), '%Y%m%d%H%M')
                output_file_name = files_list[i].split(".grd")[0]
                output_file_name += ".png"
                if levels[0] == "all":
                    if time_steps[0] == "all":
                        files_to_plot.append(file)
                        output_files.append(
                            os.path.join(
                                graphical_outputs_daily, file_specie, output_file_name
                            )
                        )
                    else:
                        for time_step in time_steps:
                            time_step_hh_mm = hour_start + int(dt / 3600) * (int(time_step))
                            time_step_datetime = simulation_start + datetime.timedelta(hours=time_step_hh_mm)
                            if time_step_datetime == file_time_step_datetime:
                                files_to_plot.append(file)
                                output_files.append(
                                    os.path.join(
                                        graphical_outputs_daily,
                                        file_specie,
                                        output_file_name,
                                    )
                                )
                    if "tavg" in file_time_step:
                        files_to_plot.append(file)
                        tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                        tavg_output_file_name = tavg_output_file_name + ".png"
                        output_files.append(
                            os.path.join(
                                graphical_outputs_daily,
                                file_specie,
                                tavg_output_file_name,
                            )
                        )
                else:
                    if time_steps[0] == "all":
                        for level in levels:
                            if file_level == processed_files_levels[int(level) - 1]:
                                files_to_plot.append(file)
                                output_files.append(
                                    os.path.join(
                                        graphical_outputs_daily,
                                        file_specie,
                                        output_file_name,
                                    )
                                )
                    else:
                        for level in levels:
                            for time_step in time_steps:
                                time_step_hh_mm = hour_start + int(dt / 3600) * (int(time_step))
                                time_step_datetime = simulation_start + datetime.timedelta(hours=time_step_hh_mm)
                                if time_step_datetime == file_time_step_datetime \
                                    and file_level == processed_files_levels[int(level) - 1]:
                                    files_to_plot.append(file)
                                    output_files.append(
                                        os.path.join(
                                            graphical_outputs_daily,
                                            file_specie,
                                            output_file_name,
                                        )
                                    )
                    for level in levels:
                        if "tavg" in file_time_step and file_level == processed_files_levels[int(level) - 1]:
                            files_to_plot.append(file)
                            tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                            tavg_output_file_name = tavg_output_file_name + ".png"
                            output_files.append(
                                os.path.join(
                                    graphical_outputs_daily,
                                    file_specie,
                                    tavg_output_file_name,
                                )
                            )
                i += 1
    if plot_ex_prob:
        for probability in exceedance_probabilities:
            try:
                os.mkdir(os.path.join(graphical_outputs_ecdf, str(probability)))
            except FileExistsError:
                print(
                    "Folder "
                    + os.path.join(graphical_outputs_ecdf, str(probability))
                    + " already exists"
                )
            for specie in species:
                try:
                    os.mkdir(
                        os.path.join(graphical_outputs_ecdf, str(probability), specie)
                    )
                except FileExistsError:
                    print(
                        "Folder "
                        + os.path.join(graphical_outputs_ecdf, str(probability), specie)
                        + " already exists"
                    )
                files_list = os.listdir(
                    os.path.join(ecdf_outputs, str(probability), specie)
                )
                for file in files_list:
                    file_path = os.path.join(
                        ecdf_outputs, str(probability), specie, file
                    )
                    file_name_splitted = file.split("_")
                    file_level = file_name_splitted[1]
                    file_time_step = file_name_splitted[2].split(".")[0]
                    output_file_name = file.split(".grd")[0]
                    output_file_name += ".png"
                    if levels[0] == "all":
                        if time_steps[0] == "all":
                            files_to_plot.append(file_path)
                            output_files.append(
                                os.path.join(
                                    graphical_outputs_ecdf,
                                    str(probability),
                                    specie,
                                    output_file_name,
                                )
                            )
                        else:
                            for time_step in time_steps:
                                if file_time_step == "{:06d}".format(int(time_step)):
                                    files_to_plot.append(file_path)
                                    output_files.append(
                                        os.path.join(
                                            graphical_outputs_ecdf,
                                            str(probability),
                                            specie,
                                            output_file_name,
                                        )
                                    )
                        if "tavg" in file_time_step:
                            files_to_plot.append(file_path)
                            tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                            tavg_output_file_name = tavg_output_file_name + ".png"
                            output_files.append(
                                os.path.join(
                                    graphical_outputs_ecdf,
                                    str(probability),
                                    specie,
                                    tavg_output_file_name,
                                )
                            )
                    else:
                        if time_steps[0] == "all":
                            for level in levels:
                                if file_level == processed_files_levels[int(level) - 1]:
                                    files_to_plot.append(file_path)
                                    output_files.append(
                                        os.path.join(
                                            graphical_outputs_ecdf,
                                            str(probability),
                                            specie,
                                            output_file_name,
                                        )
                                    )
                        else:
                            for level in levels:
                                for time_step in time_steps:
                                    if file_time_step == "{:06d}".format(
                                        int(time_step)
                                    ) and file_level == processed_files_levels[int(level) - 1]:
                                        files_to_plot.append(file_path)
                                        output_files.append(
                                            os.path.join(
                                                graphical_outputs_ecdf,
                                                str(probability),
                                                specie,
                                                output_file_name,
                                            )
                                        )
                            for level in levels:
                                if (
                                    "tavg" in file_time_step
                                    and file_level == processed_files_levels[int(level) - 1]
                                ):
                                    files_to_plot.append(file_path)
                                    tavg_output_file_name = file.split(os.sep)[
                                        -1
                                    ].split(".grd")[0]
                                    tavg_output_file_name = (
                                        tavg_output_file_name + ".png"
                                    )
                                    output_files.append(
                                        os.path.join(
                                            graphical_outputs_ecdf,
                                            str(probability),
                                            specie,
                                            tavg_output_file_name,
                                        )
                                    )
    if len(files_to_plot) == 0:
        print("No files to plot")
    else:
        if min_con == -1.0 and max_con == -1.0:
            max_con = 0
            min_con = 1000000000000000
            for file_to_plot in files_to_plot:
                ZZ = np.loadtxt(file_to_plot, skiprows=5)
                max_c = np.amax(ZZ)
                min_c = np.amin(ZZ)
                if max_c > max_con:
                    max_con = max_c
                if min_c < min_con:
                    min_con = min_c
        i = 0
        for file_to_plot in files_to_plot:
            print("plotting " + file_to_plot)
            plot_file(file_to_plot, output_files[i], dz_lines_res)
            i += 1


root = os.getcwd()

(
    plot,
    plot_ex_prob,
    time_steps,
    levels,
    days_to_plot_in,
    species,
    original_specie,
    exceedance_probabilities,
    max_number_processes,
    convert,
    persistence,
    models,
    units,
    time_av,
    min_con,
    max_con,
    plot_isolines,
    output_format,
    plot_topography_layer,
    dz_lines_res,
    plot_resolution,
    tracking_points
) = read_arguments()


(
    disgas_outputs,
    disgas_original_output_folder,
    disgas_processed_output_folder,
    ecdf_folder_name,
    disgas_ecdf,
    disgas_ecdf_tracking_points,
    disgas_persistence,
    twodee_outputs,
    twodee_original_output_folder,
    twodee_processed_output_folder,
    twodee_ecdf,
    twodee_ecdf_tracking_points,
    twodee_persistence,
    models_to_elaborate,
    twodee_output_time_step,
    model_processed_output_folder,
    ecdf_outputs,
    persistence_outputs,
    graphical_outputs,
    graphical_outputs_simulations,
    graphical_outputs_ecdf,
    graphical_outputs_ecdf_tracking_points
) = folder_structure()

days, days_to_plot = extract_days()

species_properties = gas_properties()

tavg_intervals = []
processed_files_levels = []
processed_files_steps = []
tracking_points_files = []
if tracking_points or plot_ex_prob:
    days_to_elaborate = days
else:
    days_to_elaborate = days_to_plot
for model in models_to_elaborate:
    x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, output_levels, hour_start, minute_start = domain(model)
    if tracking_points:
        stations = elaborate_tracking_points()
        # Initialize array of concentration to be used for ECDFs in the tracking points
        c = [[[[0 for i in range(0, n_time_steps + 10)] for j in range(0, len(days))]
             for k in range(0, len(stations))] for l in range(0, len(species))]
    j = 0
    for day in days_to_elaborate:
        processed_files = elaborate_day(day, model)
        if tracking_points:
            all_time_steps = extract_tracking_points(processed_files, j)
            j += 1
    if tracking_points:
        probabilistic_tracking_points()
    processed_files_levels = sort_levels(processed_files_levels)
    processed_files_steps = sorted(processed_files_steps)
    if plot_ex_prob or persistence:
        probabilistic_output(model)
    save_plots(model, min_con, max_con)