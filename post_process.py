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
from multiprocessing import Pool, cpu_count
import datetime
import linecache


def read_arguments():
    parser = argparse.ArgumentParser(description="Input data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-P",
        "--plot",
        default="False",
        help="Produce plots of the solutions and probabilistic output (if activated). True/False",)
    parser.add_argument(
        "-ECDF",
        "--calculate_ecdf",
        default="False",
        help="Calculate the Empirical Cumulative Density Function of the solution and extrapolate solutions at "
             "user-defined exceedance probabilities. True/False",
    )
    parser.add_argument(
        "-PER",
        "--persistence",
        default="False",
        help="Calculate the persistence of the gas specie, i.e. the probability to be exposed to a gas species"
             " above specified concentration thresholds for times longer than the specified exposure times for those"
             " thresholds.\n" + "Concentration thresholds and exposure times should be provided in "
                                "gas_properties.csv.\n" + " True/False")
    parser.add_argument(
        "-EX",
        "--ex_prob",
        default='',
        help="List of exceedance probabilities to be used for graphical output",
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
        "-TS", "--tracking_specie", default='', help="The original emitted specie that is tracked in the simulation"
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
    calculate_ecdf = args.calculate_ecdf
    persistence = args.persistence
    ex_prob_in = args.ex_prob
    time_steps_in = args.time_steps
    levels_in = args.levels
    days_plot_in = args.days_plot
    species_in = args.species
    original_specie = args.tracking_specie
    nproc = args.nproc
    convert = args.convert
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
    if calculate_ecdf.lower() == "true":
        calculate_ecdf = True
        if ex_prob_in == '':
            print(
                "ERROR. Please specify at least one exceedance probability to plot when --calculate_ecdf==True"
            )
            sys.exit()
    elif calculate_ecdf.lower() == "false":
        calculate_ecdf = False
    else:
        print("ERROR. Wrong value for variable -PE --plot_ex_prob")
        sys.exit()
    if plot:
        if time_steps_in == '':
            print("ERROR. Please specify at least one time step to plot")
            sys.exit()
        if levels_in == '':
            print("ERROR. Please specify at least one level to plot")
            sys.exit()
    if original_specie == '':
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
        calculate_ecdf,
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
    ecdf_folder = "output_ecdf"
    persistence_folder_name = "output_persistence"
    disgas_outputs_folder = os.path.join(root, post_processing, "disgas")
    twodee_outputs_folder = os.path.join(root, post_processing, "twodee")
    disgas_or_output_folder = os.path.join(root, original_output_folder_name, "disgas")
    twodee_or_output_folder = os.path.join(root, original_output_folder_name, "twodee")
    disgas_proc_output_folder = os.path.join(disgas_outputs_folder, processed_output_folder_name)
    twodee_proc_output_folder = os.path.join(twodee_outputs_folder, processed_output_folder_name)
    disgas_ecdf_folder = os.path.join(disgas_outputs_folder, ecdf_folder)
    twodee_ecdf_folder = os.path.join(twodee_outputs_folder, ecdf_folder)
    disgas_ecdf_tracking_points_folder = os.path.join(disgas_ecdf_folder, 'tracking_points')
    twodee_ecdf_tracking_points_folder = os.path.join(twodee_ecdf_folder, 'tracking_points')
    disgas_persistence_folder = os.path.join(disgas_outputs_folder, persistence_folder_name)
    twodee_persistence_folder = os.path.join(twodee_outputs_folder, persistence_folder_name)
    model_proc_output_folder = ''
    ecdf_outputs_folder = ''
    persistence_outputs_folder = ''
    graphical_outputs_folder = ''
    graphical_outputs_simulations_folder = ''
    graphical_outputs_ecdf_folder = ''
    graphical_outputs_ecdf_tracking_points_folder = ''
    graphical_outputs_persistence_folder = ''
    try:
        os.mkdir(post_processing)
    except FileExistsError:
        print("Folder post_processing already exists")
    if models == "disgas" or models == "all":
        try:
            os.mkdir(disgas_outputs_folder)
        except FileExistsError:
            print("Folder " + disgas_outputs_folder + " already exists")
        try:
            os.mkdir(disgas_proc_output_folder)
        except FileExistsError:
            print("Folder " + disgas_proc_output_folder + " already exists")
        try:
            os.mkdir(disgas_ecdf_folder)
        except FileExistsError:
            print("Folder " + disgas_ecdf_folder + " already exists")
        try:
            os.mkdir(disgas_ecdf_tracking_points_folder)
        except FileExistsError:
            print("Folder " + disgas_ecdf_tracking_points_folder + " already exists")
        try:
            os.mkdir(disgas_persistence_folder)
        except FileExistsError:
            print("Folder " + disgas_persistence_folder + " already exists")
    if models == "twodee" or models == "all":
        try:
            os.mkdir(twodee_outputs_folder)
        except FileExistsError:
            print("Folder " + twodee_outputs_folder + " already exists")
        try:
            os.mkdir(twodee_proc_output_folder)
        except FileExistsError:
            print("Folder " + twodee_proc_output_folder + " already exists")
        try:
            os.mkdir(twodee_ecdf_folder)
        except FileExistsError:
            print("Folder " + twodee_ecdf_folder + " already exists")
        try:
            os.mkdir(twodee_ecdf_tracking_points_folder)
        except FileExistsError:
            print("Folder " + twodee_ecdf_tracking_points_folder + " already exists")
        try:
            os.mkdir(twodee_persistence_folder)
        except FileExistsError:
            print("Folder " + twodee_persistence_folder + " already exists")
    twodee_input_file = os.path.join(root, "twodee.inp")
    twodee_output_time_step = 0
    if models == "all":
        models_to_elaborate = ["disgas", "twodee"]
    elif models == "disgas":
        models_to_elaborate = ["disgas"]
    else:
        models_to_elaborate = ["twodee"]
    for model in models_to_elaborate:
        if model == "disgas":
            model_outputs = disgas_outputs_folder
            model_proc_output_folder = disgas_proc_output_folder
            ecdf_outputs_folder = disgas_ecdf_folder
            persistence_outputs_folder = disgas_persistence_folder
        else:
            model_outputs = twodee_outputs_folder
            model_proc_output_folder = twodee_proc_output_folder
            ecdf_outputs_folder = twodee_ecdf_folder
            persistence_outputs_folder = twodee_persistence_folder
            # Read the output time interval from the twodee input file
            with open(twodee_input_file, "r") as twodee_file:
                for line in twodee_file:
                    if "OUTPUT_INTERVAL_(SEC)" in line:
                        twodee_output_time_step = float(line.split("=")[1])
            if twodee_output_time_step == 0:
                print("Unable to read the Twodee output time step")
                sys.exit()
        graphical_outputs_folder = os.path.join(model_outputs, "graphical_outputs")
        graphical_outputs_simulations_folder = os.path.join(graphical_outputs_folder, "simulations")
        graphical_outputs_ecdf_folder = os.path.join(graphical_outputs_folder, "ecdf")
        graphical_outputs_ecdf_tracking_points_folder = os.path.join(graphical_outputs_ecdf_folder, "tracking_points")
        graphical_outputs_persistence_folder = os.path.join(graphical_outputs_folder, "persistence")
        try:
            os.mkdir(graphical_outputs_folder)
        except FileExistsError:
            print("Folder " + graphical_outputs_folder + " already exists")
        try:
            os.mkdir(graphical_outputs_simulations_folder)
        except FileExistsError:
            print("Folder " + graphical_outputs_simulations_folder + " already exists")
        try:
            os.mkdir(graphical_outputs_ecdf_folder)
        except FileExistsError:
            print("Folder " + graphical_outputs_ecdf_folder + " already exists")
        try:
            os.mkdir(graphical_outputs_ecdf_tracking_points_folder)
        except FileExistsError:
            print("Folder " + graphical_outputs_ecdf_tracking_points_folder + " already exists")
        try:
            os.mkdir(graphical_outputs_persistence_folder)
        except FileExistsError:
            print("Folder " + graphical_outputs_persistence_folder + " already exists")

    return (
        disgas_outputs_folder,
        disgas_or_output_folder,
        disgas_proc_output_folder,
        ecdf_folder,
        disgas_ecdf_folder,
        disgas_ecdf_tracking_points_folder,
        disgas_persistence_folder,
        twodee_outputs_folder,
        twodee_or_output_folder,
        twodee_proc_output_folder,
        twodee_ecdf_folder,
        twodee_ecdf_tracking_points_folder,
        twodee_persistence_folder,
        models_to_elaborate,
        twodee_output_time_step,
        model_proc_output_folder,
        ecdf_outputs_folder,
        persistence_outputs_folder,
        graphical_outputs_folder,
        graphical_outputs_simulations_folder,
        graphical_outputs_ecdf_folder,
        graphical_outputs_ecdf_tracking_points_folder,
        graphical_outputs_persistence_folder
    )


def gas_properties():
    def extract_gas_properties(specie):
        global persistence
        data = pd.read_csv(gas_properties_file, on_bad_lines='skip')
        molar_ratio = None
        bg_conc = None
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
                  'persistence calculation not possible for gas specie ' + specie)
        try:
            y = np.sort(data['BG_' + specie])
            bg_conc = list(y)[0]
            if bg_conc != bg_conc:
                print('WARNING. Background concentration of ' + specie + ' not found in gas_properties.csv')
        except KeyError:
            print('WARNING. Background concentration of ' + specie + ' not found in gas_properties.csv')
        for i_exp in range(len(exp_times)):
            if float(exp_times[i_exp]) < simulation_time / n_time_steps / 3600:
                print('WARNING. For specie ' + specie + ' the exposure time ' + str(exp_times[i_exp]) +
                      ' h for the concentration threshold ' + str(conc_thresholds[i_exp]) +
                      ' ppm is less than the simulation output time step ' +
                      str(simulation_time / n_time_steps / 3600) + ' h')
                print('The persistence calculation will assume that the concentration threshold will be overcome for'
                      ' this time shorter than the time step of the simulation')
            elif float(exp_times[i_exp]) > simulation_time / 3600:
                print('WARNING. For specie ' + specie + ' the exposure time ' + str(exp_times[i_exp]) +
                      ' h for the concentration threshold ' + str(conc_thresholds[i_exp]) +
                      ' ppm is greater than the simulation duration ' + str(simulation_time / 3600) + ' h')
                exp_times[i_exp] = 'remove'
                conc_thresholds[i_exp] = 'remove'
        while len(exp_times) > 0:
            try:
                exp_times.remove('remove')
                conc_thresholds.remove('remove')
            except ValueError:
                break
        if len(exp_times) == 0:
            persistence = False
            print('WARNING. Persistence calculation not possible')
        return molar_ratio, molar_weight, conc_thresholds, exp_times, bg_conc

    gas_properties_file = os.path.join(root, "gas_properties.csv")
    try:
        with open(gas_properties_file, "r") as test_file:
            print('File ' + gas_properties_file + ' found!')
    except FileNotFoundError:
        print("File " + gas_properties_file + " not present")
        sys.exit()
    molar_ratios = []
    molar_weights = []
    concentration_thresholds = []
    exposure_times = []
    background_concentrations = []
    for specie in species:
        molar_ratio, molar_weight, concentration_thresholds_specie, exposure_times_specie, background_c = \
            extract_gas_properties(specie)
        molar_ratios.append(molar_ratio)
        molar_weights.append(molar_weight)
        concentration_thresholds.append(concentration_thresholds_specie)
        exposure_times.append(exposure_times_specie)
        background_concentrations.append(background_c)
    molar_ratio, molar_weight, concentration_thresholds_specie, exposure_times_specie, background_c = \
        extract_gas_properties(original_specie)
    molar_ratios_tracking_specie = molar_ratio
    molar_weights_tracking_specie = molar_weight
    concentration_thresholds_tracking_specie = concentration_thresholds_specie
    exposure_times_tracking_specie = exposure_times_specie
    background_concentration_tracking_specie = background_c
    species_properties = []
    for i in range(0, len(species)):
        gas_specie = {}
        gas_specie["specie_name"] = species[i]
        gas_specie["molar_ratio"] = molar_ratios[i]
        gas_specie["molar_weight"] = molar_weights[i]
        gas_specie["concentration_thresholds"] = concentration_thresholds[i]
        gas_specie["exposure_times"] = exposure_times[i]
        gas_specie["background_concentration"] = background_concentrations[i]
        species_properties.append(gas_specie)
    if original_specie not in species:
        gas_specie = {}
        gas_specie["specie_name"] = original_specie
        gas_specie["molar_ratio"] = molar_ratios_tracking_specie
        gas_specie["molar_weight"] = molar_weights_tracking_specie
        gas_specie["concentration_thresholds"] = concentration_thresholds_tracking_specie
        gas_specie["exposure_times"] = exposure_times_tracking_specie
        gas_specie["background_concentration"] = background_concentration_tracking_specie
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
    return x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, tot_time, output_levels, hour_start, minute_start


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


def probabilistic_tracking_points():
    def plot_hazard_curves(file_input, folder, min_con_tp_in, max_con_tp_in):
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
            plt.xlim(min_con_tp_in, max_con_tp_in)
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
    min_con_tp_specie = []
    max_con_tp_specie = []
    for l in range(0, len(species)):
        min_con_tp = 1000000000000000
        max_con_tp = 0
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
        min_con_tp_specie.append(min_con_tp)
        max_con_tp_specie.append(max_con_tp)
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
            plot_file_folder = os.path.join(graphical_outputs_ecdf_tracking_points, species[l])
            plot_hazard_curves(ecdf_tracking_point_file, plot_file_folder, min_con_tp_specie[l], max_con_tp_specie[l])


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
            file_time_step = file_time_step.split(".grd")[0]
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
                        try:
                            background_c = specie["background_concentration"]
                        except UnboundLocalError:
                            background_c = 0
                        try:
                            background_c = float(background_c)
                        except TypeError:
                            background_c = 0
                        if units == 'ppm':
                            species_conversion_factor = molar_ratio
                        elif units == 'kg/m3':
                            for specie in species_properties:
                                if specie["specie_name"] == specie_input:
                                    molar_weight = specie["molar_weight"]
                                if specie["specie_name"] == original_specie:
                                    molar_weight_tracking_specie = specie["molar_weight"]
                            background_c = background_c * ((molar_weight / 22.4) * (273 / t2m) * (p2m / 1013)) / 1000000
                            species_conversion_factor = molar_ratio * (molar_weight / molar_weight_tracking_specie)
                Z_converted = np.where(Z_converted < 0, 0, Z_converted)
                Z_converted += background_c
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


def calculate_persistence():
    for persistence_output_file in persistence_matrices:
        exp_time = persistence_output_file.split('_t_')[1]
        exp_time = float(exp_time.split('H')[0])
        c_thresh = persistence_output_file.split('_t_')[0]
        c_thresh = float(c_thresh.split('C_')[1])
        pers_level = persistence_output_file.split('persistence_')[1]
        pers_level = pers_level.split('.grd')[0]
        specie = persistence_output_file.split('/')[-3]
        for day_ov in overcome_matrices_all_days:
            persistence_parameters_all = overcome_matrices_all_days[day][0]
            for i_par in range(len(persistence_parameters_all)):
                persistence_parameters = persistence_parameters_all[i_par]
                if specie == persistence_parameters[0] and c_thresh == persistence_parameters[1] and \
                        exp_time == persistence_parameters[2] and pers_level == persistence_parameters[3]:
                    overcome_time_data_all = overcome_matrices_all_days[day_ov][1]
                    overcome_time_data = overcome_time_data_all[i_par]
                    for j_p in range(ny):
                        for i_p in range(nx):
                            if overcome_time_data[j_p][i_p] >= exp_time:
                                persistence_matrices[persistence_output_file][j_p][i_p] += weight_simulation


def elaborate_day(day_input):
    def prepare_persistence_calculation():
        def prepare_files(index):
            def calculate_overcome_time(output_files_to_elaborate):
                for j_overcome in range(ny):
                    for i_overcome in range(nx):
                        for file_input in output_files_to_elaborate:
                            try:
                                row = linecache.getline(file_input, j_overcome + 6)
                                if float(row.split(' ')[i_overcome]) > concentration_threshold:
                                    overcome_matrix[j_overcome][i_overcome] += simulation_time / n_time_steps / 3600
                                else:
                                    overcome_matrix[j_overcome][i_overcome] += 0
                            except IndexError:
                                print('File ' + file_input + ' not found')
                                overcome_matrix[j_overcome][i_overcome] += 0
                linecache.clearcache()

            specie_input = index[0]
            concentration_threshold = index[1]
            exposure_time = index[2]
            file_level_s = index[3]
            output_folder = os.path.join(model_processed_output_folder, day_input, specie_input)
            day_output_files = os.listdir(output_folder)
            output_files_day = []
            overcome_matrix = np.zeros((ny, nx))
            for output_file in day_output_files:
                if output_file.split('_')[1] == file_level_s:
                    output_files_day.append(os.path.join(output_folder, output_file))
            calculate_overcome_time(output_files_day)
            pers_output = os.path.join(persistence_folder, specie_input, 'C_' + str(concentration_threshold)
                                                   + '_t_' + str(exposure_time) + 'H',
                                                   'persistence_' + file_level_s + '.grd')
            return pers_output, overcome_matrix

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
                        pass
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    exposure_times = specie_dict["exposure_times"]
                    for j in range(0, len(concentration_thresholds)):
                        threshold_folder = os.path.join(specie_folder, 'C_' + str(concentration_thresholds[j]) + '_t_'
                                                        + str(exposure_times[j]) + 'H')
                        try:
                            os.mkdir(threshold_folder)
                        except FileExistsError:
                            pass

        indexes = []
        for specie in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie:
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    exposure_times = specie_dict["exposure_times"]
                    for i_threshold in range(0, len(concentration_thresholds)):
                        if exposure_times[i_threshold] * 3600 > simulation_time:
                            continue  # Skip exposure times longer than the simulation itself
                        if levels[0] == "all":
                            for j in range(0, nz):
                                indexes.append([specie, concentration_thresholds[i_threshold],
                                                exposure_times[i_threshold], processed_files_levels_elaborated[j]])
                        else:
                            all_levels = np.array(processed_files_levels_elaborated)
                            levels_indexes = [int(x) - 1 for x in levels]
                            for level in list(all_levels[levels_indexes]):
                                indexes.append([specie, concentration_thresholds[i_threshold],
                                                exposure_times[i_threshold], level])
        for i_indexes in range(0, len(indexes)):
            persistence_matrix_temp = np.zeros((ny, nx))
            persistence_output, overcome_matrix = prepare_files(indexes[i_indexes])
            persistence_matrices[persistence_output] = persistence_matrix_temp
            overcome_matrices.append(overcome_matrix)
        return indexes


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

        global tracking_points_files
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
        for l_specie in range(0, len(species)):
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
                c_interpolated_time_steps = []
                for time_step in files_time_steps:
                    files_to_use = []
                    files_levels = []
                    file_directory = ''
                    for file in files_to_interpolate:
                        file_name = file.split(os.sep)[-1]
                        file_specie = file.split(os.sep)[-2] 
                        if file_specie != species[l_specie]:
                            continue
                        file_directory = file.split(file_name)[0]
                        file_level = file_name.split('_')[1]
                        file_level = float(file_level.split('mabg')[0])
                        file_time_step = file_name.split('_')[2]
                        try:
                            file_time_step = int(file_time_step.split('.grd')[0])
                        except ValueError:
                            continue
                        if file_level in levels_for_interpolation and file_time_step == time_step and file_specie == species[l_specie]: 
                            files_to_use.append(file)
                            files_levels.append(file_level)
                    if len(files_to_use) == 0:
                        continue
                    files_levels = sorted(files_levels)
                    c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'],
                                                 files_levels,
                                                 files_to_use)
                    c_interpolated_time_steps.append(c_interpolated[0])
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
                        file_specie = file.split(os.sep)[-2] 
                        if file_specie != species[l_specie]:
                            continue
                        file_directory = file.split(file_name)[0]
                        file_level = file_name.split('_')[1]
                        file_level = float(file_level.split('mabg')[0])
                        file_time_step = file_name.split('_')[2]
                        try:
                            file_time_step = int(file_time_step.split('.grd')[0])
                        except ValueError:
                            file_time_step = file_time_step.split('.grd')[0]
                            if file_level in levels_for_interpolation and file_time_step == time_step and file_specie == species[l_specie]: 
                                files_to_use.append(file)
                                files_levels.append(file_level)
                    if len(files_to_use) == 0:
                        continue
                    files_levels = sorted(files_levels)
                    c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'],
                                                 files_levels,
                                                 files_to_use)
                    c_interpolated_time_steps.append(c_interpolated[0])
                    tracking_point_file = os.path.join(file_directory, 'TP_' + str(station['station_id']) + '.txt')
                    tracking_points_files.append(tracking_point_file)
                    with open(tracking_point_file, 'a') as tracking_point_file:
                        tracking_point_file.write(str(time_step) + '\t' + "{0:.2e}".format(c_interpolated[0]) + '\n')
                    i += 1
                c_tp[species[l_specie]][k] = {'station_id': k, 'c_tp_time_steps': c_interpolated_time_steps}
                k += 1
        return files_time_steps + files_time_averaging_steps, c_tp

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
            files_list_path.append(os.path.join(model_output_folder, file))
            models.append(model)
    for specie in species[1:]:
        files_list_path += files_list_path
        models += models
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
    time_steps_elaborated = []
    overcome_matrices = []
    all_time_steps_tp = []
    indexes_persistence = []
    processed_files_levels_elaborated = []
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
            if file_level_s not in processed_files_levels_elaborated:
                processed_files_levels_elaborated.append(file_level_s)
            time_validity = time_start + datetime.timedelta(seconds=file_time_step)
            file_validity = datetime.datetime.strftime(time_validity, '%Y%m%d%H%M')
            if time_validity not in time_steps_elaborated:
                time_steps_elaborated.append(time_validity)
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
        for file_to_check in processed_files:
            try:
                os.remove(file_to_check)
            except FileNotFoundError:
                pass
    processed_files_levels_elaborated = sorted([float(processed_files_levels_elaborated[i_lev].split('mabg')[0]) for
                                                i_lev in range(len(processed_files_levels_elaborated))])
    processed_files_levels_elaborated = ["{0:.3f}".format(processed_files_levels_elaborated[i_lev]) + 'mabg' for
                                         i_lev in range(len(processed_files_levels_elaborated))]
    n_elaborated_files = 0
    while n_elaborated_files < len(files_list_path):
        start = n_elaborated_files
        end = n_elaborated_files + max_number_processes
        if end > len(files_list_path):
            end = len(files_list_path)
        pool_files = ThreadingPool(max_number_processes)
        pool_files.map(
                converter,
                files_list_path[start:end],
                processed_files[start:end],
                species_list[start:end],
                models[start:end],
        )
        n_elaborated_files = end
        if n_elaborated_files == len(files_list_path):
            break
    if time_av is not None:
        averaged_files = []
        time_min = min(time_steps_elaborated)
        time_steps_elaborated = sorted(time_steps_elaborated)
        time_step_simulation = (time_steps_elaborated[1] - time_steps_elaborated[0]).seconds
        if model == 'twodee':
            time_min -= datetime.timedelta(seconds=time_step_simulation)
        if time_av == 0:
            time_max = max(time_steps_elaborated)
        else:
            time_max = time_min + datetime.timedelta(seconds=time_av * 3600)
            if time_max > max(time_steps_elaborated):
                time_max = max(time_steps_elaborated)
        while time_max <= max(time_steps_elaborated):
            time_diff = time_max - time_min
            if int((time_max - time_min).total_seconds()) <= time_step_simulation:
                print('Warning! Time-averaging interval smaller than or equal to the time step of the simulation')
                tavgs_intervals = []
                break
            time_max_s = datetime.datetime.strftime(time_max, '%H%M')
            if time_max_s == '0000':
                time_max_s = '2400'
            if datetime.datetime.strftime(time_min, '%H%M') + "-" + time_max_s + "-tavg" not in tavg_intervals:
                tavg_intervals.append(datetime.datetime.strftime(time_min, '%H%M') + "-" + time_max_s + "-tavg")
            for i in range(0, len(species)):
                files_to_average = []
                for level in processed_files_levels_elaborated:
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
                        +"-tavg.grd",
                    )
                    try:
                        os.remove(time_averaged_file)
                    except FileNotFoundError:
                        pass
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
                #time_min = time_max + datetime.timedelta(seconds=time_step_simulation)
                time_min = time_max
                if time_min >= max(time_steps_elaborated):
                    break
                time_max = time_min + datetime.timedelta(seconds=time_av * 3600)
                if time_max > max(time_steps_elaborated):
                    time_max = max(time_steps_elaborated)
                continue
    if persistence:
        indexes_persistence = prepare_persistence_calculation()
    if tracking_points:
        all_time_steps_tp, c_tp_time_steps = extract_tracking_points(processed_files, days.index(day_input))
    return day_input, all_time_steps_tp, processed_files_levels_elaborated, tavg_intervals, persistence_matrices, \
           indexes_persistence, overcome_matrices, c_tp


def sort_levels(input_array):
    output_array = []
    for level in input_array:
        level_float = float(level.split('mabg')[0])
        output_array.append(level_float)
    output_array = sorted(output_array)
    for i in range (0, len(output_array)):
        output_array[i] = "{0:.3f}".format(output_array[i]) + 'mabg'
    return output_array


def read_output_files_for_ecdf(ji):
    def read_file(file_input):
        try:
            row = linecache.getline(file_input, j_ecdf + 6)
            c_read = (float(row.split(' ')[i_ecdf]))
            linecache.clearcache()
            return c_read
        except IndexError:
            print('File ' + file_input + ' not found')
            return 0.0, ''

    j_ecdf = ji[0]
    i_ecdf = ji[1]
    #pool_file_read = ThreadingPool(len(files_to_process))
    #c_list_ecdf = pool_file_read.map(read_file, files_to_process)
    c_list_ecdf = []
    for file_to_read in files_to_process:
        c_list_ecdf.append(read_file(file_to_read))
    output_quantile_ecdf = np.quantile(c_list_ecdf, q=1 - probability)
    return j_ecdf, i_ecdf, output_quantile_ecdf


def write_probabilistic_file(file_input, output_to_write):
    # Create header of the processed file
    with open(file_input, "a") as processed_file:
        if output_format == "grd":
            processed_file.write("DSAA\n")
            processed_file.write(str(nx) + "  " + str(ny) + "\n")
            processed_file.write(str(x0) + "  " + str(xf) + "\n")
            processed_file.write(str(y0) + "  " + str(yf) + "\n")
            processed_file.write(
                str(np.amin(output_to_write))
                + "  "
                + str(np.amax(output_to_write))
                + "\n"
            )
        np.savetxt(processed_file, output_to_write, fmt="%.2e")


def prepare_quantile_calculation(exc_prob):
    def prepare_files(index):
        specie = index[1]
        file_level_s = index[2]
        time_step = index[3]
        ex_prob = index[0]
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
                tavg_interval_start_s = day + time_step.split('-')[0]
                tavg_interval_end_s = day + time_step.split('-')[1]
                tavg_interval_start = datetime.datetime.strptime(tavg_interval_start_s, '%Y%m%d%H%M')
                try:
                    tavg_interval_end = datetime.datetime.strptime(tavg_interval_end_s, '%Y%m%d%H%M')
                except ValueError: #in case time_step.split('-')[1] = 2400
                    if(int(time_step.split('-')[0][0:2]) + int(time_step.split('-')[1][0:2])) > 24:
                        tavg_interval_end = tavg_interval_start + datetime.timedelta(hours=int(time_step.split('-')[1][0:2]) - int(time_step.split('-')[0][0:2]))
                    else:
                        tavg_interval_end = tavg_interval_start + datetime.timedelta(hours=int(time_step.split('-')[0][0:2]) + int(time_step.split('-')[1][0:2]))
                if tavg_interval_end < tavg_interval_start:
                    tavg_interval_end = tavg_interval_start + \
                            datetime.timedelta(hours=int(time_step.split('-')[0][0:2]) +
                                int(time_step.split('-')[1][0:2]))
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
        return output_files, ecdf_output_file

    all_output_files = []
    all_ecdf_output_files = []
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
    indexes_tavg = []
    for specie in species:
        if levels[0] == "all":
            for i in range(0, nz):
                if time_steps[0] == "all":
                    for j in range(0, n_time_steps + 1):
                        indexes.append([exc_prob, specie, processed_files_levels[i], j])
                else:
                    for time_step in time_steps:
                        indexes.append([exc_prob, specie, processed_files_levels[i], time_step])
                if len(tavg_intervals) > 0:
                    for k in range(0, len(tavg_intervals)):
                        indexes_tavg.append(
                            [exc_prob, specie, processed_files_levels[i], tavg_intervals[k]]
                        )
        else:
            all_levels = np.array(processed_files_levels)
            levels_indexes = [int(x) - 1 for x in levels]
            for level in list(all_levels[levels_indexes]):
                if time_steps[0] == "all":
                    for j in range(0, n_time_steps + 1):
                        indexes.append([exc_prob, specie, level, j])
                else:
                    for time_step in time_steps:
                        indexes.append([exc_prob, specie, level, time_step])
                if len(tavg_intervals) > 0:
                    for k in range(0, len(tavg_intervals)):
                        indexes_tavg.append([exc_prob, specie, level, tavg_intervals[k]])
    for i_indexes in range(0, len(indexes)):
        a, b = prepare_files(indexes[i_indexes])
        all_output_files.append(a)
        all_ecdf_output_files.append(b)
    if len(tavg_intervals) > 0:
        for i_indexes in range(0, len(indexes_tavg)):
            a, b = prepare_files(indexes_tavg[i_indexes])
            all_output_files.append(a)
            all_ecdf_output_files.append(b)

    return all_output_files, all_ecdf_output_files


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

        try:
            os.remove(output)
        except FileNotFoundError:
            pass
        if plot_topography_layer:
            nx_top, ny_top, min_z, max_z, Z_top = resize_topography(x0, xf, y0, yf, "topography.grd")
            X_top = np.linspace(x0, xf, num=nx_top)
            Y_top = np.linspace(y0, yf, num=ny_top)
            n_levels = 100
            dz = (max_z - min_z) / n_levels
            if dz_lines_res >= max_z:
                dz_lines_res = max_z
            n_levels_lines = int((max_z - min_z) / dz_lines_res)
            dz_lines = myround((max_z - min_z) / (n_levels_lines), base=dz_lines_res)
            levels_top = np.arange(min_z + 0.0000001, max_z, dz)
            levels_top_lines = np.arange(myround(min_z, base=dz_lines_res), myround(max_z, base=dz_lines_res), dz_lines)
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        if 'persistence' in input.split(os.sep)[-1]:
            specie_name = input.split(os.sep)[-3]
            thresholds = input.split(os.sep)[-2]
            c_threshold = thresholds.split('_t_')[0]
            c_threshold = c_threshold.split('C_')[1]
            exposure_time = thresholds.split('_t_')[1]
            exposure_time = exposure_time.split('H')[0]
        else:
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
            if 'persistence' in input:
                n_levels = 100
                dp = 1. / n_levels
                levels = np.arange(0.00000000001, 1, dp)
                if output_format == "grd":
                    Z = np.loadtxt(input, skiprows=5)
        fig, ax = plt.subplots(figsize=(6, 5), dpi=plot_resolution)
        if plot_topography_layer:
            top = ax.contourf(X_top, Y_top, Z_top, levels_top, cmap="Greys", extend="max")
            top_lines = ax.contour(top, levels=levels_top_lines, colors='black', linewidths=0.05)
            top_cbar = fig.colorbar(top, orientation="horizontal", format="%.1f", shrink=0.75)
            ax.clabel(top_lines, inline=True, fontsize=2, fmt='%1.0f')
            top_cbar.ax.tick_params(labelsize=6)
            top_cbar.set_label("m a.s.l.")
        if 'persistence' in input.split(os.sep)[-1]:
            cmap = plt.get_cmap('viridis', 10)
            field = plt.contourf(X, Y, Z, levels, cmap=cmap, alpha=0.9, extend="max")
        else:
            field = plt.contourf(X, Y, Z, levels, cmap="Reds", alpha=0.9, extend="max")
        if len(plot_isolines) != 0:
            specie_isolines = ax.contour(X, Y, Z, levels=plot_isolines, colors='black', linewidths=0.2)
            ax.clabel(specie_isolines, inline=True, fontsize=4, fmt='%1.1f')
        aspect = 20
        pad_fraction = 0.5
        divider = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1.0 / aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax_c = divider.append_axes("right", size=width, pad=pad)
        if 'persistence' in input.split(os.sep)[-1]:
            cbar = fig.colorbar(field, cax=cax_c, orientation="vertical", format="%.1f")
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label("Probability")
            if units == 'ppm':
                ax.set_title(specie_name + " persistence C > " + c_threshold + ' ppm for t = ' + exposure_time
                             + ' hours')
            else:
                ax.set_title(specie_name + " persistence C > " + c_threshold + ' kg m$\mathregular{^{-3}}$ for t = '
                             + exposure_time + ' hours')
        else:
            cbar = fig.colorbar(field, cax=cax_c, orientation="vertical", format="%.1e")
            cbar.ax.tick_params(labelsize=8)
            if units == "ppm":
                cbar.set_label("C [ppm]")
                ax.set_title(specie_name + " concentration [ppm]")
            else:
                cbar.set_label("C [kg m$\mathregular{^{-3}}$]")
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
                        time_step_seconds = hour_start + dt * int(time_step)
                        time_step_datetime = simulation_start + datetime.timedelta(seconds=time_step_seconds)
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
                            time_step_seconds = hour_start + dt * int(time_step)
                            time_step_datetime = simulation_start + datetime.timedelta(seconds=time_step_seconds)
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

    if calculate_ecdf:
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
    if persistence:
        for specie in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie:
                    try:
                        os.mkdir(
                            os.path.join(graphical_outputs_persistence, specie)
                        )
                    except FileExistsError:
                        print(
                            "Folder "
                            + os.path.join(graphical_outputs_persistence, specie)
                            + " already exists"
                        )
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    exposure_times = specie_dict["exposure_times"]
                    for j in range(0, len(concentration_thresholds)):
                        try:
                            os.mkdir(
                                os.path.join(graphical_outputs_persistence, specie, 'C_' +
                                             str(concentration_thresholds[j]) + '_t_' + str(exposure_times[j]) + 'H'))
                        except FileExistsError:
                            print(
                                "Folder "
                                + os.path.join(graphical_outputs_persistence, specie, 'C_' +
                                               str(concentration_thresholds[j]) + '_t_' + str(exposure_times[j]) + 'H')
                                               + " already exists")
                        files_list = os.listdir(os.path.join(persistence_outputs, specie, 'C_' +
                                                             str(concentration_thresholds[j]) + '_t_'
                                                        + str(exposure_times[j]) + 'H'))
                        for file in files_list:
                            file_path = os.path.join(
                                persistence_outputs, specie, 'C_' + str(concentration_thresholds[j]) + '_t_'
                                                        + str(exposure_times[j]) + 'H', file)
                            files_to_plot.append(file_path)
                            persistence_plot_file_name = file.split(os.sep)[-1].split(".grd")[0]
                            persistence_plot_file_name += '.png'
                            output_files.append(os.path.join(graphical_outputs_persistence, specie, 'C_' +
                                                             str(concentration_thresholds[j]) + '_t_' +
                                                             str(exposure_times[j]) + 'H',
                                                             persistence_plot_file_name))
    if len(files_to_plot) == 0:
        print("No files to plot")
    else:
        if min_con == -1.0 and max_con == -1.0:
            for specie in species:
                files_to_plot_specie = []
                output_files_specie = []
                max_con = 0
                min_con = 1000000000000000
                i = 0
                for file_to_plot in files_to_plot:
                    if specie in file_to_plot:
                        files_to_plot_specie.append(file_to_plot)
                        output_files_specie.append(output_files[i])
                        if 'persistence' not in file_to_plot:
                            ZZ = np.loadtxt(file_to_plot, skiprows=5)
                            max_c = np.amax(ZZ)
                            min_c = np.amin(ZZ)
                            if max_c > max_con:
                                max_con = max_c
                            if min_c < min_con:
                                min_con = min_c
                    i += 1
                i = 0
                for file_to_plot in files_to_plot_specie:
                    print("plotting " + file_to_plot)
                    plot_file(file_to_plot, output_files_specie[i], dz_lines_res)
                    i += 1
        else:
            i = 0
            for file_to_plot in files_to_plot:
                print("plotting " + file_to_plot)
                plot_file(file_to_plot, output_files[i], dz_lines_res)
                i += 1


root = os.getcwd()

(
    plot,
    calculate_ecdf,
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
    graphical_outputs_ecdf_tracking_points,
    graphical_outputs_persistence
) = folder_structure()

if __name__ == '__main__':
    days, days_to_plot = extract_days()
    if max_number_processes > cpu_count():
        print('WARNING. Number of requested simulataneous processes larger than the available ' + str(cpu_count())
              + ' cores')
    tavg_intervals = []
    tracking_points_files = []
    persistence_matrices = {}
    overcome_processed_output = []
    overcome_matrices_all_days = {}
    c_tp = {}
    for specie in species:
        c_tp[specie] = {}
    for day in days:
        overcome_matrices_all_days[day] = ''
    if tracking_points or calculate_ecdf or persistence:
        days_to_elaborate = days
    else:
        days_to_elaborate = days_to_plot
    for model in models_to_elaborate:
        x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, simulation_time, output_levels, hour_start, \
            minute_start = domain(model)
        species_properties = gas_properties()
        if tracking_points:
            stations = elaborate_tracking_points()
            # Initialize array of concentration to be used for ECDFs in the tracking points
            c = [[[[0 for i in range(0, n_time_steps + 10)] for j in range(0, len(days))]
                 for k in range(0, len(stations))] for l in range(0, len(species))]
        n_completed_processes = 0
        returned_values = []
        while n_completed_processes <= len(days):
            start = n_completed_processes
            end = n_completed_processes + max_number_processes
            if end > len(days):
                end = len(days)
            pool = Pool(max_number_processes)
            returned_values_temp = pool.map(elaborate_day, days[start:end])
            for returned_value_temp in returned_values_temp:
                returned_values.append(returned_value_temp)
            pool.close()
            pool.join()
            n_completed_processes = end
            if n_completed_processes == len(days):
                break
        for returned_value in returned_values:
            all_time_steps = returned_value[1]
            processed_files_levels = returned_value[2]
            tavg_intervals = returned_value[3]
            persistence_matrices = returned_value[4]
            day_overcome_calculation = returned_value[0]
            persistence_calculation_parameters = returned_value[5]
            overcome_outputs = returned_value[6]
            c_tp = returned_value[7]
            if tracking_points:
                j_tp = days.index(day_overcome_calculation)
                for l_tp in range(len(c_tp)):
                    for k_tp in range(len(c_tp[species[l_tp]])):
                        for i_tp in range(len(c_tp[species[l_tp]][k_tp]['c_tp_time_steps'])):
                            c[l_tp][k_tp][j_tp][i_tp] = c_tp[species[l_tp]][k_tp]['c_tp_time_steps'][i_tp]
            if persistence:
                weight_simulation = 1 / len(days)
                overcome_matrices_all_days[day_overcome_calculation] = [persistence_calculation_parameters,
                                                                        overcome_outputs]
        if persistence:
            calculate_persistence()
            for persistence_output_file in persistence_matrices:
                write_probabilistic_file(persistence_output_file, persistence_matrices[persistence_output_file])
        if tracking_points:
            probabilistic_tracking_points()
        if calculate_ecdf:
            jis = [(j, i) for j in range(ny) for i in range(nx)]
            for probability in exceedance_probabilities:
                all_output_files, all_ecdf_output_files = prepare_quantile_calculation(probability)
                for ii in range(0, len(all_ecdf_output_files)):
                    output_quantile = np.zeros((ny, nx))
                    files_to_process = all_output_files[ii]
                    n_completed_processes = 0
                    while n_completed_processes <= len(jis):
                        start = n_completed_processes
                        end = n_completed_processes + max_number_processes
                        if end > len(jis):
                            end = len(jis)
                        pool_ecdf = Pool(end - start)
                        output_quantile_return = pool_ecdf.map(read_output_files_for_ecdf, jis[start:end])
                        pool_ecdf.close()
                        pool_ecdf.join()
                        for i_output in range(len(output_quantile_return)):
                            jj_ecdf = output_quantile_return[i_output][0]
                            ii_ecdf = output_quantile_return[i_output][1]
                            output_quantile[jj_ecdf, ii_ecdf] = output_quantile_return[i_output][2]
                        n_completed_processes = end
                        if n_completed_processes == len(jis):
                            break
                    try:
                        os.remove(all_ecdf_output_files[ii])
                    except FileNotFoundError:
                        pass
                    write_probabilistic_file(all_ecdf_output_files[ii], output_quantile)
        if plot:
            save_plots(model, min_con, max_con)
