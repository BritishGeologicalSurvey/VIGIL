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
    parser.add_argument("-P", "--plot", default="False",
                        help="Produce plots of the solutions and probabilistic output (if activated). True/False",)
    parser.add_argument("-ECDF", "--calculate_ecdf", default="False",
                        help="Calculate the Empirical Cumulative Density Function of the solution and extrapolate "
                             "solutions at user-defined exceedance probabilities. True/False",)
    parser.add_argument("-PER", "--persistence", default="False", help="Calculate the persistence of the gas specie, "
                        "i.e. the probability to be exposed to a gas species above specified concentration thresholds "
                        "for times longer than the specified exposure times for those thresholds.\n Concentration "
                        "thresholds and exposure times should be provided in gas_properties.csv.\n" + " True/False")
    parser.add_argument("-EX", "--ex_prob", default='', help="List of exceedance probabilities to be used for "
                                                             "graphical output",)
    parser.add_argument("-T", "--time_steps", default='', help="List of time steps to plot (integer >= 0). "
                                                               "Type all to plot all the time steps",)
    parser.add_argument("-L", "--levels", default='', help="List of vertical levels (integer >= 1) to plot. Type all "
                                                           "to plot all the levels",)
    parser.add_argument("-D", "--days_plot", default='', help="List of days to plot (YYYYMMDD). "
                                                              "Type all to plot all the days",)
    parser.add_argument("-C", "--convert", default="False", help="If True, convert output concentration into other "
                                                                 "species listed with the command -S (--species)",)
    parser.add_argument("-S", "--species", default='', help="List of gas species (e.g. CO2)")
    parser.add_argument("-TS", "--tracking_specie", default=None, help="The original emitted specie that is tracked "
                                                                       "in the simulation")
    parser.add_argument("-N", "--nproc", default=1, help="Maximum number of allowed simultaneous processes",)
    parser.add_argument("-U", "--units", default=None, help="Gas concentration units. "
                                                            "Possible options are: ppm, kg/m3",)
    parser.add_argument("-PL", "--plot_limits", default='', help="Minimum and maximum value of concentration to "
                        "display. If unspecified, they are obtained from all the outputs",)
    parser.add_argument("-PI", "--plot_isolines", default='', help="List of gas concentrations values to be used to "
                                                                   "draw isolines. Optional")
    parser.add_argument("-TA", "--time_av", default=None, help="Generate time-averaged outputs. Specify the "
                        "time-averaging interval (in hours), or 0 for averaging over the whole duration",)
    parser.add_argument("-OF", "--output_format", default="GRD", help="Select format of the processed output files. "
                                                                      "Valid options are: GRD",)
    parser.add_argument("-PT", "--plot_topography", default="False", help="Plot topography layer (True or False). "
                                                                          "Warning, it can be time-consuming!",)
    parser.add_argument("-TI", "--topography_isolines", default=100, help="Topography height contour lines spatial "
                        "resolution (in m a.s.l.). Used only if -PT True",)
    parser.add_argument("-PR", "--plot_resolution", default=600, help="Specify plot resolution in dpi")
    parser.add_argument("-TP", "--tracking_points", default="False", help="Extrapolate gas concentration at locations "
                        "specified in the file tracking_points.txt")
    args = parser.parse_args()
    plot_in = args.plot
    calculate_ecdf_in = args.calculate_ecdf
    persistence_in = args.persistence
    ex_prob_in = args.ex_prob
    time_steps_in = args.time_steps
    levels_in = args.levels
    days_plot_in = args.days_plot
    species_in = args.species
    original_specie_in = args.tracking_specie
    nproc = args.nproc
    convert_in = args.convert
    units_in = args.units
    plot_limits_in = args.plot_limits
    plot_isolines_in = args.plot_isolines
    time_av_in = args.time_av
    output_format_in = args.output_format
    plot_topography = args.plot_topography
    dz_lines_res_in = args.topography_isolines
    plot_resolution_in = args.plot_resolution
    tracking_points_in = args.tracking_points
    ex_prob = ex_prob_in.split(',')
    time_steps_in = time_steps_in.split(',')
    levels_in = levels_in.split(',')
    days_plot = days_plot_in.split(',')
    plot_limits = plot_limits_in.split(',')
    plot_isolines_s = plot_isolines_in.split(',')
    plot_isolines_in = []
    species_in = species_in.split(',')
    days_to_plot_in = []
    try:
        max_number_processes_in = int(os.environ["SLURM_NTASKS"])
    except (ValueError, KeyError):
        try:
            max_number_processes_in = int(nproc)
        except ValueError:
            print("Please provide a valid number for the maximum number of process")
            sys.exit()
    if plot_in.lower() == "true":
        plot_in = True
        if days_plot_in == '':
            print("ERROR. Please specify at least one day to plot when --plot==True")
            sys.exit()
        else:
            for day_to_plot in days_plot:
                if day_to_plot == 'all':
                    days_to_plot_in.append(day_to_plot)
                else:
                    try:
                        day_datetime = datetime.datetime.strptime(day_to_plot, '%Y%m%d')
                        days_to_plot_in.append(day_datetime.strftime('%Y%m%d'))
                    except ValueError:
                        print('ERROR. Wrong format for -D -days_plot')
                        sys.exit()
    elif plot_in.lower() == "false":
        plot_in = False
    else:
        print("ERROR. Wrong value for variable -P --plot")
        sys.exit()
    if calculate_ecdf_in.lower() == "true":
        calculate_ecdf_in = True
        if ex_prob_in == '':
            print(
                "ERROR. Please specify at least one exceedance probability to plot when --calculate_ecdf==True"
            )
            sys.exit()
    elif calculate_ecdf_in.lower() == "false":
        calculate_ecdf_in = False
    else:
        print("ERROR. Wrong value for variable -PE --plot_ex_prob")
        sys.exit()
    if plot_in:
        if time_steps_in == '':
            print("ERROR. Please specify at least one time step to plot")
            sys.exit()
        if levels_in == '':
            print("ERROR. Please specify at least one level to plot")
            sys.exit()
    if original_specie_in is None:
        print('ERROR. Please specify the name of the tracked specie')
        sys.exit()
    if species_in == '':
        print("ERROR. Please specify at least one gas specie name")
        sys.exit()
    if convert_in.lower() == "true":
        convert_in = True
    elif convert_in.lower() == "false":
        convert_in = False
    else:
        print("ERROR. Wrong value for variable -C --convert")
        sys.exit()
    if persistence_in.lower() == "true":
        persistence_in = True
    elif persistence_in.lower() == "false":
        persistence_in = False
    else:
        print("ERROR. Wrong value for variable -PER --persistence")
        sys.exit()
    exceedance_probabilities_in = []
    if ex_prob_in != '':
        for prob in ex_prob:
            exceedance_probabilities_in.append(float(prob))
    try:
        units_in = units_in.lower()
    except AttributeError:
        print("Please provide an option for -U --units")
        sys.exit()
    if units_in != "ppm" and units_in != "kg/m3":
        print("ERROR. Wrong value for variable -U --units")
        sys.exit()
    try:
        time_av_in = int(time_av_in)
    except TypeError:
        time_av_in = None
    except ValueError:
        try:
            time_av_in = int(float(time_av_in))
        except ValueError:
            print("ERROR. Please specify a valid time-averaging interval")
            sys.exit()
    min_con_in = max_con_in = -1.0
    if len(plot_limits) > 1:
        try:
            min_con_in = float(plot_limits[0])
            max_con_in = float(plot_limits[1])
        except ValueError:
            print(
                "ERROR. Please specify valid minimum and maximum concentration -PL --plot_limits"
            )
            sys.exit()
    if len(plot_isolines_s) >= 1:
        for isoline in plot_isolines_s:
            if isoline != '':
                try:
                    plot_isolines_in.append(float(isoline))
                except ValueError:
                    print("WARNING. Wrong entry for -PI --plot_isolines. Continuing discarding concentration contour "
                          "lines")
                    plot_isolines_in = []
                    break
    if output_format_in.lower() != "grd":
        print(
            "ERROR. Please specify a valid output format. Current valid options are: GRD"
        )
        sys.exit()
    else:
        output_format_in = "grd"
    if plot_topography.lower() == "true":
        plot_topography_layer_in = True
    elif plot_topography.lower() == "false":
        plot_topography_layer_in = False
    else:
        print("ERROR. Wrong value for variable -PT --plot_topography")
        sys.exit()
    if plot_topography_layer_in:
        try:
            dz_lines_res_in = float(dz_lines_res_in)
        except ValueError:
            dz_lines_res_in = 100
    try:
        plot_resolution_in = int(plot_resolution_in)
    except ValueError:
        print("ERROR. Please provide a valid number for -PR --plot_resolution")
        sys.exit()
    if tracking_points_in.lower() == 'true':
        tracking_points_in = True
        if not os.path.isfile('tracking_points.txt'):
            print("WARNING. Tracking points option activated but file tracking_points.txt not found. Continuing "
                  "without this option")
            tracking_points_in = False
    elif tracking_points_in.lower() == 'false':
        tracking_points_in = False
    else:
        print("ERROR. Wrong value for variable -TP --tracking_points")
        sys.exit()
    return (
        plot_in,
        calculate_ecdf_in,
        time_steps_in,
        levels_in,
        days_to_plot_in,
        species_in,
        original_specie_in,
        exceedance_probabilities_in,
        max_number_processes_in,
        convert_in,
        persistence_in,
        units_in,
        time_av_in,
        min_con_in,
        max_con_in,
        plot_isolines_in,
        output_format_in,
        plot_topography_layer_in,
        dz_lines_res_in,
        plot_resolution_in,
        tracking_points_in
    )


def folder_structure():
    outputs_dir = os.path.join(root, "post_processing")
    original_output_dir = os.path.join(root, "simulations", "runs")
    processed_output_dir = os.path.join(outputs_dir, "runs_processed")
    ecdf_dir = os.path.join(outputs_dir, "output_ecdf")
    ecdf_tracking_points_dir = os.path.join(ecdf_dir, 'tracking_points')
    persistence_dir = os.path.join(outputs_dir, "output_persistence")
    try:
        os.mkdir("post_processing")
    except FileExistsError:
        print("Folder post_processing already exists")
    try:
        os.mkdir(outputs_dir)
    except FileExistsError:
        print("Folder " + outputs_dir + " already exists")
    try:
        os.mkdir(processed_output_dir)
    except FileExistsError:
        print("Folder " + processed_output_dir + " already exists")
    try:
        os.mkdir(ecdf_dir)
    except FileExistsError:
        print("Folder " + ecdf_dir + " already exists")
    try:
        os.mkdir(ecdf_tracking_points_dir)
    except FileExistsError:
        print("Folder " + ecdf_tracking_points_dir + " already exists")
    try:
        os.mkdir(persistence_dir)
    except FileExistsError:
        print("Folder " + persistence_dir + " already exists")
    graphical_outputs_dir = os.path.join(outputs_dir, "graphical_outputs")
    graphical_outputs_simulations_dir = os.path.join(graphical_outputs_dir, "simulations")
    graphical_outputs_ecdf_dir = os.path.join(graphical_outputs_dir, "ecdf")
    graphical_outputs_ecdf_tracking_points_dir = os.path.join(graphical_outputs_ecdf_dir, "tracking_points")
    graphical_outputs_persistence_dir = os.path.join(graphical_outputs_dir, "persistence")
    try:
        os.mkdir(graphical_outputs_dir)
    except FileExistsError:
        print("Folder " + graphical_outputs_dir + " already exists")
    try:
        os.mkdir(graphical_outputs_simulations_dir)
    except FileExistsError:
        print("Folder " + graphical_outputs_simulations_dir + " already exists")
    try:
        os.mkdir(graphical_outputs_ecdf_dir)
    except FileExistsError:
        print("Folder " + graphical_outputs_ecdf_dir + " already exists")
    try:
        os.mkdir(graphical_outputs_ecdf_tracking_points_dir)
    except FileExistsError:
        print("Folder " + graphical_outputs_ecdf_tracking_points_dir + " already exists")
    try:
        os.mkdir(graphical_outputs_persistence_dir)
    except FileExistsError:
        print("Folder " + graphical_outputs_persistence_dir + " already exists")
    return (
        outputs_dir,
        original_output_dir,
        processed_output_dir,
        ecdf_dir,
        ecdf_tracking_points_dir,
        persistence_dir,
        graphical_outputs_dir,
        graphical_outputs_simulations_dir,
        graphical_outputs_ecdf_dir,
        graphical_outputs_ecdf_tracking_points_dir,
        graphical_outputs_persistence_dir
    )


def extract_days():
    days_formatted = []
    days_list_path = os.path.join(root, "days_list.txt")
    days_to_plot_updated = []
    with open(days_list_path, "r") as days_list:
        for line in days_list:
            day_temp = line.split(" ")[0]
            day_temp = day_temp.split("-")
            day_formatted = day_temp[0] + day_temp[1] + day_temp[2]
            days_formatted.append(day_formatted)
            try:
                if days_to_plot[0] == "all":
                    days_to_plot_updated.append(day_formatted)
                else:
                    for day_to_plot in days_to_plot:
                        if day_to_plot == day_formatted:
                            days_to_plot_updated.append(day_to_plot)
            except IndexError:
                print('No days to plot')
    return days_formatted, days_to_plot_updated


def domain():
    import glob
    simulations_folder = os.path.join(os.getcwd(), 'simulations')
    inp_files = glob.glob(simulations_folder + "/**/*.inp", recursive=True)
    disgas_inputs = False
    twodee_inputs = False
    inp_file = ''
    for input_file in inp_files:
        if 'disgas' in input_file:
            disgas_inputs = True
        if 'twodee' in input_file:
            twodee_inputs = True
    for input_file in inp_files:
        if twodee_inputs:
            if disgas_inputs:
                if 'disgas' in input_file:
                    inp_file = input_file
                    break
            else:
                if 'twodee' in input_file:
                    inp_file = input_file
                    break
        else:
            if 'disgas' in input_file:
                inp_file = input_file
                break
    output_levels_inp = []
    with open(inp_file, 'r') as input_file:
        for record in input_file:
            try:
                record_splitted = record.split("=")
                temp = record_splitted[1].split("(")
                if "SIMULATION_INTERVAL_(SEC)" in record_splitted[0]:
                    tot_time_inp = int(temp[0])
                elif "NX" in record_splitted[0]:
                    nx_inp = int(temp[0])
                elif "NY" in record_splitted[0]:
                    ny_inp = int(temp[0])
                elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                    dt_inp = int(temp[0])
                elif "DX_(M)" in record_splitted[0]:
                    dx_inp = float(temp[0])
                elif "DY_(M)" in record_splitted[0]:
                    dy_inp = float(temp[0])
                elif "X_ORIGIN_(UTM_M)" in record_splitted[0]:
                    x0_inp = float(temp[0])
                elif "Y_ORIGIN_(UTM_M)" in record_splitted[0]:
                    y0_inp = float(temp[0])
                elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                    dt_inp = float(temp[0])
                elif "HOUR" in record_splitted[0]:
                    hour_start_inp = int(temp[0])
                elif "MINUTE" in record_splitted[0]:
                    minute_start_inp = int(temp[0])
                elif "Z_LAYERS_(M)" in record_splitted[0] or "HEIGHTS_(M)" in record_splitted[0]:
                    heights = temp[0]
                    heights_list = heights.split(' ')
                    for height in heights_list:
                        try:
                            output_levels_inp.append(float(height))
                        except ValueError:
                            continue
                    output_levels_inp = sorted(output_levels_inp)
            except (IndexError, ValueError):
                continue
    yf_inp = y0_inp + (ny_inp - 1) * dy_inp
    xf_inp = x0_inp + (nx_inp - 1) * dx_inp
    n_time_steps_inp = int(tot_time_inp / dt_inp)
    nz_inp = len(output_levels_inp)
    return x0_inp, xf_inp, y0_inp, yf_inp, nx_inp, ny_inp, nz_inp, dx_inp, dy_inp, n_time_steps_inp, dt_inp, \
        tot_time_inp, output_levels_inp, hour_start_inp, minute_start_inp


def gas_properties():
    def extract_gas_properties(specie_in):
        global persistence
        data = pd.read_csv(gas_properties_file, on_bad_lines='skip')
        molar_ratio = None
        bg_conc = None
        conc_thresholds = []
        exp_times = []
        if convert:
            if specie_in == original_specie:
                molar_ratio = 1
            else:
                try:
                    x = np.sort(data[specie_in + '/' + original_specie])
                    list_x = list(x)
                    samples = random.sample(list_x, 1)
                    molar_ratio = samples[0]
                except KeyError:
                    print('ERROR. Molar ratio ' + specie_in + '/' + original_specie +
                          ' not found in gas_properties.csv')
                    exit()
        try:
            y = data['M_' + specie]
            molar_weight = list(y)[0]
            if molar_weight != molar_weight:
                print('ERROR. Molar weight of ' + specie_in + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Molar weight of ' + specie_in + ' not found in gas_properties.csv')
            sys.exit()
        try:
            y = data['CT_' + specie_in]
            conc_thresholds_temp = list(y)
            conc_thresholds = [x for x in conc_thresholds_temp if x == x]
            if len(conc_thresholds) == 0:
                print('WARNING. Concentration thresholds of ' + specie_in + ' not found in gas_properties.csv. Gas '
                'persistence calculation not possible')
        except KeyError:
            print('WARNING. Concentration thresholds of ' + specie_in + ' not found in gas_properties.csv. Gas '
            'persistence calculation not possible')
        try:
            y = data['ET_' + specie_in]
            exp_times_temp = list(y)
            exp_times = [x for x in exp_times_temp if x == x]
            if len(exp_times) == 0:
                print('WARNING. Exposure times of ' + specie_in + ' not found in gas_properties.csv. Gas '
                'persistence calculation not possible for gas specie ' + specie_in)
        except KeyError:
            print('WARNING. Exposure times of ' + specie_in + ' not found in gas_properties.csv. Gas '
            'persistence calculation not possible for gas specie ' + specie_in)
        try:
            y = np.sort(data['BG_' + specie_in])
            bg_conc = list(y)[0]
            if bg_conc != bg_conc:
                print('WARNING. Background concentration of ' + specie_in + ' not found in gas_properties.csv')
        except KeyError:
            print('WARNING. Background concentration of ' + specie_in + ' not found in gas_properties.csv')
        for i_exp in range(len(exp_times)):
            if float(exp_times[i_exp]) < simulation_time / n_time_steps / 3600:
                print('WARNING. For specie ' + specie_in + ' the exposure time ' + str(exp_times[i_exp]) +
                      ' h for the concentration threshold ' + str(conc_thresholds[i_exp]) +
                      ' ppm is less than the simulation output time step ' +
                      str(simulation_time / n_time_steps / 3600) + ' h')
                print('The persistence calculation will assume that the concentration threshold will be overcome for'
                      ' this time shorter than the time step of the simulation')
            elif float(exp_times[i_exp]) > simulation_time / 3600:
                print('WARNING. For specie ' + specie_in + ' the exposure time ' + str(exp_times[i_exp]) +
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
            print('WARNING. Persistence calculation not possible for specie ' + specie)
            persistence_sp = False
        else:
            persistence_sp = True
        return molar_ratio, molar_weight, conc_thresholds, exp_times, bg_conc, persistence_sp

    gas_properties_file = os.path.join(root, "gas_properties.csv")
    try:
        with open(gas_properties_file, "r"):
            print('File ' + gas_properties_file + ' found!')
    except FileNotFoundError:
        print("File " + gas_properties_file + " not present")
        sys.exit()
    molar_ratios = []
    molar_weights = []
    concentration_thresholds = []
    exposure_times = []
    background_concentrations = []
    persistence_calculation = []
    for specie in species:
        molar_ratio_specie, molar_weight_specie, concentration_thresholds_specie, exposure_times_specie, background_c, \
            persistence_specie = extract_gas_properties(specie)
        molar_ratios.append(molar_ratio_specie)
        molar_weights.append(molar_weight_specie)
        concentration_thresholds.append(concentration_thresholds_specie)
        exposure_times.append(exposure_times_specie)
        background_concentrations.append(background_c)
        persistence_calculation.append(persistence_specie)
    molar_ratio_specie, molar_weight_specie, concentration_thresholds_specie, exposure_times_specie, background_c, \
        persistence_tracking_specie = extract_gas_properties(original_specie)
    molar_ratios_tracking_specie = molar_ratio_specie
    molar_weights_tracking_specie = molar_weight_specie
    concentration_thresholds_tracking_specie = concentration_thresholds_specie
    exposure_times_tracking_specie = exposure_times_specie
    background_concentration_tracking_specie = background_c
    species_properties_out = []
    for i in range(0, len(species)):
        gas_specie = {"specie_name": species[i], "molar_ratio": molar_ratios[i], "molar_weight": molar_weights[i],
                      "concentration_thresholds": concentration_thresholds[i], "exposure_times": exposure_times[i],
                      "background_concentration": background_concentrations[i],
                      "persistence_calculation": persistence_calculation[i]}
        species_properties_out.append(gas_specie)
    if original_specie not in species:
        gas_specie = {"specie_name": original_specie, "molar_ratio": molar_ratios_tracking_specie,
                      "molar_weight": molar_weights_tracking_specie,
                      "concentration_thresholds": concentration_thresholds_tracking_specie,
                      "exposure_times": exposure_times_tracking_specie,
                      "background_concentration": background_concentration_tracking_specie,
                      "persistence_calculation": persistence_tracking_specie}
        species_properties_out.append(gas_specie)
    if not persistence_tracking_specie and True not in persistence_calculation:
        persistence = False
    return species_properties_out


def elaborate_tracking_points():
    import utm
    stations_out = []
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
                    print("WARNING. Invalide coordinate of the tracking point")
            else:
                station_northing = station_y
                station_easting = station_x
            if y0 <= station_northing <= yf and x0 <= station_easting <= xf and \
                    min(output_levels) <= station_z <= max(output_levels):
                station_id += 1
                stations_out.append({'station_id': station_id, 'easting': station_easting,
                                'northing': station_northing, 'elevation': station_z})
            else:
                continue
    return stations_out


def elaborate_day(day_input):
    def converter(input_file, processed_file, specie_input, model_input):
        conc = np.loadtxt(input_file, skiprows=5)
        conc[conc < 0] = 0
        molar_weight = 0
        conc_converted = np.empty_like(conc)
        if units == "ppm":
            if model_input == "disgas" or model_input == 'merged':
                file_dt = os.path.split(processed_file)[1]
                file_dt = file_dt.split("_")[2]
                file_dt = file_dt.split(".grd")[0]
                file_time_h = file_dt[-4:]
                input_file_name = input_file.split(os.sep)[-1]
                file_folder = input_file.split(input_file_name)[0]
                file_folder_daily = file_folder.split("outfiles")[0]
                if model_input == 'merged':
                    file_folder_daily = os.path.join(file_folder_daily, 'disgas')
                surface_data = os.path.join(file_folder_daily, "surface_data.txt")
                with open(surface_data) as surface_data_file:
                    for line in surface_data_file:
                        records = line.split("\t")
                        try:
                            float(records[1])
                        except IndexError:
                            continue
                        if file_time_h == records[0]:
                            t2m = float(records[2])
                            p2m = float(records[3]) / 100  # in hPa for this conversion
                            break
                for specie_property in species_properties:
                    if specie_property["specie_name"] == specie_input:
                        molar_weight = specie_property["molar_weight"]
                        conversion_factor = ((22.4 / molar_weight) * (t2m / 273) * (1013 / p2m)) * 1000000
                        conc_converted = np.multiply(conc, conversion_factor)  # convert kg/m3 to ppm
                if molar_weight == 0:
                    print('ERROR. Unable to find ' + specie_input + ' in the species properties database')
                    sys.exit()
            else:
                conc_converted = conc
        else:
            if model == "twodee":
                file_dt = os.path.split(processed_file)[1]
                file_dt = file_dt.split("_")[2]
                file_dt = file_dt.split(".grd")[0]
                file_time_h = file_dt[-4:]
                input_file_name = input_file.split(os.sep)[-1]
                file_folder = input_file.split(input_file_name)[0]
                file_folder_daily = file_folder.split("outfiles")[0]
                surface_data = os.path.join(file_folder_daily, "surface_data.txt")
                with open(surface_data) as surface_data_file:
                    for line in surface_data_file:
                        records = line.split("\t")
                        try:
                            float(records[1])
                        except IndexError:
                            continue
                        if file_time_h == records[0]:
                            t2m = float(records[2])
                            p2m = float(records[3]) / 100  # in hPa for this conversion
                            break
                for specie_property in species_properties:
                    if specie_property["specie_name"] == specie_input:
                        molar_weight = specie_property["molar_weight"]
                        conversion_factor = ((molar_weight / 22.4) * (273 / t2m) * (p2m / 1013)) / 1000000
                        conc_converted = np.multiply(conc, conversion_factor)  # convert ppm to kg/m3
                if molar_weight == 0:
                    print('ERROR. Unable to find ' + specie_input + ' in the species properties database')
                    sys.exit()
            else:
                conc_converted = conc
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
                    processed_file.write(str(np.amin(conc_converted)) + "  " + str(np.amax(conc_converted)) + "\n")
                    np.savetxt(processed_file, conc_converted, fmt="%.2e")
                else:
                    for specie_property in species_properties:
                        if specie_property["specie_name"] == specie_input:
                            molar_ratio = specie_property["molar_ratio"]
                            molar_weight = specie["molar_weight"]
                            try:
                                background_c = specie_property["background_concentration"]
                            except UnboundLocalError:
                                background_c = 0
                            try:
                                background_c = float(background_c)
                            except TypeError:
                                background_c = 0
                            if units == 'ppm':
                                species_conversion_factor = molar_ratio
                            elif units == 'kg/m3':
                                background_c = background_c * ((molar_weight / 22.4) * (273 / t2m) * (p2m / 1013)) / \
                                               1000000
                        if specie["specie_name"] == original_specie:
                            molar_weight_tracking_specie = specie["molar_weight"]
                            species_conversion_factor = molar_ratio * (molar_weight / molar_weight_tracking_specie)
                    conc_converted = np.where(conc_converted < 0, 0, conc_converted)
                    conc_converted += background_c
                    conc_converted = np.multiply(conc_converted, species_conversion_factor)
                    processed_file.write(str(np.amin(conc_converted)) + "  " + str(np.amax(conc_converted)) + "\n")
                    np.savetxt(processed_file, conc_converted, fmt="%.2e")

    def time_average(files_to_average_inp, outfile):
        conc_sum = 0
        for file_to_average in files_to_average_inp:
            if output_format == "grd":
                conc = np.loadtxt(file_to_average, skiprows=5)
                conc_sum += conc
        conc_avg = np.divide(conc_sum, len(files_to_average))
        # Create header of the processed file
        with open(outfile, "a") as processed_file:
            if output_format == "grd":
                processed_file.write("DSAA\n")
                processed_file.write(str(nx) + "  " + str(ny) + "\n")
                processed_file.write(str(x0) + "  " + str(xf) + "\n")
                processed_file.write(str(y0) + "  " + str(yf) + "\n")
                processed_file.write(
                    str(np.amin(conc_avg)) + "  " + str(np.amax(conc_avg)) + "\n"
                )
            np.savetxt(processed_file, conc_avg, fmt="%.2e")

    def prepare_persistence_calculation():
        def prepare_files(index):
            def calculate_overcome_time(output_files_to_elaborate):
                for j_overcome in range(ny):
                    for i_overcome in range(nx):
                        for file_input in output_files_to_elaborate:
                            try:
                                row = linecache.getline(file_input, j_overcome + 6)
                                if float(row.split(' ')[i_overcome]) > concentration_threshold_input:
                                    overcome_matrix_out[j_overcome][i_overcome] += simulation_time / n_time_steps / 3600
                                else:
                                    overcome_matrix_out[j_overcome][i_overcome] += 0
                            except IndexError:
                                print('File ' + file_input + ' not found')
                                overcome_matrix_out[j_overcome][i_overcome] += 0
                linecache.clearcache()

            specie_input = index[0]
            concentration_threshold_input = index[1]
            exposure_time_input = index[2]
            file_level_s_input = index[3]
            output_folder = os.path.join(processed_output_folder, day_input, specie_input)
            day_output_files = os.listdir(output_folder)
            output_files_day = []
            overcome_matrix_out = np.zeros((ny, nx))
            for output_file in day_output_files:
                if output_file.split('_')[1] == file_level_s_input:
                    output_files_day.append(os.path.join(output_folder, output_file))
            calculate_overcome_time(output_files_day)
            pers_output = os.path.join(persistence_folder, specie_input, 'C_' + str(concentration_threshold_input)
                                       + '_t_' + str(exposure_time_input) + 'H',
                                       'persistence_' + file_level_s + '.grd')
            return pers_output, overcome_matrix_out

        for specie_per in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie_per:
                    if specie_dict["persistence_calculation"]:
                        specie_folder = os.path.join(persistence_folder, specie_dict["specie_name"])
                        try:
                            os.mkdir(specie_folder)
                        except FileExistsError:
                            pass
                        concentration_thresholds = specie_dict["concentration_thresholds"]
                        exposure_times = specie_dict["exposure_times"]
                        for j in range(0, len(concentration_thresholds)):
                            threshold_folder = os.path.join(specie_folder, 'C_' + str(concentration_thresholds[j])
                                                            + '_t_' + str(exposure_times[j]) + 'H')
                            try:
                                os.mkdir(threshold_folder)
                            except FileExistsError:
                                pass

        indexes = []
        for specie_per in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie_per:
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    exposure_times = specie_dict["exposure_times"]
                    for i_threshold in range(0, len(concentration_thresholds)):
                        if float(exposure_times[i_threshold]) * 3600 > simulation_time:
                            continue  # Skip exposure times longer than the simulation itself
                        if levels[0] == "all":
                            for j in range(0, nz):
                                indexes.append([specie_per, concentration_thresholds[i_threshold],
                                                exposure_times[i_threshold], processed_files_levels_elaborated[j]])
                        else:
                            all_levels = np.array(processed_files_levels_elaborated)
                            levels_indexes = [int(x) - 1 for x in levels]
                            for level_i_threshold in list(all_levels[levels_indexes]):
                                indexes.append([specie_per, concentration_thresholds[i_threshold],
                                                exposure_times[i_threshold], level_i_threshold])
        for i_indexes in range(0, len(indexes)):
            persistence_matrix_temp = np.zeros((ny, nx))
            persistence_output, overcome_matrix = prepare_files(indexes[i_indexes])
            persistence_matrices[persistence_output] = persistence_matrix_temp
            overcome_matrices.append(overcome_matrix)
        return indexes

    def extract_tracking_points(files_to_interpolate):
        def interpolate(x, y, z, levels_interpolation, files):
            from scipy import interpolate
            x_array = np.linspace(x0, xf, nx, endpoint=True)
            y_array = np.linspace(y0, yf, ny, endpoint=True)
            conc1 = np.loadtxt(files[0], skiprows=5)
            if len(levels_interpolation) == 2:
                z_array = np.linspace(levels_interpolation[0], levels_interpolation[1], 2)
                conc2 = np.loadtxt(files[1], skiprows=5)
                conc = np.array([conc1, conc2])
                my_interpolating_function = interpolate.RegularGridInterpolator((x_array, y_array, z_array), conc.T)
                pt = np.array([x, y, z])
            else:
                conc = conc1
                my_interpolating_function = interpolate.RegularGridInterpolator((x_array, y_array), conc.T)
                pt = np.array([x, y])
            return my_interpolating_function(pt)

        global tracking_points_files
        files_time_steps = []
        files_time_averaging_steps = []
        for file_to_interpolate in files_to_interpolate:
            file_to_interpolate_name = file_to_interpolate.split(os.sep)[-1]
            file_to_interpolate_time_step = file_to_interpolate_name.split('_')[2]
            try:
                file_to_interpolate_time_step = int(file_to_interpolate_time_step.split('.grd')[0])
                if file_to_interpolate_time_step not in files_time_steps:
                    files_time_steps.append(file_to_interpolate_time_step)
            except ValueError:
                file_time_step_avg = file_to_interpolate_time_step.split('.grd')[0]
                if file_time_step_avg not in files_time_averaging_steps:
                    files_time_averaging_steps.append(file_time_step_avg)
        files_time_steps = sorted(files_time_steps)
        files_time_averaging_steps = sorted(files_time_averaging_steps)
        levels_for_interpolation = []
        for l_specie in range(0, len(species)):
            k = 0
            for station in stations:
                if min(output_levels) <= station['elevation'] <= max(output_levels):
                    for i_level in range(1, len(output_levels) + 1):
                        levels_for_interpolation = []
                        if output_levels[i_level - 1] == station['elevation']:
                            levels_for_interpolation.append(output_levels[i_level - 1])
                            break
                        elif output_levels[i_level - 1] < station['elevation'] < output_levels[i_level]:
                            levels_for_interpolation.append(output_levels[i_level - 1])
                            levels_for_interpolation.append(output_levels[i_level])
                            break
                        else:
                            continue
                c_interpolated_time_steps = []
                for time_step in files_time_steps:
                    files_to_use = []
                    files_levels = []
                    file_to_interpolate_directory = ''
                    for file_to_interpolate in files_to_interpolate:
                        file_to_interpolate_name = file_to_interpolate.split(os.sep)[-1]
                        file_to_interpolate_specie = file_to_interpolate.split(os.sep)[-2]
                        if file_to_interpolate_specie != species[l_specie]:
                            continue
                        file_to_interpolate_directory = file_to_interpolate.split(file_to_interpolate_name)[0]
                        file_to_interpolate_level = file_to_interpolate_name.split('_')[1]
                        file_to_interpolate_level = float(file_to_interpolate_level.split('mabg')[0])
                        file_to_interpolate_time_step = file_to_interpolate_name.split('_')[2]
                        try:
                            file_to_interpolate_time_step = int(file_to_interpolate_time_step.split('.grd')[0])
                        except ValueError:
                            continue
                        if file_to_interpolate_level in levels_for_interpolation and \
                                file_to_interpolate_time_step == time_step and file_to_interpolate_specie == \
                                species[l_specie]:
                            files_to_use.append(file_to_interpolate)
                            files_levels.append(file_to_interpolate_level)
                    if len(files_to_use) == 0:
                        continue
                    files_levels = sorted(files_levels)
                    c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'],
                                                 files_levels,
                                                 files_to_use)
                    c_interpolated_time_steps.append(c_interpolated[0])
                    tracking_point_file = os.path.join(file_to_interpolate_directory, 'TP_' +
                                                       str(station['station_id']) + '.txt')
                    tracking_points_files.append(tracking_point_file)
                    with open(tracking_point_file, 'a') as tracking_point_file:
                        tracking_point_file.write(str(time_step) + '\t' + "{0:.2e}".format(c_interpolated[0]) + '\n')
                for time_step in files_time_averaging_steps:
                    files_to_use = []
                    files_levels = []
                    file_to_interpolate_directory = ''
                    for file_to_interpolate in files_to_interpolate:
                        file_to_interpolate_name = file_to_interpolate.split(os.sep)[-1]
                        file_to_interpolate_specie = file_to_interpolate.split(os.sep)[-2]
                        if file_to_interpolate_specie != species[l_specie]:
                            continue
                        file_to_interpolate_directory = file_to_interpolate.split(file_name)[0]
                        file_to_interpolate_level = file_to_interpolate_name.split('_')[1]
                        file_to_interpolate_level = float(file_to_interpolate_level.split('mabg')[0])
                        file_to_interpolate_time_step = file_to_interpolate_name.split('_')[2]
                        try:
                            file_to_interpolate_time_step = int(file_to_interpolate_time_step.split('.grd')[0])
                        except ValueError:
                            file_to_interpolate_time_step = file_to_interpolate_time_step.split('.grd')[0]
                            if file_to_interpolate_level in levels_for_interpolation and \
                                    file_to_interpolate_time_step == time_step and \
                                    file_to_interpolate_specie == species[l_specie]:
                                files_to_use.append(file_to_interpolate)
                                files_levels.append(file_to_interpolate_level)
                    if len(files_to_use) == 0:
                        continue
                    files_levels = sorted(files_levels)
                    c_interpolated = interpolate(station['easting'], station['northing'], station['elevation'],
                                                 files_levels, files_to_use)
                    c_interpolated_time_steps.append(c_interpolated[0])
                    tracking_point_file = os.path.join(file_to_interpolate_directory, 'TP_' +
                                                       str(station['station_id']) + '.txt')
                    tracking_points_files.append(tracking_point_file)
                    with open(tracking_point_file, 'a') as tracking_point_file:
                        tracking_point_file.write(str(time_step) + '\t' + "{0:.2e}".format(c_interpolated[0]) + '\n')
                c_tp[species[l_specie]][k] = {'station_id': k, 'c_tp_time_steps': c_interpolated_time_steps}
                k += 1
        return files_time_steps + files_time_averaging_steps, c_tp

    run_folder = os.path.join(original_output_folder, day_input)
    run_folder_subfolders = os.listdir(run_folder)
    if 'twodee' in run_folder_subfolders:
        if 'disgas' in run_folder_subfolders:
            model = 'merged'
        else:
            model = 'disgas'
    else:
        model = 'disgas'
        if 'twodee' in run_folder_subfolders:
            model = 'merged'
        else:
            model = 'twodee'
    model_output_folder = os.path.join(run_folder, "outfiles")
    model_processed_output_folder_daily = os.path.join(processed_output_folder, day_input)
    try:
        os.mkdir(model_processed_output_folder_daily)
    except FileExistsError:
        print("Folder " + model_processed_output_folder_daily + " already exists")
    for specie in species:
        model_processed_output_folder_specie = os.path.join(model_processed_output_folder_daily, specie)
        try:
            os.mkdir(model_processed_output_folder_specie)
        except FileExistsError:
            print("Folder " + model_processed_output_folder_specie + " already exists")
        except PermissionError:  # retry
            try:
                os.mkdir(model_processed_output_folder_specie)
            except FileExistsError:
                print("Folder " + model_processed_output_folder_specie + " already exists")
    files_list_temp = os.listdir(model_output_folder)
    files_list_path = []
    files_list = []
    models = []
    for file in files_list_temp:
        if file[0:2] == "c_":
            files_list.append(file)
            files_list_path.append(os.path.join(model_output_folder, file))
            models.append(model)  # FABIO: check if "merged" does not cause any troubles here
    for specie in species[1:]:
        files_list_path += files_list_path
        models += models
    time_steps_disgas = []
    if model == 'disgas' or model == 'merged':
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
            processed_files.append(os.path.join(os.path.join(model_processed_output_folder_daily, specie),
                                                converted_file))
            processed_files_specie.append(os.path.join(os.path.join(model_processed_output_folder_daily, specie),
                                                       converted_file))
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
        start_file = n_elaborated_files
        end_file = n_elaborated_files + max_number_processes
        if end_file > len(files_list_path):
            end_file = len(files_list_path)
        pool_files = ThreadingPool(max_number_processes)
        pool_files.map(converter, files_list_path[start_file:end_file], processed_files[start_file:end_file],
                       species_list[start_file:end_file], models[start_file:end_file],)
        n_elaborated_files = end_file
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
            if int((time_max - time_min).total_seconds()) <= time_step_simulation:
                print('Warning! Time-averaging interval smaller than or equal to the time step of the simulation')
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
                        + "-tavg.grd",
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
                # time_min = time_max + datetime.timedelta(seconds=time_step_simulation)
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
        all_time_steps_tp, c_tp_time_steps = extract_tracking_points(processed_files)
    return day_input, all_time_steps_tp, processed_files_levels_elaborated, tavg_intervals, persistence_matrices, \
        indexes_persistence, overcome_matrices, c_tp


def calculate_persistence():
    for pers_output_file in persistence_matrices:
        exp_time = pers_output_file.split('_t_')[1]
        exp_time = float(exp_time.split('H')[0])
        c_thresh = pers_output_file.split('_t_')[0]
        c_thresh = float(c_thresh.split('C_')[1])
        pers_level = pers_output_file.split('persistence_')[1]
        pers_level = pers_level.split('.grd')[0]
        specie = pers_output_file.split('/')[-3]
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
                            if float(overcome_time_data[j_p][i_p]) >= exp_time:
                                persistence_matrices[pers_output_file][j_p][i_p] += weight_simulation


def probabilistic_tracking_points():
    def plot_hazard_curves(file_input, folder, min_con_tp_in, max_con_tp_in):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        c_ecdf = []
        with open(file_input, 'r') as tp_ecdf_file:
            lines = tp_ecdf_file.readlines()[1:]
            time_steps_tp = []
            for line in lines:
                entries = line.split('\t')[1:]
                time_steps_tp.append(line.split('\t')[0])
                for i in range(0, len(entries)):
                    entries[i] = float(entries[i])
                c_ecdf.append(entries)
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        tp_file = file_input.split(os.sep)[-1]
        tp_file = tp_file.split('.txt')[0]
        specie_name = file_input.split(os.sep)[-2]
        specie_name = specie_name.translate(sub)
        for i in range(0, len(c_ecdf)):
            i_tp_ecdf = 1
            if 'tavg' in str(time_steps_tp[i]):
                output_plot_file = os.path.join(folder, tp_file + '_tavg-time_interval' + str(i_tp_ecdf) + '.png')
                i_tp_ecdf += 1
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
    output_quantile_tp = [[[[0 for i in range(0, len(all_time_steps))] for m in range(0, len(quantiles))]
                       for k in range(0, len(stations))] for l_st in range(0, len(stations))]
    c_list = []
    for specie in species:
        graphical_outputs_ecdf_tracking_points_specie = os.path.join(graphical_outputs_ecdf_tracking_points_folder,
                                                                     specie)
        try:
            os.mkdir(graphical_outputs_ecdf_tracking_points_specie)
        except FileExistsError:
            print('Folder ' + graphical_outputs_ecdf_tracking_points_specie + ' already exists')
        ecdf_tracking_points_specie = os.path.join(ecdf_tracking_points_folder, specie)
        try:
            os.mkdir(ecdf_tracking_points_specie)
        except FileExistsError:
            print('Folder ' + ecdf_tracking_points_specie + ' already exists')
    min_con_tp_specie = []
    max_con_tp_specie = []
    for l_sp in range(0, len(species)):
        min_con_tp = 1000000000000000
        max_con_tp = 0
        for k in range(0, len(stations)):
            for i in range(0, len(all_time_steps)):
                for j in range(0, len(days)):
                    c_list.append(c[l_sp][k][j][i])
                for m in range(0, len(quantiles)):
                    output_quantile_tp[l_sp][k][m][i] = np.quantile(c_list, q=quantiles[m])
                    if output_quantile_tp[l_sp][k][m][i] <= min_con_tp:
                        min_con_tp = output_quantile_tp[l_sp][k][m][i]
                    if output_quantile_tp[l_sp][k][m][i] >= max_con_tp:
                        max_con_tp = output_quantile_tp[l_sp][k][m][i]
                c_list = []
        min_con_tp_specie.append(min_con_tp)
        max_con_tp_specie.append(max_con_tp)
    for l_sp in range(0, len(species)):
        for k in range(0, len(stations)):
            station = stations[k]
            ecdf_tracking_point_file = os.path.join(ecdf_tracking_points_folder, species[l_sp], 'TP_' +
                                                    str(station['station_id']) + '_ecdf.txt')
            with open(ecdf_tracking_point_file, 'w') as output_file:
                output_file.write(quantile_string + '\n')
                for i in range(0, len(all_time_steps)):
                    i_tavg = 1
                    if 'tavg' in str(all_time_steps[i]):
                        output_quantile_string = 'tavg-time_interval_' + str(i_tavg)
                        i_tavg += 1
                    else:
                        output_quantile_string = 'time_step_' + str(i + 1)
                    for m in range(0, len(quantiles)):
                        output_quantile_string += '\t' + "{0:.2e}".format(output_quantile_tp[l_sp][k][m][i])
                    output_file.write(output_quantile_string + '\n')
            plot_file_folder = os.path.join(graphical_outputs_ecdf_tracking_points_folder, species[l_sp])
            plot_hazard_curves(ecdf_tracking_point_file, plot_file_folder, min_con_tp_specie[l_sp],
                               max_con_tp_specie[l_sp])


def prepare_quantile_calculation(exc_prob):
    def prepare_files(index):
        specie_exc_prob = index[1]
        file_level_s = index[2]
        time_step_exc_prob = index[3]
        ex_prob = index[0]
        output_files = []
        time_step_s = ''
        for day_exc_prob in days:
            try:
                file_time_step = int(time_step_exc_prob)
                file_time_step = file_time_step * dt
                time_start = datetime.datetime.strptime(day_exc_prob + str(hour_start).zfill(2) +
                                                        str(minute_start).zfill(2), '%Y%m%d%H%M')
                time_validity = time_start + datetime.timedelta(seconds=file_time_step)
                file_validity = datetime.datetime.strftime(time_validity, '%Y%m%d%H%M')
                time_step_s = "{:06d}".format(int(time_step_exc_prob))
            except ValueError:
                tavg_interval_start_s = day_exc_prob + time_step_exc_prob.split('-')[0]
                tavg_interval_end_s = day_exc_prob + time_step_exc_prob.split('-')[1]
                tavg_interval_start = datetime.datetime.strptime(tavg_interval_start_s, '%Y%m%d%H%M')
                try:
                    tavg_interval_end = datetime.datetime.strptime(tavg_interval_end_s, '%Y%m%d%H%M')
                except ValueError:  # in case time_step.split('-')[1] = 2400
                    if(int(time_step_exc_prob.split('-')[0][0:2]) + int(time_step_exc_prob.split('-')[1][0:2])) > 24:
                        tavg_interval_end = tavg_interval_start + \
                                            datetime.timedelta(hours=int(time_step_exc_prob.split('-')[1][0:2]) -
                                                               int(time_step_exc_prob.split('-')[0][0:2]))
                    else:
                        tavg_interval_end = tavg_interval_start + \
                                            datetime.timedelta(hours=int(time_step_exc_prob.split('-')[0][0:2]) +
                                                               int(time_step_exc_prob.split('-')[1][0:2]))
                if tavg_interval_end < tavg_interval_start:
                    tavg_interval_end = tavg_interval_start + \
                                        datetime.timedelta(hours=int(time_step_exc_prob.split('-')[0][0:2]) +
                                                           int(time_step_exc_prob.split('-')[1][0:2]))
                tavg_interval_end_s = datetime.datetime.strftime(tavg_interval_end, '%Y%m%d%H%M')
                day_interval = datetime.datetime.strftime(tavg_interval_start, '%Y%m%d')
                if day_exc_prob != day_interval:
                    continue
                file_validity = tavg_interval_start_s + '-' + tavg_interval_end_s + '-tavg'
                time_step_s = time_step_exc_prob
            file_name = "c_" + file_level_s + "_" + file_validity + ".grd"
            output_folder = os.path.join(processed_output_folder, day_exc_prob, specie_exc_prob)
            output_files.append(os.path.join(output_folder, file_name))
        ecdf_output_file = os.path.join(ecdf_folder, str(ex_prob), specie_exc_prob, "c_" + file_level_s + "_" +
                                        time_step_s + ".grd",)
        return output_files, ecdf_output_file

    all_output_files_exc_prob = []
    all_ecdf_output_files_exc_prob = []
    for exceedance_probability in exceedance_probabilities:
        prob_folder = os.path.join(ecdf_folder, str(exceedance_probability))
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
                        indexes_tavg.append([exc_prob, specie, processed_files_levels[i], tavg_intervals[k]])
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
        all_output_files_exc_prob.append(a)
        all_ecdf_output_files_exc_prob.append(b)
    if len(tavg_intervals) > 0:
        for i_indexes in range(0, len(indexes_tavg)):
            a, b = prepare_files(indexes_tavg[i_indexes])
            all_output_files_exc_prob.append(a)
            all_ecdf_output_files_exc_prob.append(b)

    return all_output_files_exc_prob, all_ecdf_output_files_exc_prob


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
    # pool_file_read = ThreadingPool(len(files_to_process))
    # c_list_ecdf = pool_file_read.map(read_file, files_to_process)
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


def save_plots(min_con_in, max_con_in):
    import re

    def plot_file(input_file, output_file, dz_lines_res_in):
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

        def myround(z, prec=2, base=100):
            return round(base * round(float(z) / base), prec)

        def resize_topography(bottom_left_easting, top_right_easting, bottom_left_northing, top_right_northing,
                              topography):
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
            dx_top = (xf_or_top - x0_or_top) / (nx_or_top - 1)
            dy_top = (yf_or_top - y0_or_top) / (ny_or_top - 1)
            i_bottom_left = 0
            j_bottom_left = 0
            i_top_right = 0
            j_top_right = 0
            x_top_resized = x0_or_top
            y_top_resized = y0_or_top
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
            while x_top_resized <= bottom_left_easting:
                i_bottom_left += 1
                i_top_right += 1
                x_top_resized += dx_top
            while x_top_resized <= top_right_easting:
                i_top_right += 1
                x_top_resized += dx_top
            while y_top_resized <= bottom_left_northing:
                j_bottom_left += 1
                j_top_right += 1
                y_top_resized += dy_top
            while y_top_resized <= top_right_northing:
                j_top_right += 1
                y_top_resized += dy_top
            z_topography = np.loadtxt("topography.grd", skiprows=5)
            z_resized = z_topography[j_bottom_left:j_top_right, i_bottom_left:i_top_right]
            nx_resized = i_top_right - i_bottom_left
            ny_resized = j_top_right - j_bottom_left
            z_min = np.amin(z_resized)
            z_max = np.amax(z_resized)
            return nx_resized, ny_resized, z_min, z_max, z_resized

        try:
            os.remove(output_file)
        except FileNotFoundError:
            pass
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        if 'persistence' in input_file.split(os.sep)[-1]:
            specie_name = input_file.split(os.sep)[-3]
        else:
            specie_name = input_file.split(os.sep)[-2]
        specie_name = specie_name.translate(sub)
        with open(input_file) as input_file_to_plot:
            if output_format == "grd":
                conc_inp = np.loadtxt(input_file_to_plot, skiprows=5)
            x = np.linspace(x0, xf, num=nx)
            y = np.linspace(y0, yf, num=ny)
            n_levels = 10
            dc = (max_con - min_con) / n_levels
            conc_levels = np.arange(min_con + 0.0000001, max_con, dc)
            if 'persistence' in input_file:
                n_levels = 100
                dp = 1. / n_levels
                prob_levels = np.arange(0.00000000001, 1, dp)
                if output_format == "grd":
                    prob = np.loadtxt(input_file, skiprows=5)
        fig, ax = plt.subplots(figsize=(6, 5), dpi=plot_resolution)
        if plot_topography_layer:
            nx_top, ny_top, min_z, max_z, z_top = resize_topography(x0, xf, y0, yf, "topography.grd")
            x_top = np.linspace(x0, xf, num=nx_top)
            y_top = np.linspace(y0, yf, num=ny_top)
            n_levels = 100
            dz = (max_z - min_z) / n_levels
            if dz_lines_res_in >= max_z:
                dz_lines_res_in = max_z
            n_levels_lines = int((max_z - min_z) / dz_lines_res_in)
            dz_lines = myround((max_z - min_z) / n_levels_lines, base=dz_lines_res_in)
            levels_top = np.arange(min_z + 0.0000001, max_z, dz)
            levels_top_lines = np.arange(myround(min_z, base=dz_lines_res_in),
                                         myround(max_z, base=dz_lines_res_in), dz_lines)
            top = ax.contourf(x_top, y_top, z_top, levels_top, cmap="Greys", extend="max")
            top_lines = ax.contour(top, levels=levels_top_lines, colors='black', linewidths=0.05)
            top_cbar = fig.colorbar(top, orientation="horizontal", format="%.1f", shrink=0.75)
            ax.clabel(top_lines, inline=True, fontsize=2, fmt='%1.0f')
            top_cbar.ax.tick_params(labelsize=6)
            top_cbar.set_label("m a.s.l.")
        if 'persistence' in input_file.split(os.sep)[-1]:
            cmap = plt.get_cmap('viridis', 10)
            field = plt.contourf(x, y, prob, prob_levels, cmap=cmap, alpha=0.9, extend="max")
        else:
            field = plt.contourf(x, y, conc_inp, conc_levels, cmap="Reds", alpha=0.9, extend="max")
        if len(plot_isolines) != 0:
            specie_isolines = ax.contour(x, y, conc_inp, levels=plot_isolines, colors='black', linewidths=0.2)
            ax.clabel(specie_isolines, inline=True, fontsize=4, fmt='%1.1f')
        aspect = 20
        pad_fraction = 0.5
        divider = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1.0 / aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax_c = divider.append_axes("right", size=width, pad=pad)
        if 'persistence' in input_file.split(os.sep)[-1]:
            thresholds = input_file.split(os.sep)[-2]
            c_threshold = thresholds.split('_t_')[0]
            c_threshold = c_threshold.split('C_')[1]
            exposure_time = thresholds.split('_t_')[1]
            exposure_time = exposure_time.split('H')[0]
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
        plt.savefig(output_file)
        image_buffer.close()
        plt.close(fig)

    files_to_plot = []
    output_files = []

    graphical_outputs = os.path.join(outputs_folder, "graphical_outputs")
    graphical_outputs_simulations = os.path.join(graphical_outputs, "simulations")
    graphical_outputs_ecdf = os.path.join(graphical_outputs, "ecdf")
    for day_to_plot in days_to_plot:
        graphical_outputs_daily = os.path.join(graphical_outputs_simulations, day_to_plot)
        try:
            os.mkdir(graphical_outputs_daily)
        except FileExistsError:
            print("Folder " + graphical_outputs_daily + " already exists")
        model_processed_output_folder_daily = os.path.join(processed_output_folder, day_to_plot)
        model_processed_output_folder_species = []
        for specie in species:
            model_processed_output_folder_species.append(os.path.join(model_processed_output_folder_daily, specie))
        for specie in species:
            try:
                os.mkdir(os.path.join(graphical_outputs_daily, specie))
            except FileExistsError:
                print("Folder " + os.path.join(graphical_outputs_daily, specie) + " already exists")
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
            simulation_start = datetime.datetime.strptime(day_to_plot + "{:02d}".format(hour_start), '%Y%m%d%H%M')
            output_file_name = files_list[i].split(".grd")[0]
            output_file_name += ".png"
            if levels[0] == "all":
                if time_steps[0] == "all":
                    files_to_plot.append(file)
                    output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name))
                else:
                    for time_step in time_steps:
                        time_step_seconds = hour_start + dt * int(time_step)
                        time_step_datetime = simulation_start + datetime.timedelta(seconds=time_step_seconds)
                        if time_step_datetime == file_time_step_datetime:
                            files_to_plot.append(file)
                            output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name,))
                if "tavg" in file_time_step:
                    files_to_plot.append(file)
                    tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                    tavg_output_file_name = tavg_output_file_name + ".png"
                    output_files.append(os.path.join(graphical_outputs_daily, file_specie, tavg_output_file_name,))
            else:
                if time_steps[0] == "all":
                    for level in levels:
                        if file_level == processed_files_levels[int(level) - 1]:
                            files_to_plot.append(file)
                            output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name,))
                else:
                    for level in levels:
                        for time_step in time_steps:
                            time_step_seconds = hour_start + dt * int(time_step)
                            time_step_datetime = simulation_start + datetime.timedelta(seconds=time_step_seconds)
                            if time_step_datetime == file_time_step_datetime \
                                    and file_level == processed_files_levels[int(level) - 1]:
                                files_to_plot.append(file)
                                output_files.append(os.path.join(graphical_outputs_daily, file_specie,
                                                                 output_file_name,))
                for level in levels:
                    if "tavg" in file_time_step and file_level == processed_files_levels[int(level) - 1]:
                        files_to_plot.append(file)
                        tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                        tavg_output_file_name = tavg_output_file_name + ".png"
                        output_files.append(os.path.join(graphical_outputs_daily, file_specie, tavg_output_file_name,))
            i += 1

    if calculate_ecdf:
        for exceedance_probability in exceedance_probabilities:
            try:
                os.mkdir(os.path.join(graphical_outputs_ecdf, str(exceedance_probability)))
            except FileExistsError:
                print("Folder " + os.path.join(graphical_outputs_ecdf, str(exceedance_probability)) + " already exists")
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_ecdf, str(exceedance_probability), specie))
                except FileExistsError:
                    print("Folder " + os.path.join(graphical_outputs_ecdf, str(exceedance_probability), specie) +
                          " already exists")
                files_list = os.listdir(os.path.join(ecdf_folder, str(exceedance_probability), specie))
                for file in files_list:
                    file_path = os.path.join(ecdf_folder, str(exceedance_probability), specie, file)
                    file_name_splitted = file.split("_")
                    file_level = file_name_splitted[1]
                    file_time_step = file_name_splitted[2].split(".")[0]
                    output_file_name = file.split(".grd")[0]
                    output_file_name += ".png"
                    if levels[0] == "all":
                        if time_steps[0] == "all":
                            files_to_plot.append(file_path)
                            output_files.append(os.path.join(graphical_outputs_ecdf, str(exceedance_probability),
                                                             specie, output_file_name,))
                        else:
                            for time_step in time_steps:
                                if file_time_step == "{:06d}".format(int(time_step)):
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf,
                                                        str(exceedance_probability), specie, output_file_name,))
                        if "tavg" in file_time_step:
                            files_to_plot.append(file_path)
                            tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                            tavg_output_file_name = tavg_output_file_name + ".png"
                            output_files.append(os.path.join(graphical_outputs_ecdf, str(exceedance_probability),
                                                specie, tavg_output_file_name,))
                    else:
                        if time_steps[0] == "all":
                            for level in levels:
                                if file_level == processed_files_levels[int(level) - 1]:
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf,
                                                        str(exceedance_probability), specie, output_file_name,))
                        else:
                            for level in levels:
                                for time_step in time_steps:
                                    if file_time_step == "{:06d}".format(int(time_step)) and \
                                            file_level == processed_files_levels[int(level) - 1]:
                                        files_to_plot.append(file_path)
                                        output_files.append(os.path.join(graphical_outputs_ecdf,
                                                            str(exceedance_probability), specie, output_file_name,))
                            for level in levels:
                                if "tavg" in file_time_step and file_level == processed_files_levels[int(level) - 1]:
                                    files_to_plot.append(file_path)
                                    tavg_output_file_name = file.split(os.sep)[-1].split(".grd")[0]
                                    tavg_output_file_name = tavg_output_file_name + ".png"
                                    output_files.append(os.path.join(graphical_outputs_ecdf,
                                                        str(exceedance_probability), specie, tavg_output_file_name,))
    if persistence:
        for specie in species:
            for specie_dict in species_properties:
                if specie_dict["specie_name"] == specie:
                    try:
                        os.mkdir(os.path.join(graphical_outputs_persistence_folder, specie))
                    except FileExistsError:
                        print("Folder " + os.path.join(graphical_outputs_persistence_folder, specie) +
                              " already exists")
                    concentration_thresholds = specie_dict["concentration_thresholds"]
                    exposure_times = specie_dict["exposure_times"]
                    for j in range(0, len(concentration_thresholds)):
                        try:
                            os.mkdir(os.path.join(graphical_outputs_persistence_folder, specie, 'C_' +
                                     str(concentration_thresholds[j]) + '_t_' + str(exposure_times[j]) + 'H'))
                        except FileExistsError:
                            print("Folder " + os.path.join(graphical_outputs_persistence_folder, specie, 'C_' +
                                  str(concentration_thresholds[j]) + '_t_' + str(exposure_times[j]) + 'H') +
                                  " already exists")
                        files_list = os.listdir(os.path.join(persistence_folder, specie, 'C_' +
                                                str(concentration_thresholds[j]) + '_t_' +
                                                str(exposure_times[j]) + 'H'))
                        for file in files_list:
                            file_path = os.path.join(persistence_folder, specie, 'C_' + str(concentration_thresholds[j])
                                                     + '_t_' + str(exposure_times[j]) + 'H', file)
                            files_to_plot.append(file_path)
                            persistence_plot_file_name = file.split(os.sep)[-1].split(".grd")[0]
                            persistence_plot_file_name += '.png'
                            output_files.append(os.path.join(graphical_outputs_persistence_folder, specie, 'C_' +
                                                             str(concentration_thresholds[j]) + '_t_' +
                                                             str(exposure_times[j]) + 'H',
                                                             persistence_plot_file_name))
    if len(files_to_plot) == 0:
        print("No files to plot")
    else:
        if min_con_in == -1.0 and max_con_in == -1.0:
            for specie in species:
                files_to_plot_specie = []
                output_files_specie = []
                max_con_in = 0
                min_con_in = 1000000000000000
                i = 0
                for file_to_plot in files_to_plot:
                    if specie in file_to_plot:
                        files_to_plot_specie.append(file_to_plot)
                        output_files_specie.append(output_files[i])
                        if 'persistence' not in file_to_plot:
                            conc = np.loadtxt(file_to_plot, skiprows=5)
                            max_c = np.amax(conc)
                            min_c = np.amin(conc)
                            if max_c > max_con_in:
                                max_con_in = max_c
                            if min_c < min_con_in:
                                min_con_in = min_c
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
    days_to_plot,
    species,
    original_specie,
    exceedance_probabilities,
    max_number_processes,
    convert,
    persistence,
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
    outputs_folder,
    original_output_folder,
    processed_output_folder,
    ecdf_folder,
    ecdf_tracking_points_folder,
    persistence_folder,
    graphical_outputs_folder,
    graphical_outputs_simulations_folder,
    graphical_outputs_ecdf_folder,
    graphical_outputs_ecdf_tracking_points_folder,
    graphical_outputs_persistence_folder
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
    for specie_tp in species:
        c_tp[specie_tp] = {}
    for day in days:
        overcome_matrices_all_days[day] = ''
    if tracking_points or calculate_ecdf or persistence:
        days_to_elaborate = days
    else:
        days_to_elaborate = days_to_plot
    x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, simulation_time, output_levels, hour_start, \
        minute_start = domain()
    species_properties = gas_properties()
    if tracking_points:
        stations = elaborate_tracking_points()
        # Initialize array of concentration to be used for ECDFs in the tracking points
        c = [[[[0 for i in range(0, n_time_steps + 10)] for j in range(0, len(days))]
             for k in range(0, len(stations))] for l_sp in range(0, len(species))]
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
        save_plots(min_con, max_con)
