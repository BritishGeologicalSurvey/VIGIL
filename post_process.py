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
from shutil import rmtree


def read_arguments():
    parser = argparse.ArgumentParser(description="Input data")
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
        "-MO",
        "--merge_outputs",
        default="False",
        help="Merge Twodee and Disgas outputs (true or false)",
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
    models = args.models
    merge_outputs = args.merge_outputs
    units = args.units
    plot_limits_in = args.plot_limits
    time_av = args.time_av
    output_format = args.output_format
    plot_topography = args.plot_topography
    dz_lines_res = args.topography_isolines
    plot_resolution = args.plot_resolution
    ex_prob = ex_prob_in.split(',')
    time_steps = time_steps_in.split(',')
    levels = levels_in.split(',')
    days_plot = days_plot_in.split(',')
    plot_limits = plot_limits_in.split(',')
    species = species_in.split(',')
    if plot.lower() == "true":
        plot = True
        if days_plot_in == '':
            print("ERROR. Please specify at least one day to plot when --plot==True")
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
    if merge_outputs.lower() == "true":
        merge_outputs = True
    elif merge_outputs.lower() == "false":
        merge_outputs = False
    else:
        print("ERROR. Wrong value for variable -MO --merge_outputs")
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
    if len(plot_limits) > 0:
        try:
            min_con = float(plot_limits[0])
            max_con = float(plot_limits[1])
        except ValueError:
            print(
                "ERROR. Please specify valid minimum and maximum concentration -PL --plot_limits"
            )
            sys.exit()
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
    return (
        plot,
        plot_ex_prob,
        time_steps,
        levels,
        days_plot,
        species,
        original_specie,
        exceedance_probabilities,
        nproc,
        convert,
        models,
        merge_outputs,
        units,
        time_av,
        min_con,
        max_con,
        output_format,
        plot_topography_layer,
        dz_lines_res,
        plot_resolution,
    )


def folder_structure():
    original_output_folder_name = "simulations"
    post_processing = "post_processing"
    processed_output_folder_name = original_output_folder_name + "_processed"
    ecdf_folder_name = "output_ecdf"
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
    disgas_ecdf = os.path.join(disgas_outputs, ecdf_folder_name)
    twodee_processed_output_folder = os.path.join(
        twodee_outputs, processed_output_folder_name
    )
    twodee_ecdf = os.path.join(twodee_outputs, ecdf_folder_name)
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
            list_temp = os.listdir(disgas_processed_output_folder)
            for item in list_temp:
                try:
                    rmtree(
                        os.path.join(disgas_processed_output_folder, item),
                        ignore_errors=True,
                    )
                except FileNotFoundError:
                    print(
                        "Unable to remove "
                        + item
                        + " in "
                        + disgas_processed_output_folder
                    )
        try:
            os.mkdir(disgas_ecdf)
        except FileExistsError:
            print("Folder " + disgas_ecdf + " already exists")
            list_temp = os.listdir(disgas_ecdf)
            for item in list_temp:
                try:
                    rmtree(os.path.join(disgas_ecdf, item), ignore_errors=True)
                except FileNotFoundError:
                    print("Unable to remove " + item + " in " + disgas_ecdf)
    if models == "twodee" or models == "all":
        try:
            os.mkdir(twodee_outputs)
        except FileExistsError:
            print("Folder " + twodee_outputs + " already exists")
        try:
            os.mkdir(twodee_processed_output_folder)
        except FileExistsError:
            print("Folder " + twodee_processed_output_folder + " already exists")
            list_temp = os.listdir(twodee_processed_output_folder)
            for item in list_temp:
                try:
                    rmtree(
                        os.path.join(twodee_processed_output_folder, item),
                        ignore_errors=True,
                    )
                except FileNotFoundError:
                    print(
                        "Unable to remove "
                        + item
                        + " in "
                        + twodee_processed_output_folder
                    )
        try:
            os.mkdir(twodee_ecdf)
        except FileExistsError:
            print("Folder " + twodee_ecdf + " already exists")
            list_temp = os.listdir(twodee_ecdf)
            for item in list_temp:
                try:
                    rmtree(os.path.join(twodee_ecdf, item), ignore_errors=True)
                except FileNotFoundError:
                    print("Unable to remove " + item + " in " + twodee_ecdf)
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

    return (
        disgas_outputs,
        disgas_original_output_folder,
        disgas_processed_output_folder,
        ecdf_folder_name,
        disgas_ecdf,
        twodee_outputs,
        twodee_original_output_folder,
        twodee_processed_output_folder,
        twodee_ecdf,
        models_to_elaborate,
        twodee_output_time_step,
    )


def gas_properties():
    def extract_gas_properties(specie):
        data = pd.read_csv(gas_properties_file, error_bad_lines=False)
        molar_ratio = None
        molar_weight = None
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
            y = np.sort(data[specie])
            molar_weight = list(y)[0]
            if molar_weight != molar_weight:
                print('ERROR. Molar weight of ' + specie + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Molar weight of ' + specie + ' not found in gas_properties.csv')
            sys.exit()
        return molar_ratio, molar_weight

    gas_properties_file = os.path.join(root, "gas_properties.csv")
    try:
        open(gas_properties_file, "r")
    except FileNotFoundError:
        print("File " + gas_properties_file + " not present")
        sys.exit()
    molar_ratios = []
    molar_weights = []
    for specie in species:
        molar_ratio, molar_weight = extract_gas_properties(specie)
        molar_ratios.append(molar_ratio)
        molar_weights.append(molar_weight)
    molar_ratio, molar_weight = extract_gas_properties(original_specie)
    molar_ratios_tracking_specie = molar_ratio
    molar_weights_tracking_specie = molar_weight
    species_properties = []
    for i in range(0, len(species)):
        gas_specie = {}
        gas_specie["specie_name"] = species[i]
        gas_specie["molar_ratio"] = molar_ratios[i]
        gas_specie["molar_weight"] = molar_weights[i]
        species_properties.append(gas_specie)
    gas_specie = {}
    gas_specie["specie_name"] = original_specie
    gas_specie["molar_ratio"] = molar_ratios_tracking_specie
    gas_specie["molar_weight"] = molar_weights_tracking_specie
    species_properties.append(gas_specie)
    return species_properties


def domain(model):
    import re

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
                        extracted_heights = re.findall(
                            "\d+\.\d+", heights
                        )  # This extracts decimal numbers only!
                        nz = len(extracted_heights)
                        for height in extracted_heights:
                            output_levels.append(float(height))
                        output_levels = sorted(output_levels)
                    elif "OUTPUT_INTERVAL_(SEC)" in record_splitted[0]:
                        dt = float(temp[0])
                except (IndexError, ValueError):
                    continue
    yf = y0 + ny * dy
    xf = x0 + nx * dx
    n_time_steps = int(tot_time / dt)
    return x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, output_levels


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
            if days_plot[0] == "all":
                days_to_plot.append(day)
            else:
                for day_to_plot in days_plot:
                    if day_to_plot == day:
                        days_to_plot.append(day_to_plot)
    return days, days_to_plot


def converter(input_file, processed_file, specie_input, model):
    Z = np.loadtxt(input_file, skiprows=5)
    Z[Z < 0] = 0
    if units == "ppm":
        if model == "disgas":
            file_time_step = os.path.split(processed_file)[1]
            file_time_step = file_time_step.split("_")[2]
            file_time_step = int(file_time_step.split(".grd")[0])
            file_time_s = file_time_step * dt
            hours, remainder = divmod(file_time_s, 3600)
            minutes, seconds = divmod(remainder, 60)
            hours_int = int(hours)
            minutes_int = int(minutes)
            # Round to the next hours. It can be improved with a linear interpolation
            if minutes_int > 30:
                minutes_int = 0
                hours_int += 1
            if (
                hours_int > 23
            ):
                # First order approximation for the last time step. To be modified if simulations > 24 hours are to be
                # implemented
                hours_int = 23
            file_time_h = "{:02}{:02}".format(hours_int, minutes_int)
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
                (22.4 / molar_weight) * ((273 + t2m) / 273) * (1013 / p2m)
            ) * 1000000
            Z_converted = np.multiply(Z, conversion_factor)  # convert kg/m3 to ppm
        else:
            Z_converted = Z
    else:
        if model == "twodee":
            file_time_step = os.path.split(processed_file)[1]
            file_time_step = file_time_step.split("_")[2]
            file_time_step = int(file_time_step.split(".grd")[0])
            file_time_s = file_time_step * dt
            hours, remainder = divmod(file_time_s, 3600)
            minutes, seconds = divmod(remainder, 60)
            hours_int = int(hours)
            minutes_int = int(minutes)
            # Round to the next hours. It can be improved with a linear interpolation
            if minutes_int > 30:
                minutes_int = 0
                hours_int += 1
            if (
                hours_int > 23
            ):
                # First order approximation for the last time step. To be modified if simulations > 24 hours are to be
                # implemented
                hours_int = 23
            file_time_h = "{:02}{:02}".format(hours_int, minutes_int)
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
                (molar_weight / 22.4) * (273 / (273 + t2m)) * (p2m / 1013)
            ) / 1000000
            Z_converted = np.multiply(Z, conversion_factor)  # convert ppm to kg/m3
        else:
            Z_converted = Z
    # Create header of the processed file
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
    processed_file.close()


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
    converted_files = []
    processed_files = []
    species_list = []
    processed_files_species = []
    levels = []
    time_steps = []
    tavg_intervals = []
    for specie in species:
        processed_files_specie = []
        for file in files_list:
            species_list.append(specie)
            if model == "twodee":
                file_name_splitted = file.split("_")
                file_level = file_name_splitted[1]
                file_time_step = file_name_splitted[2].split(".")[0]
                file_level = float(file_level.split("cm")[0]) / 100
                try:
                    file_level_index = output_levels.index(file_level)
                except BaseException:
                    print(
                        "Cannot find the expected TWODEE output file at the level "
                        + str(file_level)
                    )
                    sys.exit()
                file_level = "{:03d}".format(file_level_index + 1)
                file_time_step = "{:06d}".format(
                    int((int(file_time_step) / twodee_output_time_step))
                )
                file = "c_" + file_level + "_" + file_time_step + ".grd"
            else:
                file_name_splitted = file.split("_")
                file_level = file_name_splitted[1]
                file_time_step = file_name_splitted[2]
                file_time_step = file_time_step.split(".")[0]
            if file_level not in levels:
                levels.append(file_level)
            if file_time_step not in time_steps:
                time_steps.append(int(file_time_step))
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
            tavg_intervals.append(str(time_min) + "-" + str(time_max) + "-tavg")
        else:
            time_max = time_min + time_av - 1
        while time_max <= max(time_steps):
            # tavg_intervals.append(str(time_min) + '-' + str(time_max) + '-tavg')
            for i in range(0, len(species)):
                files_to_average = []
                for level in levels:
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
                        + str(time_min)
                        + "-"
                        + str(time_max)
                        + "-tavg.grd",
                    )
                    for file in processed_files_species[i]:
                        file_name = os.path.split(file)[1]
                        file_level = file_name.split("_")[1]
                        file_time_step = file_name.split("_")[2]
                        file_time_step = file_time_step.split(".")[0]
                        if file_level == level:
                            if time_min <= int(file_time_step) <= time_max:
                                files_to_average.append(file)
                                averaged_files.append(file)
                    if (
                        max(time_steps) - time_max < time_av
                        and len(files_in_level) != 0
                    ):
                        for file in files_in_level:
                            if file not in averaged_files:
                                files_to_average.append(file)
                                averaged_files.append(file)
                        time_averaged_file = os.path.join(
                            os.path.join(
                                model_processed_output_folder_daily, species[i]
                            ),
                            "c_"
                            + level
                            + "_"
                            + str(time_min)
                            + "-"
                            + str(max(time_steps))
                            + "-tavg.grd",
                        )
                        tavg_intervals[-1] = (
                            str(time_min) + "-" + str(max(time_steps)) + "-tavg"
                        )
                    time_average(files_to_average, time_averaged_file)
                    files_to_average = []
            if time_av == 0:
                break
            else:
                time_min = time_max + 1
                time_max = time_min + time_av - 1
                continue
    return tavg_intervals


def probabilistic_output(model):
    def ecdf(index):
        specie = index[1]
        level = index[2]
        time_step = index[3]
        ex_prob = index[0]
        quantile = 1 - ex_prob
        output_files = []
        for day in days:
            try:
                file_name = (
                    "c_"
                    + "{:03d}".format(int(level))
                    + "_"
                    + "{:06d}".format(int(time_step))
                    + ".grd"
                )
            except ValueError:
                file_name = (
                    "c_" + "{:03d}".format(int(level)) + "_" + time_step + ".grd"
                )
            output_folder = os.path.join(model_processed_output_folder, day, specie)
            output_files.append(os.path.join(output_folder, file_name))
        try:
            ecdf_output_file = os.path.join(
                ecdf_folder,
                str(ex_prob),
                specie,
                "c_"
                + "{:03d}".format(int(level))
                + "_"
                + "{:06d}".format(int(time_step))
                + ".grd",
            )
        except ValueError:
            ecdf_output_file = os.path.join(
                ecdf_folder,
                str(ex_prob),
                specie,
                "c_" + "{:03d}".format(int(level)) + "_" + time_step + ".grd",
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
                for i in range(1, nz + 1):
                    if time_steps[0] == "all":
                        for j in range(0, n_time_steps + 1):
                            indexes.append([probability, specie, i, j])
                    else:
                        for time_step in time_steps:
                            indexes.append([probability, specie, i, time_step])
                    if len(tavg_intervals) > 0:
                        for k in range(0, len(tavg_intervals)):
                            indexes_tavg.append(
                                [probability, specie, i, tavg_intervals[k]]
                            )
            else:
                for level in levels:
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
        with open(input) as input_file:
            if output_format == "grd":
                Z = np.loadtxt(input, skiprows=5)
            X = np.linspace(x0, xf, num=nx)
            Y = np.linspace(y0, yf, num=ny)
            n_levels = 10
            dc = (max_con - min_con) / n_levels
            levels = np.arange(min_con + 0.0000001, max_con, dc)
        fig, ax = plt.subplots(figsize=(6, 5), dpi=600)
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
            ax.set_title("Gas concentration [ppm]")
        else:
            ax.set_title("Gas concentration [kg m$\mathregular{^{-3}}$]")
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
    try:
        os.mkdir(graphical_outputs)
    except FileExistsError:
        print("Folder " + graphical_outputs + " already exists")
        list_temp = os.listdir(graphical_outputs)
        for item in list_temp:
            list_temp_2 = os.listdir(os.path.join(graphical_outputs, item))
            for item_2 in list_temp_2:
                try:
                    rmtree(os.path.join(os.path.join(graphical_outputs, item), item_2))
                except FileNotFoundError:
                    print(
                        "Unable to remove "
                        + item_2
                        + " in "
                        + os.path.join(graphical_outputs, item)
                    )
    try:
        os.mkdir(graphical_outputs_simulations)
    except FileExistsError:
        print("Folder " + graphical_outputs_simulations + " already exists")
    try:
        os.mkdir(graphical_outputs_ecdf)
    except FileExistsError:
        print("Folder " + graphical_outputs_ecdf + " already exists")

    files_to_plot = []
    output_files = []
    if plot:
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
                output_file_name = files_list[i].split(".")[0]
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
                            if file_time_step == "{:06d}".format(int(time_step)):
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
                        tavg_output_file_name = file.split(os.sep)[-1].split(".")[0]
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
                            if file_level == "{:03d}".format(int(level)):
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
                                if file_time_step == "{:06d}".format(
                                    int(time_step)
                                ) and file_level == "{:03d}".format(int(level)):
                                    files_to_plot.append(file)
                                    output_files.append(
                                        os.path.join(
                                            graphical_outputs_daily,
                                            file_specie,
                                            output_file_name,
                                        )
                                    )
                    for level in levels:
                        if "tavg" in file_time_step and file_level == "{:03d}".format(
                            int(level)
                        ):
                            files_to_plot.append(file)
                            tavg_output_file_name = file.split(os.sep)[-1].split(".")[0]
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
                    output_file_name = file.split(".")[0]
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
                            tavg_output_file_name = file.split(os.sep)[-1].split(".")[0]
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
                                if file_level == "{:03d}".format(int(level)):
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
                                    ) and file_level == "{:03d}".format(int(level)):
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
                                    and file_level == "{:03d}".format(int(level))
                                ):
                                    files_to_plot.append(file_path)
                                    tavg_output_file_name = file.split(os.sep)[
                                        -1
                                    ].split(".")[0]
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
    days_plot,
    species,
    original_specie,
    exceedance_probabilities,
    nproc,
    convert,
    models,
    merge_outputs,
    units,
    time_av,
    min_con,
    max_con,
    output_format,
    plot_topography_layer,
    dz_lines_res,
    plot_resolution,
) = read_arguments()

try:
    max_number_processes = int(os.environ["SLURM_NPROCS"])
except KeyError:
    max_number_processes = int(nproc)

(
    disgas_outputs,
    disgas_original_output_folder,
    disgas_processed_output_folder,
    ecdf_folder_name,
    disgas_ecdf,
    twodee_outputs,
    twodee_original_output_folder,
    twodee_processed_output_folder,
    twodee_ecdf,
    models_to_elaborate,
    twodee_output_time_step,
) = folder_structure()

days, days_to_plot = extract_days()

species_properties = gas_properties()

for model in models_to_elaborate:
    x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps, dt, output_levels = domain(model)
    for day in days:
        tavg_intervals = elaborate_day(day, model)
    if plot_ex_prob:
        probabilistic_output(model)
    save_plots(model, min_con, max_con)
