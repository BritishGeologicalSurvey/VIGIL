#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import random
import sys, os
from pathos.multiprocessing import ThreadingPool
from io import StringIO
import argparse
from shutil import rmtree

def read_arguments():
    parser = argparse.ArgumentParser(description='Input data')
    parser.add_argument('-P', '--plot', default='False',
                        help='True: Produce plots of the solutions. False: Do not produce plots')
    parser.add_argument('-PE', '--plot_ex_prob', default='False',
                        help='True: Produce plots of the specified exceedance probabilities. False: Do not produce plots')
    parser.add_argument('-EX', '--ex_prob', nargs='+', default=[],
                        help='List of exceedence probabilities to be used for graphical output')
    parser.add_argument('-T', '--time_steps', nargs='+', default=[],
                        help='List of time steps to plot (integer >= 0). Type all to plot all the time steps')
    parser.add_argument('-L', '--levels', nargs='+', default=[],
                        help='List of vertical levels (integer >= 1) to plot. Type all to plot all the levels')
    parser.add_argument('-D', '--days_plot', nargs='+', default=[],
                        help='List of days to plot (YYYYMMDD). Type all to plot all the days')
    parser.add_argument('-C', '--convert', default='False', help='If True, convert output concentration into other species listed with the command -S (--species)')
    parser.add_argument('-S', '--species', nargs='+', default=[], help='List of gas species (e.g. CO2)')
    parser.add_argument('-N', '--nproc', default=1, help='Maximum number of allowed simultaneous processes')
    parser.add_argument('-M', '--models', default='all', help='Model outputs to post-process. Options: disgas, twodee, all')
    parser.add_argument('-MO', '--merge_outputs', default='False', help='Merge Twodee and Disgas outputs (true or false)')
    args = parser.parse_args()
    plot = args.plot
    plot_ex_prob = args.plot_ex_prob
    ex_prob = args.ex_prob
    time_steps = args.time_steps
    levels = args.levels
    days_plot = args.days_plot
    species = args.species
    nproc = args.nproc
    convert = args.convert
    models = args.models
    merge_outputs = args.merge_outputs
    if plot.lower() == 'true':
        plot = True
        if len(days_plot) == 0:
            print('ERROR. Please specify at least one day to plot when --plot==True')
            sys.exit()
    elif plot.lower() == 'false':
        plot = False
    else:
        print('ERROR. Wrong value for variable -P --plot')
        sys.exit()
    if plot_ex_prob.lower() == 'true':
        plot_ex_prob = True
        if len(ex_prob) == 0:
            print('ERROR. Please specify at least one exceedance probability to plot when --plot_ex_prob==True')
            sys.exit()
    elif plot_ex_prob.lower() == 'false':
        plot_ex_prob = False
    else:
        print('ERROR. Wrong value for variable -PE --plot_ex_prob')
        sys.exit()
    if plot or plot_ex_prob:
        if len(time_steps) == 0:
            print('ERROR. Please specify at least one time step to plot')
            sys.exit()
        if len(levels) == 0:
            print('ERROR. Please specify at least one level to plot')
            sys.exit()
    if convert.lower() == 'true':
        convert = True
        if len(species) == 0:
            print('ERROR. Please specify at least one gas specie name when --convert=True')
            sys.exit()
    elif convert.lower() == 'false':
        convert = False
        if len(species) != 0:
            del species[:]
        species.append('original_specie')
    else:
        print('ERROR. Wrong value for variable -C --convert')
        sys.exit()
    exceedance_probabilities = []
    for prob in ex_prob:
        exceedance_probabilities.append(float(prob))
    if models.lower() != 'disgas' and models.lower() != 'twodee' and models.lower() != 'all':
        print('ERROR. Wrong value for variable -M --models')
        sys.exit()
    if merge_outputs.lower() == 'true':
        merge_outputs = True
    elif merge_outputs.lower() == 'false':
        merge_outputs = False
    else:
        print('ERROR. Wrong value for variable -MO --merge_outputs')
        sys.exit()

    return plot, plot_ex_prob, time_steps, levels, days_plot, species, exceedance_probabilities, nproc, convert, models, merge_outputs

def folder_structure():
    original_output_folder_name = 'simulations'
    post_processing = 'post_processing'
    processed_output_folder_name = original_output_folder_name + '_processed'
    ecdf_folder_name = 'output_ecdf'
    disgas_outputs = os.path.join(root, post_processing, 'disgas')
    twodee_outputs = os.path.join(root, post_processing, 'twodee')
    disgas_original_output_folder = os.path.join(root, original_output_folder_name, 'disgas')
    twodee_original_output_folder = os.path.join(root, original_output_folder_name, 'twodee')
    disgas_processed_output_folder = os.path.join(disgas_outputs, processed_output_folder_name)
    disgas_ecdf = os.path.join(disgas_outputs, ecdf_folder_name)
    twodee_processed_output_folder = os.path.join(twodee_outputs, processed_output_folder_name)
    twodee_ecdf = os.path.join(twodee_outputs, ecdf_folder_name)
    os.mkdir(post_processing)
    if models == 'disgas' or models == 'all':
        try:
            os.mkdir(disgas_outputs)
        except FileExistsError:
            print('Folder ' + disgas_outputs + ' already exists')
        try:
            os.mkdir(disgas_processed_output_folder)
        except:
            print('Folder ' + disgas_processed_output_folder + ' already exists')
            list_temp = os.listdir(disgas_processed_output_folder)
            for item in list_temp:
                try:
                    rmtree(os.path.join(disgas_processed_output_folder, item), ignore_errors=True)
                except:
                    print('Unable to remove ' + item + ' in ' + disgas_processed_output_folder)
        try:
            os.mkdir(disgas_ecdf)
        except:
            print('Folder ' + disgas_ecdf + ' already exists')
            list_temp = os.listdir(disgas_ecdf)
            for item in list_temp:
                try:
                    rmtree(os.path.join(disgas_ecdf, item), ignore_errors=True)
                except:
                    print('Unable to remove ' + item + ' in ' + disgas_ecdf)
    if models == 'twodee' or models == 'all':
        try:
            os.mkdir(twodee_outputs)
        except FileExistsError:
            print('Folder ' + twodee_outputs + ' already exists')
        try:
            os.mkdir(twodee_processed_output_folder)
        except:
            print('Folder ' + twodee_processed_output_folder + ' already exists')
            list_temp = os.listdir(twodee_processed_output_folder)
            for item in list_temp:
                try:
                    rmtree(os.path.join(twodee_processed_output_folder, item), ignore_errors=True)
                except:
                    print('Unable to remove ' + item + ' in ' + twodee_processed_output_folder)
        try:
            os.mkdir(twodee_ecdf)
        except:
            print('Folder ' + twodee_ecdf + ' already exists')
            list_temp = os.listdir(twodee_ecdf)
            for item in list_temp:
                try:
                    rmtree(os.path.join(twodee_ecdf, item), ignore_errors=True)
                except:
                    print('Unable to remove ' + item + ' in ' + twodee_ecdf)
    if models == 'all':
        models_to_elaborate = ['disgas','twodee']
    elif models == 'disgas':
        models_to_elaborate = ['disgas']
    else:
        models_to_elaborate = ['twodee']
    # Read the output time interval from the twodee input file
    twodee_input_file = os.path.join(root,'twodee.inp')
    twodee_output_time_step = 0
    with open(twodee_input_file, 'r') as twodee_file:
        for line in twodee_file:
            if 'OUTPUT_INTERVAL_(SEC)' in line:
                twodee_output_time_step = float(line.split('=')[1])
    if twodee_output_time_step == 0:
        print('Unable to read the Twodee output time step')
        sys.exit()

    return disgas_outputs, disgas_original_output_folder, disgas_processed_output_folder, ecdf_folder_name, disgas_ecdf, \
           twodee_outputs, twodee_original_output_folder, twodee_processed_output_folder, twodee_ecdf, models_to_elaborate, twodee_output_time_step

def gas_properties():
    def extract_gas_properties(specie):
        data = pd.read_csv(gas_properties_file, error_bad_lines=False)
        x = np.sort(data[specie + '/H2O'])
        y = np.sort(data[specie])
        molar_weight = list(y)[0]
        list_x = list(x)
        samples = (random.sample(list_x, 1))
        molar_ratio = samples[0]
        print('The molar ratio ' + specie + '/H2O is', molar_ratio)
        return molar_ratio, molar_weight

    gas_properties_file = os.path.join(root, 'gas_properties.csv')
    try:
        open(gas_properties_file, 'r')
    except:
        print('File ' + gas_properties_file + ' not present')
        sys.exit()
    molar_ratios = []
    molar_weights = []
    for specie in species:
        molar_ratio, molar_weight = extract_gas_properties(specie)
        molar_ratios.append(molar_ratio)
        molar_weights.append(molar_weight)
    species_properties = []
    for i in range(0, len(species)):
        gas_specie = {}
        gas_specie['specie_name'] = species[i]
        gas_specie['molar_ratio'] = molar_ratios[i]
        gas_specie['molar_weight'] = molar_weights[i]
        species_properties.append(gas_specie)
    return species_properties

def domain(model):
    import re
    if model == 'disgas':
        with open(file='disgas.inp') as input_file:
            for record in input_file:
                try:
                    record_splitted = record.split('=')
                    temp = record_splitted[1].split('(')
                    if 'SIMULATION_INTERVAL_(SEC)' in record_splitted[0]:
                        tot_time = int(temp[0])
                    elif 'NX' in record_splitted[0]:
                        nx = int(temp[0])
                    elif 'NY' in record_splitted[0]:
                        ny = int(temp[0])
                    elif 'NZ' in record_splitted[0]:
                        nz = int(temp[0])
                    elif 'OUTPUT_INTERVAL_(SEC)' in record_splitted[0]:
                        dt = int(temp[0])
                    elif 'DX_(M)' in record_splitted[0]:
                        dx = float(temp[0])
                    elif 'DY_(M)' in record_splitted[0]:
                        dy = float(temp[0])
                    elif 'X_ORIGIN_(UTM_M)' in record_splitted[0]:
                        x0 = float(temp[0])
                    elif 'Y_ORIGIN_(UTM_M)' in record_splitted[0]:
                        y0 = float(temp[0])
                except:
                    continue
    else:
        with open(file='twodee.inp') as input_file:
            for record in input_file:
                try:
                    record_splitted = record.split('=')
                    temp = record_splitted[1].split('(')
                    if 'SIMULATION_INTERVAL_(SEC)' in record_splitted[0]:
                        tot_time = int(temp[0])
                    elif 'NX' in record_splitted[0]:
                        nx = int(temp[0])
                    elif 'NY' in record_splitted[0]:
                        ny = int(temp[0])
                    elif 'OUTPUT_INTERVAL_(SEC)' in record_splitted[0]:
                        dt = int(temp[0])
                    elif 'DX_(M)' in record_splitted[0]:
                        dx = float(temp[0])
                    elif 'DY_(M)' in record_splitted[0]:
                        dy = float(temp[0])
                    elif 'X_ORIGIN_(UTM_M)' in record_splitted[0]:
                        x0 = float(temp[0])
                    elif 'Y_ORIGIN_(UTM_M)' in record_splitted[0]:
                        y0 = float(temp[0])
                    elif 'HEIGHTS_(M)' in record_splitted[0]:
                        heights = temp[0]
                        extracted_heights = re.findall('\d+\.\d+', heights) # This extracts decimal numbers only!
                        nz = len(extracted_heights)
                except:
                    continue
    yf = y0 + ny * dy
    xf = x0 + nx * dx
    n_time_steps = int(tot_time / dt)
    return x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps

def extract_days():
    days = []
    days_list_path = os.path.join(root, 'days_list.txt')
    days_to_plot = []
    with open(days_list_path, 'r') as days_list:
        for line in days_list:
            day_temp = line.split(' ')[0]
            day_temp = day_temp.split('-')
            day = day_temp[0] + day_temp[1] + day_temp[2]
            days.append(day)
            if days_plot[0] == 'all':
                days_to_plot.append(day)
            else:
                for day_to_plot in days_plot:
                    if day_to_plot == day:
                        days_to_plot.append(day_to_plot)
    return days, days_to_plot

def converter(input_file, outname, specie_input, model):
    Z = np.loadtxt(input_file, skiprows=5)
    if model == 'twodee':
        Z = np.divide(Z, 1000) #convert ppm to kg/s
    if specie_input == 'original_specie':
        Z_converted = np.reshape(Z, [nx, ny])
        np.savetxt(outname, Z_converted, fmt='%.2e')
    else:
        for specie in species_properties:
            if specie['specie_name'] == specie_input:
                mol_ratio = specie['molar_ratio']
                molar_weight = specie['molar_weight']
        Z_converted = np.multiply(Z, mol_ratio)
        Z_converted = [(Z_converted / molar_weight) / (44.64 * 1000000000)]
        Z_converted = np.reshape(Z_converted, [nx, ny])
        np.savetxt(outname, Z_converted, fmt='%.2e')

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step: np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red', 'green', 'blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j, i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1],), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector
    return matplotlib.colors.LinearSegmentedColormap('colormap', cdict, 1024)

def elaborate_day(day_input, model):
    if model == 'disgas':
        model_output_folder = os.path.join(disgas_original_output_folder, day_input, 'outfiles')
        model_processed_output_folder_daily = os.path.join(disgas_processed_output_folder, day_input)
    else:
        model_output_folder = os.path.join(twodee_original_output_folder, day_input, 'outfiles')
        model_processed_output_folder_daily = os.path.join(twodee_processed_output_folder, day_input)
    try:
        os.mkdir(model_processed_output_folder_daily)
    except:
        print('Folder ' + model_processed_output_folder_daily + ' already exists')
    for specie in species:
        model_processed_output_folder_specie = os.path.join(model_processed_output_folder_daily, specie)
        try:
            os.mkdir(model_processed_output_folder_specie)
        except FileExistsError:
            print('Folder ' + model_processed_output_folder_specie + ' already exists')
        except PermissionError: #retry
            try:
                os.mkdir(model_processed_output_folder_specie)
            except FileExistsError:
                print('Folder ' + model_processed_output_folder_specie + ' already exists')
    files_list_temp = os.listdir(model_output_folder)
    files_list_path = []
    files_list = []
    models = []
    for file in files_list_temp:
        if file[0:2] == 'c_':
            files_list.append(file)
            for specie in species:
                files_list_path.append(os.path.join(model_output_folder, file))#
                models.append(model)
    converted_files = []
    outnames = []
    species_list = []
    for file in files_list:
        for specie in species:
            species_list.append(specie)
            if model == 'twodee':
                file_name_splitted = file.split('_')
                file_level = file_name_splitted[1]
                file_time_step = file_name_splitted[2].split('.')[0]
                file_level = "{:03d}".format(int(int(file_level.split('cm')[0]) / 100))
                file_time_step = "{:06d}".format(int((int(file_time_step) / twodee_output_time_step)))
                file = 'c_' + file_level + '_' + file_time_step + '.grd'
            converted_file = file
            converted_files.append(converted_file)
            outnames.append(os.path.join(os.path.join(model_processed_output_folder_daily, specie), converted_file))
    n_elaborated_files = 0
    while n_elaborated_files < len(files_list_path):
        start = n_elaborated_files
        end = n_elaborated_files + max_number_processes
        if end > len(files_list_path):
            end = len(files_list_path)
        try:
            pool_files = ThreadingPool(max_number_processes)
            pool_files.map(converter,files_list_path[start:end], outnames[start:end], species_list[start:end], models[start:end])
        except:
            print('Unable to convert files')
            sys.exit()
        n_elaborated_files = end
        if n_elaborated_files == len(files_list_path):
            break

def probabilistic_output(model):
    def ecdf(index):
        specie = index[1]
        level = index[2]
        time_step = index[3]
        quantile = index[0]
        output_files = []
        for day in days:
            file_name = 'c_' + "{:03d}".format(int(level)) + '_' + "{:06d}".format(int(time_step)) + '.grd'
            output_folder = os.path.join(model_processed_output_folder, day, specie)
            output_files.append(os.path.join(output_folder, file_name))
        ecdf_output_file = os.path.join(ecdf_folder, str(quantile), specie, 'c_' + "{:03d}".format(int(level)) + '_' + "{:06d}".format(int(time_step)) + '.grd')
        quantile = 1 - quantile
        output_quantile = np.zeros((nx, ny))
        c_arrays = []
        files_not_available = []
        for file in output_files:
            try:
                input_file = open(file)
            except:
                print('File ' + file + ' not found')
                files_not_available.append(file)
                continue
            records = []
            for line in input_file:
                records.append(line.split(' '))
            c_arrays.append(records)
        for file in files_not_available:
            output_files.remove(file)
        for i in range(0, nx):
            for j in range(0, ny):
                c_list = []
                for k in range(0, len(output_files)):
                    c_list.append(float(c_arrays[k][i][j]))
                output_quantile[i, j] = np.quantile(c_list, q=quantile)
        np.savetxt(ecdf_output_file, output_quantile, fmt='%.2e')

    if model == 'disgas':
        ecdf_folder = disgas_ecdf
        model_processed_output_folder = disgas_processed_output_folder
    else:
        ecdf_folder = twodee_ecdf
        model_processed_output_folder = twodee_processed_output_folder
    for probability in exceedance_probabilities:
        prob_folder = os.path.join(ecdf_folder, str(probability))
        try:
            os.mkdir(prob_folder)
        except:
            print('Folder ' + prob_folder + ' already exists')
        for specie in species:
            specie_folder = os.path.join(ecdf_folder, prob_folder, specie)
            try:
                os.mkdir(specie_folder)
            except:
                print('Folder ' + specie_folder + ' already exists')
    indexes = []
    pools_ecdfs = []
    n_pool = 0
    for probability in exceedance_probabilities:
        for specie in species:
            pools_ecdfs.append(n_pool)
            if levels[0] == 'all':
                for i in range(1, nz + 1):
                    if time_steps[0] == 'all':
                        for j in range(0, n_time_steps + 1):  #for j in range(1, n_time_steps):
                            indexes.append([probability, specie, i, j])
                    else:
                        for time_step in time_steps:
                            indexes.append([probability, specie, i, time_step])
            else:
                for level in levels:
                    if time_steps[0] == 'all':
                        for j in range(0, n_time_steps + 1):  #for j in range(1, n_time_steps):
                            indexes.append([probability, specie, level, j])
                    else:
                        for time_step in time_steps:
                            indexes.append([probability, specie, level, time_step])
            n_pool += 1
    n_pool = 0
    for probability in exceedance_probabilities:
        for specie in species:
            n_completed_processes = 0
            while n_completed_processes <= len(indexes):
                start = n_completed_processes
                end = n_completed_processes + max_number_processes
                if end > len(indexes):
                    end = len(indexes)
                try:
                    pools_ecdfs[n_pool] = ThreadingPool(max_number_processes)
                    pools_ecdfs[n_pool].map(ecdf, indexes[start:end])
                except:
                    print('Unable to elaborate days')
                    sys.exit()
                n_completed_processes = end
                if n_completed_processes == len(indexes):
                    break
            n_pool += 1

def save_plots(model):
    import re

    def plot_file(input,output):
        from matplotlib.ticker import FormatStrFormatter
        with open(input) as input_file:
            Z = [[float(record) for record in line.split(' ')] for line in input_file]
            fig = plt.figure(figsize=(8, 8))
            plt.title('Gas concentration [kg m$\mathregular{^{-3}}$]')
            X = np.arange(x0, xf, dx)
            Y = np.arange(y0, yf, dy)
            plt.contourf(X, Y, Z)
            plt.xlabel('X_UTM [m]')
            plt.ylabel('Y_UTM [m]')
            plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            cax = fig.add_axes([0.1, 0.03, 0.8, 0.01])
            plt.colorbar(cax = cax, cmap = dark_jet, orientation='horizontal', format='%.1e')
            image_buffer = StringIO()
            plt.savefig(output)
            image_buffer.close()
            plt.close(fig)
        input_file.close()

    if model == 'disgas':
        model_outputs = disgas_outputs
        model_processed_output_folder = disgas_processed_output_folder
        ecdf_outputs = disgas_ecdf
    else:
        model_outputs = twodee_outputs
        model_processed_output_folder = twodee_processed_output_folder
        ecdf_outputs = twodee_ecdf
    dark_jet = cmap_map(lambda x: x * 0.75, matplotlib.cm.jet)
    graphical_outputs = (os.path.join(model_outputs, 'graphical_outputs'))
    graphical_outputs_simulations = (os.path.join(graphical_outputs, 'simulations'))
    graphical_outputs_ecdf = (os.path.join(graphical_outputs, 'ecdf'))
    try:
        os.mkdir(graphical_outputs)
    except FileExistsError:
        print('Folder ' + graphical_outputs + ' already exists')
        list_temp = os.listdir(graphical_outputs)
        for item in list_temp:
            list_temp_2 = os.listdir(os.path.join(graphical_outputs, item))
            for item_2 in list_temp_2:
                try:
                    rmtree(os.path.join(os.path.join(graphical_outputs,item),item_2))
                except:
                    print('Unable to remove ' + item_2 + ' in ' + os.path.join(graphical_outputs,item))
    try:
        os.mkdir(graphical_outputs_simulations)
    except FileExistsError:
        print('Folder ' + graphical_outputs_simulations + ' already exists')
    try:
        os.mkdir(graphical_outputs_ecdf)
    except FileExistsError:
        print('Folder ' + graphical_outputs_ecdf + ' already exists')

    files_to_plot = []
    output_files = []
    if plot:
        for day in days_to_plot:
            graphical_outputs_daily = os.path.join(graphical_outputs_simulations, day)
            try:
                os.mkdir(graphical_outputs_daily)
            except FileExistsError:
                print('Folder ' + graphical_outputs_daily + ' already exists')
            model_processed_output_folder_daily = os.path.join(model_processed_output_folder, day)
            model_processed_output_folder_species = []
            for specie in species:
                model_processed_output_folder_species.append(os.path.join(model_processed_output_folder_daily, specie))
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_daily, specie))
                except FileExistsError:
                    print('Folder ' + os.path.join(graphical_outputs_daily, specie) + ' already exists')
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
                file_specie = re.sub('\W+', '', file_specie[0])
                file_name_splitted = files_list[i].split('_')
                file_level = file_name_splitted[1]
                file_time_step = file_name_splitted[2].split('.')[0]
                output_file_name = files_list[i].split('.')[0]
                output_file_name += '.png'
                if levels[0] == 'all':
                    if time_steps[0] == 'all':
                        files_to_plot.append(file)
                        output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name))
                    else:
                        for time_step in time_steps:
                            if file_time_step == "{:06d}".format(int(time_step)):
                                files_to_plot.append(file)
                                output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name))
                else:
                    if time_steps[0] == 'all':
                        for level in levels:
                            if file_level == "{:03d}".format(int(level)):
                                files_to_plot.append(file)
                                output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name))
                    else:
                        for level in levels:
                            for time_step in time_steps:
                                if file_time_step == "{:06d}".format(int(time_step)) and file_level == "{:03d}".format(int(level)):
                                    files_to_plot.append(file)
                                    output_files.append(os.path.join(graphical_outputs_daily, file_specie, output_file_name))
                i += 1
    if plot_ex_prob:
        for probability in exceedance_probabilities:
            try:
                os.mkdir(os.path.join(graphical_outputs_ecdf,str(probability)))
            except FileExistsError:
                print('Folder ' + os.path.join(graphical_outputs_ecdf, str(probability)) + ' already exists')
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_ecdf, str(probability), specie))
                except FileExistsError:
                    print(
                        'Folder ' + os.path.join(graphical_outputs_ecdf, str(probability), specie) + ' already exists')
                files_list = os.listdir(os.path.join(ecdf_outputs, str(probability), specie))
                for file in files_list:
                    file_path = os.path.join(ecdf_outputs, str(probability), specie, file)
                    file_name_splitted = file.split('_')
                    file_level = file_name_splitted[1]
                    file_time_step = file_name_splitted[2].split('.')[0]
                    output_file_name = file.split('.')[0]
                    output_file_name += '.png'
                    if levels[0] == 'all':
                        if time_steps[0] == 'all':
                            files_to_plot.append(file_path)
                            output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), specie, output_file_name))
                        else:
                            for time_step in time_steps:
                                if file_time_step == "{:06d}".format(int(time_step)):
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), specie, output_file_name))
                    else:
                        if time_steps[0] == 'all':
                            for level in levels:
                                if file_level == "{:03d}".format(int(level)):
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), specie, output_file_name))
                        else:
                            for level in levels:
                                for time_step in time_steps:
                                    if file_time_step == "{:06d}".format(int(time_step)) and file_level == "{:03d}".format(int(level)):
                                        files_to_plot.append(file_path)
                                        output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), specie, output_file_name))
    if len(files_to_plot) == 0:
        print('No files to plot')
    else:
        i = 0
        for file_to_plot in files_to_plot:
            print('plotting ' + file_to_plot)
            plot_file(file_to_plot,output_files[i])
            i += 1

root = os.getcwd()

plot, plot_ex_prob, time_steps, levels, days_plot, species, exceedance_probabilities, nproc, convert, models, merge_outputs = read_arguments()

try:
    max_number_processes = int(os.environ["SLURM_NPROCS"])
except:
    max_number_processes = int(nproc)

disgas_outputs, disgas_original_output_folder, disgas_processed_output_folder, ecdf_folder_name, disgas_ecdf, twodee_outputs, \
twodee_original_output_folder, twodee_processed_output_folder, twodee_ecdf, models_to_elaborate, twodee_output_time_step = folder_structure()

days, days_to_plot = extract_days()

if convert:
    species_properties = gas_properties()
for model in models_to_elaborate:
    x0, xf, y0, yf, nx, ny, nz, dx, dy, n_time_steps = domain(model)
    for day in days:
        elaborate_day(day, model)
    probabilistic_output(model)
    save_plots(model)
