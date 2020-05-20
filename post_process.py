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
    if plot == 'True':
        plot = True
        if len(days_plot) == 0:
            print('ERROR. Please specify at least one day to plot when --plot==True')
            sys.exit()
    elif plot == 'False':
        plot = False
    else:
        print('ERROR. Wrong value for variable -P --plot')
        sys.exit()
    if plot_ex_prob == 'True':
        plot_ex_prob = True
        if len(ex_prob) == 0:
            print('ERROR. Please specify at least one exceedance probability to plot when --plot_ex_prob==True')
            sys.exit()
    elif plot_ex_prob == 'False':
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
    if convert == 'True':
        convert = True
        if len(species) == 0:
            print('ERROR. Please specify at least one gas specie name when --convert=True')
            sys.exit()
    elif convert == 'False':
        convert = False
        if len(species) != 0:
            species = []
    else:
        print('ERROR. Wrong value for variable -C --convert')
        sys.exit()
    exceedance_probabilities = []
    for prob in ex_prob:
        exceedance_probabilities.append(float(prob))
    return plot, plot_ex_prob, time_steps, levels, days_plot, species, exceedance_probabilities, nproc, convert

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

    gas_properties_file = os.path.join(cwd, 'gas_properties.csv')
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

def domain():
    with open(file='disgas.inp') as input_file:
        for record in input_file:
            command = record.split(' = ')
            if command[0] == '  SIMULATION_INTERVAL_(SEC)':
                tot_time = int(command[1])
            elif command[0] == '  NX    ':
                nx = int(command[1])
            elif command[0] == '  NY    ':
                ny = int(command[1])
            elif command[0] == '  NZ    ':
                nz = int(command[1])
            elif command[0] == '  OUTPUT_INTERVAL_(SEC)':
                dt = command[1]
                dt = int(dt.split(' ')[0])
            elif command[0] == '  DX_(M)':
                dx = float(command[1])
            elif command[0] == '  DY_(M)':
                dy = float(command[1])
            elif command[0] == '  X_ORIGIN_(UTM_M)':
                x0 = float(command[1])
            elif command[0] == '  Y_ORIGIN_(UTM_M)':
                y0 = float(command[1])
    yf = y0 + (ny - 1) * dy
    xf = x0 + (nx - 1) * dx
    n_time_steps = int(tot_time / dt)
    return x0, xf, y0, yf, nx, ny, nz, n_time_steps

def extract_days():
    days = []
    days_list_path = os.path.join(cwd, 'days_list.txt')
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

def converter(input_file, outname, specie_input):
    Z = np.loadtxt(input_file, skiprows=5)
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

def elaborate_day(day_input):
    disgas_output_folder = os.path.join(disgas_original_output_folder, day_input, 'outfiles')
    disgas_converted_output_folder_daily = os.path.join(disgas_converted_output_folder, day_input)
    try:
        os.mkdir(disgas_converted_output_folder_daily)
    except:
        print('Folder ' + disgas_converted_output_folder_daily + ' already exists')
    for specie in species:
        disgas_converted_output_folder_specie = os.path.join(disgas_converted_output_folder_daily, specie)
        try:
            os.mkdir(disgas_converted_output_folder_specie)
        except FileExistsError:
            print('Folder ' + disgas_converted_output_folder_specie + ' already exists')
        except PermissionError: #retry
            try:
                os.mkdir(disgas_converted_output_folder_specie)
            except FileExistsError:
                print('Folder ' + disgas_converted_output_folder_specie + ' already exists')
    files_list_temp = os.listdir(disgas_output_folder)
    files_list_path = []
    files_list = []
    for file in files_list_temp:
        if file[0:2] == 'c_':
            files_list.append(file)
            for specie in species:
                files_list_path.append(os.path.join(disgas_output_folder, file))
    converted_files = []
    outnames = []
    species_list = []
    for file in files_list:
        for specie in species:
            species_list.append(specie)
            converted_file = specie + '_' + file
            converted_files.append(converted_file)
            outnames.append(os.path.join(os.path.join(disgas_converted_output_folder_daily, specie), converted_file))
    n_elaborated_files = 0
    while n_elaborated_files < len(files_list_path):
        start = n_elaborated_files
        end = n_elaborated_files + max_number_processes
        if end > len(files_list_path):
            end = len(files_list_path)
            try:
                pool_files = ThreadingPool(max_number_processes)
                pool_files.map(converter,files_list_path[start:end], outnames[start:end], species_list[start:end])
            except:
                print('Unable to convert files')
                sys.exit()
        n_elaborated_files = end
        if n_elaborated_files == len(files_list_path):
            break

def ecdf(index):
    specie = index[1]
    level = index[2]
    time_step = index[3]
    quantile = index[0]
    output_files = []
    for day in days:
        file_name = specie + '_c_' + "{:03d}".format(int(level)) + '_' + "{:06d}".format(int(time_step)) + '.grd'
        output_folder = os.path.join(disgas_converted_output_folder, day, specie)
        output_files.append(os.path.join(output_folder, file_name))
    ecdf_output_file = os.path.join(ecdf_folder_name, str(quantile), specie,
                                    specie + '_c_' + "{:03d}".format(int(level)) + '_' + "{:06d}".format(int(time_step)) + '.grd')
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

def save_plots():
    dark_jet = cmap_map(lambda x: x * 0.75, matplotlib.cm.jet)
    graphical_outputs = (os.path.join(cwd, 'graphical_outputs'))
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

    def plot_file(input,output):
        with open(input) as input_file:
            Z = [[float(record) for record in line.split(' ')] for line in input_file]
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_subplot(111)
            ax1.imshow(Z, extent=[x0, xf, y0, yf], cmap=dark_jet,aspect='auto')
            image_buffer = StringIO()
            fig.savefig(output)
            image_buffer.close()
            plt.close(fig)
        input_file.close()

    files_to_plot = []
    output_files = []
    if plot:
        for day in days_to_plot:
            disgas_converted_output_folder_daily = os.path.join(disgas_converted_output_folder, day)
            disgas_converted_output_folder_species = []
            for specie in species:
                disgas_converted_output_folder_species.append(os.path.join(disgas_converted_output_folder_daily, specie))
            graphical_outputs_daily = os.path.join(graphical_outputs_simulations, day)
            try:
                os.mkdir(graphical_outputs_daily)
            except FileExistsError:
                print('Folder ' + graphical_outputs_daily + ' already exists')
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_daily,specie))
                except FileExistsError:
                    print('Folder ' + os.path.join(graphical_outputs_daily,specie) + ' already exists')
            files_list_path = []
            files_list = []
            for folder in disgas_converted_output_folder_species:
                files_list_temp = os.listdir(folder)
                for file in files_list_temp:
                    files_list.append(file)
                    files_list_path.append(os.path.join(folder,file))
            i = 0
            for file in files_list_path:
                file_name_splitted = files_list[i].split('_')
                file_specie = file_name_splitted[0]
                file_level = file_name_splitted[2]
                file_time_step = file_name_splitted[3].split('.')[0]
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
                print('Folder ' + os.path.join(disgas_ecdf, str(probability)) + ' already exists')
            for specie in species:
                try:
                    os.mkdir(os.path.join(graphical_outputs_ecdf,str(probability),specie))
                except FileExistsError:
                    print('Folder ' + os.path.join(graphical_outputs_ecdf,str(probability),specie) + ' already exists')
                files_list = os.listdir(os.path.join(disgas_ecdf, str(probability), specie))
                for file in files_list:
                    file_path = os.path.join(disgas_ecdf, str(probability), specie,file)
                    file_name_splitted = file.split('_')
                    file_specie = file_name_splitted[0]
                    file_level = file_name_splitted[2]
                    file_time_step = file_name_splitted[3].split('.')[0]
                    output_file_name = file.split('.')[0]
                    output_file_name += '.png'
                    if levels[0] == 'all':
                        if time_steps[0] == 'all':
                            files_to_plot.append(file_path)
                            output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), file_specie, output_file_name))
                        else:
                            for time_step in time_steps:
                                if file_time_step == "{:06d}".format(int(time_step)):
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), file_specie, output_file_name))
                    else:
                        if time_steps[0] == 'all':
                            for level in levels:
                                if file_level == "{:03d}".format(int(level)):
                                    files_to_plot.append(file_path)
                                    output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), file_specie, output_file_name))
                        else:
                            for level in levels:
                                for time_step in time_steps:
                                    if file_time_step == "{:06d}".format(int(time_step)) and file_level == "{:03d}".format(int(level)):
                                        files_to_plot.append(file_path)
                                        output_files.append(os.path.join(graphical_outputs_ecdf, str(probability), file_specie, output_file_name))
    if len(files_to_plot) == 0:
        print('No files to plot')
    else:
        i = 0
        for file_to_plot in files_to_plot:
            print('plotting ' + file_to_plot)
            plot_file(file_to_plot,output_files[i])
            i += 1

plot, plot_ex_prob, time_steps, levels, days_plot, species, exceedance_probabilities, nproc, convert = read_arguments()

try:
    max_number_processes = int(os.environ["SLURM_NPROCS"])
except:
    max_number_processes = int(nproc)

cwd = os.getcwd()
original_output_folder_name = 'simulations'
converted_output_folder_name = original_output_folder_name + '_converted'
ecdf_folder_name = 'output_ecdf'
disgas_original_output_folder = os.path.join(cwd,original_output_folder_name)
disgas_converted_output_folder= os.path.join(cwd,converted_output_folder_name)
disgas_ecdf = os.path.join(cwd,ecdf_folder_name)
try:
    os.mkdir(disgas_converted_output_folder) 
except:
    print('Folder ' + disgas_converted_output_folder + ' already exists')
    list_temp = os.listdir(disgas_converted_output_folder)
    for item in list_temp:
        try:
            rmtree(os.path.join(disgas_converted_output_folder,item),ignore_errors=True)
        except:
            print('Unable to remove ' + item + ' in ' + disgas_converted_output_folder)

x0, xf, y0, yf, nx, ny, nz, n_time_steps = domain()

days, days_to_plot = extract_days()

species_properties = gas_properties()

if convert:
    for day in days:
        elaborate_day(day)

try:
    os.mkdir(disgas_ecdf)
except:
    print('Folder ' + disgas_ecdf + ' already exists')
    list_temp = os.listdir(disgas_ecdf)
    for item in list_temp:
        try:
            rmtree(os.path.join(disgas_ecdf,item), ignore_errors=True)
        except:
            print('Unable to remove ' + item + ' in ' + disgas_ecdf)

for probability in exceedance_probabilities:
    prob_folder = os.path.join(disgas_ecdf, str(probability))
    try:
        os.mkdir(prob_folder)
    except:
        print('Folder ' + prob_folder + ' already exists')
    for specie in species:
        specie_folder = os.path.join(disgas_ecdf, prob_folder,specie)
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
                    for j in range(1,n_time_steps):
                        indexes.append([probability, specie, i, j])
                else:
                    for time_step in time_steps:
                        indexes.append([probability, specie, i, time_step])
        else:
            for level in levels:
                if time_steps[0] == 'all':
                    for j in range(1,n_time_steps):
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
            ps = []
            start = n_completed_processes
            end = n_completed_processes + max_number_processes
            if end > len(indexes):
                end = len(indexes)
            try:
                pools_ecdfs[n_pool] = ThreadingPool(max_number_processes)
                pools_ecdfs[n_pool].map(ecdf,indexes[start:end])
            except:
                print('Unable to elaborate days')
                sys.exit()
            n_completed_processes = end
            if n_completed_processes == len(indexes):
                break
        n_pool += 1

save_plots()
