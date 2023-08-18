from random import sample
import os
import shutil
import subprocess
import argparse
import numpy as np
import sys
import utm


def read_arguments():
    parser = argparse.ArgumentParser(description='Input data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-N', '--nproc', default=1, help='Maximum number of allowed simultaneous processes',)
    parser.add_argument('-RT', '--run_type', default='new', help='Specify if the simulation is a new one or a restart.'
                        'Possible options are: new, restart')
    parser.add_argument('-CS', '--continuous_simulation', default='False', help='Specify if the simulation is '
                        'continuous between the specified start and end dates. Possible options are True or False')
    parser.add_argument('-RS', '--random_sources', default='off', help='on: randomly select NS locations from a '
                        'probability map. off: fixed source locations from file sources_input.txt',)
    parser.add_argument('-NS', '--nsources', default='random', help='Specify a number for a fixed number of sources. '
                        'If random, then randomly select the number of sources from an interval',)
    parser.add_argument('-SINT', '--sources_interval', default='', help='Type the minimum and maximum number of '
                        'sources',)
    parser.add_argument('-SLOC', '--source_location', default='', help='Coordinate type (UTM/GEO), latitude/northing, '
                        'longitude/easting, elevation (above ground in m) of 1 fixed source',)
    parser.add_argument('-SDX', '--source_dx', default=999999, help='Extension [m] along the X direction of 1 single '
                        'source. Option valid for Twodee only',)
    parser.add_argument('-SDY', '--source_dy', default=999999, help='Extension [m] along the Y direction of 1 single '
                        'source. Option valid for Twodee only',)
    parser.add_argument('-SDUR', '--source_dur', default=0, help='Emission duration [s] of 1 single source. Option '
                        'valid for Twodee only',)
    parser.add_argument('-D', '--domain', default='', help='Coordinates type (UTM/GEO), coordinates (latitude/northing,'
                        ' longitude/easting) of the bottom left corner and top right corner of the domain',)
    parser.add_argument('-NX', '--nx', default=-1, help='Number of grid cells along the x-direction. If not provided, '
                        'the grid spacing along the x-direction must be provided',)
    parser.add_argument('-NY', '--ny', default=-1, help='Number of grid cells along the y-direction. If not provided, '
                        'the grid spacing along the y-direction must be provided',)
    parser.add_argument('-DX', '--dx', default=-1, help='Grid spacing (in m) along the x-direction. If not provided, '
                        'the number of grid cells along the x-direction must be provided',)
    parser.add_argument('-DY', '--dy', default=-1, help='Grid spacing (in m) along the y-direction. If not provided, '
                        'the number of grid cells along the y-direction must be provided',)
    parser.add_argument('-SEM', '--source_emission', default='999', help='Source emission rate [kg/s]. If specified, '
                        'it is assigned to all the sources in the domain',)
    parser.add_argument('-RER', '--random_emission', default='off', help='on: randomly assign emission rate for each '
                        'source in the domain sampled from a flux.txt file. off: use specified emission rate',)
    parser.add_argument('-RD', '--run_duration', default=24, help='Run duration (hours). Currently fractions of hours '
                        'or duration > 24 hours are not allowed')
    parser.add_argument('-OI', '--output_interval', default=1, help='Output interval (hours). Currently fractions of '
                        'hours are not allowed')
    parser.add_argument('-OH', '--output_heights', default='', help='List of output heights (comma separated) in m '
                        'above the ground')
    parser.add_argument('-DI', '--diagno', default='on', help='on or off, to run Diagno. Turn it off only if Diagno has'
                        ' already been run')
    parser.add_argument('-DM', '--dispersion_model', default='off', help='Twodee, Disgas or Automatic')
    parser.add_argument('-US', '--use_slurm', default='False', help='True or False, to use SLURM Workload Manager')
    parser.add_argument('-SP', '--slurm_partition', default='', help='Name of the cluster partition to run the Slurm '
                        'jobs')
    parser.add_argument('-TS', '--tracking_specie', default='', help='The original emitted specie that is tracked in '
                        'the simulation')
    args = parser.parse_args()
    nproc = args.nproc
    run_type_in = args.run_type
    continuous_simulation_in = args.continuous_simulation
    random_sources_in = args.random_sources
    nsources_in = args.nsources
    source_location_in = args.source_location
    sources_interval_in = args.sources_interval
    domain_in = args.domain
    nx_in = args.nx
    ny_in = args.ny
    dx_in = args.dx
    dy_in = args.dy
    source_emission_in = args.source_emission
    random_emission_in = args.random_emission
    tracking_specie_in = args.tracking_specie
    try:
        max_number_processes_in = int(nproc)
    except ValueError:
        print('Please provide a valid number for the maximum number of process')
        sys.exit()
    model = args.dispersion_model
    diagno = args.diagno
    use_slurm = args.use_slurm
    slurm_partition = args.slurm_partition
    source_dx_in = args.source_dx
    source_dy_in = args.source_dy
    source_dur_in = args.source_dur
    source_easting_in = source_northing_in = source_el_in = 0
    sources_interval_in = sources_interval_in.split(',')
    source_location = source_location_in.split(',')
    domain = domain_in.split(',')
    run_type_in = run_type_in.lower()
    run_duration_in = args.run_duration
    output_interval_in = args.output_interval
    output_heights_in = args.output_heights
    if run_type_in != 'new' and run_type_in != 'restart':
        print('ERROR. Please provide a valid entry for -RT --run_type')
        sys.exit()
    if continuous_simulation_in.lower() == 'true':
        continuous_simulation_in = True
    elif continuous_simulation_in.lower() == 'false':
        continuous_simulation_in = False
    else:
        print('ERROR. Wrong value for variable -CS --continuous_simulation')
        sys.exit()
    try:
        source_emission_in = float(source_emission_in)
    except ValueError:
        print('Please provide a valid number for the emission rate of the source')
        sys.exit()
    if len(domain) != 5:
        print('ERROR. Please provide valid entries for -D --domain')
        sys.exit()
    else:
        coordinates_type = domain[0]
        if coordinates_type == 'GEO':
            bottom_left_1 = float(domain[1])
            bottom_left_2 = float(domain[2])
            top_right_1 = float(domain[3])
            top_right_2 = float(domain[4])
            if (-90 <= bottom_left_1 <= 90 and -180 <= bottom_left_2 <= 180) and \
                    (-90 <= top_right_1 <= 90 and -180 <= top_right_2 <= 180):  # identify valid geographic coordinates
                try:
                    out_utm = utm.from_latlon(bottom_left_1, bottom_left_2)
                    bottom_left_easting_in = float(out_utm[0])
                    bottom_left_northing_in = float(out_utm[1])
                except ValueError:
                    print('ERROR. Please provide valid coordinates for the bottom left corner of the domain')
                    sys.exit()
                try:
                    out_utm = utm.from_latlon(top_right_1, top_right_2)
                    top_right_easting_in = float(out_utm[0])
                    top_right_northing_in = float(out_utm[1])
                except ValueError:
                    print('ERROR. Please provide valid coordinates for the top right corner of the domain')
                    sys.exit()
            else:
                print('ERROR. Please provide valid coordinates')
                sys.exit()
        elif coordinates_type == 'UTM':
            bottom_left_northing_in = float(domain[1])
            bottom_left_easting_in = float(domain[2])
            top_right_northing_in = float(domain[3])
            top_right_easting_in = float(domain[4])
        else:
            print('ERROR. Please provide a valide type of coordinates (UTM or GEO)')
            sys.exit()
        if bottom_left_northing_in == top_right_northing_in or bottom_left_easting_in == top_right_easting_in:
            print('ERROR. Coordinates of the corners cannot coincide')
            sys.exit()
        if bottom_left_northing_in > top_right_northing_in and bottom_left_easting_in > top_right_easting_in:
            # Check coordinates are in the proper order, otherwise swap
            temp = bottom_left_northing_in
            bottom_left_northing_in = top_right_northing_in
            top_right_northing_in = temp
            temp = bottom_left_easting_in
            bottom_left_easting_in = top_right_easting_in
            top_right_easting_in = temp
        with open(topography, 'r') as topography_file:
            for pos, line in enumerate(topography_file):
                if pos == 2:
                    x_left_top = float(line.split(' ')[0])
                    x_right_top = float(line.split(' ')[1])
                elif pos == 3:
                    y_bottom_top = float(line.split(' ')[0])
                    y_top_top = float(line.split(' ')[1])
            if bottom_left_northing_in < y_bottom_top or bottom_left_easting_in < x_left_top or \
                    top_right_northing_in > y_top_top or top_right_easting_in > x_right_top:
                print('ERROR. Defined computational domain not consistent with topography.grd')
                print('Topography domain extent')
                print('X_west (m) = ' + str(x_left_top) + ' X_east (m) = ' + str(x_right_top))
                print('Y_south (m) = ' + str(y_bottom_top) + ' Y_north (m) = ' + str(y_top_top))
                print('Computational domain extent')
                print('X_west (m) = ' + str(bottom_left_easting_in) + ' X_east (m) = ' + str(top_right_easting_in))
                print('Y_south (m) = ' + str(bottom_left_northing_in) + ' Y_north (m) = ' + str(top_right_northing_in))
                sys.exit()
    try:
        nx_in = int(nx_in)
    except ValueError:
        print('Please provide a valid number for -NX --nx')
        sys.exit()
    try:
        ny_in = int(ny_in)
    except ValueError:
        print('Please provide a valid number for -NY --ny')
        sys.exit()
    try:
        dx_in = float(dx_in)
    except ValueError:
        print('Please provide a valid number for -DX --dx')
        sys.exit()
    try:
        dy_in = float(dy_in)
    except ValueError:
        print('Please provide a valid number for -DY --dy')
        sys.exit()
    if nx_in == -1 or ny_in == -1:
        if dx_in == -1 or dy_in == -1:
            print('ERROR. Either (NX, NY) or (DX, DY) must be specified')
            sys.exit()
        elif dx_in != -1 and dy_in != -1:
            if dx_in < 0 or dx_in > (top_right_easting_in - bottom_left_easting_in) or dy_in < 0 \
                    or dy_in > (top_right_northing_in - bottom_left_northing_in):
                print('ERROR. Please provide a valid number for (DX, DY)')
                sys.exit()
            else:
                nx_in = int((top_right_easting_in - bottom_left_easting_in) / dx_in)
                ny_in = int((top_right_northing_in - bottom_left_northing_in) / dy_in)
    elif nx_in != -1 and ny_in != -1:
        if nx_in < 0 or ny_in < 0:
            print('ERROR. Please provide a valid number for (NX, NY)')
            sys.exit()
        else:
            dx_in = (top_right_easting_in - bottom_left_easting_in) / float(nx_in)
            dy_in = (top_right_northing_in - bottom_left_northing_in) / float(ny_in)
    # Check provided nx, ny or dx, dy match the provided domain extent, otherwise correct
    if bottom_left_easting_in + dx_in * nx_in != top_right_easting_in:
        top_right_easting_in = bottom_left_easting_in + dx_in * nx_in
    if bottom_left_northing_in + dy_in * ny_in != top_right_northing_in:
        top_right_northing_in = bottom_left_northing_in + dy_in * ny_in
    if random_sources_in == 'on':
        try:
            np.loadtxt('probability_map.grd', skiprows=5)
        except ValueError:
            print('Please provide a valid probability_map.grd file when random_sources option is on')
            sys.exit()
        if nsources_in == 'random':
            if len(sources_interval_in) != 2:
                print('ERROR. Please specify the minimum and maximum number of sources with -SINT --sources_interval')
                sys.exit()
        else:
            try:
                nsources_in = int(nsources_in)
            except ValueError:
                print('Please provide a valid integer for -NS --nsources')
                sys.exit()
        if random_emission_in == 'off' and source_emission_in == 999:
            print('ERROR. random_sources set to on requires either random_emission set to on or a specified '
                  'source_emission')
            sys.exit()
    else:
        if random_sources_in != 'off':
            print('ERROR. Valid options for -RS --random_sources are on and off')
            sys.exit()
        else:
            try:
                sources_file = open('sources_input.txt', 'r')
                sources_file.close()
            except FileNotFoundError:
                print('File sources_input.txt not found. Using one source from input data')
                if len(source_location) != 4:
                    print('ERROR. Please provide valid entries for -SLOC --sources_location')
                    sys.exit()
                else:
                    coordinates_type = source_location[0]
                    if coordinates_type == 'GEO':
                        if -90 <= float(source_location[1]) <= 90 and -180 <= float(source_location[2]) <= 180:
                            # identify geographic coordinates
                            try:
                                out_utm = utm.from_latlon(float(source_location[1]), float(source_location[2]))
                                source_easting_in = float(out_utm[0])
                                source_northing_in = float(out_utm[1])
                            except ValueError:
                                print('ERROR. Please provide valid coordinates of the source location')
                                sys.exit()
                    elif coordinates_type == 'UTM':
                        source_easting_in = float(source_location[1])
                        source_northing_in = float(source_location[0])
                        if not bottom_left_easting_in <= source_easting_in <= top_right_easting_in or \
                                not bottom_left_northing_in <= source_northing_in <= top_right_northing_in:
                            print('ERROR. Location not within the domain')
                            sys.exit()
                    else:
                        print('ERROR. Please provide a valide type of coordinates (UTM or GEO)')
                        sys.exit()
                    if float(source_location[2]) < 0:
                        print('ERROR. Please provide a valid value for the source elevation in m above ground (>= 0 m)')
                        sys.exit()
                    else:
                        source_el_in = float(source_location[2])
    if random_emission_in == 'on':
        try:
            sources_file = open('flux.txt', 'r')
            sources_file.close()
        except FileNotFoundError:
            print('ERROR. File flux.txt not found')
            sys.exit()
    elif random_emission_in != 'off':
        print('ERROR. Valid options for -RER --random_sources are on and off')
        sys.exit()
    if model.lower() == 'twodee':
        twodee = True
        disgas = False
    elif model.lower() == 'disgas':
        twodee = False
        disgas = True
    elif model.lower() == 'automatic':
        twodee = True
        disgas = True
        if continuous_simulation_in:
            print('ERROR. Option Automatic for -DM --dispersion_model is currently not compatible with the continuous '
                  'simulation mode')
            sys.exit()
    else:
        print('ERROR. Please provide a valid entry for the variable -DM --dispersion_model')
        sys.exit()
    if diagno.lower() == 'on':
        diagno = True
    elif diagno.lower() == 'off':
        diagno = False
    else:
        print('ERROR. Please provide a valid entry for the variable -DI --diagno')
        sys.exit()
    if source_dx_in != 999999:
        try:
            source_dx_in = float(source_dx_in)
        except ValueError:
            print('ERROR. Please provide a valid number for the variable -SDX --source_dx')
            sys.exit()
    if source_dy_in != 999999:
        try:
            source_dy_in = float(source_dy_in)
        except ValueError:
            print('ERROR. Please provide a valid number for the variable -SDY --source_dy')
            sys.exit()
    if source_dur_in != 0:
        try:
            source_dur_in = float(source_dur_in)
        except ValueError:
            print('ERROR. Please provide a valid number for the variable -SDUR --source_dur')
            sys.exit()
    if use_slurm.lower() == 'true':
        use_slurm = True
        if slurm_partition == '':
            print('ERROR. Cluster partition must be declared if use_slurm is on')
    elif use_slurm.lower() == 'false':
        use_slurm = False
    else:
        print('ERROR. Please provide a valid entry for the variable -US --use_slurm')
        sys.exit()
    if tracking_specie_in == '':
        print('ERROR. Please specify the name of the tracked specie -TS --tracking_specie')
        sys.exit()
    slurm_partition = slurm_partition.lower()
    try:
        run_duration_in = int(run_duration_in)
        if run_duration_in < 1:
            run_duration_in = 1
        elif run_duration_in > 24:
            run_duration_in = 24
    except ValueError:
        print('ERROR. Please provide a valid entry for the variable -RD --run_duration')
        sys.exit()
    try:
        output_interval_in = int(output_interval_in)
        if output_interval_in < 1:
            output_interval_in = 1
    except ValueError:
        print('ERROR. Please provide a valid entry for the variable -OR --output_interval')
        sys.exit()
    try:
        output_heights_in = output_heights_in.split(',')
        output_heights_in_str = ''
        output_heights_in = [float(output_height) for output_height in output_heights_in]
        if 0.0 not in output_heights_in:  # 0 m level is necessary for DISGAS is the source height is set to 0 m abg
            output_heights_in.append(0.0)
        if 10.0 > max(output_heights_in):
            output_heights_in.append(10.0)  # 10.0 is the reference level for wind measurements in diagno.inp. This
            # can be generalized, at the moment ZSWIND in diagno.inp is not changed in VIGIL
        nz_in = len(output_heights_in)
        for height in sorted(output_heights_in):
            output_heights_in_str += str(height) + ' '
    except ValueError:
        print('ERROR. Please provide a valid entry for the variable -OH --output_heights')
        sys.exit()
    except IndexError:
        print('ERROR. Please provide a valid entry for the variable -OH --output_heights')
        sys.exit()
    return (
        run_type_in,
        continuous_simulation_in,
        max_number_processes_in,
        random_sources_in,
        nsources_in,
        sources_interval_in,
        source_easting_in,
        source_northing_in,
        source_el_in,
        source_emission_in,
        random_emission_in,
        bottom_left_northing_in,
        bottom_left_easting_in,
        top_right_northing_in,
        top_right_easting_in,
        nx_in,
        ny_in,
        nz_in,
        dx_in,
        dy_in,
        source_dx_in,
        source_dy_in,
        source_dur_in,
        diagno,
        twodee,
        disgas,
        use_slurm,
        slurm_partition,
        tracking_specie_in,
        run_duration_in,
        output_interval_in,
        output_heights_in_str,
    )


def prepare_days():
    days_prepared = []
    try:
        raw_days = []  # store the days as originally formatted
        with open(os.path.join(root, 'days_list.txt'), 'r', encoding='utf-8', errors='surrogateescape') as \
                days_list_file:
            for line in days_list_file:
                raw_days.append(line)
    except FileNotFoundError:
        print('ERROR. Unable to find days_list.txt file. Please restart activating DIAGNO')
        sys.exit()
    for k in range(0, len(raw_days)):
        temp = raw_days[k].split(' ')
        temp = temp[0].split('-')
        days_prepared.append(temp[0] + temp[1] + temp[2])
    return sorted(days_prepared)


def pre_process(run_mode):
    def extract_gas_properties():
        import pandas as pd
        data = pd.read_csv('gas_properties.csv', on_bad_lines='skip')
        try:
            y = np.sort(data['R_' + tracking_specie])
            r_gas = list(y)[0]
            if r_gas != r_gas:
                print('ERROR. Specific gas constant of ' + tracking_specie + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Specific gas constant of ' + tracking_specie + ' not found in gas_properties.csv')
            sys.exit()
        try:
            y = np.sort(data['T_' + tracking_specie])
            t_gas = list(y)[0]
            if t_gas != t_gas:
                print('ERROR. Gas temperature of ' + tracking_specie + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Gas temperature of ' + tracking_specie + ' not found in gas_properties.csv')
            sys.exit()
        try:
            y = np.sort(data['M_' + tracking_specie])
            m_gas = list(y)[0]
            if m_gas != m_gas:
                print('ERROR. Gas molar weight of ' + tracking_specie + ' not found in gas_properties.csv')
                sys.exit()
        except KeyError:
            print('ERROR. Gas molar weight of ' + tracking_specie + ' not found in gas_properties.csv')
            sys.exit()
        try:
            y = np.sort(data['BG_' + tracking_specie])
            bg_conc = list(y)[0]
            if bg_conc != bg_conc:
                print('WARNING. Background concentration of ' + tracking_specie + ' not found in gas_properties.csv')
        except KeyError:
            print('ERROR. Background concentration of ' + tracking_specie + ' not found in gas_properties.csv')
            sys.exit()
        rho_g = 101325 / (r_gas * t_gas)
        rho_g_20 = 101325 / (r_gas * 293)  # Time-averaged gas density for the twodee.inp file
        return r_gas, t_gas, rho_g, rho_g_20, m_gas, bg_conc

    def sample_random_sources(
        n_sources_inp, input_file, dur_min_inp, dur_max_inp, source_size_min_inp, source_size_max_inp
    ):
        from random import choices

        with open(input_file) as probability_file:
            i = 1
            for line_probability in probability_file:
                if i == 2:
                    nx_prob = int(line_probability.split(' ')[0])
                    ny_prob = int(line_probability.split(' ')[1])
                    i += 1
                    continue
                elif i == 3:
                    xmin = float(line_probability.split(' ')[0])
                    xmax = float(line_probability.split(' ')[1])
                    i += 1
                    continue
                elif i == 4:
                    ymin = float(line_probability.split(' ')[0])
                    ymax = float(line_probability.split(' ')[1])
                    i += 1
                    continue
                else:
                    i += 1
                    continue
        probabilities_input = np.loadtxt(input_file, skiprows=5)
        location_cum_indexes = []
        location_indexes = []
        sampled_probabilities = []
        dxs = []
        dys = []
        durs = []
        k = 0
        for i in range(0, nx_prob):
            for j_prob in range(0, ny_prob):
                location_cum_indexes.append(k)
                location_indexes.append([i, j_prob])
                k += 1
        x = np.linspace(xmin, xmax, nx_prob)
        y = np.linspace(ymin, ymax, ny_prob)
        for i in range(0, nx_prob):
            for j_prob in range(0, ny_prob):
                sampled_probabilities.append(probabilities_input[i, j_prob])
        selected_locations = choices(location_cum_indexes, sampled_probabilities, k=n_sources_inp)
        source_size = source_size_min_inp
        while source_size <= source_size_max_inp:
            dxs.append(source_size)
            dys.append(source_size)
            source_size += 0.1
        source_duration = dur_min_inp
        while source_duration <= dur_max_inp:
            durs.append(source_duration)
            source_duration += 60
        for location in selected_locations:
            row = location_indexes[location][0]
            column = location_indexes[location][1]
            xpr = x[row]
            ypr = y[column]
            probability = sampled_probabilities[location]
            random_eastings.append(xpr)
            random_northings.append(ypr)
            random_elevations.append(0.0)
            random_probabilities.append(probability)
            random_fluxes.append(99999999)
            random_dx.append(choices(dxs, k=1)[0])
            random_dy.append(choices(dys, k=1)[0])
            random_dur.append(choices(durs, k=1)[0])
            random_temperatures.append(gas_temperature)
        return (
            random_eastings,
            random_northings,
            random_elevations,
            random_probabilities,
            random_fluxes,
            random_dx,
            random_dy,
            random_dur,
            random_temperatures,
        )

    def fluxes():
        fluxes_in = []
        with open('flux.txt') as flux_file:
            for line_flux in flux_file:
                fluxes_in.append(float(line_flux))
        flux_file.close()
        x = np.sort(fluxes_in)
        list_x = list(x)
        sampled_flux = sample(list_x, 1)
        return sampled_flux

    # Extract temperature and gas constant of the tracking specie
    gas_constant, gas_temperature, gas_density, gas_density_20, gas_molar_weight, gas_background_concentration = \
        extract_gas_properties()
    # Set source files for DISGAS and TWODEE. TWODEE also needs source start and stop time, to be generalized
    random_eastings = []
    random_northings = []
    random_elevations = []
    random_probabilities = []
    random_fluxes = []
    random_dx = []
    random_dy = []
    random_dur = []
    gas_fluxes = []
    random_temperatures = []
    dur_min = 1
    dur_max = 86400
    source_size_min = 0.1
    source_size_max = 100
    if random_sources == 'on':
        if nsources == 'random':
            n_sources = [*range(int(sources_interval[0]), int(sources_interval[1]) + 1)]
        else:
            n_sources = [nsources]
        n_random_sources = sample(n_sources, 1)[0]
        (
            random_eastings,
            random_northings,
            random_elevations,
            random_probabilities,
            random_fluxes,
            random_dx,
            random_dy,
            random_dur,
            random_temperatures,
        ) = sample_random_sources(
            n_random_sources,
            'probability_map.grd',
            dur_min,
            dur_max,
            source_size_min,
            source_size_max,
        )
    easting = random_eastings
    northing = random_northings
    elevations = random_elevations
    probabilities = random_probabilities
    fluxes_input = random_fluxes
    # dx_sources_twodee = random_dx
    # dy_sources_twodee = random_dy
    # FABIO: new lines
    dx_src = random_dx
    dy_src = random_dy
    dur = random_dur
    source_temperatures = random_temperatures
    n_sources = 0
    try:
        with open('sources_input.txt', 'r', encoding='utf-8', errors='surrogateescape') as locations_file:
            for line in locations_file:
                try:
                    records = line.split('\t')
                    easting.append(float(records[0]))
                    northing.append(float(records[1]))
                    elevations.append(float(records[2]))
                    probabilities.append(float(records[3]))
                    fluxes_input.append(float(records[4]))
                    source_temperatures.append(float(records[5]))
                    # FABIO: new lines
                    dx_src.append(float(records[6]))
                    dy_src.append(float(records[7]))
                    dur.append(float(records[8]))
                    ####################
                    # if twodee_on:
                    #     try:
                    #         dx_sources_twodee.append(float(records[6]))
                    #         dy_sources_twodee.append(float(records[7]))
                    #         dur.append(float(records[8]))
                    #     except IndexError:
                    #         dx_sources_twodee.append(1)
                    #         dy_sources_twodee.append(1)
                    #         dur.append(86400)
                    n_sources += 1
                except ValueError:
                    continue
                except IndexError:
                    print('ERROR. Badly formatted sources_input.txt file (probably one or more data are missing')
                    sys.exit()
    except FileNotFoundError:
        if random_sources != 'on':
            easting.append(source_easting)
            northing.append(source_northing)
            elevations.append(source_el)
            probabilities.append(1.0)
            fluxes_input.append(source_emission)
            source_temperatures.append(gas_temperature)  # If problems with reading sources_input.txt, assign
            # the same temperauture to all sources
            # FABIO: new lines
            dx_src.append(source_dx)
            dy_src.append(source_dy)
            ######################
            # dx_sources_twodee.append(source_dx)
            # dy_sources_twodee.append(source_dy)
            dur.append(source_dur)
        n_sources = len(easting)

    for j_source in range(0, n_sources):
        if source_emission != 999:
            gas_fluxes.append(source_emission)
        else:
            if fluxes_input[j_source] == 99999999 and random_emission == 'on':
                gas_fluxes.append(fluxes()[0])
            else:
                gas_fluxes.append(fluxes_input[j_source])
    simulations = os.path.join(root, 'simulations')
    diagno = os.path.join(simulations, 'diagno')
    runs_folder = os.path.join(simulations, 'runs')
    # Set runs folder
    try:
        os.mkdir(simulations)
    except FileExistsError:
        print('Folder ' + simulations + ' already exists')
    try:
        os.mkdir(runs_folder)
    except FileExistsError:
        print('Folder ' + runs_folder + ' already exists')
        if not diagno_on:
            shutil.rmtree(runs_folder)
            os.mkdir(runs_folder)
    for i_day in range(0, len(days)):
        day = days[i_day]
        diagno_daily = os.path.join(diagno, str(day))
        if diagno_on:
            met_data = os.path.join(simulations, str(day))
        else:
            met_data = os.path.join(simulations, 'diagno', str(day))
        run = os.path.join(runs_folder, str(day))
        try:
            os.mkdir(run)
        except FileExistsError:
            print('Folder ' + run + ' already exists')
        if continuous_simulation:
            if i_day == 0:
                run_mode = 'new'
            else:
                run_mode = 'restart'
        # Set DISGAS folder
        if disgas_on:
            disgas_daily = os.path.join(run, 'disgas')
            runs_disgas.append(disgas_daily)
            outfiles_disgas = os.path.join(disgas_daily, 'outfiles')
            try:
                os.mkdir(disgas_daily)
            except FileExistsError:
                print('Folder ' + disgas_daily + ' already exists')
            disgas_input = os.path.join(disgas_daily, 'disgas.inp')
            with open(os.path.join(disgas_daily, 'source.dat'), 'w', encoding='utf-8', errors='surrogateescape',) as \
                    source_file:
                for j_source in range(0, n_sources):
                    # source_file.write('{0:7.3f}'.format(easting[j_source]) + ' ' +
                    #                   '{0:7.3f}'.format(northing[j_source]) + ' ' +
                    #                   '{0:7.2f}'.format(elevations[j_source]) + ' ' + str(gas_fluxes[j_source]) +
                    #                   '\n')
                    # FABIO: new lines
                    if dx_src[j_source] * dy_src[j_source] > dx * dy:
                        # Split the source across the computational cells, since currently DISGAS does only allow point
                        # sources
                        nx_source = int(dx_src[j_source] / dx) + 1
                        ny_source = int(dy_src[j_source] / dy) + 1
                        x_or_source = easting[j_source] - dx_src[j_source] / 2
                        y_or_source = northing[j_source] - dy_src[j_source] / 2
                        splitted_flux = gas_fluxes[j_source] / (nx_source * ny_source)
                        x_source = x_or_source
                        y_source = y_or_source
                        for jj_source in range(0, ny_source):
                            for ii_source in range(0, nx_source):
                                source_file.write('{0:7.3f}'.format(x_source) + ' ' + '{0:7.3f}'.format(y_source) + ' '
                                                  + '{0:7.2f}'.format(elevations[j_source]) + ' ' + str(splitted_flux) +
                                                  '\n')
                                x_source += dx
                            y_source += dy
                    else:
                        source_file.write('{0:7.3f}'.format(easting[j_source]) + ' ' +
                                          '{0:7.3f}'.format(northing[j_source]) + ' ' +
                                          '{0:7.2f}'.format(elevations[j_source]) + ' ' + str(gas_fluxes[j_source]) +
                                          '\n')
                    #################################
            source_file.close()
            roughness_file_exist = True
            try:
                shutil.copyfile(os.path.join(root, 'roughness.grd'), os.path.join(disgas_daily, 'roughness.grd'),)
            except (FileExistsError, FileNotFoundError):
                roughness_file_exist = False
            try:
                shutil.copyfile(os.path.join(met_data, 'surface_data.txt'), os.path.join(disgas_daily,
                                                                                         'surface_data.txt'),)
            except (FileExistsError, FileNotFoundError):
                print('ERROR with surface_data.txt')
            # read and memorize disgas.inp file
            disgas_input_records = []
            with open(disgas_original, 'r', encoding='utf-8', errors='surrogateescape') as disgas_or_input:
                for line in disgas_or_input:
                    disgas_input_records.append(line)
            with open(disgas_input, 'w', encoding='utf-8', errors='surrogateescape') as disgas_input_file:
                for record in disgas_input_records:
                    if 'ROUGHNESS_MODEL' in record:
                        roughness_command = record.split('=')[1]
                        roughness_command = roughness_command.split('(')[0]
                        if 'MATRIX' in roughness_command:
                            if not roughness_file_exist:
                                print('Warning! ROUGHNESS_MODEL set to MATRIX in disgas.inp and roughness file '
                                    'does not exist')
                                print('Setting ROUGHNESS_MODEL to UNIFORM')
                    if 'YEAR' in record:
                        disgas_input_file.write('  YEAR   = ' + day[0:4] + '\n')
                    elif 'MONTH' in record:
                        disgas_input_file.write('  MONTH  = ' + day[4:6] + '\n')
                    elif 'DAY' in record:
                        disgas_input_file.write('  DAY    = ' + day[6:8] + '\n')
                    elif 'HOUR' in record:
                        try:
                            hour_start = float(record.split('=')[1])
                        except ValueError:
                            hour_start = record.split('=')[1].strip()
                            hour_start = float(hour_start.split(' ')[0])
                        if continuous_simulation:
                            if i_day > 0:
                                hour_start = 0
                        disgas_input_file.write('  HOUR   = ' + '{0:2.0f}'.format(hour_start) + '\n')
                    elif 'SIMULATION_INTERVAL_(SEC)' in record:
                        simulation_interval = run_duration * 3600
                        disgas_input_file.write('  SIMULATION_INTERVAL_(SEC) = ' +
                                                '{0:7.0f}'.format(simulation_interval) + '\n')
                    elif 'NX' in record:
                        disgas_input_file.write('  NX     = ' + str(nx) + '\n')
                    elif 'NY' in record:
                        disgas_input_file.write('  NY     = ' + str(ny) + '\n')
                    elif 'NZ' in record:
                        disgas_input_file.write('  NZ     = ' + str(nz) + '\n')
                    elif 'Z_LAYERS' in record:
                        disgas_input_file.write('  Z_LAYERS_(M) = ' + output_heights + '\n')
                    elif 'DX_(M)' in record:
                        disgas_input_file.write('  DX_(M) = ' + str(dx) + '\n')
                    elif 'DY_(M)' in record:
                        disgas_input_file.write('  DY_(M) = ' + str(dy) + '\n')
                    elif 'X_ORIGIN_(UTM_M)' in record:
                        disgas_input_file.write('  X_ORIGIN_(UTM_M) = ' + str(bottom_left_easting) + '\n')
                    elif 'Y_ORIGIN_(UTM_M)' in record:
                        disgas_input_file.write('  Y_ORIGIN_(UTM_M) = ' + str(bottom_left_northing) + '\n')
                    elif 'RESTART_RUN' in record:
                        if run_mode == 'restart':
                            disgas_input_file.write('  RESTART_RUN = YES\n')
                        else:
                            disgas_input_file.write('  RESTART_RUN = NO\n')
                    elif 'RESET_TIME' in record:
                        if run_mode == 'restart':
                            disgas_input_file.write('  RESET_TIME = YES\n')
                        else:
                            disgas_input_file.write('  RESET_TIME = NO\n')
                    elif 'TOPOGRAPHY_FILE_PATH' in record:
                        disgas_input_file.write('   TOPOGRAPHY_FILE_PATH   = ' +
                                                os.path.join(diagno_daily, 'topography.grd') + ' \n')
                    elif 'ROUGHNESS_FILE_PATH' in record:
                        disgas_input_file.write('   ROUGHNESS_FILE_PATH   = ' +
                                                os.path.join(disgas_daily, 'roughness.grd') + ' \n')
                    elif 'RESTART_FILE_PATH' in record:
                        disgas_input_file.write('   RESTART_FILE_PATH   = ' +
                                                os.path.join(disgas_daily, 'restart.dat') + ' \n')
                    elif 'SOURCE_FILE_PATH' in record:
                        disgas_input_file.write('   SOURCE_FILE_PATH   = ' +
                                                os.path.join(disgas_daily, 'source.dat') + ' \n')
                    elif 'WIND_FILE_PATH' in record:
                        disgas_input_file.write('   WIND_FILE_PATH   = ' +
                                                os.path.join(disgas_daily, 'winds.dat') + ' \n')
                    elif 'DIAGNO_FILE_PATH' in record:
                        disgas_input_file.write('   DIAGNO_FILE_PATH   = ' +
                                                os.path.join(diagno_daily, 'diagno.out') + ' \n')
                    elif 'OUTPUT_DIRECTORY' in record:
                        disgas_input_file.write('   OUTPUT_DIRECTORY    = ' + outfiles_disgas + ' \n')
                    elif 'OUTPUT_INTERVAL' in record:
                        output_interval_sec = output_interval * 3600
                        disgas_input_file.write('  OUTPUT_INTERVAL_(SEC) = ' + str(output_interval_sec) + '\n')
                    else:
                        disgas_input_file.write(record)
            shutil.copy(disgas_input, disgas_original)
        if twodee_on:
            twodee_daily = os.path.join(run, 'twodee')
            runs_twodee.append(twodee_daily)
            outfiles_twodee = os.path.join(twodee_daily, 'outfiles')
            try:
                os.mkdir(twodee_daily)
            except FileExistsError:
                print('Folder ' + twodee_daily + ' already exists')
            twodee_input = os.path.join(twodee_daily, 'twodee.inp')
            try:
                shutil.copyfile(os.path.join(root, 'roughness.grd'), os.path.join(twodee_daily, 'roughness.grd'),)
            except FileNotFoundError:
                print('Unable to find a valid roughness file for TWODEE')
                sys.exit()
            try:
                shutil.copyfile(os.path.join(met_data, 'surface_data.txt'), os.path.join(twodee_daily,
                                'surface_data.txt'),)
            except (FileExistsError, FileNotFoundError):
                print('ERROR with surface_data.txt')
            with open(os.path.join(twodee_daily, 'source.dat'), 'w', encoding='utf-8', errors='surrogateescape',) as \
                    source_file:
                for j in range(0, n_sources):
                    # FABIO: new lines
                    source_file.write('{0:7.3f}'.format(easting[j]) + ' ' + '{0:7.3f}'.format(northing[j]) + ' '
                                      + '{0:7.3f}'.format(gas_fluxes[j]) + ' ' + '{0:7.2f}'.format(dx_src[j]) + ' '
                                      + '{0:7.2f}'.format(dy_src[j]) + ' KG_SEC 0 ' + '{0:7.3f}'.format(dur[j]) + '\n')
                    ##############################
                    # source_file.write('{0:7.3f}'.format(easting[j]) + ' ' + '{0:7.3f}'.format(northing[j]) + ' '
                    #     + '{0:7.3f}'.format(gas_fluxes[j]) + ' ' + '{0:7.2f}'.format(dx_sources_twodee[j]) + ' '
                    #    + '{0:7.2f}'.format(dy_sources_twodee[j]) + ' KG_SEC 0 ' + '{0:7.3f}'.format(dur[j]) + '\n')
            source_file.close()
            # read and memorize twodee.inp file. First read and memorize the domain properties of diagno, since
            # twodee's domain must coincide with it
            with open(os.path.join(met_data, 'diagno.inp'), 'r') as diagno_input_data:
                for line in diagno_input_data:
                    if 'NX\n' in line:
                        nx_diagno = int(line.split('NX')[0])
                    elif 'NY\n' in line:
                        ny_diagno = int(line.split('NY')[0])
                    elif 'UTMXOR' in line:
                        xor_diagno = float(line.split('UTMXOR')[0]) * 1000
                    elif 'UTMYOR' in line:
                        yor_diagno = float(line.split('UTMYOR')[0]) * 1000
            twodee_input_records = []
            with open(twodee_original, 'r', encoding='utf-8', errors='surrogateescape') as twodee_or_input:
                for line in twodee_or_input:
                    twodee_input_records.append(line)
            with open(twodee_input, 'w', encoding='utf-8', errors='surrogateescape') as twodee_input_file:
                for record in twodee_input_records:
                    if 'YEAR' in record:
                        twodee_input_file.write('  YEAR   = ' + day[0:4] + '\n')
                    elif 'MONTH' in record:
                        twodee_input_file.write('  MONTH  = ' + day[4:6] + '\n')
                    elif 'DAY' in record:
                        twodee_input_file.write('  DAY    = ' + day[6:8] + '\n')
                    elif 'HOUR' in record and 'NC_VAR_HOUR_NAME' not in record:
                        try:
                            hour_start = float(record.split('=')[1])
                        except ValueError:
                            hour_start = record.split('=')[1].strip()
                            hour_start = float(hour_start.split(' ')[0])
                        if continuous_simulation:
                            if i_day > 0:
                                hour_start = 0
                        twodee_input_file.write('  HOUR   = ' + '{0:2.0f}'.format(hour_start) + '\n')
                    elif 'SIMULATION_INTERVAL_(SEC)' in record:
                        # for the moment this is managed exactly like DISGAS, but would require the RESET_TIME option to
                        # be implemented in TWODEE as well
                        simulation_interval = run_duration * 3600
                        twodee_input_file.write('  SIMULATION_INTERVAL_(SEC) = ' +
                                                '{0:7.0f}'.format(simulation_interval) + '\n')
                    elif 'NX' in record:
                        twodee_input_file.write('  NX     = ' + str(nx_diagno) + '\n')
                    elif 'NY' in record:
                        twodee_input_file.write('  NY     = ' + str(ny_diagno) + '\n')
                    elif 'DX_(M)' in record:
                        twodee_input_file.write('  DX_(M) = ' + str(dx) + '\n')
                    elif 'DY_(M)' in record:
                        twodee_input_file.write('  DY_(M) = ' + str(dy) + '\n')
                    elif 'X_ORIGIN_(UTM_M)' in record:
                        twodee_input_file.write('  X_ORIGIN_(UTM_M) = ' + str(xor_diagno) + '\n')
                    elif 'Y_ORIGIN_(UTM_M)' in record:
                        twodee_input_file.write('  Y_ORIGIN_(UTM_M) = ' + str(yor_diagno) + '\n')
                    elif 'AMBIENT_GAS_DENSITY_20C_(KG/M3)' in record:
                        twodee_input_file.write('  AMBIENT_GAS_DENSITY_20C_(KG/M3) = 1.204\n')
                    elif 'DENSE_GAS_DENSITY_20C_(KG/M3)' in record:
                        twodee_input_file.write('  DENSE_GAS_DENSITY_20C_(KG/M3) = ' +
                                                '{0:6.3f}'.format(gas_density_20) + '\n')
                    elif 'AVERAGED_TEMPERATURE_(C)' in record:
                        twodee_input_file.write('  AVERAGED_TEMPERATURE_(C) = ' + str(gas_temperature) + '\n')
                        # Average temperature initialized with the temperature specified in gas_properties.csv.
                        # If T for each source is provided in sources_input.txt, this will be reassigned
                    elif 'RESTART_RUN' in record:
                        if run_mode == 'restart':
                            twodee_input_file.write('  RESTART_RUN = YES\n')
                        else:
                            twodee_input_file.write('  RESTART_RUN = NO\n')
                    elif 'OUTPUT_INTERVAL' in record:
                        output_interval_sec = output_interval * 3600
                        twodee_input_file.write('  OUTPUT_INTERVAL_(SEC) = ' + str(output_interval_sec) + '\n')
                    elif 'OUTPUT_DIRECTORY' in record:
                        twodee_input_file.write('   OUTPUT_DIRECTORY   = ' + outfiles_twodee + ' \n')
                    elif 'HEIGHTS_(M)' in record:
                        twodee_input_file.write('     HEIGHTS_(M)          = ' + output_heights + '\n')
                    elif 'CONCENTRATION_BG' in record:
                        twodee_input_file.write('     CONCENTRATION_BG     = ' + str(gas_background_concentration) +
                                                '\n')
                    elif 'TOPOGRAPHY_FILE' in record and 'TOPOGRAPHY_FILE_FORMAT' not in record:
                        twodee_input_file.write('   TOPOGRAPHY_FILE   = ' +
                                                os.path.join(diagno_daily, 'topography.grd') + ' \n')
                    elif 'ROUGHNESS_FILE' in record and 'ROUGHNESS_FILE_FORMAT' not in record:
                        twodee_input_file.write('   ROUGHNESS_FILE   = ' +
                                                os.path.join(twodee_daily, 'roughness.grd') + ' \n')
                    elif 'SOURCE_FILE' in record:
                        twodee_input_file.write('   SOURCE_FILE   = ' + os.path.join(twodee_daily, 'source.dat') +
                                                ' \n')
                    elif 'SURF_DATA_FILE' in record:
                        twodee_input_file.write('   SURF_DATA_FILE   = ' +
                                                os.path.join(diagno_daily, 'surface_data.txt') + ' \n')
                    elif 'DIAGNO_FILE' in record:
                        twodee_input_file.write('   DIAGNO_FILE   = ' + os.path.join(diagno_daily, 'diagno.out') +
                                                ' \n')
                    elif 'RESTART_FILE' in record:
                        twodee_input_file.write('   RESTART_FILE   = ' + os.path.join(twodee_daily, 'restart.dat') +
                                                ' \n')
                    else:
                        twodee_input_file.write(record)
            shutil.copy(twodee_input, twodee_original)
    # FABIO: new lines
    return easting, northing, elevations, dx_src, dy_src, dur, gas_fluxes, source_temperatures, gas_density, \
           gas_constant, gas_molar_weight
    ######################################
    # return easting, northing, elevations, dx_sources_twodee, dy_sources_twodee, dur, gas_fluxes, source_temperatures, \
    #    gas_density, gas_constant, gas_molar_weight


def run_diagno(max_np):
    def prepare_diagno():
        diagno = os.path.join(root, 'simulations', 'diagno')
        try:
            os.mkdir(diagno)
        except FileExistsError:
            print('Folder ' + diagno + ' already exists')
        for i_day in range(0, len(days)):
            day_diagno = days[i_day]
            diagno_daily = os.path.join(diagno, str(day_diagno))
            path = os.path.join(root, 'simulations', str(day_diagno))
            files = os.listdir(path)
            try:
                os.mkdir(diagno_daily)
            except FileExistsError:
                print('Folder ' + diagno_daily + ' already exists')
            for f in files:
                path_f = os.path.join(path, f)
                try:
                    shutil.move(path_f, diagno_daily)
                except FileExistsError:
                    print('File ' + f + ' already present in ' + diagno)
            shutil.copy(topography, os.path.join(diagno_daily, 'topography.grd'))
            shutil.rmtree(path)

    prepare_diagno()
    n_elaborated_days = 0
    n_node = 0
    max_height = 0
    node = ''
    while n_elaborated_days <= len(days):
        ps = []
        start_day = n_elaborated_days
        end_day = n_elaborated_days + max_np
        if end_day > len(days):
            end_day = len(days)
        try:
            for day in days[start_day:end_day]:
                if len(nodes_list) > 0:
                    try:
                        node = nodes_list[n_node]
                    except IndexError:
                        node = ''
                diagno_folder = os.path.join(root, 'simulations', 'diagno', day)
                # read and memorize diagno.inp file
                diagno_input_records = []
                diagno_input = os.path.join(diagno_folder, 'diagno.inp')
                with open(diagno_input, 'r', encoding='utf-8', errors='surrogateescape') as diagno_or_input:
                    for line in diagno_or_input:
                        diagno_input_records.append(line)
                new_record_string = ''
                old_record = ''
                heights = []
                for record in diagno_input_records:
                    if 'CELLZB' in record:
                        old_record = record
                        new_record_string = ''
                        str_heights = record.split('CELLZB(1:NZ+1)')[0].split(' ')
                        for height in str_heights:
                            try:
                                heights.append(float(height))
                            except ValueError:
                                continue
                        if 10. not in heights:
                            heights.append(10.)
                        heights = sorted(heights)
                        max_height = max(heights)
                        for height in heights:
                            new_record_string += str(height) + ' '
                        new_record_string += ' CELLZB(1:NZ+1) (m) ' + '\n'
                if new_record_string != '':
                    diagno_input_records = [new_record_string if item == old_record else item for item in
                                            diagno_input_records]
                with open(diagno_input, 'w', encoding='utf-8', errors='surrogateescape') as diagno_input_file:
                    for record in diagno_input_records:
                        if 'NX' in record:
                            if disgas_on:
                                diagno_input_file.write(str(nx + 2) + '          NX\n')
                            else:
                                diagno_input_file.write(str(nx) + '          NX\n')
                        elif 'NY' in record:
                            if disgas_on:
                                diagno_input_file.write(str(ny + 2) + '          NY\n')
                            else:
                                diagno_input_file.write(str(ny) + '          NY\n')
                        elif 'NZ' in record and ('CELLZB' not in record and 'NZPRNT' not in record):
                            diagno_input_file.write(str(len(heights) - 1) + '          NZ\n')
                        elif 'DXK' in record:
                            diagno_input_file.write('{0:7.3f}'.format(dx / 1000) + '          DXK (km)\n')
                        elif 'DYK' in record:
                            diagno_input_file.write('{0:7.3f}'.format(dy / 1000) + '          DYK (km)\n')
                        elif 'UTMXOR' in record:
                            if disgas_on:
                                diagno_input_file.write('{0:7.3f}'.format((bottom_left_easting - dx) / 1000) + '      '
                                                        'UTMXOR (km)\n')
                            else:
                                diagno_input_file.write('{0:7.3f}'.format(bottom_left_easting / 1000) + '      '
                                                        'UTMXOR (km)\n')
                        elif 'UTMYOR' in record:
                            if disgas_on:
                                diagno_input_file.write('{0:7.3f}'.format((bottom_left_northing - dy) / 1000) + '      '
                                                        'UTMYOR  (km)\n')
                            else:
                                diagno_input_file.write('{0:7.3f}'.format(bottom_left_northing / 1000) + '      '
                                                        'UTMYOR  (km)\n')
                        else:
                            diagno_input_file.write(record)
                os.chdir(diagno_folder)
                if slurm:
                    try:
                        p = subprocess.Popen(['srun', '-n', '1', '--partition=' + partition, 'presfc', '&'])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' presfc')
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(['presfc'])
                    except FileNotFoundError:
                        print('Unable to run presfc for DIAGNO')
                        sys.exit()
                p.wait()
                ps.append(p)
                if shutil.which('preupr') is not None and \
                        (not os.path.exists(os.path.join(root, 'weather_stations_list.txt'))):
                    if slurm:
                        try:
                            p = subprocess.Popen(['srun', '-n', '1', '--partition=' + partition, 'preupr', '&'])
                        except FileNotFoundError:
                            print('Unable to run srun -n 1 --partition=' + partition + ' preupr')
                            sys.exit()
                    else:
                        try:
                            p = subprocess.Popen(['preupr'])
                        except FileNotFoundError:
                            print('Unable to run preupr for DIAGNO')
                            sys.exit()
                    p.wait()
                    ps.append(p)
                if slurm:
                    try:
                        p = subprocess.Popen(['srun', '-n', '1', '--partition=' + partition, '--nodelist=' + node,
                                              'diagno', '&'])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' --nodelist=' + node + ' diagno')
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(['diagno'])
                    except FileNotFoundError:
                        print('Unable to run DIAGNO')
                        sys.exit()
                ps.append(p)
                if len(nodes_list) > 0:
                    n_node += 1
                    if n_node >= len(nodes_list):
                        n_node = 0
        except BaseException:
            print('Unable to process weather data with Diagno')
            sys.exit()
        n_elaborated_days = end_day
        if n_elaborated_days == len(days):
            break
        for p in ps:
            p.wait()
    for p in ps:
        p.wait()
    print('All weather data have been successfully processed with Diagno')
    os.chdir(root)
    return max_height


def read_diagno_outputs():
    from scipy.io import FortranFile

    def prepare_temporary_simulations(day_simulation):
        try:
            os.remove(os.path.join(runs_folder, day_simulation, 'source.dat'))
        except FileNotFoundError:
            pass
        disgas_folder = os.path.join(runs_folder, day_simulation, 'disgas')
        twodee_folder = os.path.join(runs_folder, day_simulation, 'twodee')
        runs_disgas.append(disgas_folder)
        runs_twodee.append(twodee_folder)
        disgas_inp = os.path.join(disgas_folder, 'disgas.inp')
        twodee_inp = os.path.join(twodee_folder, 'twodee.inp')
        disgas_inp_new = os.path.join(disgas_folder, 'disgas.inp_new')
        twodee_inp_new = os.path.join(twodee_folder, 'twodee.inp_new')
        disgas_input_records = []
        twodee_input_records = []
        try:
            os.mkdir(os.path.join(disgas_folder, 'outfiles'))
        except FileExistsError:
            print('Folder ' + os.path.join(disgas_folder, 'outfiles') + ' already exists')
        # Copy the necessary files into the temporary twodee folder
        try:
            os.mkdir(os.path.join(twodee_folder, 'outfiles'))
        except FileExistsError:
            print('Folder ' + os.path.join(twodee_folder, 'outfiles') + ' already exists')
        with open(os.path.join(disgas_folder, 'source.dat'), 'w', encoding='utf-8', errors='surrogateescape',) \
                as source_file:
            for j_source in index_sources_disgas:
                source_file.write(
                    '{0:7.3f}'.format(easting_sources[j_source])
                    + ' '
                    + '{0:7.3f}'.format(northing_sources[j_source])
                    + ' '
                    + '{0:7.2f}'.format(elevation_sources[j_source])
                    + ' '
                    + str(gas_fluxes_sources[j_source])
                    + '\n')
        with open(os.path.join(twodee_folder, 'source.dat'), 'w', encoding='utf-8', errors='surrogateescape',) \
                as source_file:
            for j_source in index_sources_twodee:
                source_file.write(
                    '{0:7.3f}'.format(easting_sources[j_source])
                    + ' '
                    + '{0:7.3f}'.format(northing_sources[j_source])
                    + ' '
                    + '{0:7.3f}'.format(gas_fluxes_sources[j_source])
                    + ' '
                    + '{0:7.2f}'.format(dx_sources[j_source])
                    + ' '
                    + '{0:7.2f}'.format(dy_sources[j_source])
                    + ' KG_SEC 0 '
                    + '{0:7.3f}'.format(dur_sources[j_source])
                    + '\n')
        # Update the input files where necessary
        with open(disgas_inp, 'r', encoding='utf-8', errors='surrogateescape') as disgas_or_input:
            for line_input in disgas_or_input:
                disgas_input_records.append(line_input)
        with open(twodee_inp, 'r', encoding='utf-8', errors='surrogateescape') as twodee_or_input:
            for line_input in twodee_or_input:
                twodee_input_records.append(line_input)
        with open(disgas_inp_new, 'w', encoding='utf-8', errors='surrogateescape') as disgas_input_file:
            for record in disgas_input_records:
                if 'ROUGHNESS_FILE_PATH' in record:
                    disgas_input_file.write('   ROUGHNESS_FILE_PATH   = ' +
                                            os.path.join(disgas_folder, 'roughness.grd') + ' \n')
                elif 'RESTART_FILE_PATH' in record:
                    disgas_input_file.write('   RESTART_FILE_PATH   = ' + os.path.join(disgas_folder, 'restart.dat') +
                                            ' \n')
                elif 'SOURCE_FILE_PATH' in record:
                    disgas_input_file.write('   SOURCE_FILE_PATH   = ' + os.path.join(disgas_folder, 'source.dat')
                        + ' \n')
                elif 'WIND_FILE_PATH' in record:
                    disgas_input_file.write('   WIND_FILE_PATH   = ' + os.path.join(disgas_folder, 'winds.dat') + ' \n')
                elif 'OUTPUT_DIRECTORY' in record:
                    disgas_input_file.write('   OUTPUT_DIRECTORY    = ' + os.path.join(disgas_folder, 'outfiles') +
                                            ' \n')
                else:
                    disgas_input_file.write(record)
        with open(twodee_inp_new, 'w', encoding='utf-8', errors='surrogateescape') as twodee_input_file:
            for record in twodee_input_records:
                if 'NX' in record:
                    twodee_input_file.write('  NX     = ' + str(nx_diagno) + '\n')
                elif 'NY' in record:
                    twodee_input_file.write('  NY     = ' + str(ny_diagno) + '\n')
                elif 'X_ORIGIN_(UTM_M)' in record:
                    twodee_input_file.write('  X_ORIGIN_(UTM_M) = ' + str(xor_diagno) + '\n')
                elif 'Y_ORIGIN_(UTM_M)' in record:
                    twodee_input_file.write('  Y_ORIGIN_(UTM_M) = ' + str(yor_diagno) + '\n')
                elif 'OUTPUT_DIRECTORY' in record:
                    twodee_input_file.write('   OUTPUT_DIRECTORY   = ' + os.path.join(twodee_folder, 'outfiles') +
                                            ' \n')
                elif 'ROUGHNESS_FILE' in record and 'ROUGHNESS_FILE_FORMAT' not in record:
                    twodee_input_file.write('   ROUGHNESS_FILE   = ' + os.path.join(twodee_folder, 'roughness.grd') +
                                            ' \n')
                elif 'SOURCE_FILE' in record:
                    twodee_input_file.write('   SOURCE_FILE   = ' + os.path.join(twodee_folder, 'source.dat') + ' \n')
                elif 'RESTART_FILE' in record:
                    twodee_input_file.write('   RESTART_FILE   = ' + os.path.join(twodee_folder, 'restart.dat') + ' \n')
                else:
                    twodee_input_file.write(record)
        os.rename(disgas_inp_new, disgas_inp)
        os.rename(twodee_inp_new, twodee_inp)

    runs_folder = os.path.join(root, 'simulations', 'runs')
    tz0 = 298
    p_air = 101325
    r_air = 287
    for day in days:
        g_prime_sources = []
        rho_gas_sources = []
        twodee_sources = []
        disgas_sources = []
        index_sources_twodee = []
        index_sources_disgas = []
        diagno_output_file = os.path.join(root, 'simulations', 'diagno', day, 'diagno.out')
        diagno_input_file = os.path.join(root, 'simulations', 'diagno', day, 'diagno.inp')
        with open(diagno_input_file, 'r') as diagno_input_data:
            for line in diagno_input_data:
                if 'NX\n' in line:
                    nx_diagno = int(line.split('NX')[0])
                if 'NY\n' in line:
                    ny_diagno = int(line.split('NY')[0])
                elif 'NZ\n' in line:
                    nz_diagno = int(line.split('NZ')[0])
                elif 'TINF\n' in line:
                    tz0 = float(line.split('TINF')[0])
        rho_air = p_air / (r_air * tz0)
        for i_source in range(0, len(gas_fluxes_sources)):
            t_source = temperatures_sources[i_source]
            rho_gas_sources.append(p_air / (r_tracking_specie * t_source))
            g_prime_sources.append((9.81 * (rho_gas_sources[-1] - rho_air)) / rho_air)
        f = FortranFile(diagno_output_file, 'r')
        vx_time_average_sources = [0 for flux in gas_fluxes_sources]
        vy_time_average_sources = [0 for flux in gas_fluxes_sources]
        n_times_source_counted = [0 for flux in gas_fluxes_sources]
        n_steps = 0
        while True:
            try:
                vx_average_sources = [0 for flux in gas_fluxes_sources]
                vy_average_sources = [0 for flux in gas_fluxes_sources]
                for nline in range(0, 2):
                    f.read_reals(dtype='float32')
                n_steps += 1
                xor_diagno, yor_diagno = f.read_reals(dtype='float32')
                f.read_reals(dtype='float32')
                dx_diagno, dy_diagno = f.read_reals(dtype='float32')
                z_diagno = f.read_reals(dtype='float32')
                min_diff = 100000000
                k_z = 0
                for k in range(0, len(z_diagno)):
                    z_diff = abs(10 - z_diagno[k])
                    if z_diff <= min_diff:
                        min_diff = z_diff
                        k_z = k
                if k_z == 0:
                    print('WARNING. Unable to find level z = 10 m in diagno.out. Reading the wind output from the '
                          'first level above the ground')
                for nline in range(0, 3):
                    f.read_reals(dtype='float32')
                for k in range(0, nz_diagno):
                    for j in range(0, ny_diagno):
                        y_domain = yor_diagno + j * dy_diagno
                        vx_vector = f.read_reals(dtype='float32')
                        for i_v in range(0, len(vx_vector)):
                            x_domain = xor_diagno + i_v * dx_diagno
                            for i_source in range(0, len(gas_fluxes_sources)):
                                x_source = easting_sources[i_source]
                                y_source = northing_sources[i_source]
                                if k == k_z and x_domain - dx_diagno < x_source < x_domain + dx_diagno and \
                                        y_domain - dy_diagno < y_source < y_domain + dy_diagno:
                                    vx_average_sources[i_source] += vx_vector[i_v]
                                    n_times_source_counted[i_source] += 1
                vx_time_average_sources = [vx_time_average_sources[i_source] + vx_average_sources[i_source] /
                                           n_times_source_counted[i_source]
                                           for i_source in range(0, len(gas_fluxes_sources))]
                n_times_source_counted = [0 for flux in gas_fluxes_sources]
                for k in range(0, nz_diagno):
                    for j in range(0, ny_diagno):
                        y_domain = yor_diagno + j * dy_diagno
                        vy_vector = f.read_reals(dtype='float32')
                        for i_v in range(0, len(vy_vector)):
                            x_domain = xor_diagno + i_v * dx_diagno
                            for i_source in range(0, len(gas_fluxes_sources)):
                                x_source = easting_sources[i_source]
                                y_source = northing_sources[i_source]
                                if k == k_z and x_domain - dx_diagno < x_source < x_domain + dx_diagno and \
                                        y_domain - dy_diagno < y_source < y_domain + dy_diagno:
                                    vy_average_sources[i_source] += vy_vector[i_v]
                                    n_times_source_counted[i_source] += 1
                vy_time_average_sources = [vy_time_average_sources[i_source] + vy_average_sources[i_source] /
                                           n_times_source_counted[i_source] for i_source in
                                           range(0, len(gas_fluxes_sources))]
                n_times_source_counted = [0 for flux in gas_fluxes_sources]
                # Skip vz
                for k in range(0, nz_diagno):
                    for j in range(0, ny_diagno):
                        f.read_reals(dtype='float32')
            except BaseException:
                break
        vx_time_average_sources = [vx_time_average_sources[i_source] / n_steps
                                   for i_source in range(0, len(vx_time_average_sources))]
        vy_time_average_sources = [vy_time_average_sources[i_source] / n_steps
                                   for i_source in range(0, len(vy_time_average_sources))]
        wind_time_average_sources = [(vx_time_average_sources[i_source] ** 2 +
                                      vy_time_average_sources[i_source] ** 2) ** 0.5
                                     for i_source in range(0, len(vx_time_average_sources))]
        q_average_sources = [gas_fluxes_sources[i_source] / rho_gas_sources[i_source]
                             for i_source in range(0, len(gas_fluxes_sources))]
        for i_source in range(0, len(wind_time_average_sources)):
            if twodee_on:
                area_source = dx_sources[i_source] * dy_sources[i_source]
            else:
                area_source = dx * dy
            r_source = (area_source / 3.1416) ** 0.5
            try:
                ri_source = (1 / wind_time_average_sources[i_source] ** 2) * ((g_prime_sources[i_source] *
                             q_average_sources[i_source]) / r_source) ** (2 / 3)
            except ValueError:  # When gprime <0, for the moment we set Ri = 0.01 (small number)
                ri_source = 0.01
            if ri_source > 0.25: # FABIO: less restrictive condition
                twodee_sources.append(ri_source)
                index_sources_twodee.append(i_source)
            else:
                disgas_sources.append(ri_source)
                index_sources_disgas.append(i_source)
        if len(index_sources_twodee) != 0 and len(index_sources_disgas) != 0:  # The run needs to be split
            for run in runs_disgas:
                if day in run:
                    runs_disgas.remove(run)
                    break
            for run in runs_twodee:
                if day in run:
                    runs_twodee.remove(run)
                    break
            prepare_temporary_simulations(day)
        elif len(index_sources_twodee) == 0 and len(index_sources_disgas) != 0:
            for run in runs_twodee:
                if day in run:
                    runs_twodee.remove(run)
                    break
        else:
            for run in runs_disgas:
                if day in run:
                    runs_disgas.remove(run)
                    break


def run_simulations(max_np):
    import datetime
    n_elaborated_runs = 0
    n_node = 0
    node = ''
    if continuous_simulation:
        max_np = 1
    while n_elaborated_runs <= len(runs):
        ps = []
        start_run = n_elaborated_runs
        end_run = n_elaborated_runs + max_np
        if end_run > len(runs):
            end_run = len(runs)
        try:
            for run in runs[start_run:end_run + 1]:
                if 'disgas' in run:
                    solver = 'disgas'
                elif 'twodee' in run:
                    solver = 'twodee'
                else:
                    print('ERROR. Unable to identify the solver to use in run_simulations')
                    sys.exit()
                day = run.split(os.sep)[-1]
                if len(nodes_list) > 0:
                    try:
                        node = nodes_list[n_node]
                    except IndexError:
                        node = ''
                input_file = os.path.join(run, solver + '.inp')
                log_file = os.path.join(run, solver + '_log.txt')
                if continuous_simulation and day != days[0]:
                    current_day = datetime.datetime.strptime(day, '%Y%m%d')
                    previous_day = current_day - datetime.timedelta(days=1)
                    previous_day = previous_day.strftime('%Y%m%d')
                    previous_day_folder = os.path.join(root, 'simulations', 'runs', str(previous_day))
                    try:
                        shutil.copyfile(os.path.join(previous_day_folder, 'restart.dat'),
                                        os.path.join(run, 'restart.dat'))
                    except FileNotFoundError:
                        print('ERROR. Restart run requested but file ' +
                              os.path.join(previous_day_folder, 'restart.dat') + ' not found')
                        sys.exit()
                if slurm:
                    try:
                        p = subprocess.Popen(['srun', '-n', '1', '--partition=' + partition, '--nodelist=' + node,
                                              solver, input_file, log_file])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' --nodelist=' + node + ' ' + solver
                              + ' ' + input_file + ' ' + log_file)
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen([solver, input_file, log_file])
                    except FileNotFoundError:
                        print('Unable to run ' + solver.upper())
                        sys.exit()
                ps.append(p)
                if len(nodes_list) > 0:
                    n_node += 1
                    if n_node >= len(nodes_list):
                        n_node = 0
        except BaseException:
            print('Unable to run the simulations')
            sys.exit()
        n_elaborated_runs = end_run
        if n_elaborated_runs == len(days):
            break
        for p in ps:
            p.wait()
    for p in ps:
        p.wait()


def converter(run_in):
    def convert_to_ppm(input_file, output_file):
        import datetime
        splitted_path = run_in.split(os.sep)
        day = splitted_path[-2]
        met_data = run_in.split('runs')[0]
        met_data = os.path.join(met_data, 'diagno', day)
        first_records = []
        with open(input_file, 'r') as input_file_read:
            i_line = 1
            for line in input_file_read:
                first_records.append(line)
                if i_line > 3:
                    break
                else:
                    i_line += 1
        c_kgm3 = np.loadtxt(input_file, skiprows=5)
        c_kgm3[c_kgm3 < 0] = 0
        file_time_step = os.path.split(output_file)[1]
        file_time_step = file_time_step.split('_')[2]
        file_time_step = file_time_step.split('.grd')[0]
        file_time_step = int(file_time_step[-4:])
        # To be generalized if time step < 1 h or if duration > 1 day
        file_validity = output_interval * file_time_step
        if file_validity > 23:
            file_validity = 23  # Because diagno outputs do not go beyond 23 hours
        file_validity = datetime.timedelta(hours=file_validity)         
        surface_data = os.path.join(met_data, 'surface_data.txt')
        with open(surface_data) as surface_data_file:
            for line in surface_data_file:
                try:
                    records = line.split('\t')
                    if file_validity == datetime.timedelta(hours=int(records[0][0:2])):  # To be generalized if
                        # time step < 1h
                        t2m = float(records[2])
                        p2m = float(records[3]) / 100  # in hPa for this conversion
                        break
                except ValueError:
                    continue
        conversion_factor = ((22.4 / m_tracking_specie) * (t2m / 273) * (1013 / p2m)) * 1000000
        c_ppm = c_kgm3 * conversion_factor  # convert kg/m3 to ppm
        with open(output_file, 'w') as converted_output_file:
            for record in first_records:
                converted_output_file.write(record)
            converted_output_file.write(str(np.amin(c_ppm)) + '  ' + str(np.amax(c_ppm)) + '\n')
            np.savetxt(converted_output_file, c_ppm, fmt='%.5e')

    outfiles_folder = os.path.join(run_in, 'outfiles')
    temp_folder = os.path.join(outfiles_folder, 'temp')
    try:
        os.mkdir(temp_folder)
    except FileExistsError:
        shutil.rmtree(temp_folder)
        os.mkdir(temp_folder)
    files_to_convert = []
    converted_files = []
    for file in os.listdir(outfiles_folder):
        if 'c_' in file and '000000' not in file:
            shutil.move(os.path.join(outfiles_folder, file), os.path.join(temp_folder, file))
            files_to_convert.append(os.path.join(temp_folder, file))
            converted_files.append(os.path.join(outfiles_folder, file))
    for i in range(0, len(converted_files)):
        convert_to_ppm(files_to_convert[i], converted_files[i])
    shutil.rmtree(temp_folder)


def match_grid(run_in):
    from scipy.interpolate import RegularGridInterpolator
    x_disgas = np.linspace(bottom_left_easting + dx / 2, top_right_easting - dx / 2, nx)
    y_disgas = np.linspace(bottom_left_northing + dy / 2, top_right_northing - dy / 2, ny)
    x_grd_disgas, y_grd_disgas = np.meshgrid(x_disgas, y_disgas)
    x_twodee = np.linspace(bottom_left_easting - dx, top_right_easting + dx, nx)
    y_twodee = np.linspace(bottom_left_northing - dy, top_right_northing + dy, ny)
    splitted_run_in = run_in.split(os.path.sep)
    try:
        splitted_run_in.remove('')
        splitted_run_in.remove('disgas')
        splitted_run_in.remove('twodee')
    except ValueError:
        pass
    outfiles_folder = os.path.sep
    for element in splitted_run_in:
        outfiles_folder += element + os.path.sep
    outfiles_folder = os.path.join(outfiles_folder, 'outfiles')
    temp_folder = os.path.join(outfiles_folder, 'temp')
    try:
        os.mkdir(temp_folder)
    except FileExistsError:
        shutil.rmtree(temp_folder)
        os.mkdir(temp_folder)
    files_to_convert = []
    converted_files = []
    for file in os.listdir(outfiles_folder):
        if 'c_' in file:
            shutil.move(os.path.join(outfiles_folder, file), os.path.join(temp_folder, file))
            files_to_convert.append(os.path.join(temp_folder, file))
            converted_files.append(os.path.join(outfiles_folder, file))
    for i in range(0, len(converted_files)):
        twodee_output = np.transpose(np.loadtxt(files_to_convert[i], skiprows=5))
        interp = RegularGridInterpolator((x_twodee, y_twodee), twodee_output)
        twodee_new_grid = interp((x_grd_disgas, y_grd_disgas))
        with open(converted_files[i], 'w') as twodee_interpolated_output:
            twodee_interpolated_output.write('DSAA\n')
            twodee_interpolated_output.write(str(nx) + '  ' + str(ny) + '\n')
            twodee_interpolated_output.write(str(bottom_left_easting + dx / 2) + '  ' + str(top_right_easting - dx / 2)
                                             + '\n')
            twodee_interpolated_output.write(str(bottom_left_northing + dy / 2) + '  ' + str(top_right_northing - dy /
                                             2) + '\n')
            twodee_interpolated_output.write(str(np.min(twodee_new_grid)) + '   ' + str(np.max(twodee_new_grid)) + '\n')
            np.savetxt(twodee_interpolated_output, twodee_new_grid, fmt='%.4e')
    shutil.rmtree(temp_folder)


def elaborate_outputs():
    for run in runs_disgas:
        splitted_run_in = run.split(os.path.sep)
        try:
            splitted_run_in.remove('')
            splitted_run_in.remove('disgas')
        except ValueError:
            pass
        outfiles_folder = os.path.sep
        for element in splitted_run_in:
            outfiles_folder += element + os.path.sep
        twodee_folder = os.path.join(outfiles_folder, 'twodee')
        if os.path.exists(twodee_folder):
            twodee_subfolders = os.listdir(twodee_folder)
            if 'outfiles' not in twodee_subfolders:
                shutil.rmtree(twodee_folder)
        outfiles_folder = os.path.join(outfiles_folder, 'outfiles')
        disgas_outfiles_folder = os.path.join(run, 'outfiles')
        shutil.copytree(disgas_outfiles_folder, outfiles_folder)
        if len(runs_twodee) == 0:
            shutil.rmtree(disgas_outfiles_folder)
    for run in runs_twodee:
        splitted_run_in = run.split(os.path.sep)
        try:
            splitted_run_in.remove('')
            splitted_run_in.remove('twodee')
        except ValueError:
            pass
        outfiles_folder = os.path.sep
        for element in splitted_run_in:
            outfiles_folder += element + os.path.sep
        disgas_folder = os.path.join(outfiles_folder, 'disgas')
        if os.path.exists(disgas_folder):
            disgas_subfolders = os.listdir(disgas_folder)
            if 'outfiles' not in disgas_subfolders:
                shutil.rmtree(disgas_folder)
        outfiles_folder = os.path.join(outfiles_folder, 'outfiles')
        twodee_outfiles_folder = os.path.join(run, 'outfiles')
        try:
            shutil.copytree(twodee_outfiles_folder, outfiles_folder)
            for output_file in os.listdir(outfiles_folder):
                if 'c_' not in output_file: # remove all other twodee output files
                    os.remove(os.path.join(outfiles_folder, output_file))
        except FileExistsError:
            # We are in a split simulation
            disgas_outfiles_folder = os.path.join(disgas_folder, 'outfiles')
            output_heights_v_str = output_heights.split(' ')
            output_heights_v = []
            for output_height in output_heights_v_str:
                try:
                    output_heights_v.append(float(output_height))
                except ValueError:
                    continue
            for i_level in range(1, len(output_heights_v) + 1):
                level_disgas = str(i_level).zfill(3)
                level_twodee = str(int(output_heights_v[i_level - 1] * 100)).zfill(4) + 'cm'
                for i_time in range(1, 25):
                    time_disgas = str(i_time).zfill(6)
                    time_twodee = str(i_time * output_interval * 3600).zfill(6)
                    original_disgas_file = os.path.join(disgas_outfiles_folder, 'c_' + level_disgas + '_' + time_disgas
                                                        + '.grd')
                    original_twodee_file = os.path.join(twodee_outfiles_folder, 'c_' + level_twodee + '_' + time_twodee
                                                        + '.grd')
                    # For the moment the merged file has the same name formatting of disgas outputs
                    merged_output_file = os.path.join(outfiles_folder, 'c_' + level_disgas + '_' + time_disgas + '.grd')
                    first_records = []
                    with open(original_disgas_file, 'r') as input_file_read:
                        i_line = 1
                        for line in input_file_read:
                            first_records.append(line)
                            if i_line > 3:
                                break
                            else:
                                i_line += 1
                    c_disgas = np.loadtxt(original_disgas_file, skiprows=5)
                    c_twodee = np.loadtxt(original_twodee_file, skiprows=5)
                    c_merged = np.add(c_disgas, c_twodee)
                    with open(merged_output_file, 'w') as new_merged_output_file:
                        for record in first_records:
                            new_merged_output_file.write(record)
                        new_merged_output_file.write(str(np.amin(c_merged)) + '  ' + str(np.amax(c_merged)) + '\n')
                        np.savetxt(new_merged_output_file, c_merged, fmt='%.5e')
        shutil.rmtree(twodee_outfiles_folder)


root = os.getcwd()
disgas_original = os.path.join(root, 'disgas.inp')
twodee_original = os.path.join(root, 'twodee.inp')
topography = os.path.join(root, 'topography.grd')

(
    run_type,
    continuous_simulation,
    max_number_processes,
    random_sources,
    nsources,
    sources_interval,
    source_easting,
    source_northing,
    source_el,
    source_emission,
    random_emission,
    bottom_left_northing,
    bottom_left_easting,
    top_right_northing,
    top_right_easting,
    nx,
    ny,
    nz,
    dx,
    dy,
    source_dx,
    source_dy,
    source_dur,
    diagno_on,
    twodee_on,
    disgas_on,
    slurm,
    partition,
    tracking_specie,
    run_duration,
    output_interval,
    output_heights,
) = read_arguments()


if disgas_on == 'off' and twodee_on == 'off':
    print('Both DISGAS and TWODEE are turned off')
    sys.exit()

if slurm:
    list_available_nodes = []
    result = subprocess.run(['sinfo', '-o=%N'], stdout=subprocess.PIPE)
    sinfo_node_output = result.stdout.decode('utf-8')
    node_key = sinfo_node_output.split('=')[2]
    node_key = node_key.split('[')[0]
    result = subprocess.run(['sinfo'], stdout=subprocess.PIPE)
    sinfo_output = result.stdout.decode('utf-8').split('\n')
    nodes_ranges_temp = []
    nodes_ranges = []
    for output in sinfo_output:
        if partition in output and 'idle' in output:
            try:
                nodes_list_formatted = output.split('idle ')[1]
            except IndexError:
                nodes_list_formatted = output.split('idle~ ')[1]
            nodes_list = nodes_list_formatted.split(node_key)[1]
            if '[' in nodes_list:
                nodes_list = nodes_list.split('[')[1]
                nodes_list = nodes_list.split(']')[0]
                nodes_list = nodes_list.split(',')
                for nodes_range in nodes_list:
                    nodes_ranges.append(nodes_range)
            else:
                nodes_ranges.append(nodes_list)
    for node_range in nodes_ranges:
        node_range = node_range.split('-')
        start = int(node_range[0])
        try:
            end = int(node_range[1])
        except IndexError:
            end = start
        for i_node in range(start, end + 1):
            list_available_nodes.append(node_key + '{:0>2}'.format(i_node))
    result = subprocess.run(['sinfo', '-o=%c'], stdout=subprocess.PIPE)
    sinfo_core_output = result.stdout.decode('utf-8')
    ncpus_per_node = int(sinfo_core_output.split('=')[-1])
    n_nodes = - (-max_number_processes // ncpus_per_node)
    nodes_list = list_available_nodes[0:n_nodes]
else:
    nodes_list = []

gas_properties_file = os.path.join(root, 'gas_properties.csv')
try:
    open(gas_properties_file)
except FileNotFoundError:
    print('ERROR. File gas_properties.csv not found')
    sys.exit()

days = prepare_days()

runs_disgas = []
runs_twodee = []

easting_sources, northing_sources, elevation_sources, dx_sources, dy_sources, dur_sources, gas_fluxes_sources, \
    temperatures_sources, rho_tracking_specie, r_tracking_specie, m_tracking_specie = pre_process(run_type)

if diagno_on:
    max_height_diagno = run_diagno(max_number_processes)

if disgas_on and twodee_on:  # We are in automatic scenario detection mode, hence we need to read the meteorological
    # data at the source locations to calculate the Richardson number at the source and assign the right solver to each
    # simulation, for this reason we re-initialize the runs vectors
    read_diagno_outputs()

if len(runs_disgas) == 0:
    disgas_on = False
elif len(runs_twodee) == 0:
    twodee_on = False

runs = runs_disgas + runs_twodee
run_simulations(max_number_processes)

if disgas_on:
    for simulation in runs_disgas:
        converter(simulation)  # Convert disgas output in ppm to be consistent with twodee outputs in case of merging
        # and/or automatic scenario

if twodee_on:
    if disgas_on:
        for simulation in runs_twodee:
            match_grid(simulation)  # Interpolate twodee output onto disgas grid, this is needed for the following
            # interpolation

elaborate_outputs()  # To copy each outfiles folder into the general outfiles folder and to merge simulation outputs
# when runs have been split
