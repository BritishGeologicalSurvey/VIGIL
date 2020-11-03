from random import sample
import os
import shutil
import subprocess
import argparse
import numpy as np
import sys
import utm

def read_arguments():
    parser = argparse.ArgumentParser(description='Input data')
    parser.add_argument('-N', '--nproc', default=1, help='Maximum number of allowed simultaneous processes')
    parser.add_argument('-RS','--random_sources', default='off', help='on: randomly select NS locations from a probability map. off: fixed source locations from file sources_input.txt')
    parser.add_argument('-NS','--nsources',default='random', help='Specify a number for a fixed number of sources. If random, then randomly select the number of sources from an interval')
    parser.add_argument('-SINT','--sources_interval', nargs='+', default=[], help='Type the minimum and maximum number of sources')
    parser.add_argument('-SLOC','--source_location',nargs='+', default=[], help='Coordinate type (UTM/GEO), latitude/northing, longitude/easting, elevation (above ground in m) of 1 fixed source')
    parser.add_argument('-SDX','--source_dx', default=999999, help = 'Extension [m] along the X direction of 1 single source. Option valid for Twodee only')
    parser.add_argument('-SDY', '--source_dy', default=999999, help='Extension [m] along the Y direction of 1 single source. Option valid for Twodee only')
    parser.add_argument('-SDUR', '--source_dur', default=0, help='Emission duration [s] of 1 single source. Option valid for Twodee only')
    parser.add_argument('-D','--domain',nargs='+', default=[], help='Coordinates type (UTM/GEO), coordinates (latitude/northing, longitude/easting) of the bottom left corner and top right corner of the domain')
    parser.add_argument('-SEM','--source_emission',default='999',help='Source emission rate [kg/s]. If specified, it is assigned to all the sources in the domain')
    parser.add_argument('-RER','--random_emission',default='off',help='on: randomly assign emission rate for each source in the domain sampled from a flux.csv file. off: use specified emission rate')
    parser.add_argument('-TD', '--twodee', default='off',help='on or off, to run Twodee')
    parser.add_argument('-DG', '--disgas', default='off', help='on or off, to run Disgas')
    args = parser.parse_args()
    nproc = args.nproc
    random_sources = args.random_sources
    nsources = args.nsources
    source_location = args.source_location
    sources_interval = args.sources_interval
    domain = args.domain
    source_emission = args.source_emission
    random_emission = args.random_emission
    max_number_processes = int(nproc)
    twodee = args.twodee
    disgas = args.disgas
    source_dx = args.source_dx
    source_dy = args.source_dy
    source_dur = args.source_dur
    source_easting = source_northing = source_el = 0
    try:
        source_emission = float(source_emission)
    except:
        print('Pleae provide a valid number for the emission rate of the source')
        sys.exit()
    if len(domain) == 0 or len(domain) > 5:
        print('ERROR. Please provide valid entries for -D --domain')
        sys.exit()
    else:
        coordinates_type = domain[0]
        if coordinates_type == 'GEO':
            bottom_left_1 = float(domain[1])
            bottom_left_2 = float(domain[2])
            top_right_1 = float(domain[3])
            top_right_2 = float(domain[4])
            if (-90 <= bottom_left_1 <= 90 and -180 <= bottom_left_2 <= 180) and (
                    -90 <= top_right_1 <= 90 and -180 <= top_right_2 <= 180):  # identify valid geographic coordinates
                try:
                    out_utm = utm.from_latlon(bottom_left_1, bottom_left_2)
                    bottom_left_easting = float(out_utm[0])
                    bottom_left_northing = float(out_utm[1])
                except:
                    print('ERROR. Please provide valid coordinates for the bottom left corner of the domain')
                    sys.exit()
                try:
                    out_utm = utm.from_latlon(top_right_1, top_right_2)
                    top_right_easting = float(out_utm[0])
                    top_right_northing = float(out_utm[1])
                except:
                    print('ERROR. Please provide valid coordinates for the top right corner of the domain')
                    sys.exit()
            else:
                print('ERROR. Please provide valid coordinates')
                sys.exit()
        elif coordinates_type == 'UTM':
            bottom_left_northing = float(domain[1])
            bottom_left_easting = float(domain[2])
            top_right_northing = float(domain[3])
            top_right_easting = float(domain[4])
        else:
            print('ERROR. Please provide a valide type of coordinates (UTM or GEO)')
            sys.exit()
        if bottom_left_northing == top_right_northing or bottom_left_easting == top_right_easting:
            print('ERROR. Coordinates of the corners cannot coincide')
            sys.exit()
        if bottom_left_northing > top_right_northing and bottom_left_easting > top_right_easting:  # Check coordinates are in the proper order, otherwise swap
            temp = bottom_left_northing
            bottom_left_northing = top_right_northing
            top_right_northing = temp
            temp = bottom_left_easting
            bottom_left_easting = top_right_easting
            top_right_easting = temp
    if random_sources == 'on':
        try:
            np.loadtxt('probability_map.txt')
        except:
            print('Please provide a valid probability_map.txt file when random_sources option is on')
            sys.exit()
        if nsources == 'random':
            if len(sources_interval) == 0 or len(sources_interval) > 2:
                print('ERROR. Please specify the minimum and maximum number of sources with -SINT --sources_interval')
                sys.exit()
        else:
            try:
                nsources = int(nsources)
            except:
                print('Please provide a valid integer for -NS --nsources')
                sys.exit()
        if random_emission == 'off' and source_emission == 999:
            print('ERROR. random_sources set to on requires either random_emission set to on or a specified source_emission')
            sys.exit()
    else:
        if random_sources != 'off':
            print('Valid options for -RS --random_sources are on and off')
            sys.exit()
        else:
            try:
                sources_file =  open('sources_input.txt','r')
                sources_file.close()
            except:
                print('File sources_input.txt not found. Using one source from input data')
                if len(source_location) == 0 or len(source_location) > 4:
                    print('ERROR. Please provide valid entries for -SLOC --sources_location')
                    sys.exit()
                else:
                    coordinates_type = source_location[0]
                    if coordinates_type == 'GEO':
                        if -90 <= float(source_location[1]) <= 90 and -180 <= float(source_location[2]) <= 180: #identify geographic coordinates
                            try:
                                out_utm = utm.from_latlon(float(source_location[1]), float(source_location[2]))
                                source_easting = float(out_utm[0])
                                source_northing = float(out_utm[1])
                            except:
                                print('Please provide valid coordinates of the source location')
                                sys.exit()
                    elif coordinates_type == 'UTM':
                        source_easting = float(source_location[1])
                        source_northing = float(source_location[0])
                        if not bottom_left_easting <= source_easting <= top_right_easting or not bottom_left_northing <= source_northing <= top_right_northing:
                            print('Location not within the domain')
                            sys.exit()
                    else:
                        print('ERROR. Please provide a valide type of coordinates (UTM or GEO)')
                        sys.exit()
                    if float(source_location[2]) < 0:
                        print('Please provide a valid value for the source elevation in m above ground (>= 0 m)')
                        sys.exit()
                    else:
                        source_el = float(source_location[2])
    if random_emission == 'on':
        try:
            sources_file = open('flux.csv', 'r')
            sources_file.close()
        except:
            print('ERROR. File flux.csv not found')
            sys.exit()
    elif random_emission != 'off':
        print('Valid options for -RER --random_sources are on and off')
        sys.exit()
    if twodee.lower() == 'on':
        twodee_on = True
    elif twodee.lower() == 'off':
        twodee_on = False
    else:
        print('Please provide a valid entry for the variable -TD --twodee')
        sys.exit()
    if disgas.lower() == 'on':
        disgas_on = True
    elif disgas.lower() == 'off':
        disgas_on = False
    else:
        print('Please provide a valid entry for the variable -DG --disgas')
        sys.exit()
    if source_dx != 999999:
        try:
            source_dx = float(source_dx)
        except:
            print('Please provide a valid number for the variable -SDX --source_dx')
            sys.exit()
    if source_dy != 999999:
        try:
            source_dy = float(source_dy)
        except:
            print('Please provide a valid number for the variable -SDY --source_dy')
            sys.exit()
    if source_dur != 0:
        try:
            source_dur = float(source_dur)
        except:
            print('Please provide a valid number for the variable -SDUR --source_dur')
            sys.exit()
    return max_number_processes, random_sources, nsources, sources_interval, source_easting, source_northing, source_el, \
           source_emission, random_emission, bottom_left_northing, bottom_left_easting, top_right_northing, top_right_easting, \
           source_dx, source_dy, source_dur, twodee_on, disgas_on

def pre_process():
    def sample_random_sources(n_sources, input_file, xmin, xmax, ymin, ymax, dur_min, dur_max, source_size_min, source_size_max):
        from random import choices

        probabilities_input = np.loadtxt(input_file)
        nx = probabilities_input.shape[0]
        ny = probabilities_input.shape[1]
        location_cum_indexes = []
        location_indexes = []
        probabilities = []
        dxs = []
        dys = []
        durs = []
        k = 0
        for i in range(0, nx):
            for j in range(0, ny):
                location_cum_indexes.append(k)
                location_indexes.append([i, j])
                k += 1
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        for i in range(0, nx):
            for j in range(0, ny):
                probabilities.append(probabilities_input[i, j])
        selected_locations = choices(location_cum_indexes, probabilities, k=n_sources)
        source_size = source_size_min
        while source_size <= source_size_max:
            dxs.append(source_size)
            dys.append(source_size)
            source_size += 0.1
        source_duration = dur_min
        while source_duration <= dur_max:
            durs.append(source_duration)
            source_duration += 60
        for location in selected_locations:
            row = location_indexes[location][0]
            column = location_indexes[location][1]
            xpr = x[row]
            ypr = y[column]
            probability = probabilities[location]
            random_eastings.append(xpr)
            random_northings.append(ypr)
            random_elevations.append(0.0)
            random_probabilities.append(probability)
            random_fluxes.append(99999999)
            random_dx.append(choices(dxs, k = 1)[0])
            random_dy.append(choices(dys, k = 1)[0])
            random_dur.append(choices(durs, k = 1)[0])

        return random_eastings, random_northings, random_elevations, random_probabilities, random_elevations, random_fluxes, random_dx, random_dy, random_dur

    def fluxes():
        import numpy as np
        import pandas as pd

        data = pd.read_csv('flux.csv', error_bad_lines=False)
        x = np.sort(data['flux'])
        y = np.arange(1, len(x) + 1) / len(x)
        list_x = list(x)
        sampled_flux = (sample(list_x, 1))
        return sampled_flux

    def topography_converter(original_topography_file):
        # To convert the DIAGNO-DISGAS topography file into a Twodee-readable topography file
        with open(original_topography_file, 'r', encoding='utf-8') as diagno_grid:
            line_number = 1
            lines = []
            for line in diagno_grid:
                lines.append(line)
        NX = int(lines[1].split(' ')[0])
        NY = int(lines[1].split(' ')[1].split('\n')[0])
        X0 = float(lines[2].split(' ')[0])
        XF = float(lines[2].split(' ')[1].split('\n')[0])
        Y0 = float(lines[3].split(' ')[0])
        DX = round(((XF - X0) / NX), 1)
        with open(os.path.join(infiles_twodee,'topography.grd'), 'w', encoding='utf-8') as twodee_grid:
            twodee_grid.write('nx ' + str(NX) + '\n')
            twodee_grid.write('ny ' + str(NY) + '\n')
            twodee_grid.write('xllcorner ' + str(X0) + '\n')
            twodee_grid.write('yllcorner ' + str(Y0) + '\n')
            twodee_grid.write('dx ' + str(DX) + '\n')
            twodee_grid.write('no 123\n')
            for line in range(5, len(lines)):
                twodee_grid.write(lines[line])

    # Set source files for DISGAS and TWODEE. TWODEE also needs source start and stop time, to be generalized
    random_eastings = []
    random_northings = []
    random_elevations = []
    random_probabilities = []
    random_fluxes = []
    random_dx = []
    random_dy = []
    random_dur = []
    dur_min = 1
    dur_max = 86400
    source_size_min = 0.1
    source_size_max = 100
    if random_sources == 'on':
        if nsources == 'random':
            Nsources = [*range(sources_interval[0],sources_interval[1] + 1)]
        else:
            Nsources = [nsources]
        n_random_sources = sample(Nsources, 1)[0]
        random_eastings, random_northings, random_elevations, random_probabilities, random_elevations, random_fluxes, \
        random_dx, random_dy, random_dur = sample_random_sources(n_random_sources, 'probability_map.txt', bottom_left_easting, top_right_easting,
                                                                 bottom_left_northing, top_right_northing, dur_min, dur_max, source_size_min, source_size_max)
    easting = random_eastings
    northing = random_northings
    elevations = random_elevations
    probabilities = random_probabilities
    fluxes_input = random_fluxes
    dx = random_dx
    dy = random_dy
    dur = random_dur
    n_sources = 0
    try:
        with open('sources_input.txt', 'r', encoding="utf-8", errors="surrogateescape") as locations_file:
            for line in locations_file:
                try:
                    records = line.split('\t')
                    easting.append(float(records[0]))
                    northing.append(float(records[1]))
                    elevations.append(float(records[2]))
                    probabilities.append(float(records[3]))
                    fluxes_input.append(float(records[4]))
                    dx.append(float(records[5]))
                    dy.append(float(records[6]))
                    dur.append(float(records[7]))
                    n_sources += 1
                except:
                    continue
    except:
        easting.append(source_easting)
        northing.append(source_northing)
        elevations.append(source_el)
        probabilities.append(1.0)
        fluxes_input.append(source_emission)
        dx.append(source_dx)
        dy.append(source_dy)
        dur.append(source_dur)
    n_sources = len(easting)

    raw_days = [] # store the days as originally formatted
    days = [] #store days in format YYYYMMDD as per folder name
    # Set DIAGNO folder
    diagno = os.path.join(root, 'simulations', 'diagno')
    try:
        os.mkdir(diagno)
    except:
        print('Folder ' + diagno + ' already exists')
    if disgas_on:
        # Set DISGAS folder
        disgas = os.path.join(root, 'simulations', 'disgas')
        try:
            os.mkdir(disgas)
        except:
            print('Folder ' + disgas + ' already exists')
    if twodee_on:
        # Set DISGAS folder
        twodee = os.path.join(root, 'simulations', 'twodee')
        try:
            os.mkdir(twodee)
        except:
            print('Folder ' + twodee + ' already exists')
    # read days_list file
    with open('days_list.txt','r',encoding="utf-8", errors="surrogateescape") as days_list_file:
        for line in days_list_file:
            raw_days.append(line)
    days_list_file.close()
    i=0
    for day in raw_days:
        temp = raw_days[i].split(' ')
        temp = temp[0].split('-')
        days.append(temp[0]+temp[1]+temp[2])
        i+=1
    for day in days:
        path = os.path.join(root,'simulations', str(day))
        files = os.listdir(path)
        diagno_daily = os.path.join(diagno, str(day))
        try:
            os.mkdir(diagno_daily)
        except:
            print('Folder ' + diagno_daily + ' already exists')
        for f in files:
            path_f = os.path.join(path,f)
            try:
                shutil.move(path_f, diagno_daily)
            except:
                print('File ' + f + ' already present in ' + diagno)
        shutil.copy(topography_original, os.path.join(diagno_daily, 'topography.grd'))
        # Set DISGAS folder
        if disgas_on:
            disgas_daily = os.path.join(disgas, str(day))
            infiles = os.path.join(disgas_daily, 'infiles')
            outfiles = os.path.join(disgas_daily, 'outfiles')
            try:
                os.mkdir(disgas_daily)
            except:
                print('Folder ' + disgas_daily + ' already exists')
            try:
                os.mkdir(infiles)
            except:
                print('Folder infiles already exists in ' + str(disgas_daily))
            try:
                os.mkdir(outfiles)
                if not outfiles.endswith(os.path.sep):
                    outfiles += os.path.sep
            except:
                print('Folder outfiles already exists in ' + str(disgas_daily))
            disgas_input = os.path.join(infiles,'disgas.inp')
            with open(os.path.join(infiles,'source.dat'), 'w', encoding="utf-8", errors="surrogateescape") as source_file:
                for i in range(0, n_sources):
                    if source_emission != 999:
                        gas_flux = source_emission
                    else:
                        if fluxes_input[i] == 99999999 and random_emission == 'on':
                            gas_flux = fluxes()[0]
                        else:
                            gas_flux = fluxes_input[i]
                    source_file.write('{0:7.3f}'.format(easting[i]) + ' ' + '{0:7.3f}'.format(northing[i]) + ' '
                                      + '{0:7.2f}'.format(elevations[i]) + ' ' + '{0:7.3f}'.format(gas_flux) + '\n')
            source_file.close()
            roughness_file_exist = True
            try:
                shutil.copyfile(os.path.join(root,'roughness_disgas.grd'),os.path.join(infiles,'roughness.grd'))
            except:
                roughness_file_exist = False
            # read and memorize disgas.inp file
            disgas_input_records = []
            with open(disgas_original, 'r', encoding="utf-8", errors="surrogateescape") as disgas_or_input:
                for line in disgas_or_input:
                    disgas_input_records.append(line)
            disgas_or_input.close()
            with open(disgas_input,'w', encoding="utf-8", errors="surrogateescape") as disgas_input:
                for record in disgas_input_records:
                    if 'ROUGHNESS_MODEL' in record:
                        roughness_command = record.split('=')[1]
                        roughness_command = roughness_command.split('(')[0]
                        if 'MATRIX' in roughness_command:
                            if not roughness_file_exist:
                                print('Warning! ROUGHNESS_MODEL set to MATRIX in disgas.inp and roughness file does not exist')
                                print('Setting ROUGHNESS_MODEL to UNIFORM')
                    if 'YEAR' in record:
                        disgas_input.write('  YEAR   = ' + day[0:4] + '\n')
                    elif 'MONTH' in record:
                        disgas_input.write('  MONTH  = ' + day[4:6] + '\n')
                    elif 'DAY' in record:
                        disgas_input.write('  DAY    = ' + day[6:8] + '\n')
                    elif 'TOPOGRAPHY_FILE_PATH' in record:
                        disgas_input.write('   TOPOGRAPHY_FILE_PATH   = ' + os.path.join(diagno_daily, 'topography.grd') + ' \n')
                    elif 'ROUGHNESS_FILE_PATH' in record:
                        disgas_input.write('   ROUGHNESS_FILE_PATH   = ' + os.path.join(infiles, 'roughness.grd') + ' \n')
                    elif 'RESTART_FILE_PATH' in record:
                        disgas_input.write('   RESTART_FILE_PATH   = ' + os.path.join(infiles, 'restart.dat') + ' \n')
                    elif 'SOURCE_FILE_PATH' in record:
                        disgas_input.write('   SOURCE_FILE_PATH   = ' + os.path.join(infiles, 'source.dat') + ' \n')
                    elif 'WIND_FILE_PATH' in record:
                        disgas_input.write('   WIND_FILE_PATH   = ' + os.path.join(infiles, 'winds.dat') + ' \n')
                    elif 'DIAGNO_FILE_PATH' in record:
                        disgas_input.write('   DIAGNO_FILE_PATH   = ' + os.path.join(diagno_daily, 'diagno.out') + ' \n')
                    elif 'OUTPUT_DIRECTORY' in record:
                        disgas_input.write('   OUTPUT_DIRECTORY    = ' + outfiles + ' \n')
                    else:
                        disgas_input.write(record)
            disgas_input.close()
        if twodee_on:
            twodee_daily = os.path.join(twodee, str(day))
            infiles_twodee = os.path.join(twodee_daily, 'infiles')
            outfiles_twodee = os.path.join(twodee_daily, 'outfiles')
            try:
                os.mkdir(twodee_daily)
            except:
                print('Folder ' + twodee_daily + ' already exists')
            try:
                os.mkdir(infiles_twodee)
            except:
                print('Folder infiles already exists in ' + str(twodee_daily))
            try:
                os.mkdir(outfiles_twodee)
                if not outfiles_twodee.endswith(os.path.sep):
                    outfiles_twodee += os.path.sep
            except:
                print('Folder outfiles already exists in ' + str(twodee_daily))
            twodee_input = os.path.join(twodee_daily, 'twodee.inp')
            topography_converter(topography_original)
            shutil.move(os.path.join(diagno_daily,'surface_data.txt'),os.path.join(infiles_twodee,'surface_data.txt'))
            shutil.copyfile(os.path.join(root,'roughness_twodee.grd'),os.path.join(infiles_twodee,'roughness.grd'))
            with open(os.path.join(infiles_twodee,'source.dat'), 'w', encoding="utf-8", errors="surrogateescape") as source_file:
                for i in range(0, n_sources):
                    if source_emission != 999:
                        gas_flux = source_emission
                    else:
                        if fluxes_input[i] == 99999999 and random_emission == 'on':
                            gas_flux = fluxes()[0]
                        else:
                            gas_flux = fluxes_input[i]
                    source_file.write('{0:7.3f}'.format(easting[i]) + ' ' + '{0:7.3f}'.format(northing[i]) + ' ' + '{0:7.3f}'.format(gas_flux) + ' '
                                      + '{0:7.2f}'.format(dx[i]) + ' ' + '{0:7.2f}'.format(dy[i]) + ' KG_SEC 0 ' + '{0:7.3f}'.format(dur[i]) + '\n')
            source_file.close()
            # read and memorize twodee.inp file
            twodee_input_records = []
            with open(twodee_original, 'r', encoding="utf-8", errors="surrogateescape") as twodee_or_input:
                for line in twodee_or_input:
                    twodee_input_records.append(line)
            twodee_or_input.close()
            with open(twodee_input, 'w', encoding="utf-8", errors="surrogateescape") as twodee_input:
                for record in twodee_input_records:
                    if 'YEAR' in record:
                        twodee_input.write('  YEAR   = ' + day[0:4] + '\n')
                    elif 'MONTH' in record:
                        twodee_input.write('  MONTH  = ' + day[4:6] + '\n')
                    elif 'DAY' in record:
                        twodee_input.write('  DAY    = ' + day[6:8] + '\n')
                    elif 'INPUT_DIRECTORY' in record:
                        twodee_input.write('   INPUT_DIRECTORY   = ' + infiles_twodee + ' \n')
                    elif 'OUTPUT_DIRECTORY' in record:
                        twodee_input.write('   OUTPUT_DIRECTORY   = ' + outfiles_twodee + ' \n')
                    else:
                        twodee_input.write(record)
            twodee_input.close()
        shutil.rmtree(path)
    return days

def run_diagno():
    n_elaborated_days = 0
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                diagno_folder = os.path.join(root, 'simulations', 'diagno', day)
                os.chdir(diagno_folder)
                try:
                    p = subprocess.Popen(['srun', '-n', '1', 'presfc','&'])
                except:
                    p = subprocess.Popen(['presfc'])
                p.wait()
                ps.append(p)
                try:
                    p = subprocess.Popen(['srun', '-n', '1', 'preupr','&'])
                except:
                    p = subprocess.Popen(['preupr'])
                p.wait()
                ps.append(p)
                try:
                    p = subprocess.Popen(['srun', '-n', '1', 'diagno','&'])
                except:
                    p = subprocess.Popen(['diagno'])
                ps.append(p)
            for p in ps:
                p.wait()
        except:
            print('Unable to process weather data with Diagno')
            sys.exit()
        print('DIAGNO successfully processed days ' + str(days[start:end]))
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break
    print('All weather data have been successfully processed with Diagno')
    os.chdir(root)

def run_disgas():
    n_elaborated_days = 0
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                disgas_folder = os.path.join(root, 'simulations', 'disgas', day, 'infiles')
                disgas_input_file = os.path.join(disgas_folder, 'disgas.inp')
                try:
                    p = subprocess.Popen(['srun', '-n', '1', 'disgas', disgas_input_file,'&'])
                except:
                    p = subprocess.Popen(['disgas', disgas_input_file])
                ps.append(p)
            for p in ps:
                p.wait()
        except:
            print('Unable to run DISGAS')
            sys.exit()
        print('DISGAS successfully processed days ' + str(days[start:end]))
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break

def run_twodee():
    n_elaborated_days = 0
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                diagno = os.path.join(root, 'simulations', 'diagno', day)
                twodee_folder = os.path.join(root, 'simulations', 'twodee', day)
                shutil.copyfile(os.path.join(diagno,'diagno.out'),os.path.join(twodee_folder,'infiles','diagno.out'))
                twodee_input_file = os.path.join(twodee_folder, 'twodee.inp')
                twodee_log_file = os.path.join(twodee_folder, 'twodee.log')
                try:
                    p = subprocess.Popen(['srun', '-n', '1', 'Twodee.2.2.x', twodee_input_file, twodee_log_file,'&'])
                except:
                    p = subprocess.Popen(['Twodee.2.2.x', twodee_input_file, twodee_log_file])
                ps.append(p)
            for p in ps:
                p.wait()
        except:
            print('Unable to run TWODEE')
            sys.exit()
        print('TWODEE successfully processed days ' + str(days[start:end]))
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break

root = os.getcwd()
disgas_original = os.path.join(root,'disgas.inp')
twodee_original = os.path.join(root,'twodee.inp')
topography_original = os.path.join(root,'topography.grd')

max_number_processes, random_sources, nsources, sources_interval, source_easting, source_northing, source_el, \
source_emission, random_emission, bottom_left_northing, bottom_left_easting, top_right_northing, top_right_easting, \
source_dx, source_dy, source_dur, twodee_on, disgas_on = read_arguments()

if disgas_on == 'off' and twodee_on == 'off':
    print('Both DISGAS and TWODEE are turned off')
    sys.exit()

days = pre_process()

run_diagno()

if disgas_on:
    run_disgas()

if twodee_on:
    run_twodee()

