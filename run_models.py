from random import sample
import os
import shutil
import subprocess
import argparse
import numpy as np
import sys
import utm


def read_arguments():
    parser = argparse.ArgumentParser(description="Input data", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-N",
        "--nproc",
        default=1,
        help="Maximum number of allowed simultaneous processes",
    )
    parser.add_argument('-RT', '--run_type', default='new', help='Specify if the simulation is a new one or a restart.'
                                                                 'Possible options are: new, restart')
    parser.add_argument('-CS', '--continuous_simulation', default='False', help='Specify if the simulation is '
                                                                                'continuous between the specified '
                                                                                'start and end dates. Possible options '
                                                                                'are True or False')
    parser.add_argument(
        "-RS",
        "--random_sources",
        default="off",
        help="on: randomly select NS locations from a probability map. off: fixed source locations from file "
             "sources_input.txt",
    )
    parser.add_argument(
        "-NS",
        "--nsources",
        default="random",
        help="Specify a number for a fixed number of sources. If random, then randomly select the number of sources "
             "from an interval",
    )
    parser.add_argument(
        "-SINT",
        "--sources_interval",
        default='',
        help="Type the minimum and maximum number of sources",
    )
    parser.add_argument(
        "-SLOC",
        "--source_location",
        default='',
        help="Coordinate type (UTM/GEO), latitude/northing, longitude/easting, elevation (above ground in m) of 1 "
             "fixed source",
    )
    parser.add_argument(
        "-SDX",
        "--source_dx",
        default=999999,
        help="Extension [m] along the X direction of 1 single source. Option valid for Twodee only",
    )
    parser.add_argument(
        "-SDY",
        "--source_dy",
        default=999999,
        help="Extension [m] along the Y direction of 1 single source. Option valid for Twodee only",
    )
    parser.add_argument(
        "-SDUR",
        "--source_dur",
        default=0,
        help="Emission duration [s] of 1 single source. Option valid for Twodee only",
    )
    parser.add_argument(
        "-D",
        "--domain",
        default='',
        help="Coordinates type (UTM/GEO), coordinates (latitude/northing, longitude/easting) of the bottom left corner "
             "and top right corner of the domain",
    )
    parser.add_argument(
        "-NX",
        "--nx",
        default=-1,
        help="Number of grid cells along the x-direction. If not provided, the grid spacing along the x-direction "
             "must be provided",
    )
    parser.add_argument(
        "-NY",
        "--ny",
        default=-1,
        help="Number of grid cells along the y-direction. If not provided, the grid spacing along the y-direction "
             "must be provided",
    )
    parser.add_argument(
        "-DX",
        "--dx",
        default=-1,
        help="Grid spacing (in m) along the x-direction. If not provided, the number of grid cells along the "
             "x-direction must be provided",
    )
    parser.add_argument(
        "-DY",
        "--dy",
        default=-1,
        help="Grid spacing (in m) along the y-direction. If not provided, the number of grid cells along the "
             "y-direction must be provided",
    )
    parser.add_argument(
        "-SEM",
        "--source_emission",
        default="999",
        help="Source emission rate [kg/s]. If specified, it is assigned to all the sources in the domain",
    )
    parser.add_argument(
        "-RER",
        "--random_emission",
        default="off",
        help="on: randomly assign emission rate for each source in the domain sampled from a flux.txt file. off: use "
             "specified emission rate",
    )
    parser.add_argument(
        "-DI", "--diagno", default="on", help="on or off, to run Diagno. Turn it off only if Diagno has already "
                                              "been run"
    )
    parser.add_argument(
        "-DM", "--dispersion_model", default="off", help="Twodee, Disgas or Automatic"
    )
    parser.add_argument(
        "-US", "--use_slurm", default="False", help="True or False, to use SLURM Workload Manager"
    )
    parser.add_argument(
        "-SP", "--slurm_partition", default="", help="Name of the cluster partition to run the Slurm jobs"
    )
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
    try:
        max_number_processes_in = int(nproc)
    except ValueError:
        print("Please provide a valid number for the maximum number of process")
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
    if run_type_in != 'new' and run_type_in != 'restart':
        print('ERROR. Please provide a valid entry for -RT --run_type')
        sys.exit()
    if continuous_simulation_in.lower() == "true":
        continuous_simulation_in = True
    elif continuous_simulation_in.lower() == "false":
        continuous_simulation_in = False
    else:
        print("ERROR. Wrong value for variable -CS --continuous_simulation")
        sys.exit()
    try:
        source_emission_in = float(source_emission_in)
    except ValueError:
        print("Please provide a valid number for the emission rate of the source")
        sys.exit()
    if len(domain) != 5:
        print("ERROR. Please provide valid entries for -D --domain")
        sys.exit()
    else:
        coordinates_type = domain[0]
        if coordinates_type == "GEO":
            bottom_left_1 = float(domain[1])
            bottom_left_2 = float(domain[2])
            top_right_1 = float(domain[3])
            top_right_2 = float(domain[4])
            if (-90 <= bottom_left_1 <= 90 and -180 <= bottom_left_2 <= 180) and (
                -90 <= top_right_1 <= 90 and -180 <= top_right_2 <= 180
            ):  # identify valid geographic coordinates
                try:
                    out_utm = utm.from_latlon(bottom_left_1, bottom_left_2)
                    bottom_left_easting_in = float(out_utm[0])
                    bottom_left_northing_in = float(out_utm[1])
                except ValueError:
                    print(
                        "ERROR. Please provide valid coordinates for the bottom left corner of the domain"
                    )
                    sys.exit()
                try:
                    out_utm = utm.from_latlon(top_right_1, top_right_2)
                    top_right_easting_in = float(out_utm[0])
                    top_right_northing_in = float(out_utm[1])
                except ValueError:
                    print(
                        "ERROR. Please provide valid coordinates for the top right corner of the domain"
                    )
                    sys.exit()
            else:
                print("ERROR. Please provide valid coordinates")
                sys.exit()
        elif coordinates_type == "UTM":
            bottom_left_northing_in = float(domain[1])
            bottom_left_easting_in = float(domain[2])
            top_right_northing_in = float(domain[3])
            top_right_easting_in = float(domain[4])
        else:
            print("ERROR. Please provide a valide type of coordinates (UTM or GEO)")
            sys.exit()
        if (
            bottom_left_northing_in == top_right_northing_in
            or bottom_left_easting_in == top_right_easting_in
        ):
            print("ERROR. Coordinates of the corners cannot coincide")
            sys.exit()
        if (
            bottom_left_northing_in > top_right_northing_in
            and bottom_left_easting_in > top_right_easting_in
        ):  # Check coordinates are in the proper order, otherwise swap
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
        print("Please provide a valid number for -NX --nx")
        sys.exit()
    try:
        ny_in = int(ny_in)
    except ValueError:
        print("Please provide a valid number for -NY --ny")
        sys.exit()
    try:
        dx_in = float(dx_in)
    except ValueError:
        print("Please provide a valid number for -DX --dx")
        sys.exit()
    try:
        dy_in = float(dy_in)
    except ValueError:
        print("Please provide a valid number for -DY --dy")
        sys.exit()
    if nx_in == -1 or ny_in == -1:
        if dx_in == -1 or dy_in == -1:
            print("ERROR. Either (NX, NY) or (DX, DY) must be specified")
            sys.exit()
        elif dx_in != -1 and dy_in != -1:
            if dx_in < 0 or dx_in > (top_right_easting_in - bottom_left_easting_in) or dy_in < 0 \
                    or dy_in > (top_right_northing_in - bottom_left_northing_in):
                print("ERROR. Please provide a valid number for (DX, DY)")
                sys.exit()
            else:
                nx_in = int((top_right_easting_in - bottom_left_easting_in) / dx_in)
                ny_in = int((top_right_northing_in - bottom_left_northing_in) / dy_in)
    elif nx_in != -1 and ny_in != -1:
        if nx_in < 0 or ny_in < 0:
            print("ERROR. Please provide a valid number for (NX, NY)")
            sys.exit()
        else:
            dx_in = (top_right_easting_in - bottom_left_easting_in) / float(nx_in)
            dy_in = (top_right_northing_in - bottom_left_northing_in) / float(ny_in)
    # Check provided nx, ny or dx, dy match the provided domain extent, otherwise correct
    if bottom_left_easting_in + dx_in * nx_in != top_right_easting_in:
        top_right_easting_in = bottom_left_easting_in + dx_in * nx_in
    if bottom_left_northing_in + dy_in * ny_in != top_right_northing_in:
        top_right_northing_in = bottom_left_northing_in + dy_in * ny_in

    if random_sources_in == "on":
        try:
            np.loadtxt("probability_map.grd", skiprows=5)
        except ValueError:
            print(
                "Please provide a valid probability_map.grd file when random_sources option is on"
            )
            sys.exit()
        if nsources_in == "random":
            if len(sources_interval_in) != 2:
                print(
                    "ERROR. Please specify the minimum and maximum number of sources with -SINT --sources_interval"
                )
                sys.exit()
        else:
            try:
                nsources_in = int(nsources_in)
            except ValueError:
                print("Please provide a valid integer for -NS --nsources")
                sys.exit()
        if random_emission_in == "off" and source_emission_in == 999:
            print(
                "ERROR. random_sources set to on requires either random_emission set to on or a "
                "specified source_emission"
            )
            sys.exit()
    else:
        if random_sources_in != "off":
            print("ERROR. Valid options for -RS --random_sources are on and off")
            sys.exit()
        else:
            try:
                sources_file = open("sources_input.txt", "r")
                sources_file.close()
            except FileNotFoundError:
                print(
                    "File sources_input.txt not found. Using one source from input data"
                )
                if len(source_location) != 4:
                    print(
                        "ERROR. Please provide valid entries for -SLOC --sources_location"
                    )
                    sys.exit()
                else:
                    coordinates_type = source_location[0]
                    if coordinates_type == "GEO":
                        if (
                            -90 <= float(source_location[1]) <= 90
                            and -180 <= float(source_location[2]) <= 180
                        ):  # identify geographic coordinates
                            try:
                                out_utm = utm.from_latlon(
                                    float(source_location[1]), float(source_location[2])
                                )
                                source_easting_in = float(out_utm[0])
                                source_northing_in = float(out_utm[1])
                            except ValueError:
                                print(
                                    "ERROR. Please provide valid coordinates of the source location"
                                )
                                sys.exit()
                    elif coordinates_type == "UTM":
                        source_easting_in = float(source_location[1])
                        source_northing_in = float(source_location[0])
                        if (
                            not bottom_left_easting_in
                            <= source_easting_in
                            <= top_right_easting_in
                            or not bottom_left_northing_in
                            <= source_northing_in
                            <= top_right_northing_in
                        ):
                            print("ERROR. Location not within the domain")
                            sys.exit()
                    else:
                        print(
                            "ERROR. Please provide a valide type of coordinates (UTM or GEO)"
                        )
                        sys.exit()
                    if float(source_location[2]) < 0:
                        print(
                            "ERROR. Please provide a valid value for the source elevation in m above ground (>= 0 m)"
                        )
                        sys.exit()
                    else:
                        source_el_in = float(source_location[2])
    if random_emission_in == "on":
        try:
            sources_file = open("flux.txt", "r")
            sources_file.close()
        except FileNotFoundError:
            print("ERROR. File flux.txt not found")
            sys.exit()
    elif random_emission_in != "off":
        print("ERROR. Valid options for -RER --random_sources are on and off")
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
    else:
        print("ERROR. Please provide a valid entry for the variable -DM --dispersion_model")
        sys.exit()
    if diagno.lower() == "on":
        diagno = True
    elif diagno.lower() == "off":
        diagno = False
    else:
        print("ERROR. Please provide a valid entry for the variable -DI --diagno")
        sys.exit()
    if source_dx_in != 999999:
        try:
            source_dx_in = float(source_dx_in)
        except ValueError:
            print("ERROR. Please provide a valid number for the variable -SDX --source_dx")
            sys.exit()
    if source_dy_in != 999999:
        try:
            source_dy_in = float(source_dy_in)
        except ValueError:
            print("ERROR. Please provide a valid number for the variable -SDY --source_dy")
            sys.exit()
    if source_dur_in != 0:
        try:
            source_dur_in = float(source_dur_in)
        except ValueError:
            print("ERROR. Please provide a valid number for the variable -SDUR --source_dur")
            sys.exit()
    if use_slurm.lower() == "true":
        use_slurm = True
        if slurm_partition == '':
            print('ERROR. Cluster partition must be declared if use_slurm is on')
    elif use_slurm.lower() == "false":
        use_slurm = False
    else:
        print("ERROR. Please provide a valid entry for the variable -US --use_slurm")
        sys.exit()
    slurm_partition = slurm_partition.lower()
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
    )


def prepare_days():
    try:
        raw_days = []  # store the days as originally formatted
        formatted_days = []  # store days in format YYYYMMDD as per folder name
        # read days_list file
        with open(
            os.path.join(root, "days_list.txt"), "r", encoding="utf-8", errors="surrogateescape"
        ) as days_list_file:
            for line in days_list_file:
                raw_days.append(line)
    except FileNotFoundError:
        print('ERROR. Unable to find days_list.txt file. Please restart activating DIAGNO')
        sys.exit()
    for k in range(0, len(raw_days)):
        temp = raw_days[k].split(" ")
        temp = temp[0].split("-")
        days.append(temp[0] + temp[1] + temp[2])
    return sorted(formatted_days)


def pre_process(run_type):
    def sample_random_sources(
        n_sources, input_file, dur_min, dur_max, source_size_min, source_size_max
    ):
        from random import choices

        with open(input_file) as probability_file:
            i = 1
            for line in probability_file:
                if i == 2:
                    nx = int(line.split(" ")[0])
                    ny = int(line.split(" ")[1])
                    i += 1
                    continue
                elif i == 3:
                    xmin = float(line.split(" ")[0])
                    xmax = float(line.split(" ")[1])
                    i += 1
                    continue
                elif i == 4:
                    ymin = float(line.split(" ")[0])
                    ymax = float(line.split(" ")[1])
                    i += 1
                    continue
                else:
                    i += 1
                    continue
        probabilities_input = np.loadtxt(input_file, skiprows=5)
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
            random_dx.append(choices(dxs, k=1)[0])
            random_dy.append(choices(dys, k=1)[0])
            random_dur.append(choices(durs, k=1)[0])

        return (
            random_eastings,
            random_northings,
            random_elevations,
            random_probabilities,
            random_fluxes,
            random_dx,
            random_dy,
            random_dur,
        )

    def fluxes():
        import numpy as np

        fluxes_in = []
        with open("flux.txt") as flux_file:
            for line in flux_file:
                fluxes_in.append(float(line))
        flux_file.close()
        x = np.sort(fluxes_in)
        list_x = list(x)
        sampled_flux = sample(list_x, 1)
        return sampled_flux

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
    if random_sources == "on":
        if nsources == "random":
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
        ) = sample_random_sources(
            n_random_sources,
            "probability_map.grd",
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
    dx_sources = random_dx
    dy_sources = random_dy
    dur = random_dur
    n_sources = 0
    try:
        with open(
            "sources_input.txt", "r", encoding="utf-8", errors="surrogateescape"
        ) as locations_file:
            for line in locations_file:
                try:
                    records = line.split("\t")
                    easting.append(float(records[0]))
                    northing.append(float(records[1]))
                    elevations.append(float(records[2]))
                    probabilities.append(float(records[3]))
                    fluxes_input.append(float(records[4]))
                    if twodee_on:
                        dx_sources.append(float(records[5]))
                        dy_sources.append(float(records[6]))
                        dur.append(float(records[7]))
                    n_sources += 1
                except ValueError:
                    continue
    except BaseException:
        if random_sources != "on":
            easting.append(source_easting)
            northing.append(source_northing)
            elevations.append(source_el)
            probabilities.append(1.0)
            fluxes_input.append(source_emission)
            dx_sources.append(source_dx)
            dy_sources.append(source_dy)
            dur.append(source_dur)
        n_sources = len(easting)

    disgas = os.path.join(root, "simulations", "disgas")
    twodee = os.path.join(root, "simulations", "twodee")
    diagno = os.path.join(root, "simulations", "diagno")
    if diagno_on:
        # Set DIAGNO folder
        try:
            os.mkdir(diagno)
        except FileExistsError:
            print("Folder " + diagno + " already exists")
    if disgas_on:
        # Set DISGAS folder
        try:
            os.mkdir(disgas)
        except FileExistsError:
            print("Folder " + disgas + " already exists")
    if twodee_on:
        # Set TWODEE folder
        try:
            os.mkdir(twodee)
        except FileExistsError:
            print("Folder " + twodee + " already exists")
    days = prepare_days()
    for i in range(0, len(days)):
        day = days[i]
        if continuous_simulation:
            if i == 0:
                run_type = 'new'
            else:
                run_type = 'restart'
        diagno_daily = os.path.join(diagno, str(day))
        path = os.path.join(root, "simulations", str(day))
        if diagno_on:
            files = os.listdir(path)
            try:
                os.mkdir(diagno_daily)
            except FileExistsError:
                print("Folder " + diagno_daily + " already exists")
            for f in files:
                path_f = os.path.join(path, f)
                try:
                    shutil.move(path_f, diagno_daily)
                except FileExistsError:
                    print("File " + f + " already present in " + diagno)
            shutil.copy(topography, os.path.join(diagno_daily, "topography.grd"))

        # Set DISGAS folder
        if disgas_on:
            disgas_daily = os.path.join(disgas, str(day))
            outfiles = os.path.join(disgas_daily, "outfiles")
            try:
                os.mkdir(disgas_daily)
            except FileExistsError:
                print("Folder " + disgas_daily + " already exists")
            try:
                os.mkdir(outfiles)
                if not outfiles.endswith(os.path.sep):
                    outfiles += os.path.sep
            except FileExistsError:
                print("Folder outfiles already exists in " + str(disgas_daily))
            disgas_input = os.path.join(disgas_daily, "disgas.inp")
            with open(
                os.path.join(disgas_daily, "source.dat"),
                "w",
                encoding="utf-8",
                errors="surrogateescape",
            ) as source_file:
                for j in range(0, n_sources):
                    if source_emission != 999:
                        gas_flux = source_emission
                    else:
                        if fluxes_input[j] == 99999999 and random_emission == "on":
                            gas_flux = fluxes()[0]
                        else:
                            gas_flux = fluxes_input[j]
                    source_file.write(
                        "{0:7.3f}".format(easting[j])
                        + " "
                        + "{0:7.3f}".format(northing[j])
                        + " "
                        + "{0:7.2f}".format(elevations[j])
                        + " "
                        + str(gas_flux)
                        + "\n"
                    )
            source_file.close()
            roughness_file_exist = True
            try:
                shutil.copyfile(
                    os.path.join(root, "roughness.grd"),
                    os.path.join(disgas_daily, "roughness.grd"),
                )
            except (FileExistsError, FileNotFoundError):
                roughness_file_exist = False
            try:
                shutil.move(
                    os.path.join(diagno_daily, "surface_data.txt"),
                    os.path.join(disgas_daily, "surface_data.txt"),
                )
            except (FileExistsError, FileNotFoundError):
                print('ERROR with surface_data.txt')
            # read and memorize disgas.inp file
            disgas_input_records = []
            with open(
                disgas_original, "r", encoding="utf-8", errors="surrogateescape"
            ) as disgas_or_input:
                for line in disgas_or_input:
                    disgas_input_records.append(line)
            with open(
                disgas_input, "w", encoding="utf-8", errors="surrogateescape"
            ) as disgas_input_file:
                for record in disgas_input_records:
                    if "ROUGHNESS_MODEL" in record:
                        roughness_command = record.split("=")[1]
                        roughness_command = roughness_command.split("(")[0]
                        if "MATRIX" in roughness_command:
                            if not roughness_file_exist:
                                print(
                                    "Warning! ROUGHNESS_MODEL set to MATRIX in disgas.inp and roughness file "
                                    "does not exist"
                                )
                                print("Setting ROUGHNESS_MODEL to UNIFORM")
                    if "YEAR" in record:
                        disgas_input_file.write("  YEAR   = " + day[0:4] + "\n")
                    elif "MONTH" in record:
                        disgas_input_file.write("  MONTH  = " + day[4:6] + "\n")
                    elif "DAY" in record:
                        disgas_input_file.write("  DAY    = " + day[6:8] + "\n")
                    elif "HOUR" in record:
                        try:
                            hour_start = float(record.split('=')[1])
                        except ValueError:
                            hour_start = record.split('=')[1].strip()
                            hour_start = float(hour_start.split(' ')[0])
                        if continuous_simulation:
                            if i > 0:
                                hour_start = 0
                        disgas_input_file.write("  HOUR   = " + "{0:2.0f}".format(hour_start) + "\n")
                    elif 'SIMULATION_INTERVAL_(SEC)' in record:
                        try:
                            simulation_interval = float(record.split('=')[1])
                        except ValueError:
                            simulation_interval = record.split('=')[1].strip()
                            simulation_interval = float(simulation_interval.split(' ')[0])
                        if simulation_interval + hour_start * 3600 > 86400:
                            simulation_interval = 86400
                        disgas_input_file.write("  SIMULATION_INTERVAL_(SEC) = " +
                                                "{0:7.0f}".format(simulation_interval) + "\n")
                    elif 'NX' in record:
                        disgas_input_file.write("  NX     = " + str(nx) + "\n")
                    elif 'NY' in record:
                        disgas_input_file.write("  NY     = " + str(ny) + "\n")
                    elif 'DX_(M)' in record:
                        disgas_input_file.write("  DX_(M) = " + str(dx) + "\n")
                    elif 'DY_(M)' in record:
                        disgas_input_file.write("  DY_(M) = " + str(dy) + "\n")
                    elif 'X_ORIGIN_(UTM_M)' in record:
                        disgas_input_file.write("  X_ORIGIN_(UTM_M) = " + str(bottom_left_easting) + "\n")
                    elif 'Y_ORIGIN_(UTM_M)' in record:
                        disgas_input_file.write("  Y_ORIGIN_(UTM_M) = " + str(bottom_left_northing) + "\n")
                    elif 'RESTART_RUN' in record:
                        if run_type == 'restart':
                            disgas_input_file.write("  RESTART_RUN = YES\n")
                        else:
                            disgas_input_file.write("  RESTART_RUN = NO\n")
                    elif 'RESET_TIME' in record:
                        if run_type == 'restart':
                            disgas_input_file.write("  RESET_TIME = YES\n")
                        else:
                            disgas_input_file.write("  RESET_TIME = NO\n")
                    elif "TOPOGRAPHY_FILE_PATH" in record:
                        disgas_input_file.write(
                            "   TOPOGRAPHY_FILE_PATH   = "
                            + os.path.join(diagno_daily, "topography.grd")
                            + " \n"
                        )
                    elif "ROUGHNESS_FILE_PATH" in record:
                        disgas_input_file.write(
                            "   ROUGHNESS_FILE_PATH   = "
                            + os.path.join(disgas_daily, "roughness.grd")
                            + " \n"
                        )
                    elif "RESTART_FILE_PATH" in record:
                        disgas_input_file.write(
                                "   RESTART_FILE_PATH   = "
                                + os.path.join(disgas_daily, "restart.dat")
                                + " \n"
                            )
                    elif "SOURCE_FILE_PATH" in record:
                        disgas_input_file.write(
                            "   SOURCE_FILE_PATH   = "
                            + os.path.join(disgas_daily, "source.dat")
                            + " \n"
                        )
                    elif "WIND_FILE_PATH" in record:
                        disgas_input_file.write(
                            "   WIND_FILE_PATH   = "
                            + os.path.join(disgas_daily, "winds.dat")
                            + " \n"
                        )
                    elif "DIAGNO_FILE_PATH" in record:
                        disgas_input_file.write(
                            "   DIAGNO_FILE_PATH   = "
                            + os.path.join(diagno_daily, "diagno.out")
                            + " \n"
                        )
                    elif "OUTPUT_DIRECTORY" in record:
                        disgas_input_file.write(
                            "   OUTPUT_DIRECTORY    = " + outfiles + " \n"
                        )
                    else:
                        disgas_input_file.write(record)
            shutil.copy(disgas_input, disgas_original)
        if twodee_on:
            twodee_daily = os.path.join(twodee, str(day))
            outfiles_twodee = os.path.join(twodee_daily, "outfiles")
            try:
                os.mkdir(twodee_daily)
            except FileExistsError:
                print("Folder " + twodee_daily + " already exists")
            try:
                os.mkdir(outfiles_twodee)
                if not outfiles_twodee.endswith(os.path.sep):
                    outfiles_twodee += os.path.sep
            except FileExistsError:
                print("Folder outfiles already exists in " + str(twodee_daily))
            twodee_input = os.path.join(twodee_daily, "twodee.inp")
            try:
                shutil.copyfile(
                    os.path.join(root, "roughness.grd"),
                    os.path.join(twodee_daily, "roughness.grd"),
                )
            except FileNotFoundError:
                print("Unable to find a valid roughness file for TWODEE")
                sys.exit()
            try:
                shutil.copyfile(
                    os.path.join(diagno_daily, "surface_data.txt"),
                    os.path.join(twodee_daily, "surface_data.txt"),
                )
            except (FileExistsError, FileNotFoundError):
                print('ERROR with surface_data.txt')
            with open(
                os.path.join(twodee_daily, "source.dat"),
                "w",
                encoding="utf-8",
                errors="surrogateescape",
            ) as source_file:
                for j in range(0, n_sources):
                    if source_emission != 999:
                        gas_flux = source_emission
                    else:
                        if fluxes_input[j] == 99999999 and random_emission == "on":
                            gas_flux = fluxes()[0]
                        else:
                            gas_flux = fluxes_input[j]
                    source_file.write(
                        "{0:7.3f}".format(easting[j])
                        + " "
                        + "{0:7.3f}".format(northing[j])
                        + " "
                        + "{0:7.3f}".format(gas_flux)
                        + " "
                        + "{0:7.2f}".format(dx_sources[j])
                        + " "
                        + "{0:7.2f}".format(dy_sources[j])
                        + " KG_SEC 0 "
                        + "{0:7.3f}".format(dur[j])
                        + "\n"
                    )
            source_file.close()
            # read and memorize twodee.inp file
            twodee_input_records = []
            with open(
                twodee_original, "r", encoding="utf-8", errors="surrogateescape"
            ) as twodee_or_input:
                for line in twodee_or_input:
                    twodee_input_records.append(line)
            with open(
                twodee_input, "w", encoding="utf-8", errors="surrogateescape"
            ) as twodee_input_file:
                for record in twodee_input_records:
                    if "YEAR" in record:
                        twodee_input_file.write("  YEAR   = " + day[0:4] + "\n")
                    elif "MONTH" in record:
                        twodee_input_file.write("  MONTH  = " + day[4:6] + "\n")
                    elif "DAY" in record:
                        twodee_input_file.write("  DAY    = " + day[6:8] + "\n")
                    elif "HOUR" in record and "NC_VAR_HOUR_NAME" not in record:
                        try:
                            hour_start = float(record.split('=')[1])
                        except ValueError:
                            hour_start = record.split('=')[1].strip()
                            hour_start = float(hour_start.split(' ')[0])
                        if continuous_simulation:
                            if i > 0:
                                hour_start = 0
                        twodee_input_file.write("  HOUR   = " + "{0:2.0f}".format(hour_start) + "\n")
                    elif 'SIMULATION_INTERVAL_(SEC)' in record:
                        # for the moment this is managed exactly like DISGAS, but would require the RESET_TIME option to
                        # be implemented in TWODEE as well
                        try:
                            simulation_interval = float(record.split('=')[1])
                        except ValueError:
                            simulation_interval = record.split('=')[1].strip()
                            simulation_interval = float(simulation_interval.split(' ')[0])
                        if simulation_interval + hour_start * 3600 > 86400:
                            simulation_interval -= hour_start * 3600
                        twodee_input_file.write("  SIMULATION_INTERVAL_(SEC) = " +
                                                "{0:7.0f}".format(simulation_interval) + "\n")
                    elif 'NX' in record:
                        twodee_input_file.write("  NX     = " + str(nx) + "\n")
                    elif 'NY' in record:
                        twodee_input_file.write("  NY     = " + str(ny) + "\n")
                    elif 'DX_(M)' in record:
                        twodee_input_file.write("  DX_(M) = " + str(dx) + "\n")
                    elif 'DY_(M)' in record:
                        twodee_input_file.write("  DY_(M) = " + str(dy) + "\n")
                    elif 'X_ORIGIN_(UTM_M)' in record:
                        twodee_input_file.write("  X_ORIGIN_(UTM_M) = " + str(bottom_left_easting) + "\n")
                    elif 'Y_ORIGIN_(UTM_M)' in record:
                        twodee_input_file.write("  Y_ORIGIN_(UTM_M) = " + str(bottom_left_northing) + "\n")
                    elif 'RESTART_RUN' in record:
                        if run_type == 'restart':
                            twodee_input_file.write("  RESTART_RUN = YES\n")
                        else:
                            twodee_input_file.write("  RESTART_RUN = NO\n")
                    elif "OUTPUT_DIRECTORY" in record:
                        twodee_input_file.write(
                            "   OUTPUT_DIRECTORY   = " + outfiles_twodee + " \n"
                        )
                    elif (
                        "TOPOGRAPHY_FILE" in record
                        and "TOPOGRAPHY_FILE_FORMAT" not in record
                    ):
                        twodee_input_file.write(
                            "   TOPOGRAPHY_FILE   = "
                            + os.path.join(diagno_daily, "topography.grd")
                            + " \n"
                        )
                    elif (
                        "ROUGHNESS_FILE" in record
                        and "ROUGHNESS_FILE_FORMAT" not in record
                    ):
                        twodee_input_file.write(
                            "   ROUGHNESS_FILE   = "
                            + os.path.join(twodee_daily, "roughness.grd")
                            + " \n"
                        )
                    elif "SOURCE_FILE" in record:
                        twodee_input_file.write(
                            "   SOURCE_FILE   = "
                            + os.path.join(twodee_daily, "source.dat")
                            + " \n"
                        )
                    elif "SURF_DATA_FILE" in record:
                        twodee_input_file.write(
                            "   SURF_DATA_FILE   = "
                            + os.path.join(diagno_daily, "surface_data.txt")
                            + " \n"
                        )
                    elif "DIAGNO_FILE" in record:
                        twodee_input_file.write(
                            "   DIAGNO_FILE   = "
                            + os.path.join(diagno_daily, "diagno.out")
                            + " \n"
                        )
                    elif "RESTART_FILE" in record:
                        twodee_input_file.write(
                                "   RESTART_FILE   = "
                                + os.path.join(twodee_daily, "restart.dat")
                                + " \n"
                            )
                    else:
                        twodee_input_file.write(record)
            shutil.copy(twodee_input, twodee_original)
        if diagno_on:
            shutil.rmtree(path)
    return days


def run_diagno(max_number_processes):
    n_elaborated_days = 0
    n_node = 0
    node = ''
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                if len(nodes_list) > 0:
                    try:
                        node = nodes_list[n_node]
                    except IndexError:
                        node = ''
                diagno_folder = os.path.join(root, "simulations", "diagno", day)
                # read and memorize diagno.inp file
                diagno_input_records = []
                diagno_input = os.path.join(diagno_folder, 'diagno.inp')
                with open(
                        diagno_input, "r", encoding="utf-8", errors="surrogateescape"
                ) as diagno_or_input:
                    for line in diagno_or_input:
                        diagno_input_records.append(line)
                with open(
                        diagno_input, "w", encoding="utf-8", errors="surrogateescape"
                ) as diagno_input_file:
                    for record in diagno_input_records:
                        if 'NX' in record:
                            if disgas_on:
                                diagno_input_file.write(str(nx + 2) + '          NX\n')
                            elif twodee_on:
                                diagno_input_file.write(str(nx) + '          NX\n')
                        elif 'NY' in record:
                            if disgas_on:
                                diagno_input_file.write(str(ny + 2) + '          NY\n')
                            elif twodee_on:
                                diagno_input_file.write(str(ny) + '          NY\n')
                        elif 'DXK' in record:
                            diagno_input_file.write("{0:7.3f}".format(dx / 1000) + '          DXK (km)\n')
                        elif 'DYK' in record:
                            diagno_input_file.write("{0:7.3f}".format(dy / 1000) + '          DYK (km)\n')
                        elif 'UTMXOR' in record:
                            if disgas_on:
                                diagno_input_file.write("{0:7.3f}".format((bottom_left_easting - dx) / 1000) + '      '
                                                        'UTMXOR (km)\n')
                            elif twodee_on:
                                diagno_input_file.write("{0:7.3f}".format(bottom_left_easting / 1000) + '      '
                                                        'UTMXOR (km)\n')
                        elif 'UTMYOR' in record:
                            if disgas_on:
                                diagno_input_file.write("{0:7.3f}".format((bottom_left_northing - dy) / 1000) + '      '
                                                        'UTMYOR  (km)\n')
                            elif twodee_on:
                                diagno_input_file.write("{0:7.3f}".format(bottom_left_northing / 1000) + '      '
                                                        'UTMYOR  (km)\n')
                        else:
                            diagno_input_file.write(record)
                os.chdir(diagno_folder)
                if slurm:
                    try:
                        p = subprocess.Popen(["srun", "-n", "1", '--partition=' + partition, "presfc", "&"])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' presfc')
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(["presfc"])
                    except FileNotFoundError:
                        print('Unable to run presfc for DIAGNO')
                        sys.exit()
                p.wait()
                ps.append(p)
                if shutil.which('preupr') is not None and \
                        (not os.path.exists(os.path.join(root, 'weather_stations_list.txt'))):
                    if slurm:
                        try:
                            p = subprocess.Popen(["srun", "-n", "1", '--partition=' + partition, "preupr", "&"])
                        except FileNotFoundError:
                            print('Unable to run srun -n 1 --partition=' + partition + ' preupr')
                            sys.exit()
                    else:
                        try:
                            p = subprocess.Popen(["preupr"])
                        except FileNotFoundError:
                            print('Unable to run preupr for DIAGNO')
                            sys.exit()
                    p.wait()
                    ps.append(p)
                if slurm:
                    try:
                        p = subprocess.Popen(["srun", "-n", "1", '--partition=' + partition, '--nodelist=' + node,
                                              "diagno", "&"])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' --nodelist=' + node + ' diagno')
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(["diagno"])
                    except FileNotFoundError:
                        print('Unable to run DIAGNO')
                        sys.exit()
                ps.append(p)
                if len(nodes_list) > 0:
                    n_node += 1
                    if n_node >= len(nodes_list):
                        n_node = 0
        except BaseException:
            print("Unable to process weather data with Diagno")
            sys.exit()
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break
        for p in ps:
            p.wait()
    for p in ps:
        p.wait()
    print("All weather data have been successfully processed with Diagno")
    os.chdir(root)


def read_diagno_outputs():
    from scipy.io import FortranFile

    nx_source = [i_source for i_source in range(43, 44)]
    ny_source = [j_source for j_source in range(56, 57)]
    f = FortranFile('diagno.out', 'r')
    vx = []
    vy = []
    wind = []
    while True:
        try:
            for nline in range(0, 3):
                f.read_reals(dtype='float32')
            nx, ny, nz = f.read_reals(dtype='int')
            for nline in range(0, 5):
                f.read_reals(dtype='float32')
            npoints = 0
            vx_average = 0
            for j in range(0, ny):
                vx_vector = f.read_reals(dtype='float32')
                if j in ny_source:
                    for i_v in range(0, len(vx_vector)):
                        if i_v in nx_source:
                            vx_average += vx_vector[i]
                            npoints += 1
            vx.append(vx_average / npoints)
            for k in range(1, nz):
                for j in range(0, ny):
                    f.read_reals(dtype='float32')
            npoints = 0
            vy_average = 0
            for j in range(0, ny):
                vy_vector = f.read_reals(dtype='float32')
                if j in ny_source:
                    for i_v in range(0, len(vy_vector)):
                        if i_v in nx_source:
                            vy_average += vy_vector[i]
                            npoints += 1
            vy.append(vy_average / npoints)
            wind.append((vx[-1] ** 2 + vy[-1] ** 2) ** 0.5)
            for k in range(1, nz):
                for j in range(0, ny):
                    f.read_reals(dtype='float32')
            # Skip vz
            for k in range(0, nz):
                for j in range(0, ny):
                    f.read_reals(dtype='float32')
        except BaseException:
            break


def run_disgas(max_number_processes):
    import datetime
    disgas = os.path.join(root, "simulations", "disgas")
    n_elaborated_days = 0
    n_node = 0
    node = ''
    if continuous_simulation:
        max_number_processes = 1
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                if len(nodes_list) > 0:
                    try:
                        node = nodes_list[n_node]
                    except IndexError:
                        node = ''
                disgas_folder = os.path.join(disgas, day)
                disgas_input_file = os.path.join(disgas_folder, "disgas.inp")
                disgas_log_file = os.path.join(
                    disgas_folder, "disgas_log_" + day + ".txt"
                )
                if continuous_simulation and day != days[0]:
                    current_day = datetime.datetime.strptime(day, '%Y%m%d')
                    previous_day = current_day - datetime.timedelta(days=1)
                    previous_day = previous_day.strftime('%Y%m%d')
                    disgas_previous_day = os.path.join(disgas, str(previous_day))
                    try:
                        shutil.copyfile(os.path.join(disgas_previous_day, "restart.dat"),
                                        os.path.join(disgas_folder, "restart.dat"))
                    except FileNotFoundError:
                        print('ERROR. Restart run requested but file ' +
                              os.path.join(disgas_previous_day, "restart.dat") + ' not found')
                        sys.exit()
                if not diagno_on:
                    disgas_output_folder = os.path.join(disgas_folder, 'outfiles')
                    try:
                        os.mkdir(disgas_output_folder)
                    except FileExistsError:
                        shutil.rmtree(disgas_output_folder)
                        os.mkdir(disgas_output_folder)
                if slurm:
                    try:
                        p = subprocess.Popen(["srun", "-n", "1", '--partition=' + partition, '--nodelist=' + node,
                                              "disgas", disgas_input_file, disgas_log_file])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' --nodelist=' + node + ' disgas ' +
                              disgas_input_file + ' ' + disgas_log_file)
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(["disgas", disgas_input_file, disgas_log_file])
                    except FileNotFoundError:
                        print('Unable to run DISGAS')
                        sys.exit()
                ps.append(p)
                if len(nodes_list) > 0:
                    n_node += 1
                    if n_node >= len(nodes_list):
                        n_node = 0
        except BaseException:
            print("Unable to run DISGAS")
            sys.exit()
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break
        for p in ps:
            p.wait()
    for p in ps:
        p.wait()


def run_twodee(max_number_processes):
    import datetime
    twodee = os.path.join(root, "simulations", "twodee")
    n_elaborated_days = 0
    n_node = 0
    node = ''
    if continuous_simulation:
        max_number_processes = 1
    while n_elaborated_days <= len(days):
        ps = []
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > len(days):
            end = len(days)
        try:
            for day in days[start:end]:
                if len(nodes_list) > 0:
                    try:
                        node = nodes_list[n_node]
                    except IndexError:
                        node = ''
                diagno = os.path.join(root, "simulations", "diagno", day)
                twodee_folder = os.path.join(twodee, day)
                shutil.copyfile(
                    os.path.join(diagno, "diagno.out"),
                    os.path.join(twodee_folder, "diagno.out"),
                )
                twodee_input_file = os.path.join(twodee_folder, "twodee.inp")
                twodee_log_file = os.path.join(
                    twodee_folder, "twodee_log_" + day + ".txt"
                )
                if continuous_simulation and day != days[0]:
                    current_day = datetime.datetime.strptime(day, '%Y%m%d')
                    previous_day = current_day - datetime.timedelta(days=1)
                    previous_day = previous_day.strftime('%Y%m%d')
                    twodee_previous_day = os.path.join(twodee, str(previous_day))
                    try:
                        shutil.copyfile(os.path.join(twodee_previous_day, "restart.dat"),
                                        os.path.join(twodee_folder, "restart.dat"))
                    except FileNotFoundError:
                        print('ERROR. Restart run requested but file ' +
                              os.path.join(twodee_previous_day, "restart.dat") + ' not found')
                        sys.exit()
                if not diagno_on:
                    twodee_output_folder = os.path.join(twodee_folder, 'outfiles')
                    try:
                        os.mkdir(twodee_output_folder)
                    except FileExistsError:
                        shutil.rmtree(twodee_output_folder)
                        os.mkdir(twodee_output_folder)
                if slurm:
                    try:
                        p = subprocess.Popen(["srun", "-n", "1", '--partition=' + partition,  '--nodelist=' + node,
                                              "twodee", twodee_input_file, twodee_log_file, ])
                    except FileNotFoundError:
                        print('Unable to run srun -n 1 --partition=' + partition + ' --nodelist=' + node + ' ' +
                              twodee_input_file + ' ' + twodee_log_file)
                        sys.exit()
                else:
                    try:
                        p = subprocess.Popen(["twodee", twodee_input_file, twodee_log_file])
                    except FileNotFoundError:
                        print('Unable to run TWODEE')
                        sys.exit()
                ps.append(p)
                if len(nodes_list) > 0:
                    n_node += 1
                    if n_node >= len(nodes_list):
                        n_node = 0
        except BaseException:
            print("Unable to run TWODEE")
            sys.exit()
        n_elaborated_days = end
        if n_elaborated_days == len(days):
            break
        for p in ps:
            p.wait()
    for p in ps:
        p.wait()


root = os.getcwd()
disgas_original = os.path.join(root, "disgas.inp")
twodee_original = os.path.join(root, "twodee.inp")
topography = os.path.join(root, "topography.grd")

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
) = read_arguments()

if disgas_on == "off" and twodee_on == "off":
    print("Both DISGAS and TWODEE are turned off")
    sys.exit()

if slurm:
    list_available_nodes = []
    result = subprocess.run(['sinfo', '-o=%N'], stdout=subprocess.PIPE)
    sinfo_node_output = result.stdout.decode("utf-8")
    node_key = sinfo_node_output.split('=')[2]
    node_key = node_key.split('[')[0]
    result = subprocess.run(['sinfo'], stdout=subprocess.PIPE)
    sinfo_output = result.stdout.decode("utf-8").split('\n')
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
        for i in range(start, end + 1):
            list_available_nodes.append(node_key + '{:0>2}'.format(i))
    result = subprocess.run(['sinfo', '-o=%c'], stdout=subprocess.PIPE)
    sinfo_core_output = result.stdout.decode("utf-8")
    ncpus_per_node = int(sinfo_core_output.split('=')[-1])
    n_nodes = - (-max_number_processes // ncpus_per_node)
    nodes_list = list_available_nodes[0:n_nodes]
else:
    nodes_list = []

days = pre_process(run_type)

if diagno_on:
    run_diagno(max_number_processes)

if disgas_on and twodee_on: #We are in automatic scenario detection mode, hence we need to read the meteorological data
# at the source locations
    read_diagno_outputs()

if disgas_on:
    run_disgas(max_number_processes)

if twodee_on:
    run_twodee(max_number_processes)
