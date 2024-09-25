# VIGIL - automatic probabilistic VolcanIc Gas dIspersion modeLling

> A collection of Python scripts for automating workflow of gas dispersion
> modelling (weather data download and processing; running DISGAS and TWODEE
> models; post-process and plot results) for use in gas pollution modelling.


## Package content and description

The package contains the following files:

- weather.py

```bash
Python script that prepares the weather data needed by Diagno to compute the wind field for DISGAS. 
The weather data can be either retrieved from ECMWF ERA5 database or from time series of data from weather stations that are stored in files in the working folder.
	+ ERA5 Reanalysis data
	In this mode, the script is designed to randomly sample N days from a time interval defined by the user. The number of days N is also specified by the user in input. It is possible to retrieve date from one single day by setting the end date equal to the start date and the number of samples to 1
	+ Weather station data
	In this mode, the script is design to extract weather data in the time interval specified by the user from selected weather data file. The script reads the number of files, file location and name from weather_stations.txt; the file data should be stored in the folder weather_stations.
The following flags control the execution of weather.py:
usage: weather.py [-h] [-M MODE] [-RT RUN_TYPE] [-CS CONTINUOUS_SIMULATION] [-S START_DATE] [-E END_DATE]
                  [-SY SAMPLED_YEARS] [-SM SAMPLED_MONTHS] [-SD SAMPLED_DAYS] [-V VOLC] [-LAT LAT] [-LON LON] [-EL ELEV]
                  [-NS SAMPLES] [-ERA5 ERA5] [-WST STATION] [-N NPROC] [-TD TWODEE] [-DG DISGAS]
Input data optional arguments:
  -h, --help            show this help message and exit
  -M MODE, --mode MODE  Possible options: reanalysis, forecast. If reanalysis, either ERA5 or WST options should be on. If
                        forecast, GFS data will be downloaded and processed (default: reanalysis)
  -RT RUN_TYPE, --run_type RUN_TYPE
                        Specify if the simulation is a new one or a restart. Possible options are: new, restart (default:
                        new)
  -CS CONTINUOUS_SIMULATION, --continuous_simulation CONTINUOUS_SIMULATION
                        Specify if the simulation is continuous between the specified start and end dates. Possible options
                        are True or False (default: False)
  -S START_DATE, --start_date START_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY (default: 999)
  -E END_DATE, --end_date END_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY (default: 999)
  -SY SAMPLED_YEARS, --sampled_years SAMPLED_YEARS
                        Specify years to sample from the time interval (default: )
  -SM SAMPLED_MONTHS, --sampled_months SAMPLED_MONTHS
                        Specify months to sample from the time interval (default: )
  -SD SAMPLED_DAYS, --sampled_days SAMPLED_DAYS
                        Specify days to sample from the time interval (default: )
  -V VOLC, --volc VOLC  This is the volcano ID based on the Smithsonian Institute IDs (default: 999)
  -LAT LAT, --lat LAT   Volcano latitude (default: 999)
  -LON LON, --lon LON   Volcano longitude (default: 999)
  -EL ELEV, --elev ELEV
                        Volcano elevation (default: 999)
  -NS SAMPLES, --samples SAMPLES
                        Number of days to sample (default: 1)
  -ERA5 ERA5, --era5 ERA5
                        True: Use ERA5 reanalysis. False: Do not use ERA5 reanalysis (default: False)
  -WST STATION, --station STATION
                        True: Use weather station data. False: Do not use weather station data (default: False)
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes (default: 1)
  -TD TWODEE, --twodee TWODEE
                        on or off, to prepare additional weather data files for Twodee. (default: off)
  -DG DISGAS, --disgas DISGAS
                        on or off, to run Disgas (default: off)

```

- run_models.py

```bash
Python script to run Diagno and DISGAS for the days sampled with weather.py. 
The following flags control the execution of run_models.py:
usage: run_models.py [-h] [-N NPROC] [-RT RUN_TYPE] [-CS CONTINUOUS_SIMULATION] [-RS RANDOM_SOURCES] [-NS NSOURCES]
                     [-SINT SOURCES_INTERVAL] [-SLOC SOURCE_LOCATION] [-SDX SOURCE_DX] [-SDY SOURCE_DY]
                     [-SDUR SOURCE_DUR] [-D DOMAIN] [-NX NX] [-NY NY] [-DX DX] [-DY DY] [-SEM SOURCE_EMISSION]
                     [-RER RANDOM_EMISSION] [-PDEM PROB_DISTR_EMISSION] [-PDPAR PROB_DISTR_PARAMS] [-RD RUN_DURATION]
                     [-OI OUTPUT_INTERVAL] [-OH OUTPUT_HEIGHTS] [-DI DIAGNO] [-DM DISPERSION_MODEL] [-US USE_SLURM]
                     [-SP SLURM_PARTITION] [-TS TRACKING_SPECIE]

Input data

optional arguments:
  -h, --help            show this help message and exit
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes (default: 1)
  -RT RUN_TYPE, --run_type RUN_TYPE
                        Specify if the simulation is a new one or a restart.Possible options are: new, restart
                        (default: new)
  -CS CONTINUOUS_SIMULATION, --continuous_simulation CONTINUOUS_SIMULATION
                        Specify if the simulation is continuous between the specified start and end dates. Possible
                        options are True or False (default: False)
  -RS RANDOM_SOURCES, --random_sources RANDOM_SOURCES
                        on: randomly select NS locations from a probability map. off: fixed source locations from file
                        sources_input.txt (default: off)
  -NS NSOURCES, --nsources NSOURCES
                        Specify a number for a fixed number of sources. If random, then randomly select the number of
                        sources from an interval (default: random)
  -SINT SOURCES_INTERVAL, --sources_interval SOURCES_INTERVAL
                        Type the minimum and maximum number of sources (default: )
  -SLOC SOURCE_LOCATION, --source_location SOURCE_LOCATION
                        Coordinate type (UTM/GEO), latitude/northing, longitude/easting, elevation (above ground in m)
                        of 1 fixed source (default: )
  -SDX SOURCE_DX, --source_dx SOURCE_DX
                        Extension [m] along the X direction of 1 single source. Option valid for Twodee only (default:
                        999999)
  -SDY SOURCE_DY, --source_dy SOURCE_DY
                        Extension [m] along the Y direction of 1 single source. Option valid for Twodee only (default:
                        999999)
  -SDUR SOURCE_DUR, --source_dur SOURCE_DUR
                        Emission duration [s] of 1 single source. Option valid for Twodee only (default: 0)
  -D DOMAIN, --domain DOMAIN
                        Coordinates type (UTM/GEO), coordinates (latitude/northing, longitude/easting) of the bottom
                        left corner and top right corner of the domain (default: )
  -NX NX, --nx NX       Number of grid cells along the x-direction. If not provided, the grid spacing along the
                        x-direction must be provided (default: -1)
  -NY NY, --ny NY       Number of grid cells along the y-direction. If not provided, the grid spacing along the
                        y-direction must be provided (default: -1)
  -DX DX, --dx DX       Grid spacing (in m) along the x-direction. If not provided, the number of grid cells along the
                        x-direction must be provided (default: -1)
  -DY DY, --dy DY       Grid spacing (in m) along the y-direction. If not provided, the number of grid cells along the
                        y-direction must be provided (default: -1)
  -SEM SOURCE_EMISSION, --source_emission SOURCE_EMISSION
                        Source emission rate [kg/s]. If specified, it is assigned to all the sources in the domain
                        (default: )
  -RER RANDOM_EMISSION, --random_emission RANDOM_EMISSION
                        on: randomly assign emission rate for each source in the domain. off: use specified emission
                        rate (default: off)
  -PDEM PROB_DISTR_EMISSION, --prob_distr_emission PROB_DISTR_EMISSION
                        Probability distribution function to randomly sample the emission rate. Options: uniform,
                        normal, ecdf (default: )
  -PDPAR PROB_DISTR_PARAMS, --prob_distr_params PROB_DISTR_PARAMS
                        If -PDEM=uniform: minimum, maximum. If -PDEM=gaussian: median, standard deviation (default: )
  -RD RUN_DURATION, --run_duration RUN_DURATION
                        Run duration (hours). Currently fractions of hours or duration > 24 hours are not allowed
                        (default: 24)
  -OI OUTPUT_INTERVAL, --output_interval OUTPUT_INTERVAL
                        Output interval (hours). Currently fractions of hours are not allowed (default: 1)
  -OH OUTPUT_HEIGHTS, --output_heights OUTPUT_HEIGHTS
                        List of output heights (comma separated) in m above the ground (default: )
  -DI DIAGNO, --diagno DIAGNO
                        on or off, to run Diagno. Turn it off only if Diagno has already been run (default: on)
  -DM DISPERSION_MODEL, --dispersion_model DISPERSION_MODEL
                        Twodee, Disgas, Automatic, None (default: off)
  -US USE_SLURM, --use_slurm USE_SLURM
                        True or False, to use SLURM Workload Manager (default: False)
  -SP SLURM_PARTITION, --slurm_partition SLURM_PARTITION
                        Name of the cluster partition to run the Slurm jobs (default: )
  -TS TRACKING_SPECIE, --tracking_specie TRACKING_SPECIE
                        The original emitted specie that is tracked in the simulation (default: )

```

- post_process.py

```bash
Python script to process the outputs of run_models.py to produce plots, statistical outputs and tracking points. 
The following flags control the execution of post_process.py:
usage: post_process.py [-h] [-P PLOT] [-ECDF CALCULATE_ECDF] [-PER PERSISTENCE] [-EX EX_PROB] [-T TIME_STEPS] [-L LEVELS]
                       [-D DAYS_PLOT] [-C CONVERT] [-S SPECIES] [-TS TRACKING_SPECIE] [-N NPROC] [-U UNITS]
                       [-PL PLOT_LIMITS] [-PI PLOT_ISOLINES] [-TA TIME_AV] [-OF OUTPUT_FORMAT] [-PT PLOT_TOPOGRAPHY]
                       [-TI TOPOGRAPHY_ISOLINES] [-PR PLOT_RESOLUTION] [-TP TRACKING_POINTS]
Input data optional arguments:
  -h, --help            show this help message and exit
  -P PLOT, --plot PLOT  Produce plots of the solutions and probabilistic output (if activated). True/False (default: False)
  -ECDF CALCULATE_ECDF, --calculate_ecdf CALCULATE_ECDF
                        Calculate the Empirical Cumulative Density Function of the solution and extrapolate solutions at
                        user-defined exceedance probabilities. True/False (default: False)
  -PER PERSISTENCE, --persistence PERSISTENCE
                        Calculate the persistence of the gas specie, i.e. the probability to be exposed to a gas species
                        above specified concentration thresholds for times longer than the specified exposure times for
                        those thresholds. Concentration thresholds and exposure times should be provided in
                        gas_properties.csv. True/False (default: False)
  -EX EX_PROB, --ex_prob EX_PROB
                        List of exceedance probabilities to be used for graphical output (default: )
  -T TIME_STEPS, --time_steps TIME_STEPS
                        List of time steps to plot (integer >= 0). Type all to plot all the time steps (default: )
  -L LEVELS, --levels LEVELS
                        List of vertical levels (integer >= 1) to plot. Type all to plot all the levels (default: )
  -D DAYS_PLOT, --days_plot DAYS_PLOT
                        List of days to plot (YYYYMMDD). Type all to plot all the days (default: )
  -C CONVERT, --convert CONVERT
                        If True, convert output concentration into other species listed with the command -S (--species)
                        (default: False)
  -S SPECIES, --species SPECIES
                        List of gas species (e.g. CO2) (default: )
  -TS TRACKING_SPECIE, --tracking_specie TRACKING_SPECIE
                        The original emitted specie that is tracked in the simulation (default: None)
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes (default: 1)
  -U UNITS, --units UNITS
                        Gas concentration units. Possible options are: ppm, kg/m3 (default: None)
  -PL PLOT_LIMITS, --plot_limits PLOT_LIMITS
                        Minimum and maximum value of concentration to display. If unspecified, they are obtained from all
                        the outputs (default: )
  -PI PLOT_ISOLINES, --plot_isolines PLOT_ISOLINES
                        List of gas concentrations values to be used to draw isolines. Optional (default: )
  -TA TIME_AV, --time_av TIME_AV
                        Generate time-averaged outputs. Specify the time-averaging interval (in hours), or 0 for averaging
                        over the whole duration (default: None)
  -OF OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        Select format of the processed output files. Valid options are: GRD (default: GRD)
  -PT PLOT_TOPOGRAPHY, --plot_topography PLOT_TOPOGRAPHY
                        Plot topography layer (True or False). Warning, it can be time-consuming! (default: False)
  -TI TOPOGRAPHY_ISOLINES, --topography_isolines TOPOGRAPHY_ISOLINES
                        Topography height contour lines spatial resolution (in m a.s.l.). Used only if -PT True (default:
                        100)
  -PR PLOT_RESOLUTION, --plot_resolution PLOT_RESOLUTION
                        Specify plot resolution in dpi (default: 600)
  -TP TRACKING_POINTS, --tracking_points TRACKING_POINTS
                        Extrapolate gas concentration at locations specified in the file tracking_points.txt (default:
                        False)


```

### Dependencies and installation instructions

#### System dependencies

The following software are required (and should be available on the system
$PATH):

Gas modelling software:

- DIAGNO v1.2.3: the diagnostic wind model (Douglas et al., 1990)
  Link: http://datasim.ov.ingv.it/models/diagno.html
- DISGAS v2.5.3: the dilute gas dispersion simulation tool (Costa et al., 2005; Costa and Macedonio, 2016)
  Link: http://datasim.ov.ingv.it/models/disgas.html
- TWODEE v2.6: the dense gas dispersion simulation tool (Hankin and Britter, 1999; Folch et al., 2009)
  Link: http://datasim.ov.ingv.it/models/twodee.html

Weather data software:

- wgrib2
  Link: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/ 
  The script assumes the executable is in the system PATH
- grib-tools
  Linux: Distributed as part of `eccodes` (required for ERA5 data and reanalysis mode)


- CDSAPI client key (required for ERA5 data and reanalysis mode)
  The user needs to register to: https://cds-beta.climate.copernicus.eu/
  Once the registration is approved, to get the data follow the instructions here: https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+documentation
  The user needs to install the personal key in a .cdsapirc file, to save in different locations depending on the OS. Please read the instructions.


#### Python dependencies

With Conda, it is possible to set a virtual environmnent with all the required dependencies specific for VIGIL. This simplifies the 
installation of the different packages and the management of the Python installation in the system.

Instructions for setting the Conda environment:

1) create the environment with all the needed additional packages:
	`conda env create -f environment.yml`
2) activate the environment with:
	`conda activate vigil`
3) to exit from the environment:
	`conda deactivate`

If required, the name of the environment can be changed by editing
`environment.yml`.

### Running an example case

Example 1 runs DIAGNO against forecast data and can be run without an API key.
Change to the `example_1` directory and run the commands in `commands.txt` in
order.

Ensure that the `simulations` directory is empty before running `weather.py`.

At the end, plots can be found in `post_processing` directory.

## Developers

`VIGIL` was created by and is maintained by:

+ Fabio Dioguardi. British Geological Survey, The Lyell Centre, Edinburgh, United Kingdom. Email: fabiod@bgs.ac.uk
+ Silvia Massaro. Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Bologna, Bologna, Italy. Email: silvia.massaro@ingv.it
+ John A Stevenson (@volcan01010). British Geological Survey, The Lyell Centre, Edinburgh, United Kingdom.

## Licence

VIGIL is distributed under the GPL v3.0 licence. Copyright: © BGS / UKRI 2021

## References

Costa, A., Macedonio, G., Chiodini, G., 2005. Numerical model of gas dispersion emitted from volcanic sources. Annals of Geophysics, vol. 48, n.4/5.
Costa, A. and Macedonio, G., 2016. DISGAS-2.0: a model for passive DISpersion of GAS. Rapporti tecnici INGV. Istituto Nazionale Di Geofisica e Vulcanologia, Italy, 332, 2039e7941.
Dioguardi, F., Massaro, S., Chiodini, G., Costa, A., Folch, A., Macedonio, G., Sandri, L., Selva, J., Tamburello, G., 2022. VIGIL: A Python tool for automatized probabilistic VolcanIc Gas dIspersion modeLling. Annals of Geophysics, 65(1), DM107, https://doi.org/10.4401/ag-8796. 
Douglas, S.G., Kessler, R.C., and Carr, E.L., 1990. User's guide for the Urban Airshed Model. Volume 3. User's manual for the Diagnostic Wind Model (No. PB-91-131243/XAB). Systems Applications, Inc., San Rafael, CA (USA).
Folch, A., Costa, A., and Hankin, R.K., 2009. TWODEE-2: a shallow layer model for dense gas dispersion on complex topography. Computers & Geosciences, 35(3), 667-674.
Hankin, R.K.S., and Britter, R.E., 1999. Twodee: the Health and Safety Laboratory's shallow layer model for heavy gas dispersion Part 3: Experimental validation (Thorney Island). Journal of hazardous materials, 66(3), 239-261.
