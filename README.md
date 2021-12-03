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
usage: weather.py [-h] [-M MODE] [-RT RUN_TYPE] [-CS CONTINUOUS_SIMULATION] [-S START_DATE] [-E END_DATE] [-SY SAMPLED_YEARS [SAMPLED_YEARS ...]] [-SM SAMPLED_MONTHS [SAMPLED_MONTHS ...]] [-SD SAMPLED_DAYS [SAMPLED_DAYS ...]] [-V VOLC] [-LAT LAT] [-LON LON]
                  [-EL ELEV] [-NS SAMPLES] [-ERA5 ERA5] [-WST STATION] [-N NPROC] [-TD TWODEE] [-DG DISGAS]
  -h, --help            show this help message and exit
  -M MODE, --mode MODE  Possible options: reanalysis, forecast. If reanalysis, either ERA5 or WST options should be
                        on. If forecast, GFS data will be downloaded and processed
  -RT RUN_TYPE, --run_type RUN_TYPE
                        Specify if the simulation is a new one or a restart. Possible options are: new, restart
  -CS CONTINUOUS_SIMULATION, --continuous_simulation CONTINUOUS_SIMULATION
                        Specify if the simulation is continuous between the specified start and end dates. Possible options are True or False
  -S START_DATE, --start_date START_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY
  -E END_DATE, --end_date END_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY
  -SY SAMPLED_YEARS [SAMPLED_YEARS ...], --sampled_years SAMPLED_YEARS [SAMPLED_YEARS ...]
                        Specify years to sample from the time interval (comma separated values)
  -SM SAMPLED_MONTHS [SAMPLED_MONTHS ...], --sampled_months SAMPLED_MONTHS [SAMPLED_MONTHS ...]
                        Specify months to sample from the time interval (comma separated values)
  -SD SAMPLED_DAYS [SAMPLED_DAYS ...], --sampled_days SAMPLED_DAYS [SAMPLED_DAYS ...]
                        Specify days to sample from the time interval (comma separated values)
  -V VOLC, --volc VOLC  This is the volcano ID based on the Smithsonian
                        Institute IDs
  -LAT LAT, --lat LAT   Volcano latitude
  -LON LON, --lon LON   Volcano longitude
  -EL ELEV, --elev ELEV
                        Volcano elevation
  -NS SAMPLES, --samples SAMPLES
                        Number of days to sample
  -ERA5 ERA5, --ERA5 ERA5
                        True: Use ERA5 reanalysis. False: Do not use ERA5
                        reanalysis
  -WST STATION, --station STATION
                        True: Use weather station data. False: Do not use
                        weather station data
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes
  -TD TWODEE, --twodee TWODEE
                        on or off, to prepare additional weather data files for Twodee.
  -DG DISGAS, --disgas DISGAS
                        on or off, to run Disgas
```

- run_models.py

```bash
Python script to run Diagno and DISGAS for the days sampled with weather.py. 
The following flags control the execution of hazard_fumaroles.py:
usage: run_models.py [-h] [-N NPROC] [-RT RUN_TYPE] [-CS CONTINUOUS_SIMULATION] [-RS RANDOM_SOURCES] [-NS NSOURCES]
                     [-SINT SOURCES_INTERVAL [SOURCES_INTERVAL ...]] [-SLOC SOURCE_LOCATION [SOURCE_LOCATION ...]]
                     [-SDX SOURCE_DX] [-SDY SOURCE_DY] [-SDUR SOURCE_DUR] [-D DOMAIN [DOMAIN ...]]
                     [-SEM SOURCE_EMISSION] [-RER RANDOM_EMISSION] [-TD TWODEE] [-DG DISGAS]
  -h, --help            show this help message and exit
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes
  -RT RUN_TYPE, --run_type RUN_TYPE
                        Specify if the simulation is a new one or a restart. Possible options are: new, restart
  -CS CONTINUOUS_SIMULATION, --continuous_simulation CONTINUOUS_SIMULATION
                        Specify if the simulation is continuous between the specified start and end dates. Possible options are True or False
  -RS RANDOM_SOURCES, --random_sources RANDOM_SOURCES
                        on: randomly select NS locations from a probability
                        map. off: fixed source locations
  -NS NSOURCES, --nsources NSOURCES
                        Specify a number for a fixed number of sources. If
                        random, then randomly select the number of sources
                        from an interval
  -SINT SOURCES_INTERVAL [SOURCES_INTERVAL ...], --sources_interval SOURCES_INTERVAL [SOURCES_INTERVAL ...]
                        Type the minimum and maximum number of sources  (comma separated values)
  -SLOC SOURCE_LOCATION [SOURCE_LOCATION ...], --source_location SOURCE_LOCATION [SOURCE_LOCATION ...]
                        Coordinate type (UTM/GEO), latitude/northing,
                        longitude/easting, elevation (above ground in m) of 1
                        fixed source (comma separated values)
  -SDX SOURCE_DX, --source_dx SOURCE_DX
                        Extension [m] along the X direction of 1 single source. Option valid for Twodee only
  -SDY SOURCE_DY, --source_dy SOURCE_DY
                        Extension [m] along the Y direction of 1 single source. Option valid for Twodee only
  -SDUR SOURCE_DUR, --source_dur SOURCE_DUR
                        Emission duration [s] of 1 single source. Option valid for Twodee only
  -D DOMAIN [DOMAIN ...], --domain DOMAIN [DOMAIN ...]
                        Coordinates type (UTM/GEO), coordinates
                        (latitude/northing, longitude/easting) of the bottom
                        left corner and top right corner of the domain (comma separated values)
  -SEM SOURCE_EMISSION, --source_emission SOURCE_EMISSION
                        Source emission rate [kg/s]. If specified, it is
                        assigned to all the sources in the domain
  -RER RANDOM_EMISSION, --random_emission RANDOM_EMISSION
                        on: randomly assign emission rate for each source in
                        the domain sampled from a flux.csv file. off: use
                        specified emission rate
  -DI DIAGNO, --diagno DIAGNO
                        on or off, to run Diagno. Turn it off only if Diagno has already been run
  -TD TWODEE, --twodee TWODEE
                        on or off, to run Twodee
  -DG DISGAS, --disgas DISGAS
                        on or off, to run Disgas
```

- post_process.py

```bash
Python script that:
	+ reads the DISGAS outputs produced by hazard_fumaroles.py and based on the list of days simulated, which is read from the file days_list.txt generated by hazard_fumaroles.txt
	+ converts the DISGAS outputs (currently in H2O concentration) in concentration of other gas species specified by the user in input and based on the gas species properties made available by the user in the file gas_properties.csv. The converted concentrations are stored in the folder simulation_converted
	+ calculates the converted outputs at user's specified exceedance probabilities, time steps and vertical levels; these are stored in the folder output_ecdf
	+ plot the converted outputs and those at selected exceedance probabilities at user's selectd time steps and vertical levels; the plots are stored in the folder graphical_outputs.
The following flags control the execution of post_process.py:
usage: post_process.py [-h] [-P PLOT] [-PE PLOT_EX_PROB] [-EX EX_PROB [EX_PROB ...]] [-T TIME_STEPS [TIME_STEPS ...]]
                       [-L LEVELS [LEVELS ...]] [-D DAYS_PLOT [DAYS_PLOT ...]] [-C CONVERT] [-S SPECIES [SPECIES ...]]
                       [-N NPROC] [-M MODELS] [-MO MERGE_OUTPUTS] [-U UNITS] [-PL PLOT_LIMITS [PLOT_LIMITS ...]] [-TA TIME_AV] 
		       [-OF OUTPUT_FORMAT] [-PT PLOT_TOPOGRAPHY] [-PR PLOT_RESOLUTION] [-TP TRACKING_POINTS]

  -h, --help            show this help message and exit
  -P PLOT, --plot PLOT  True: Produce plots of the solutions. False: Do not
                        produce plots
  -PE PLOT_EX_PROB, --plot_ex_prob PLOT_EX_PROB
                        True: Produce plots of the specified exceedance
                        probabilities. False: Do not produce plots
  -EX EX_PROB [EX_PROB ...], --ex_prob EX_PROB [EX_PROB ...]
                        List of exceedance probabilities to be used for
                        graphical output (comma separated values)
  -T TIME_STEPS [TIME_STEPS ...], --time_steps TIME_STEPS [TIME_STEPS ...]
                        List of time steps to plot (integer >= 0). Type all to
                        plot all the time steps (comma separated values)
  -L LEVELS [LEVELS ...], --levels LEVELS [LEVELS ...]
                        List of vertical levels (integer >= 1) to plot. Type
                        all to plot all the levels (comma separated values)
  -D DAYS_PLOT [DAYS_PLOT ...], --days_plot DAYS_PLOT [DAYS_PLOT ...]
                        List of days to plot (YYYYMMDD). Type all to plot all
                        the days (comma separated values)
  -C CONVERT, --convert CONVERT
                        If True, convert output concentration into other species listed with the command -S (--species)
  -S SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]
                        List of gas species (e.g. CO2) (comma separated values)
  -TS TRACKING_SPECIE, --tracking_specie TRACKING_SPECIE
                        The original emitted specie that is tracked in the simulation
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes
  -M MODELS, --models MODELS
                        Model outputs to post-process. Options: disgas, twodee, all
  -MO MERGE_OUTPUTS, --merge_outputs MERGE_OUTPUTS
                        Merge Twodee and Disgas outputs (true or false)
  -U UNITS, --units UNITS
                        Gas concentration units. Possible options are: ppm, kg/m3
  -PL PLOT_LIMITS [PLOT_LIMITS ...], --plot_limits PLOT_LIMITS [PLOT_LIMITS ...]
                        Minimum and maximum value of concentration to display. If unspecified, they are obtained from all the outputs (comma separated values)
  -PI PLOT_ISOLINES, --plot_isolines PLOT_ISOLINES
                        List of gas concentrations values to be used to draw isolines. Optional
  -TA TIME_AV, --time_av TIME_AV
                        Generate time-averaged outputs. Specify the time-averaging interval (in hours), or 0 for averaging over the whole duration
  -OF OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        Select format of the processed output files. Valid options are: GRD
  -PT PLOT_TOPOGRAPHY, --plot_topography PLOT_TOPOGRAPHY
                        Plot topography layer (True or False). Warning, it can be time-consuming!
  -TI TOPOGRAPHY_ISOLINES, --topography_isolines TOPOGRAPHY_ISOLINES
                        Topography height a.s.l. contour lines spatial resolution (in m). Used only if -PT True
  -PR PLOT_RESOLUTION, --plot_resolution PLOT_RESOLUTION
                        Specify plot resolution in dpi
  -TP TRACKING_POINTS, --tracking_points TRACKING_POINTS
                        Extrapolate gas concentration at locations specified in the file tracking_points.txt
```

### Dependencies and installation instructions

#### System dependencies

The following software are required (and should be available on the system
$PATH):

Gas modelling software:

- DIAGNO v1.1.6: the diagnostic wind model (Douglas et al., 1990)
  Link: http://datasim.ov.ingv.it/models/diagno.html
- DISGAS v2.2.2: the dilute gas dispersion simulation tool (Costa et al., 2005; Costa and Macedonio, 2016)
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
  The user needs to register to: https://cds.climate.copernicus.eu/cdsapp#!/home
  Once the registration is approved, to get the data follow the instructions here: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5
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
Douglas, S.G., Kessler, R.C., and Carr, E.L., 1990. User's guide for the Urban Airshed Model. Volume 3. User's manual for the Diagnostic Wind Model (No. PB-91-131243/XAB). Systems Applications, Inc., San Rafael, CA (USA).
Folch, A., Costa, A., and Hankin, R.K., 2009. TWODEE-2: a shallow layer model for dense gas dispersion on complex topography. Computers & Geosciences, 35(3), 667-674.
Hankin, R.K.S., and Britter, R.E., 1999. Twodee: the Health and Safety Laboratory's shallow layer model for heavy gas dispersion Part 3: Experimental validation (Thorney Island). Journal of hazardous materials, 66(3), 239-261.
