# APVGDM - Automatic Probabilistic Volcanic Gas Dispersion Modelling
Fabio Dioguardi. British Geological Survey, The Lyell Centre, Edinburgh, United Kingdom. Email: fabiod@bgs.ac.uk
Silvia Massaro. Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Bologna, Bologna, Italy. Email: silvia.massaro@ingv.it

### PACKAGE CONTENT AND DESCRIPTION ###
The package contains the following files:
- weather.py
Python script that prepares the weather data needed by Diagno to compute the wind field for DISGAS. 
The weather data can be either retrieved from ECMWF ERA5 database or from time series of data from weather stations that are stored in files in the working folder.
	+ ERA5 Reanalysis data
	In this mode, the script is designed to randomly sample N days from a time interval defined by the user. The number of days N is also specified by the user in input. It is possible to retrieve date from one single day by setting the end date equal to the start date and the number of samples to 1
	+ Weather station data
	In this mode, the script is design to extract weather data in the time interval specified by the user from selected weather data file. The script reads the number of files, file location and name from weather_stations.txt; the file data should be stored in the folder weather_stations. Currently, the following format is assumed:
	yyyy mm dd HH MM SS H.R.(%) T-Air(°C) Irrad.(W/m2) D-Vent(°N) V-Vent(km/h) Rain(mm) P-Air(hPa) Batt.(V) T-Air2(°C).
	Only the wind direction (D-Vent) and speed (V-Vent) are used, but it is important to provide the other data fields (placing random numbers or NaN). This will be improved in the future. [GENERALIZATION NEEDED]
The following flags control the execution of weather.py:
usage: weather.py [-h] [-S START_DATE] [-E END_DATE] [-V VOLC] [-LAT LAT]
                  [-LON LON] [-EL ELEV] [-NS SAMPLES] [-ERA5 ERA5]
                  [-WST STATION] [-N NPROC]  [-TD TWODEE] [-DG DISGAS]
  -h, --help            show this help message and exit
  -S START_DATE, --start_date START_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY
  -E END_DATE, --end_date END_DATE
                        Start date of the sampling period. Format: DD/MM/YYYY
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

- run_models.py
Python script to run Diagno and DISGAS for the days sampled with weather.py. 
The following flags control the execution of hazard_fumaroles.py:
usage: run_models.py [-h] [-N NPROC] [-RS RANDOM_SOURCES] [-NS NSOURCES]
                     [-SINT SOURCES_INTERVAL [SOURCES_INTERVAL ...]] [-SLOC SOURCE_LOCATION [SOURCE_LOCATION ...]]
                     [-SDX SOURCE_DX] [-SDY SOURCE_DY] [-SDUR SOURCE_DUR] [-D DOMAIN [DOMAIN ...]]
                     [-SEM SOURCE_EMISSION] [-RER RANDOM_EMISSION] [-TD TWODEE] [-DG DISGAS]
  -h, --help            show this help message and exit
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes
  -RS RANDOM_SOURCES, --random_sources RANDOM_SOURCES
                        on: randomly select NS locations from a probability
                        map. off: fixed source locations
  -NS NSOURCES, --nsources NSOURCES
                        Specify a number for a fixed number of sources. If
                        random, then randomly select the number of sources
                        from an interval
  -SINT SOURCES_INTERVAL [SOURCES_INTERVAL ...], --sources_interval SOURCES_INTERVAL [SOURCES_INTERVAL ...]
                        Type the minimum and maximum number of sources
  -SLOC SOURCE_LOCATION [SOURCE_LOCATION ...], --source_location SOURCE_LOCATION [SOURCE_LOCATION ...]
                        Coordinate type (UTM/GEO), latitude/northing,
                        longitude/easting, elevation (above ground in m) of 1
                        fixed source
  -SDX SOURCE_DX, --source_dx SOURCE_DX
                        Extension [m] along the X direction of 1 single source. Option valid for Twodee only
  -SDY SOURCE_DY, --source_dy SOURCE_DY
                        Extension [m] along the Y direction of 1 single source. Option valid for Twodee only
  -SDUR SOURCE_DUR, --source_dur SOURCE_DUR
                        Emission duration [s] of 1 single source. Option valid for Twodee on
  -D DOMAIN [DOMAIN ...], --domain DOMAIN [DOMAIN ...]
                        Coordinates type (UTM/GEO), coordinates
                        (latitude/northing, longitude/easting) of the bottom
                        left corner and top right corner of the domain
  -SEM SOURCE_EMISSION, --source_emission SOURCE_EMISSION
                        Source emission rate [kg/s]. If specified, it is
                        assigned to all the sources in the domain
  -RER RANDOM_EMISSION, --random_emission RANDOM_EMISSION
                        on: randomly assign emission rate for each source in
                        the domain sampled from a flux.csv file. off: use
                        specified emission rate
  -TD TWODEE, --twodee TWODEE
                        on or off, to run Twodee
  -DG DISGAS, --disgas DISGAS
                        on or off, to run Disgas

- post_process.py
Python script that:
	+ reads the DISGAS outputs produced by hazard_fumaroles.py and based on the list of days simulated, which is read from the file days_list.txt generated by hazard_fumaroles.txt
	+ converts the DISGAS outputs (currently in H2O concentration) in concentration of other gas species specified by the user in input and based on the gas species properties made available by the user in the file gas_properties.csv. The converted concentrations are stored in the folder simulation_converted
	+ calculates the converted outputs at user's specified exceedance probabilities, time steps and vertical levels; these are stored in the folder output_ecdf
	+ plot the converted outputs and those at selected exceedance probabilities at user's selectd time steps and vertical levels; the plots are stored in the folder graphical_outputs.
The following flags control the execution of post_process.py:
usage: post_process.py [-h] [-P PLOT] [-PE PLOT_EX_PROB] [-EX EX_PROB [EX_PROB ...]] [-T TIME_STEPS [TIME_STEPS ...]]
                       [-L LEVELS [LEVELS ...]] [-D DAYS_PLOT [DAYS_PLOT ...]] [-C CONVERT] [-S SPECIES [SPECIES ...]]
                       [-N NPROC] [-M MODELS] [-MO MERGE_OUTPUTS] [-U UNITS] [-TA TIME_AV]

  -h, --help            show this help message and exit
  -P PLOT, --plot PLOT  True: Produce plots of the solutions. False: Do not
                        produce plots
  -PE PLOT_EX_PROB, --plot_ex_prob PLOT_EX_PROB
                        True: Produce plots of the specified exceedance
                        probabilities. False: Do not produce plots
  -EX EX_PROB [EX_PROB ...], --ex_prob EX_PROB [EX_PROB ...]
                        List of exceedence probabilities to be used for
                        graphical output
  -T TIME_STEPS [TIME_STEPS ...], --time_steps TIME_STEPS [TIME_STEPS ...]
                        List of time steps to plot (integer >= 0). Type all to
                        plot all the time steps
  -L LEVELS [LEVELS ...], --levels LEVELS [LEVELS ...]
                        List of vertical levels (integer >= 1) to plot. Type
                        all to plot all the levels
  -D DAYS_PLOT [DAYS_PLOT ...], --days_plot DAYS_PLOT [DAYS_PLOT ...]
                        List of days to plot (YYYYMMDD). Type all to plot all
                        the days
  -C CONVERT, --convert CONVERT
                        If True, convert output concentration into other species listed with the command -S (--species)
  -S SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]
                        List of gas species (e.g. CO2)
  -N NPROC, --nproc NPROC
                        Maximum number of allowed simultaneous processes
  -M MODELS, --models MODELS
                        Model outputs to post-process. Options: disgas, twodee, all
  -MO MERGE_OUTPUTS, --merge_outputs MERGE_OUTPUTS
                        Merge Twodee and Disgas outputs (true or false)
  -U UNITS, --units UNITS
                        Gas concentration units. Possible options are: ppm, kg/m3
  -TA TIME_AV, --time_av TIME_AV
                        Generate time-averaged outputs. Specify the time-averaging interval (in hours), or 0 for averaging over the whole duration

### DEPENDENCIES AND INSTALLATION INSTRUCTIONS ###
The following software are required:
- wgrib2
  Link: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/ 
  The script assumes the executable is in the system PATH
- grib-tools
  Windows: use chocolatey to install it. https://chocolatey.org/packages/grib-tools
  Linux: Install eccodes
  In both cases, add the binaries folder to the system PATH
- Additional python packages needed
  utm, cdsapi, pandas, xlrd
- CDSAPI client key
  The user needs to register to: https://cds.climate.copernicus.eu/cdsapp#!/home
  Once the registration is approved, to get the data follow the instructions here: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5
  The user needs to install the personal key in a .cdsapirc file, to save in different locations depending on the OS. Please read the instructions.

With Conda, it is possible to set a virtual environmnent with all the required dependencies specific for REFIR. This simplifies the 
installation of the different packages and the management of the Python installation in the system.
Instructions for setting the Conda environment:
1) create the environment with all the needed additional packages:
	conda create --name name_of_environment python=3.7 utm cdsapi pandas xlrd
2) activate the environment with:
	conda activate name_of_environment
3) to exit from the environment:
	conda deactivate 