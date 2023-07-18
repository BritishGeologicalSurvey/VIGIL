import datetime
from random import sample
from pathos.multiprocessing import ThreadingPool
from shutil import copy, rmtree
import os
import utm
import argparse
import pandas as pd
import sys


def read_arguments():
    parser = argparse.ArgumentParser(description="Input data", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-M", "--mode", default="reanalysis", help="Possible options: reanalysis, forecast. If "
                        "reanalysis, either ERA5 or WST options should be on. \n If forecast, GFS data will be "
                        "downloaded and processed")
    parser.add_argument('-RT', '--run_type', default='new', help='Specify if the simulation is a new one or a restart. '
                        'Possible options are: new, restart')
    parser.add_argument('-CS', '--continuous_simulation', default='False', help='Specify if the simulation is '
                        'continuous between the specified start and end dates. Possible options are True or False')
    parser.add_argument("-S", "--start_date", default=999, help="Start date of the sampling period. Format: DD/MM/YYYY")
    parser.add_argument("-E", "--end_date", default=999, help="Start date of the sampling period. Format: DD/MM/YYYY")
    parser.add_argument("-SY", "--sampled_years", default='', help="Specify years to sample from the time interval")
    parser.add_argument("-SM", "--sampled_months", default='', help="Specify months to sample from the time interval")
    parser.add_argument("-SD", "--sampled_days", default='', help="Specify days to sample from the time interval")
    parser.add_argument("-V", "--volc", default=999, help="This is the volcano ID based on the Smithsonian Institute "
                        "IDs")
    parser.add_argument("-LAT", "--lat", default=999, help="Volcano latitude")
    parser.add_argument("-LON", "--lon", default=999, help="Volcano longitude")
    parser.add_argument("-EL", "--elev", default=999, help="Volcano elevation")
    parser.add_argument("-NS", "--samples", default=1, help="Number of days to sample")
    parser.add_argument("-ERA5", "--era5", default="False", help="True: Use ERA5 reanalysis. False: Do not use ERA5 "
                        "reanalysis")
    parser.add_argument("-WST", "--station", default="False", help="True: Use weather station data. False: Do not use "
                        "weather station data",)
    parser.add_argument("-N", "--nproc", default=1, help="Maximum number of allowed simultaneous processes")
    parser.add_argument("-TD", "--twodee", default="off", help="on or off, to prepare additional weather data files "
                        "for Twodee.")
    parser.add_argument("-DG", "--disgas", default="off", help="on or off, to run Disgas")
    args = parser.parse_args()
    mode_in = args.mode
    run_type_in = args.run_type
    continuous_simulation_in = args.continuous_simulation
    start_date = args.start_date
    end_date = args.end_date
    sampled_years_input = args.sampled_years
    sampled_months_input = args.sampled_months
    sampled_days_input = args.sampled_days
    volc_id = int(args.volc)
    volc_lat_in = float(args.lat)
    volc_lon_in = float(args.lon)
    elevation_in = float(args.elev)
    era5_on_in = args.era5
    weather_station_on_in = args.station
    nproc = args.nproc
    twodee = args.twodee
    disgas = args.disgas
    mode_in = mode_in.lower()
    if mode_in != "reanalysis" and mode_in != "forecast":
        print("ERROR. Wrong value for variable -M --mode")
        sys.exit()
    run_type_in = run_type_in.lower()
    if run_type_in != 'new' and run_type_in != 'restart':
        print('ERROR. Please provide a valid entry for -RT --run_type')
    if continuous_simulation_in.lower() == "true":
        continuous_simulation_in = True
    elif continuous_simulation_in.lower() == "false":
        continuous_simulation_in = False
    else:
        print("ERROR. Wrong value for variable -CS --continuous_simulation")
        sys.exit()
    if era5_on_in.lower() == "true":
        era5_on_in = True
    elif era5_on_in.lower() == "false":
        era5_on_in = False
    else:
        print("ERROR. Wrong value for variable -ERA5 --ERA5")
        sys.exit()
    if weather_station_on_in.lower() == "true":
        weather_station_on_in = True
    elif weather_station_on_in.lower() == "false":
        weather_station_on_in = False
    else:
        print("ERROR. Wrong value for variable --station")
        sys.exit()
    if mode_in == "reanalysis":
        if not era5_on_in and not weather_station_on_in:
            print("ERROR. Either ERA5 or weather station data should be activated in reanalysis mode")
            sys.exit()
    elif mode_in == "forecast":
        if era5_on_in:
            print("WARNING. ERA5 data cannot be used in forecast mode. Turning ERA5 off")
            era5_on_in = False
        if weather_station_on_in:
            print("WARNING. Weather station data cannot be used in forecast mode. Turning weather stations off")
            weather_station_on_in = False
    if weather_station_on_in and era5_on_in:
        print("ERROR. It is currently not possible to use both reanalysis and weather station data")
        sys.exit()
    if volc_lat_in == 999 and volc_lon_in == 999 and elevation_in == 999:
        try:
            database = pd.read_excel("https://webapps.bgs.ac.uk/research/volcanoes/esp/volcanoExport.xlsx",
                                     sheet_name="volcanoes",)
            nrows = database.shape[0]
            row = 0
            while True:
                if database["SMITHSONIAN_ID"][row] == volc_id:
                    elevation_in = database["ELEVATION_m"][row]
                    volc_lat_in = database["LATITUDE"][row]
                    volc_lon_in = database["LONGITUDE"][row]
                    break
                else:
                    row += 1
                    if row >= nrows:
                        raise NameError('Volcano ID not found')
        except NameError:
            print("Unable to retrieve volcano information from the ESPs database and no location information have been "
                  "provided")
            raise
    sampled_days_in = []
    sampled_months_in = []
    sampled_years_in = []
    if not continuous_simulation_in:
        if sampled_days_input != '':
            sampled_days_in = sampled_days_input.split(',')
            for sampled_day_in in sampled_days_in:
                try:
                    datetime.datetime.strptime(sampled_day_in, '%d')
                except ValueError:
                    print('ERROR. Please provide a valid entry for -SD --sampled_days')
                    sys.exit()
        if sampled_months_input != '':
            sampled_months_in = sampled_months_input.split(',')
            for sampled_month_in in sampled_months_in:
                try:
                    datetime.datetime.strptime(sampled_month_in, '%m')
                except ValueError:
                    print('ERROR. Please provide a valid entry for -SM --sampled_months')
                    sys.exit()
        if sampled_years_input != '':
            sampled_years_in = sampled_years_input.split(',')
            for sampled_year_in in sampled_years_in:
                try:
                    datetime.datetime.strptime(sampled_year_in, '%Y')
                except ValueError:
                    print('ERROR. Please provide a valid entry for -SY --sampled_years')
                    sys.exit()
    try:
        max_number_processes_in = int(nproc)
    except ValueError:
        print("Please provide a valid number for the maximum number of process")
        sys.exit()
    out_utm = utm.from_latlon(volc_lat_in, volc_lon_in)
    easting_in = out_utm[0]
    northing_in = out_utm[1]
    # Check location is inside the topography
    topography_file = os.path.join(root, 'topography.grd')
    with open(topography_file, 'r') as topography:
        topography_content = topography.readlines()
    x0 = float(topography_content[2].split(' ')[0])
    xf = float(topography_content[2].split(' ')[1])
    y0 = float(topography_content[3].split(' ')[0])
    yf = float(topography_content[3].split(' ')[1])
    if easting_in <= x0 or easting_in >= xf or northing_in <= y0 or northing_in >= yf:
        print('Warning. Weather data location outside the domain of the topography file. '
              'Moving the location to the centre of the domain')
        easting_in = (x0 + xf) / 2
        northing_in = (y0 + yf) / 2
    try:
        start_date = start_date.split("/")
        stop_date = end_date.split("/")
        d_start_s = start_date[0]
        mo_start_s = start_date[1]
        y_start_s = start_date[2]
        d_stop_s = stop_date[0]
        mo_stop_s = stop_date[1]
        y_stop_s = stop_date[2]
        analysis_start_in = y_start_s + mo_start_s + d_start_s
        analysis_stop_in = y_stop_s + mo_stop_s + d_stop_s
        time_start_in = datetime.datetime(int(y_start_s), int(mo_start_s), int(d_start_s))
        time_stop_in = datetime.datetime(int(y_stop_s), int(mo_stop_s), int(d_stop_s))
    except IndexError:
        print("Unable to process start and end date")
        print("Please ensure the dates are provided in the format DD/MM/YYYY (e.g. 23/06/2018)")
        sys.exit()
    if not continuous_simulation_in:
        try:
            nsamples_in = int(args.samples)
        except ValueError:
            print('ERROR. Wrong entry for variable -NS --samples')
            sys.exit()
    else:
        nsamples_in = (time_stop_in - time_start_in).days + 1
        run_type_in = 'restart'
    if twodee.lower() == "on":
        twodee_on_in = True
    elif twodee.lower() == "off":
        twodee_on_in = False
    else:
        print("Please provide a valid entry for the variable -TD --twodee")
        sys.exit()
    if disgas.lower() == "on":
        disgas_on_in = True
    elif disgas.lower() == "off":
        disgas_on_in = False
    else:
        print("Please provide a valid entry for the variable -DG --disgas")
        sys.exit()
    return (
        mode_in,
        run_type_in,
        continuous_simulation_in,
        nsamples_in,
        time_start_in,
        time_stop_in,
        analysis_start_in,
        analysis_stop_in,
        era5_on_in,
        weather_station_on_in,
        elevation_in,
        volc_lat_in,
        volc_lon_in,
        easting_in,
        northing_in,
        max_number_processes_in,
        twodee_on_in,
        disgas_on_in,
        sampled_years_in,
        sampled_months_in,
        sampled_days_in,
    )


def automatic_weather(analysis_start_in):
    def era5_retrieve(lon_source, lat_source, retrieved_day):
        def era5_request_pressure(folder_grib, year_pl, month_pl, day_pl):
            import cdsapi

            check_pl_era5 = 1
            grib_file = os.path.join(folder_grib, "pressure_levels.grib")
            print("Downloading file from ERA5 database")
            c = cdsapi.Client()
            try:
                print("Retrieving pressure-levels data")
                c.retrieve(
                    "reanalysis-era5-pressure-levels",
                    {
                        "pressure_level": [
                            "1",
                            "2",
                            "3",
                            "5",
                            "7",
                            "10",
                            "20",
                            "30",
                            "50",
                            "70",
                            "100",
                            "125",
                            "150",
                            "175",
                            "200",
                            "225",
                            "250",
                            "300",
                            "350",
                            "400",
                            "450",
                            "500",
                            "550",
                            "600",
                            "650",
                            "700",
                            "750",
                            "775",
                            "800",
                            "825",
                            "850",
                            "875",
                            "900",
                            "925",
                            "950",
                            "975",
                            "1000",
                        ],
                        "variable": [
                            "geopotential",
                            "u_component_of_wind",
                            "v_component_of_wind",
                            "temperature",
                        ],
                        "time": [
                            "00:00",
                            "01:00",
                            "02:00",
                            "03:00",
                            "04:00",
                            "05:00",
                            "06:00",
                            "07:00",
                            "08:00",
                            "09:00",
                            "10:00",
                            "11:00",
                            "12:00",
                            "13:00",
                            "14:00",
                            "15:00",
                            "16:00",
                            "17:00",
                            "18:00",
                            "19:00",
                            "20:00",
                            "21:00",
                            "22:00",
                            "23:00",
                        ],
                        "product_type": "reanalysis",
                        "year": year_pl,
                        "day": day_pl,
                        "month": month_pl,
                        "area": area,
                        "format": "grib",
                    },
                    grib_file,
                )
            except (Exception, ConnectionError):
                print("Unable to retrieve ERA5 pressure level data")
                check_pl_era5 = 0
            return check_pl_era5

        def era5_request_single(folder_grib, year_sl, month_sl, day_sl):
            import cdsapi

            check_sl_era5 = 1
            print("Downloading file from ERA5 database")
            grib_file = os.path.join(folder_grib, "surface.grib")
            c = cdsapi.Client()
            try:
                print("Retrieving single levels data")
                c.retrieve(
                    "reanalysis-era5-single-levels",
                    {
                        "variable": [
                            "10m_u_component_of_wind",
                            "10m_v_component_of_wind",
                            "2m_temperature",
                            "soil_temperature_level_1",
                            "surface_pressure",
                        ],
                        "time": [
                            "00:00",
                            "01:00",
                            "02:00",
                            "03:00",
                            "04:00",
                            "05:00",
                            "06:00",
                            "07:00",
                            "08:00",
                            "09:00",
                            "10:00",
                            "11:00",
                            "12:00",
                            "13:00",
                            "14:00",
                            "15:00",
                            "16:00",
                            "17:00",
                            "18:00",
                            "19:00",
                            "20:00",
                            "21:00",
                            "22:00",
                            "23:00",
                        ],
                        "product_type": "reanalysis",
                        "year": year_sl,
                        "day": day_sl,
                        "month": month_sl,
                        "area": area,
                        "format": "grib",
                    },
                    grib_file,
                )
            except (Exception, ConnectionError):
                print("Unable to retrieve ERA5 single level data")
                check_sl_era5 = 0
            return check_sl_era5

        retrieved_day_s = str(retrieved_day)
        try:
            year_retday = retrieved_day_s[0:4]
            month_retday = retrieved_day_s[5:7]
            day_retday = retrieved_day_s[8:10]
            data_folder_era5 = os.path.join(simulations, year_retday + month_retday + day_retday)
            slon_source = str(lon_source)
            slat_source = str(lat_source)
            date_bis = year_retday + month_retday + day_retday
            wtfile = os.path.join(data_folder, "weather_data_" + date_bis)
            wtfile_sl = os.path.join(data_folder, "weather_data_sl_" + date_bis)
            lat_n = int(lat_source) + 2
            lat_s = int(lat_source) - 2
            lon_w = int(lon_source) - 2
            lon_e = int(lon_source) + 2
            area = [lat_n, lon_w, lat_s, lon_e]
            # Retrieve files
            check_pl = era5_request_pressure(data_folder_era5, year_retday, month_retday, day_retday)

            check_sl = era5_request_single(data_folder_era5, year_retday, month_retday, day_retday)

            if check_pl == 0:
                with open("log_weather.txt", "a+", encoding="utf-8", errors="surrogateescape") as logger_pl:
                    logger_pl.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                                 ' Unable to download pressure level data for day ' + retrieved_day_s + '\n')
                return
            if check_sl == 0:
                with open("log_weather.txt", "a+", encoding="utf-8", errors="surrogateescape") as logger_sl:
                    logger_sl.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                                 ' Unable to download surface data for day ' + retrieved_day_s + '\n')
                return
            # Convert grib1 to grib2 with the NOAA Perl script. To make it more portable and avoiding the need to set up
            # many paths, I have included in the package also the required files and scripts that are originally
            # available in the grib2 installation folder
            print("Converting grib1 data to grib2")
            pl_grib_file = os.path.join(data_folder_era5, "pressure_levels.grib ")
            sl_grib_file = os.path.join(data_folder_era5, "surface.grib ")
            os_cmd = os.system("srun -n 1 grib_set -s edition=2 " + pl_grib_file + wtfile)
            if os_cmd != 0:
                print('Command srun not found in the system')
                os.system("grib_set -s edition=2 " + pl_grib_file + wtfile)
            os_cmd = os.system("srun -n 1 grib_set -s edition=2 " + sl_grib_file + wtfile_sl)
            if os_cmd != 0:
                print('Command srun not found in the system')
                os.system("grib_set -s edition=2 " + sl_grib_file + wtfile_sl)
            wtfile_prof = os.path.join(data_folder_era5, "profile_" + date_bis + ".txt")
            wtfile_sl_location = os.path.join(data_folder_era5, "data_location_" + date_bis + ".txt")
            print("Saving weather data along the vertical at the vent location")
            os_cmd = os.system("srun -n 1 wgrib2 " + wtfile + " -s -lon " + slon_source + " " + slat_source + "  >"
                               + wtfile_prof)
            if os_cmd != 0:
                print('Command srun not found in the system')
                os.system("wgrib2 " + wtfile + " -s -lon " + slon_source + " " + slat_source + "  >" + wtfile_prof)
            os_cmd = os.system("srun -n 1 wgrib2 " + wtfile_sl + " -s -lon " + slon_source + " " + slat_source + "  >"
                               + wtfile_sl_location)
            if os_cmd != 0:
                print('Command srun not found in the system')
                os.system("wgrib2 " + wtfile_sl + " -s -lon " + slon_source + " " + slat_source + "  >" +
                          wtfile_sl_location)
            # Split wtfile_prof into multiple file, each one for a specific time step
            splitlen = 148
            outputbase = os.path.join(data_folder_era5, "profile_")
            input_prof = open(wtfile_prof, "r", encoding="utf-8", errors="surrogateescape")
            count = 0
            dest = None
            steps = []
            for line_prof in input_prof:
                if count % splitlen == 0:
                    if dest:
                        dest.close()
                    first_line = line_prof.split(":")
                    val = first_line[2].split("d=")
                    dest = open(outputbase + val[1] + ".txt", "w", encoding="utf-8", errors="surrogateescape",)
                    steps.append(val[1])
                dest.write(line_prof)
                count += 1
            input_prof.close()
            dest.close()
            # Split data_location into multiple file, each one for a specific time step
            splitlen = 5
            outputbase = os.path.join(data_folder_era5, "data_location_")
            input_sl = open(wtfile_sl_location, "r", encoding="utf-8", errors="surrogateescape")
            count = 0
            dest = None
            steps = []
            for line_sl in input_sl:
                if count % splitlen == 0:
                    if dest:
                        dest.close()
                    first_line = line_sl.split(":")
                    val = first_line[2].split("d=")
                    dest = open(outputbase + val[1] + ".txt", "w", encoding="utf-8", errors="surrogateescape",)
                    steps.append(val[1])
                dest.write(line_sl)
                count += 1
            input_sl.close()
            dest.close()
        except BaseException:
            with open("log_weather.txt", "a+", encoding="utf-8", errors="surrogateescape") as logger_era5:
                logger_era5.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                             ' Unable to download weather data for day ' + retrieved_day_s + '\n')

    def gfs_retrieve(lon_source, lat_source, nfcst, time_in):
        import urllib.request
        import urllib.error
        from datetime import datetime, timedelta
        import shutil
        import time

        def wtfile_download(url_to_download, wtfile_to_download):
            print("Downloading forecast file " + url)
            try:
                urllib.request.urlretrieve(url_to_download, wtfile_to_download)
            except (urllib.error.HTTPError, urllib.error.URLError):
                try:
                    time.sleep(60)
                    urllib.request.urlretrieve(url_to_download, wtfile_to_download)
                except (urllib.error.HTTPError, urllib.error.URLError):
                    try:
                        time.sleep(60)
                        urllib.request.urlretrieve(url_to_download, wtfile_to_download)
                    except (urllib.error.HTTPError, urllib.error.URLError):
                        print('Unable to retrieve URL ' + url_to_download)
                        sys.exit()

        def interpolate(slon, slat, corner_lon, corner_lat, zoom_to_location, wtfile_in, wtfile_interpolated,
                        wtfile_profile):
            from shutil import copyfile
            if zoom_to_location:
                os_cmd = os.system("srun -n 1 wgrib2 " + wtfile_in + " -set_grib_type same -new_grid_winds earth "
                                   "-new_grid latlon " +
                                   corner_lon + ":400:0.01 " + corner_lat + ":400:0.01 " + wtfile_interpolated)
                if os_cmd != 0:
                    print('Command srun not found in the system')
                    os.system("wgrib2 " + wtfile_in + " -set_grib_type same -new_grid_winds earth -new_grid latlon " +
                          corner_lon + ":400:0.01 " + corner_lat + ":400:0.01 " + wtfile_interpolated)
            else:
                copyfile(wtfile_in, wtfile_interpolated)
            print("Saving weather data along the vertical at the vent location")
            os_cmd = os.system("srun -n 1 wgrib2 " + wtfile_in + " -s -lon " + slon + " " + slat + " >" +
                               wtfile_profile)
            if os_cmd != 0:
                print('Command srun not found in the system')
                os.system("wgrib2 " + wtfile_interpolated + " -s -lon " + slon + " " + slat + " >" + wtfile_profile)

        slon_source_left = str(lon_source - 2)
        slon_source_right = str(lon_source + 2)
        slat_source_bottom = str(lat_source - 2)
        slat_source_top = str(lat_source + 2)
        if lon_source < 0:
            lon_source = 360 + lon_source
        slon_source = str(lon_source)
        slat_source = str(lat_source)
        lon_corner = str(int(lon_source - 2))
        lat_corner = str(int(lat_source - 2))
        urls = []
        wtfiles = []
        wtfiles_int = []
        wtfiles_prof = []
        zooms = []
        lon_corners = []
        lat_corners = []
        slon_sources = []
        slat_sources = []
        if time_in != 999:
            now = time_in
            now_str = str(time_in)
        else:
            now = datetime.utcnow()
            now_str = str(datetime.utcnow())
        day_before = str(time_in - timedelta(1))
        year_now = now_str[0:4]
        month_now = now_str[5:7]
        day_now = now_str[8:10]
        hour = now_str[11:13]
        year_yst = day_before[0:4]
        month_yst = day_before[5:7]
        day_yst = day_before[8:10]
        data_folder_anl = os.path.join(simulations, year_now + month_now + day_now)
        # Find last GFS analysis
        ihour = int(hour)
        if 0 <= ihour < 6:
            ianl = 0
        elif 6 <= ihour < 12:
            ianl = 6
        elif 12 <= ihour < 18:
            ianl = 12
        else:
            ianl = 18
        anl = "{:02d}".format(ianl)
        year_anl = year_now
        month_anl = month_now
        day_anl = day_now
        url = ("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + year_anl + month_anl + day_anl + "/" + anl)
        while True:
            try:
                print(url)
                urllib.request.urlopen(url)
                break
            except urllib.error.HTTPError:
                ianl = ianl - 6
                print("Analysis file at " + anl + "z not yet available. Retrieving the latest available")
            except urllib.error.URLError:
                ianl = ianl - 6
                print("Analysis file at " + anl + "z not yet available. Retrieving the latest available")
            if ianl < 0:  # this is in case the analysis at 00z is not available; in this case, ianl = -6 from above,
                # hence must be corrected. Additionally, the variable ianl will be updated later otherwise it would
                # affect ifcst
                anl = "18"
                ianl = 18
                year_anl = year_yst
                month_anl = month_yst
                day_anl = day_yst
                now -= timedelta(1)
                day_before = str(now - timedelta(1))
                year_yst = day_before[0:4]
                month_yst = day_before[5:7]
                day_yst = day_before[8:10]
            elif 0 <= ianl < 10:
                anl = "0" + str(ianl)
            else:
                anl = str(ianl)
            url = ("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + year_anl + month_anl + day_anl + "/" +
                   anl)
        print("Most up to date GFS analysis: " + url)

        # Retrieve weather data that best matches current time
        time_anl = datetime.strptime(year_anl + month_anl + day_anl + anl, "%Y%m%d%H")
        ifcst = int((time_in - time_anl).total_seconds() / 3600)
        if ianl < 0:
            ianl = 18
            ifcst = ihour + 6
        # Check all forecast files are available
        max_ifcst = ifcst + nfcst
        time_profile = time_in
        while ifcst < max_ifcst:
            fcst = "f" + "{:03d}".format(ifcst)
            wtfile_dwnl = "gfs.t" + anl + "z.pgrb2.0p25." + fcst
            url = ("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + year_anl + month_anl + day_anl + "/" +
                   anl + "/atmos/" + wtfile_dwnl)
            try:
                urllib.request.urlopen(url)
            except urllib.error.HTTPError:
                ianl = ianl - 6
                ifcst = ifcst + 6
                print("Forecast file " + wtfile_dwnl + " not yet available. Retrieving the equivalent from the "
                      "previous forecast")
                urls = []
                wtfiles = []
                max_ifcst = ifcst + nfcst
            except urllib.error.URLError:
                ianl = ianl - 6
                ifcst = ifcst + 6
                print("Forecast file " + wtfile_dwnl + " not yet available. Retrieving the equivalent from the "
                      "previous forecast")
                urls = []
                wtfiles = []
                max_ifcst = ifcst + nfcst
            if ianl < 0:
                ianl = 18
                year_anl = year_yst
                month_anl = month_yst
                day_anl = day_yst
                anl = str(ianl)
            elif 0 <= ianl < 10:
                anl = "0" + str(ianl)
            else:
                anl = str(ianl)
            fcst = "f" + "{:03d}".format(ifcst)
            wtfile_dwnl = "gfs.t" + anl + "z.pgrb2.0p25." + fcst
            year_profile = str(time_profile)[0:4]
            month_profile = str(time_profile)[5:7]
            day_profile = str(time_profile)[8:10]
            hour_profile = str(time_profile)[11:13]
            abs_validity = year_profile + month_profile + day_profile + hour_profile
            wtfile_int = os.path.join(data_folder_anl, "weather_data_interpolated_" + year_anl + month_anl + day_anl +
                                      anl + "_" + fcst,)
            wtfile = os.path.join(data_folder_anl, "weather_data_" + year_anl + month_anl + day_anl + anl + "_" + fcst,)
            wtfile_prof = os.path.join(data_folder_anl, "profile_" + abs_validity + ".txt")
            try:
                url = ("https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=" + wtfile_dwnl
                       + "&all_lev=on&var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on&var_PRES=on&subregion=&leftlon="
                       + slon_source_left + "&rightlon=" + slon_source_right + "&toplat=" + slat_source_top
                       + "&bottomlat=" + slat_source_bottom + "&dir=%2Fgfs." + year_anl + month_anl + day_anl + "%2F"
                       + anl + "%2Fatmos")
                urllib.request.urlopen(url)
                zoom = False
            except (urllib.error.HTTPError, urllib.error.URLError):
                time.sleep(60)
                try:
                    url = ("https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=" + wtfile_dwnl
                           + "&all_lev=on&var_HGT=on&var_TMP=on&var_UGRD=on&var_VGRD=on&var_PRES=on&subregion=&leftlon="
                           + slon_source_left + "&rightlon=" + slon_source_right + "&toplat=" + slat_source_top
                           + "&bottomlat=" + slat_source_bottom + "&dir=%2Fgfs." + year_anl + month_anl + day_anl
                           + "%2F" + anl + "%2Fatmos")
                    urllib.request.urlopen(url)
                    zoom = False
                except (urllib.error.HTTPError, urllib.error.URLError):
                    url = ("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + year_anl + month_anl + day_anl
                           + "/" + anl + "/atmos/" + wtfile_dwnl)
                    urllib.request.urlopen(url)
                    zoom = True
            urls.append(url)
            ifcst += 1
            wtfiles.append(wtfile)
            wtfiles_prof.append(wtfile_prof)
            wtfiles_int.append(wtfile_int)
            zooms.append(zoom)
            lon_corners.append(lon_corner)
            lat_corners.append(lat_corner)
            slon_sources.append(slon_source)
            slat_sources.append(slat_source)
            time_profile += timedelta(hours=1)
        n_downloaded_days = 0
        pools_dwnl = []
        n_pool_dwnl = 0
        while n_downloaded_days <= nfcst:
            end_day = n_downloaded_days + 24
            if end_day > nfcst:
                end_day = nfcst
            n_downloaded_days = end_day
            pools_dwnl.append(n_pool_dwnl)
            n_pool_dwnl += 1
            if n_downloaded_days == nfcst:
                break
        n_downloaded_days = 0
        n_pool_dwnl = 0
        while n_downloaded_days <= nfcst:
            start_day = n_downloaded_days
            end_day = n_downloaded_days + 24
            if end_day > nfcst:
                end_day = nfcst
            try:
                pools_dwnl[n_pool_dwnl] = ThreadingPool(24)
                pools_dwnl[n_pool_dwnl].map(wtfile_download, urls[start_day:end_day], wtfiles[start_day:end_day])
            except BaseException:
                print("Unable to doanload weather data")
                sys.exit()
            n_downloaded_days = end_day
            n_pool_dwnl += 1
            if n_downloaded_days == nfcst:
                break
        try:
            pool_dwnl = ThreadingPool(nfcst)
            pool_dwnl.map(
                interpolate,
                slon_sources,
                slat_sources,
                lon_corners,
                lat_corners,
                zooms,
                wtfiles,
                wtfiles_int,
                wtfiles_prof,
            )
        except BaseException:
            print("Error in weather data interpolation")
        if nfcst > 24:
            simulation_day_current = time_in
            simulation_day = simulation_day_current
            for i_file in range(0, len(wtfiles)):
                simulation_day_s = str(simulation_day)
                year_simday = simulation_day_s[0:4]
                month_simday = simulation_day_s[5:7]
                day_simday = simulation_day_s[8:10]
                data_folder_new = os.path.join(simulations, year_simday + month_simday + day_simday)
                try:
                    os.mkdir(data_folder_new)
                except FileExistsError:
                    print('Folder ' + data_folder_new + ' already exists')
                try:
                    shutil.move(wtfiles[i_file], data_folder_new)
                except shutil.Error:
                    print('File ' + wtfiles[i_file] + ' already present in ' + data_folder_new)
                try:
                    shutil.move(wtfiles_int[i_file], data_folder_new)
                except shutil.Error:
                    print('File ' + wtfiles_int[i_file] + ' already present in ' + data_folder_new)
                try:
                    shutil.move(wtfiles_prof[i_file], data_folder_new)
                except shutil.Error:
                    print('File ' + wtfiles_prof[i_file] + ' already present in ' + data_folder_new)
                simulation_day += timedelta(hours=1)

    def prepare_diagno_files(data_folder_diagno, year_diagno, month_diagno, day_diagno):
        def extract_grib_data(folder_profile, validity_profile, wtfile_profile):
            from math import atan2, pi

            file_profile = open(wtfile_profile, "r", encoding="utf-8", errors="surrogateescape")
            records1 = []
            records2 = []
            nrecords = 0
            for line_profile in file_profile:
                nrecords += 1
                records1.append(line_profile.split(":"))
                records2.append(line_profile.split("val="))
            u_tmp = []
            v_tmp = []
            hgt_tmp = []
            t_tmp = []
            u_profile = []
            v_profile = []
            t_profile = []
            wind_profile = []
            direction_profile = []
            hgt_profile = []
            i_record = 0
            while i_record < nrecords:
                if "mb" in records1[i_record][4] and " mb " not in records1[i_record][4]:
                    if records1[i_record][3] == "UGRD":
                        u_tmp.append(float(records2[i_record][1]))
                    elif records1[i_record][3] == "VGRD":
                        v_tmp.append(float(records2[i_record][1]))
                    elif records1[i_record][3] == "GP":
                        hgt_tmp.append(float(records2[i_record][1]) / 9.8066)
                    elif records1[i_record][3] == "HGT":
                        hgt_tmp.append(float(records2[i_record][1]))
                    elif records1[i_record][3] == "TMP":
                        t_tmp.append(float(records2[i_record][1]))
                i_record += 1
            j = 0
            for i_profile in range(len(u_tmp) - 1, -1, -1):
                if hgt_tmp[i_profile] > elevation:
                    hgt_profile.append(hgt_tmp[i_profile] - elevation)
                    u_profile.append(u_tmp[i_profile])
                    v_profile.append(v_tmp[i_profile])
                    t_profile.append(t_tmp[i_profile])
                    wind_profile.append((u_profile[j] ** 2 + v_profile[j] ** 2) ** 0.5)
                    wind_dir_degrees = atan2(u_profile[j], v_profile[j]) * 180 / pi
                    direction_profile.append(wind_dir_degrees + 180)
                    j += 1
                    if j > 19:
                        break
                else:
                    continue
            prof_file = os.path.join(folder_profile, "profile_data_" + validity_profile + ".txt")
            wt_output = open(prof_file, "w", encoding="utf-8", errors="surrogateescape")
            wt_output.write("HGT[m abg] U[m/s]     V[m/s]  WIND[m/s]  WIND_DIR[deg]      T[K]\n")
            for i_profile in range(j - 1, -1, -1):
                wt_output.write("%8.2f %8.2f %10.2f %10.2f %10.2f %10.2f\n" % (hgt_profile[i_profile],
                                                                               u_profile[i_profile],
                                                                               v_profile[i_profile],
                                                                               wind_profile[i_profile],
                                                                               direction_profile[i_profile],
                                                                               t_profile[i_profile]))
            wt_output.close()
            gamma_pl_profile = -(t_profile[-1] - t_profile[0]) / ((hgt_profile[-1] - hgt_profile[0]) / 1000)
            return wind_profile, direction_profile, hgt_profile, gamma_pl_profile

        def extract_grib_data_sl(folder_sl, validity_sl, wtfile_sl):
            from math import atan2, pi

            file_sl = open(wtfile_sl, "r", encoding="utf-8", errors="surrogateescape")
            records1 = []
            records2 = []
            u_sl = None
            v_sl = None
            t2m_sl = None
            pz0_sl = None
            tz0_sl = None
            nrecords = 0
            for line_sl in file_sl:
                nrecords += 1
                records1.append(line_sl.split(":"))
                records2.append(line_sl.split("val="))
            i_sl = 0
            if mode == "reanalysis":
                while i_sl < nrecords:
                    if records1[i_sl][3] == "UGRD":
                        u_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "VGRD":
                        v_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "TMP":
                        t2m_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "PRES":
                        pz0_sl = float(records2[i_sl][1])
                    else:
                        tz0_sl = float(records2[i_sl][1])
                    i_sl += 1
            else:
                while i_sl < nrecords:
                    if records1[i_sl][3] == "UGRD" and records1[i_sl][4] == "10 m above ground":
                        u_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "VGRD" and records1[i_sl][4] == "10 m above ground":
                        v_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "TMP" and records1[i_sl][4] == "2 m above ground":
                        t2m_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "PRES" and records1[i_sl][4] == "surface":
                        pz0_sl = float(records2[i_sl][1])
                    elif records1[i_sl][3] == "TMP" and records1[i_sl][4] == "surface":
                        tz0_sl = float(records2[i_sl][1])
                    i_sl += 1
            if u_sl is None or v_sl is None or t2m_sl is None or pz0_sl is None or tz0_sl is None:
                print('Unable to read all necessary meteorological parameters from ' + wtfile_sl_location_step)
                sys.exit()
            wind_sl_sl = (u_sl ** 2 + v_sl ** 2) ** 0.5
            wind_dir_degrees = atan2(u_sl, v_sl) * 180 / pi
            direction_sl_sl = wind_dir_degrees + 180
            prof_file = os.path.join(folder_sl, "data_location_data_" + validity_sl + ".txt")
            wt_output = open(prof_file, "w", encoding="utf-8", errors="surrogateescape")
            wt_output.write("    U[m/s]     V[m/s]  WIND[m/s]  WIND_DIR[deg]    T2m[K]    Tz0[K]   Pz0[Pa]\n")
            wt_output.write("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n" % (u_sl, v_sl, wind_sl_sl,
                                                                                    direction_sl_sl, t2m_sl, tz0_sl,
                                                                                    pz0_sl))
            wt_output.close()
            gamma_sl_sl = (t2m - tz0) / 2.0
            return u_sl, v_sl, t2m_sl, wind_sl, direction_sl, tz0_sl, gamma_sl_sl, pz0_sl

        files_list = os.listdir(data_folder_diagno)
        path = os.path.normpath(data_folder_diagno)
        splitted_path = path.split(os.sep)
        retrieved_day_s = splitted_path[-1]
        profile_data_files = []
        surface_data_files = []
        for file in files_list:
            if "profile_" in file and "profile_" + year_diagno + month_diagno + day_diagno + ".txt" != file:
                profile_data_files.append(file)
            if "data_location_" in file and "data_location_" + year_diagno + month_diagno + day_diagno + ".txt" != file:
                surface_data_files.append(file)
        profile_data_files = sorted(profile_data_files)
        surface_data_files = sorted(surface_data_files)
        wind_direction_string = ""
        wind_speed_string = ""
        heights_string = ""
        gamma_string_1 = ""
        gamma_string_2 = ""
        gamma_string_3 = ""
        tref_vector = []
        tsoil_vector = []
        press_vector = []
        try:
            diagno_preupr = open(os.path.join(data_folder_diagno, "preupr.dat"), "w", encoding="utf-8",
                                 errors="surrogateescape",)
            diagno_preupr.write(str(n_stations) + "        NSTA\n")
            diagno_preupr.write("20         LEVELS\n")
            diagno_preupr.write("0          NSTRHR\n")
            diagno_preupr.write("23         NENDHR\n")
            diagno_preupr.write("0.5          TDIF\n")
            diagno_preupr.write("13          NCELL\n")
            diagno_preupr.write("0. 1. 2. 4. 8. 16. 24. 32. 40. 60. 80. 100. 250. 500.  CELLZB\n")
            diagno_preupr.write(year_diagno[2:4] + "       KYEAR\n")
            diagno_preupr.write(str(int(month_diagno)) + "       KMONTH\n")
            diagno_preupr.write(str(int(day_diagno)) + "       KDAY\n")
            diagno_preupr.write("2          IOPT\n")
            diagno_preupr.write("ST01" + " " + "{0:7.1f}".format(easting / 1000) + "{0:7.1f}".format(northing / 1000) +
                                "{0:7.1f}".format(elevation) + "\n")
            for file in profile_data_files:
                validity = file.split("profile_")[1]
                validity = validity.split(".txt")[0]
                year_diagno = validity[0:4]
                month_diagno = validity[4:6]
                day_diagno = validity[6:8]
                hour = validity[8:10]
                wtfile_prof_step = os.path.join(data_folder_diagno, file)
                # Extract and elaborate weather data
                wind, direction, height, gamma_pl = extract_grib_data(data_folder_diagno, validity, wtfile_prof_step)
                for i_wind in range(0, len(wind) - 1):
                    try:
                        heights_string += "{:>5}".format(str(int(round(height[i]))))
                    except TypeError:
                        heights_string += "   -1"
                    try:
                        if wind[i_wind] > 50:  # This is the maximum speed limit accepted by DIAGNO
                            wind_speed_string += "   -1"
                        else:
                            wind_speed_string += "{:>5}".format(str(int(round(wind[i_wind] * 10))))
                    except TypeError:
                        wind_speed_string += "   -1"
                    try:
                        wind_direction_string += "{:>5}".format(str(int(round(direction[i_wind]))))
                    except TypeError:
                        wind_direction_string += "   -1"
                diagno_preupr.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "ST01" + " " + hour + "00" + " " + "MET" +
                                    heights_string + "\n")
                diagno_preupr.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "ST01" + " " + hour + "00" + " " + "DEG" +
                                    wind_direction_string + "\n")
                diagno_preupr.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "ST01" + " " + hour + "00" + " " + "MPS" +
                                    wind_speed_string + "\n")
                wind_direction_string = ""
                wind_speed_string = ""
                heights_string = ""
            wind_direction_sl_string = ""
            wind_speed_sl_string = ""
            um_string_1 = ""
            um_string_2 = ""
            um_string_3 = ""
            vm_string_1 = ""
            vm_string_2 = ""
            vm_string_3 = ""
            tinf = 0
            if mode == "reanalysis":
                files_to_iterate = surface_data_files
                file_keyword = "data_location_"
            else:
                files_to_iterate = profile_data_files
                file_keyword = "profile_"
            for file in files_to_iterate:
                validity = file.split(file_keyword)[1]
                validity = validity.split(".txt")[0]
                year_diagno = validity[0:4]
                month_diagno = validity[4:6]
                day_diagno = validity[6:8]
                # hour = validity[8:10]
                wtfile_sl_location_step = os.path.join(data_folder_diagno, file)
                u, v, t2m, wind_sl, direction_sl, tz0, gamma_sl, pz0 = extract_grib_data_sl(data_folder_diagno,
                                                                                            validity,
                                                                                            wtfile_sl_location_step)
                tinf += 0.5 * (tz0 + t2m)
                tref_vector.append(t2m)
                tsoil_vector.append(tz0)
                press_vector.append(pz0)
                # Extract and elaborate weather data
                try:
                    wind_speed_sl_string += "{:>3}".format(str(int(round(wind_sl * 10))))
                except TypeError:
                    wind_speed_sl_string += " -1"
                try:
                    wind_direction_sl_string += "{:>3}".format(str(int(round(direction_sl))))
                except TypeError:
                    wind_direction_sl_string += " -1"
                if len(um_string_1) < 48:
                    um_string_1 += "{:<6.1f}".format(u)
                else:
                    if len(um_string_2) < 48:
                        um_string_2 += "{:<6.1f}".format(u)
                    else:
                        um_string_3 += "{:<6.1f}".format(u)
                if len(vm_string_1) < 48:
                    vm_string_1 += "{:<6.1f}".format(v)
                else:
                    if len(vm_string_2) < 48:
                        vm_string_2 += "{:<6.1f}".format(v)
                    else:
                        vm_string_3 += "{:<6.1f}".format(v)
                if len(gamma_string_1) < 56:
                    gamma_string_1 += "{:<7.1f}".format(gamma_sl)
                else:
                    if len(gamma_string_2) < 56:
                        gamma_string_2 += "{:<7.1f}".format(gamma_sl)
                    else:
                        gamma_string_3 += "{:<7.1f}".format(gamma_sl)
            tinf = tinf / float(len(profile_data_files))
            str_tinf = "{:<4.1f}".format(tinf)
            str_tinf += "       TINF\n"
            um_string_1 += "     UM\n"
            um_string_2 += "     UM\n"
            um_string_3 += "     UM\n"
            vm_string_1 += "     VM\n"
            vm_string_2 += "     VM\n"
            vm_string_3 += "     VM\n"
            gamma_string_1 += "         GAMMA (K/km)\n"
            gamma_string_2 += "         GAMMA (K/km)\n"
            gamma_string_3 += "         GAMMA (K/km)\n"
            # Memorize the needed records in the original presfc.dat
            try:
                presfc_file_records = []
                diagno_presfc = open(os.path.join(data_folder_diagno, "presfc.dat"), "r", encoding="utf-8",
                                     errors="surrogateescape",)
                for line_presfc in diagno_presfc:
                    presfc_file_records.append(line_presfc)
                diagno_presfc.close()
                os.remove("presfc.dat")
                diagno_presfc = open(os.path.join(data_folder_diagno, "presfc.dat"), "w", encoding="utf-8",
                                     errors="surrogateescape",)
                diagno_presfc.write(str(n_stations) + "        NSTA\n")
                diagno_presfc.write("0        NSTRHR\n")
                diagno_presfc.write("23       NENDHR\n")
                diagno_presfc.write("0.5      TDIF\n")
                diagno_presfc.write(year_diagno[2:4] + "       KYEAR\n")
                diagno_presfc.write(str(int(month_diagno)) + "       KMONTH\n")
                diagno_presfc.write(str(int(day_diagno)) + "       KDAY\n")
                diagno_presfc.write(presfc_file_records[7])
                diagno_presfc.write("ST02" + "  " + "{0:7.1f}".format(easting / 1000) +
                                    "{0:7.1f}".format(northing / 1000) + "\n")
                diagno_presfc.write(presfc_file_records[8])
                diagno_presfc.write(presfc_file_records[9])
                diagno_presfc.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "  ST02" + "  WD" + "  DEG  " +
                                    wind_direction_sl_string + "\n")
                diagno_presfc.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "  ST02" + "  WS" + "  MPS  " +
                                    wind_speed_sl_string + "\n")
            except BaseException:
                diagno_presfc = open(os.path.join(data_folder_diagno, "presfc.dat"), "w", encoding="utf-8",
                                     errors="surrogateescape",)
                diagno_presfc.write(str(n_stations) + "        NSTA\n")
                diagno_presfc.write("0        NSTRHR\n")
                diagno_presfc.write("23       NENDHR\n")
                diagno_presfc.write("0.5      TDIF\n")
                diagno_presfc.write(year_diagno[2:4] + "       KYEAR\n")
                diagno_presfc.write(str(int(month_diagno)) + "       KMONTH\n")
                diagno_presfc.write(str(int(day_diagno)) + "       KDAY\n")
                diagno_presfc.write("ST02" + "  " + "{0:7.1f}".format(easting / 1000) +
                                    "{0:7.1f}".format(northing / 1000) + "\n")
                diagno_presfc.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "  ST02" + "  WD" + "  DEG  " +
                                    wind_direction_sl_string + "\n")
                diagno_presfc.write("{:2}".format(year_diagno[2:4]) + "{:>2}".format(str(int(month_diagno))) +
                                    "{:>2}".format(str(int(day_diagno))) + "  ST02" + "  WS" + "  MPS  " +
                                    wind_speed_sl_string + "\n")
            try:
                diagno_records = []
                diagno = open(os.path.join(data_folder_diagno, "diagno.inp"), "r", encoding="utf-8",
                              errors="surrogateescape",)
                for line_diagno in diagno:
                    diagno_records.append(line_diagno)
                diagno_records[12] = '    2         NWIND\n'
                diagno_records[42] = gamma_string_1
                diagno_records[43] = gamma_string_2
                diagno_records[44] = gamma_string_3
                diagno_records[47] = str_tinf
                diagno_records[51] = um_string_1
                diagno_records[52] = um_string_2
                diagno_records[53] = um_string_3
                diagno_records[54] = vm_string_1
                diagno_records[55] = vm_string_2
                diagno_records[56] = vm_string_3
                with open(os.path.join(data_folder_diagno, "diagno.inp"), "w", encoding="utf-8",
                          errors="surrogateescape",) as diagno:
                    diagno.writelines(diagno_records)
            except BaseException:
                print("Unable to process diagno.inp")
        except BaseException:
            with open("log_weather.txt", "a+", encoding="utf-8", errors="surrogateescape") as logger_diagno:
                logger_diagno.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                             ' Unable to prepare DIAGNO input file for day ' + retrieved_day_s + '\n')
        new_file_list = os.listdir(data_folder_diagno)
        for file in new_file_list:
            if file.startswith('data_') or file.startswith('profile_') or file.startswith('weather_') or \
                    file.endswith('.grb'):
                os.remove(os.path.join(data_folder_diagno, file))
        return tref_vector, tsoil_vector, press_vector

    def extract_station_data(station_data_files_in, eastings_in, northings_in, zst_in, data_folder_in, time_in):
        global n_stations
        import math

        n_wt_stations = len(station_data_files_in)
        wind_direction_station_strings = []
        wind_speed_station_strings = []
        stations_id = []
        um_string_1 = ""
        um_string_2 = ""
        um_string_3 = ""
        vm_string_1 = ""
        vm_string_2 = ""
        vm_string_3 = ""
        gamma_string_1 = ""
        gamma_string_2 = ""
        gamma_string_3 = ""
        um_avg = []
        vm_avg = []
        gamma_avg = []
        tref_vector = []
        tsoil_vector = []
        press_vector = []
        um = [[0 for x in range(24)] for y in range(n_wt_stations)]
        vm = [[0 for x in range(24)] for y in range(n_wt_stations)]
        gamma = [[0 for x in range(24)] for y in range(n_wt_stations)]
        year_start = str(time_in)[0:4]
        month_start = str(int(str(time_in)[5:7]))
        day_start = str(int(str(time_in)[8:10]))
        ii = 0
        tinf = 0
        n_tinf = 0
        for station_file in station_data_files_in:
            file = open(station_file, encoding="utf-8", errors="surrogateescape")
            records_wst = []
            wind_direction_station_string = ""
            wind_speed_station_string = ""
            n_line = 1
            for line_wst in file:
                if n_line != 1:  # Skip header line
                    records_wst.append(line_wst.split(","))
                n_line += 1
            jj = 0
            for i_wst in range(0, len(records_wst)):
                try:
                    time_record = records_wst[i_wst][0][0:8]
                    day_record = datetime.datetime.strptime(time_record, "%Y%m%d")
                except ValueError:
                    continue
                if day_record == time_in:
                    try:
                        wind_speed = float(records_wst[i_wst][3]) * 0.28
                        wind_station = int(round(float(records_wst[i_wst][3]) * 0.28 * 10))
                        wind_speed_station_string += "{:>3}".format(str(wind_station))
                    except TypeError:
                        wind_speed = 0
                        wind_speed_station_string += " -1"
                    try:
                        wind_direction = math.radians(float(records_wst[i_wst][2]))
                        wind_direction_station = int(round(float(records_wst[i_wst][2])))
                        wind_direction_station_string += "{:>3}".format(str(wind_direction_station))
                    except TypeError:
                        wind_direction = 0
                        wind_direction_station_string += " -1"
                    try:
                        pressure = float(records_wst[i_wst][4]) * 100
                    except TypeError:
                        pressure = 101300
                    press_vector.append(pressure)
                    try:
                        tref_wst = float(records_wst[i_wst][1]) + 273.15
                    except ValueError:
                        tref_wst = 273.15
                    tref_vector.append(tref_wst)
                    try:
                        tsoil_wst = float(records_wst[i_wst][5]) + 273.15
                    except ValueError:
                        tsoil_wst = tref_wst
                    if math.isnan(tsoil_wst):
                        tsoil_wst = tref_wst
                    tsoil_vector.append(tsoil_wst)
                    try:
                        um[ii][jj] = wind_speed * math.cos(wind_direction)
                        vm[ii][jj] = wind_speed * math.sin(wind_direction)
                        gamma[ii][jj] = (tref_wst - tsoil_wst) / zst_in[ii]
                    except ValueError:
                        um[ii][jj] = 0
                        vm[ii][jj] = 0
                        gamma[ii][jj] = 0
                    jj += 1
                    tinf += 0.5 * (tref_wst + tsoil_wst)
                    n_tinf += 1
                else:
                    # Skip line
                    continue
            wind_direction_station_strings.append(wind_direction_station_string)
            wind_speed_station_strings.append(wind_speed_station_string)
            ii += 1
        for j in range(0, 24):
            um_average = 0
            vm_average = 0
            gamma_average = 0
            for i_st in range(0, n_wt_stations):
                um_average += um[i_st][j]
                vm_average += vm[i_st][j]
                gamma_average += gamma[i_st][j]
            um_average = um_average / n_wt_stations
            vm_average = vm_average / n_wt_stations
            gamma_average = gamma_average / n_wt_stations
            um_avg.append(um_average)
            vm_avg.append(vm_average)
            gamma_avg.append(gamma_average)
        for j in range(0, 24):
            if len(um_string_1) < 48:
                um_string_1 += "{:<6.1f}".format(um_avg[j])
            else:
                if len(um_string_2) < 48:
                    um_string_2 += "{:<6.1f}".format(um_avg[j])
                else:
                    um_string_3 += "{:<6.1f}".format(um_avg[j])
            if len(vm_string_1) < 48:
                vm_string_1 += "{:<6.1f}".format(vm_avg[j])
            else:
                if len(vm_string_2) < 48:
                    vm_string_2 += "{:<6.1f}".format(vm_avg[j])
                else:
                    vm_string_3 += "{:<6.1f}".format(vm_avg[j])
            if len(gamma_string_1) < 56:
                gamma_string_1 += "{:<7.1f}".format(gamma_avg[j])
            else:
                if len(gamma_string_2) < 56:
                    gamma_string_2 += "{:<7.1f}".format(gamma_avg[j])
                else:
                    gamma_string_3 += "{:<7.1f}".format(gamma_avg[j])
        um_string_1 += "     UM\n"
        um_string_2 += "     UM\n"
        um_string_3 += "     UM\n"
        vm_string_1 += "     VM\n"
        vm_string_2 += "     VM\n"
        vm_string_3 += "     VM\n"
        gamma_string_1 += "         GAMMA (K/km)\n"
        gamma_string_2 += "         GAMMA (K/km)\n"
        gamma_string_3 += "         GAMMA (K/km)\n"
        tinf = tinf / n_tinf
        str_tinf = "{:<4.1f}".format(tinf)
        str_tinf += "       TINF\n"
        n_stations = n_wt_stations
        for n_station in range(0, n_stations):
            stations_id.append("ST" + "{:02d}".format(n_station + 1))
            if (wind_speed_station_strings[n_station]
                and wind_direction_station_strings[n_station]
                    == " -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1"):
                print("No valid weather data available for the specified time interval in " +
                      station_data_files_in[n_station])
                stations_id.pop(-1)
                wind_speed_station_strings.pop(n_station)
                wind_direction_station_strings.pop(n_station)
                n_stations -= 1
        diagno_presfc = open(os.path.join(data_folder_in, "presfc.dat"), "w", encoding="utf-8",
                             errors="surrogateescape",)
        diagno_presfc.write(str(n_stations) + "        NSTA\n")
        diagno_presfc.write("0        NSTRHR\n")
        diagno_presfc.write("23       NENDHR\n")
        diagno_presfc.write("0.5      TDIF\n")
        diagno_presfc.write(year_start[2:4] + "       KYEAR\n")
        diagno_presfc.write(month_start + "       KMONTH\n")
        diagno_presfc.write(day_start + "       KDAY\n")
        for n_station in range(0, n_stations):
            diagno_presfc.write(stations_id[n_station] + "  " + "{0:7.1f}".format(eastings_in[n_station]) +
                                "{0:7.1f}".format(northings_in[n_station]) + "\n")
        for n_station in range(0, n_stations):
            station_id = "  ST" + "{:02d}".format(n_station + 1)
            diagno_presfc.write("{:2}".format(year_start[2:4]) + "{:>2}".format(month_start) +
                                "{:>2}".format(day_start) + station_id + "  WD" + "  DEG  " +
                                wind_direction_station_strings[n_station] + "\n")
            diagno_presfc.write("{:2}".format(year_start[2:4]) + "{:>2}".format(month_start) + "{:>2}".format(day_start)
                                + station_id + "  WS" + "  MPS  " + wind_speed_station_strings[n_station] + "\n")
        try:
            diagno_records = []
            diagno = open(os.path.join(data_folder_in, "diagno.inp"), "r", encoding="utf-8", errors="surrogateescape",)
            for line_diagno in diagno:
                diagno_records.append(line_diagno)
            diagno_records[12] = '    ' + str(n_stations) + '         NWIND\n'
            diagno_records[42] = gamma_string_1
            diagno_records[43] = gamma_string_2
            diagno_records[44] = gamma_string_3
            diagno_records[47] = str_tinf
            diagno_records[51] = um_string_1
            diagno_records[52] = um_string_2
            diagno_records[53] = um_string_3
            diagno_records[54] = vm_string_1
            diagno_records[55] = vm_string_2
            diagno_records[56] = vm_string_3
            with open(os.path.join(data_folder_in, "diagno.inp"), "w", encoding="utf-8", errors="surrogateescape",) as \
                    diagno:
                diagno.writelines(diagno_records)
        except BaseException:
            print("Unable to process diagno.inp")
        return tref_vector, tsoil_vector, press_vector

    def save_surface_data(tref_surf, tsoil_surf, press_surf, data_folder_surf):
        with open(os.path.join(data_folder_surf, "surface_data.txt"), "w", encoding="UTF-8") as surface_data:
            surface_data.write("Time [HHMM] Tref[K]  Tz0[K]   P[Pa]\n")
            for i_surf in range(0, len(tref_surf)):
                surface_data.write("{:02d}".format(i_surf) + "00" + "\t" + str(tref_surf[i_surf]) + "\t" +
                                   str(tsoil_surf[i_surf]) + "\t" + str(press_surf[i_surf]) + "\n")

    tref = None
    tsoil = None
    press = None
    analysis_start_s = str(analysis_start_in)
    year = analysis_start_s[0:4]
    month = analysis_start_s[5:7]
    day = analysis_start_s[8:10]
    data_folder = os.path.join(root, "simulations", year + month + day)
    try:
        os.mkdir(data_folder)
    except FileExistsError:
        print("Folder " + data_folder + " already exists")
    try:
        copy("diagno.inp", os.path.join(data_folder, "diagno.inp"))
    except FileNotFoundError:
        print("File diagno.inp not found")
    if mode == "forecast":
        if analysis_start_in == time_start:
            print("Retrieving GFS data for day " + str(analysis_start_in)[0:10])
            gfs_retrieve(volc_lon, volc_lat, ((time_stop - time_start).days + 1) * 24, analysis_start_in)
    if era5_on:
        print("Retrieving ERA5 data for day " + str(analysis_start_in)[0:10])
        era5_retrieve(volc_lon, volc_lat, analysis_start_in)
    if mode == "forecast" or era5_on:
        tref, tsoil, press = prepare_diagno_files(data_folder, year, month, day)
    if weather_station_on:
        stations_input = open("weather_stations_list.txt", "r", encoding="utf-8-sig", errors="surrogateescape",)
        print("Analysing weather station data for day " + str(analysis_start_in)[0:10])
        records = []
        for line in stations_input:
            records.append(line.split("\n")[0])
        n_weather_stations = int(records[0])
        i = 1
        eastings = []
        northings = []
        station_data_files = []
        zst = []
        while i <= n_weather_stations:
            station_lat = float(records[i].split("\t")[0])
            station_lon = float(records[i].split("\t")[1])
            out_utm = utm.from_latlon(station_lat, station_lon)
            eastings.append(out_utm[0] / 1000)
            northings.append(out_utm[1] / 1000)
            zst.append(float(records[i].split("\t")[2]))
            station_data_file = records[i].split("\t")[3]
            try:
                test = open(os.path.join(root, "weather_stations", station_data_file), "r")
                test.close()
                station_data_files.append(os.path.join(root, "weather_stations", station_data_file))
                i += 1
            except FileNotFoundError:
                i += 1
                continue
        tref, tsoil, press = extract_station_data(station_data_files, eastings, northings, zst, data_folder,
                                                  analysis_start_in)
    if tref is None or tsoil is None or press is None:
        print('Unable to save surface data because some data is missing')
        sys.exit()
    save_surface_data(tref, tsoil, press, data_folder)
    if not os.path.exists(os.path.join(data_folder, "upper.dat")):
        with open(os.path.join(data_folder, "upper.dat"), "w", encoding="utf-8") as fake_upper:  # To not let DIAGNO
            # crash when using just presfc
            fake_upper.write("Fake upper.dat file")
        fake_upper.close()


root = os.getcwd()
simulations = os.path.join(root, "simulations")
(
    mode,
    run_type,
    continuous_simulation,
    nsamples,
    time_start,
    time_stop,
    analysis_start,
    analysis_stop,
    era5_on,
    weather_station_on,
    elevation,
    volc_lat,
    volc_lon,
    easting,
    northing,
    max_number_processes,
    twodee_on,
    disgas_on,
    sampled_years,
    sampled_months,
    sampled_days,
) = read_arguments()


n_stations = 1  # FABIO: to generalize to allow more weather data point when using NWP data

if not disgas_on and not twodee_on:
    print("Both DISGAS and TWODEE are turned off")
    sys.exit()

delta = time_stop - time_start
days_list = []
sampled_periods_end = []
# Generate a list of days in the time interval defined above
sampled_years_int = [int(i) for i in sampled_years]
sampled_months_int = [float(i) for i in sampled_months]
sampled_days_int = [float(i) for i in sampled_days]
for i in range(delta.days + 1):
    day_datetime = time_start + datetime.timedelta(days=i)
    if day_datetime.year in sampled_years_int or len(sampled_years_int) == 0:
        if day_datetime.month in sampled_months_int or len(sampled_months_int) == 0:
            if day_datetime.day in sampled_days_int or len(sampled_days_int) == 0:
                days_list.append(day_datetime)
# Create a new list of nsamples randomly-sampled periods
try:
    sampled_period_days = sample(days_list, nsamples)
except ValueError:
    print("ERROR. The sample is larger than the population of days")
    sys.exit()
try:
    os.mkdir(simulations)
except FileExistsError:
    print("Folder " + simulations + " already exists")
    if run_type == 'new':
        for filename in os.listdir(simulations):
            file_path = os.path.join(simulations, filename)
            rmtree(file_path)
if os.path.exists("log_weather.txt"):
    os.remove("log_weather.txt")
if os.path.exists("days_list.txt"):
    os.remove("days_list.txt")
logger = open("log_weather.txt", "w", encoding="utf-8", errors="surrogateescape")
logger.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
             " Unable to complete weather data processing for the following days\n")
# Create a text file with the list of sampled days
days_list_file = open("days_list.txt", "a+", encoding="utf-8", errors="surrogateescape")
for i in range(nsamples):
    days_list_file.write(str(sampled_period_days[i]) + "\n")
days_list_file.close()
if mode == 'forecast' and continuous_simulation:
    for sampled_day in sorted(sampled_period_days):
        automatic_weather(sampled_day)
else:
    n_elaborated_days = 0
    pools = []
    n_pool = 0
    while n_elaborated_days <= nsamples:
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > nsamples:
            end = nsamples
        n_elaborated_days = end
        pools.append(n_pool)
        n_pool += 1
        if n_elaborated_days == nsamples:
            break
    n_elaborated_days = 0
    n_pool = 0
    while n_elaborated_days <= nsamples:
        start = n_elaborated_days
        end = n_elaborated_days + max_number_processes
        if end > nsamples:
            end = nsamples
        try:
            pools[n_pool] = ThreadingPool(max_number_processes)
            pools[n_pool].map(automatic_weather, sampled_period_days[start:end])
        except BaseException:
            print("Unable to process weather data")
            sys.exit()
        n_elaborated_days = end
        n_pool += 1
        if n_elaborated_days == nsamples:
            break

    attempt = 1
    while attempt < 10:
        days_to_reelaborate = []
        simulation_folders = os.listdir(simulations)
        if mode == "reanalysis" and weather_station_on and not era5_on:
            n_expected_weather_files = 4
        else:
            n_expected_weather_files = 5
        for folder in simulation_folders:
            if len(os.listdir(os.path.join(simulations, folder))) < n_expected_weather_files:
                days_to_reelaborate.append(datetime.datetime.strptime(folder, "%Y%m%d"))
        if len(days_to_reelaborate) == 0:
            with open("log_weather.txt", "a", encoding="utf-8", errors="surrogateescape") as log_file:
                log_file.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                               ' There are no days to re-elaborate\n')
            break
        else:
            with open("log_weather.txt", "a", encoding="utf-8", errors="surrogateescape") as log_file:
                log_file.write(datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d-%H:%M:%S") +
                               ' Re-elaboratoring the following days\n')
                for day_to_reelaborate in days_to_reelaborate:
                    log_file.write(datetime.datetime.strftime(day_to_reelaborate, "%Y-%m-%d %H:%M:%S") + '\n')
            n_elaborated_days = 0
            pools = []
            n_pool = 0
            while n_elaborated_days <= len(days_to_reelaborate):
                start = n_elaborated_days
                end = n_elaborated_days + max_number_processes
                if end > len(days_to_reelaborate):
                    end = len(days_to_reelaborate)
                n_elaborated_days = end
                pools.append(n_pool)
                n_pool += 1
                if n_elaborated_days == len(days_to_reelaborate):
                    break
            n_elaborated_days = 0
            n_pool = 0
            while n_elaborated_days <= len(days_to_reelaborate):
                start = n_elaborated_days
                end = n_elaborated_days + max_number_processes
                if end > len(days_to_reelaborate):
                    end = len(days_to_reelaborate)
                try:
                    pools[n_pool] = ThreadingPool(max_number_processes)
                    pools[n_pool].map(automatic_weather, days_to_reelaborate[start:end])
                except BaseException:
                    print("Unable to process weather data")
                    sys.exit()
                n_elaborated_days = end
                n_pool += 1
                if n_elaborated_days == len(days_to_reelaborate):
                    break
            attempt += 1
