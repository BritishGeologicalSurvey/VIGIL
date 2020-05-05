from random import sample,randrange
import os
import shutil
from pathos.multiprocessing import Pool, ThreadingPool
import multiprocessing

def pre_process():
    def random_sources(n_fumaroles):
        import numpy as np
        from matplotlib import pyplot as plt

        #############################################
        # CREAZIONE MATRICE PESI
        #############################################
        n = 98
        heatmap_data = np.ones([n, n])  # valore di default 1 (np.ones crea una matrice con tutti valori 1)

        # present-day ACTIVE FUMAROLES
        heatmap_data[48:50, 46:50] += 5

        # FRACTURES
        # settore NS
        heatmap_data[39:50, 46:47] += 3  # righe da 40 a 47 e colonna 46: assegno valore 3
        heatmap_data[41:42, 44:46] += 4
        # settore E
        heatmap_data[50:51, 47:54] += 3
        heatmap_data[51:52, 52:54] += 3
        # settore NO fract Faujas
        heatmap_data[42:44, 41:43] += 3
        heatmap_data[42:43, 39:41] += 4
        heatmap_data[43:44, 38:40] += 4
        heatmap_data[44:47, 37:38] += 4
        # settore NE fract du NordEst
        heatmap_data[42:44, 49:50] += 3
        heatmap_data[40:42, 50:51] += 3
        # settore S past activity
        heatmap_data[65:66, 46:48] += 4
        heatmap_data[69:71, 47:48] += 4
        # Ty fault
        heatmap_data[62:64, 50:52] += 3
        heatmap_data[64:68, 51:53] += 3
        heatmap_data[60:62, 49:51] += 3
        heatmap_data[68:70, 52:54] += 3

        # settore SE fract del crater sud ???
        heatmap_data[51:57, 48:50] += 3
        # settore NO fract du NordOvest
        heatmap_data[44:46, 44:45] += 3
        heatmap_data[45:47, 45:46] += 3

        # PAST fumarolic activity
        # Nord 1976-77
        heatmap_data[37:38, 45:48] += 4
        # NE 1976-77
        heatmap_data[44:47, 55:56] += 4
        # Est 1976-77
        heatmap_data[50:53, 55:59] += 4
        heatmap_data[55:58, 54:57] += 4
        # SSE 1976-77 Morue Mitan
        heatmap_data[60:62, 50:55] += 4
        # SE 1956
        heatmap_data[53:54, 50:51] += 4
        # Sud 2017
        heatmap_data[61:62, 49:50] += 4

        # inside DOME AREA
        y, x = np.ogrid[-49:n - 49, -45:n - 45]
        mask = x * x + y * y <= 11 * 11
        heatmap_data[mask] += 2

        # area a peso 0 aldisopra delle strutture di collasso calderico passato
        bordo = [60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 58, 57, 54, 51, 51, 49, 47, 45, 41, 39, 39, 36,
                 34, 34, 33, 33, 31, 31, 28, 28, 27, 27, 27, 26, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
                 28, 28, 31, 30, 27, 27, 27, 27, 27, 27, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 31, 31, 31, 31, 31,
                 33, 33, 33, 33, 35, 35, 35, 35, 35, 35, 35, 40, 40, 40, 40, 43, 43, 43, 45, 45, 45, 45, 45, 45]
        for c, r in enumerate(bordo):
            heatmap_data[:r, c] = 0

        ##############################
        # Calcolo probabilità best guess
        ##############################

        heatmap_data = heatmap_data / heatmap_data.sum()
        np.savetxt('result_probabilities.txt', heatmap_data, fmt='%.2e')  # save a text file of the probabilities

        #
        # --The x and y coordinates from the grd file
        #

        xmin = 641525
        xmax = 645475
        x = np.linspace(xmin, xmax, n)

        ymin = 1772530
        ymax = 1776480
        y = np.linspace(ymin, ymax, n)

        ######################
        # Randomly selecting probabilities and their rows and columns from the heatmap matrix in order to take at
        # the same row and column the x and y
        ######################

        file = open('source.txt', 'w')

        random_index = randrange(len(heatmap_data))
        for i in range(0, n_fumaroles-1):
            row = randrange(0, 97)
            column = randrange(0, 97)
            print(row, column)
            probability = heatmap_data[row, column]
            xpr = x[row]
            ypr = y[column]
            print("at x", xpr, "at y", ypr, "Randomly selected probability - ", probability)
            file.write(str(xpr))
            file.write("|")
            file.write(str(ypr))
            file.write("|")
            file.write(str(probability))
            file.write('\n')
        file.close()

    def fluxes():
        import numpy as np
        from matplotlib import pyplot as plt
        import pandas as pd

        data = pd.read_csv('flux.csv', error_bad_lines=False)
        x = np.sort(data['flux'])
        y = np.arange(1, len(x) + 1) / len(x)
        # plt.plot (x, y, marker= '.', markersize=3, linestyle='-', color='orange')
        # plt.xlabel('H2O flux (kg/s)')
        # plt.ylabel ('ECDF')
        # plt.margins(0.02)
        # plt.show()

        # Sample randomly 100 data
        with open('random_list_flux.txt', 'w') as outfile:
            list_x = list(x)
            sampled_flux = (sample(list_x, 1))
            # strFormat = len(samples) * '{:10f} '
            # formattedList = strFormat.format(*samples)
        return sampled_flux
    #read and memorize disgas.inp file
    disgas_input_records = []
    with open('disgas.inp','r',encoding="utf-8", errors="surrogateescape") as disgas_or_input:
        for line in disgas_or_input:
            disgas_input_records.append(line)

    Nf=[2,3,4,5] # fumaroles number vector to be sampled
    raw_days = [] # store the days as originally formatted
    days = [] #store days in format YYYYMMDD as per folder name
    # read days_list file
    with open('days_list.txt','r',encoding="utf-8", errors="surrogateescape") as days_list_file:
        for line in days_list_file:
            raw_days.append(line)

    i=0
    for day in raw_days:
        temp = raw_days[i].split(' ')
        temp = temp[0].split('-')
        days.append(temp[0]+temp[1]+temp[2])
        i+=1

    for day in days:
        path = os.path.join(os.getcwd(),'simulations',str(day))  # To modify accordingly
        rawdata = os.path.join(path,'raw_data')
        infiles = os.path.join(path, 'infiles')
        outfiles = os.path.join(path, 'outfiles')
        presfc = os.path.join(rawdata, 'presfc.dat')
        preupr = os.path.join(rawdata, 'preupr.dat')
        diagno = os.path.join(rawdata, 'diagno.inp')
        disgas_input = os.path.join(path,'infiles','disgas.inp')
        try:
            os.mkdir(rawdata)
        except:
            print('Folder raw_data already exists in '+str(path))
        files = os.listdir(path)
        for f in files:
            path_f = os.path.join(path,f)
            if f != 'raw_data' and f != 'infiles' and f != 'outfiles':
                shutil.move(path_f,rawdata)
                #shutil.copy(path_f, raw_data)
        try:
            os.mkdir(infiles)
        except:
            print('Folder infiles already exists in ' + str(path))
        try:
            os.mkdir(outfiles)
        except:
            print('Folder outfiles already exists in ' + str(path))
        n_fumaroles = sample(Nf,1)[0]
        random_sources(n_fumaroles)
        records=[]
        easting=[]
        northing=[]
        with open('source.txt','r',encoding="utf-8", errors="surrogateescape") as locations_file:
            i=0
            for line in locations_file:
                records.append(line.split('|'))
                easting.append(records[i][0])
                northing.append(records[i][1])
                i+=1
        with open('source.dat','w', encoding="utf-8", errors="surrogateescape") as source_file:
            for i in range(0,n_fumaroles-1):
                h2o_flux = fluxes()
                source_file.write(str(easting[i])+' '+str(northing[i]) + ' 0. ' + str(h2o_flux[0]) + '\n')
        try:
            shutil.move('source.txt',os.path.join(rawdata,'source.txt'))
            shutil.move('source.dat',os.path.join(infiles,'source.dat'))
            shutil.move(presfc,os.path.join(infiles,'presfc.dat'))
            shutil.move(preupr, os.path.join(infiles, 'preupr.dat'))
            shutil.move(diagno, os.path.join(infiles, 'diagno.inp'))
        except:
             print('Files already there')
        with open(disgas_input,'w', encoding="utf-8", errors="surrogateescape") as disgas:
            for i in range(0,len(disgas_input_records)):
                if i == 3:
                    disgas.write('  YEAR   = ' + day[0:4] + '\n')
                    #disgas_input[3] = '  YEAR   = ' + day[0:4]
                elif i == 4:
                    disgas.write('  MONTH  = ' + day[4:6] + '\n')
                    #disgas_input[4] = '  MONTH  = ' + day[4:6]
                elif i == 5:
                    disgas.write('  DAY    = ' + day[6:8] + '\n')
                    #disgas_input[5] = '  DAY    = ' + day[6:8]
                else:
                    disgas.write(disgas_input_records[i])
        shutil.copy('topography.grd',os.path.join(infiles,'topography.grd'))
    #Move all infiles into the working folder
    for day in days:
        path = os.path.join(local_path, 'simulations', str(day))
        infiles_local = os.path.join(path, 'infiles')
        infiles_scratch = os.path.join(scratch_path, 'simulations', str(day), 'infiles/')
        outfiles = os.path.join(scratch_path, 'simulations', str(day), 'outfiles/')
        shutil.move(infiles_local,infiles_scratch)
        os.makedirs(outfiles)
        os.chdir(infiles_scratch)
        disgas_input_records = []
        with open('disgas.inp', 'r', encoding="utf-8", errors="surrogateescape") as disgas_or_input:
            for line in disgas_or_input:
                disgas_input_records.append(line)
        with open('disgas.inp', 'w', encoding="utf-8", errors="surrogateescape") as disgas:
            for i in range(0, len(disgas_input_records)):
                if i == 45:
                    disgas.write('   TOPOGRAPHY_FILE_PATH   = ' + infiles_scratch + 'topography.grd ' + '\n')
                elif i == 46:
                    disgas.write('   ROUGHNESS_FILE_PATH   = ' + infiles_scratch + 'roughness.grd ' + '\n')
                elif i == 47:
                    disgas.write('   RESTART_FILE_PATH   = ' + infiles_scratch + 'restart.dat ' + '\n')
                elif i == 48:
                    disgas.write('   SOURCE_FILE_PATH   = ' + infiles_scratch + 'source.dat ' + '\n')
                elif i == 49:
                    disgas.write('   WIND_FILE_PATH   = ' + infiles_scratch + 'winds.dat ' + '\n')
                elif i == 50:
                    disgas.write('   DIAGNO_FILE_PATH   = ' + infiles_scratch + 'diagno.out ' + '\n')
                elif i == 51:
                    disgas.write('   OUTPUT_DIRECTORY    = ' + outfiles + '\n')
                else:
                    disgas.write(disgas_input_records[i])

    return days

def run_model(day):
    infiles = os.path.join(scratch_path, 'simulations', str(day),'infiles')  # To modify accordingly
    os.chdir(infiles)
    os.system('presfc && preupr && diagno && disgas disgas.inp &')

scratch_path = '/scratch/silviamassaro'
local_path = os.getcwd()
days = pre_process()

path = os.getcwd()
for day in days:
    run_model(day)



