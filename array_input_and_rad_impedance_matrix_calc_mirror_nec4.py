# Python3
Software_version = '2021.02.03'
Software_name = 'Array Input and Radiation Impedance Matrix Calculator for NEC4 results'

"""
            array_input_and_rad_impedance_matrix_calc
   ***  Antenna Array Input and Radiation Impedance Matrices Calculator   ***
Program reads the result files of NEC modeling (.nec) where each (of non mirrored) dipole
of the antenna array is excited one by one (in each file), finds frequencies of analysis,
excitation sources, loads and their values, currents in sources and loads, antenna patterns.
Then the program calculates mutual impedances, radiation resistances, and store them to
text files (.txt).

The initial version was created for AA of 25 GURT dipoles of 1 polarization of
incoming waves, number of segments per dipole is 325 (Multi frequencies) then it was
extended for other possibilities, but still you need to analyze carefully the results
* no_of_wires_per_dipole
    98 for UTR-2 6x5 array 574 segments
    182  for UTR-2 6x5 array 1018 segments
"""

# *************************************************************
#                         PARAMETERS                          *
# *************************************************************
path_to_data = 'DATA/'

no_of_wires_per_dipole = 182      # Number of wires per dipole (not segments!) needed to find loads  # 434
array_input_num = 30              # Array inputs number (total number of dipoles in array or inputs of dipoles)
num_of_freq = 25                  # Maximal possible number of frequencies analyzed
print_or_not = 1

# Rectangular UTR-2 antenna array of 6 * 5 = 30 dipoles
no_of_dip = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 19, 20, 21, 25, 26]  # dipoles being excited
mirror_from = [ 1,  2,  3,  4,  7,  8,  9, 13, 14, 19]	       # Active dipoles under the diagonal
mirror_into = [30, 29, 28, 27, 24, 23, 22, 18, 17, 12]         # Passive dipoles above the diagonal

# Linear antenna array 1 * 6 = 6 dipoles
# no_of_dip = [1,  2,  3]          # 3 dipoles being excited
# mirror_from = [1, 2, 3]	       # Active dipoles under the diagonal
# mirror_into = [6, 5, 4]         # Passive dipoles above the diagonal

# Square antenna array 5 * 5 = 25 dipoles
# no_of_dip =    [1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 12, 13, 16, 17, 21]  # 15 dipoles being excited
# mirror_from = [1,  2,  3,  4,  6,  7,  8,  11, 12, 16]	                  # Active dipoles under the diagonal
# mirror_into = [25, 24, 23, 22, 20, 19, 18, 15, 14, 10]                      # Passive dipoles above the diagonal

# Linear antenna array 1 * 25 dipoles
# no_of_dip =    [1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13]  # Dipols being excited (numbers in NEC out files names)
# mirror_from = [1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12]	   # Active dipoles under the diagonal
# mirror_into = [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14]       # Passive dipoles above the diagonal

pi = 3.141593
# *************************************************************
#                    IMPORT LIBRARIES                         *
# *************************************************************
import os
import sys
import time
# from scipy.integrate import simps
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from package_common_modules.find_line_in_text_file import find_line_in_text_file
from package_common_modules.intergration_trap_2d import intergration_trap_2d
from package_common_modules.figure_color_map import figure_color_map


# *************************************************************
#                       MAIN PROGRAM                          *
# *************************************************************

print('\n\n\n\n\n    ---------------------------------------------------------------------------------')
print('    ***                                                                          ***')
print('    *** ', Software_name, ' ***')
print('    ***                            v.', Software_version, '                                ***')
print('    ***       Program analyzes NEC output files with it turn excitation of       ***')
print('    ***         array elements (one source other terminals are loaded),          ***')
print('    ***             reads and calculates self and mutual impedances,             ***')
print('    ***                radiation resistances of array elements.                  ***')
print('    ***               Works with frequency sweep results of NEC4                 ***')
print('    ***                              (c) YeS 2020                                ***')
print('    --------------------------------------------------------------------------------\n\n\n')

# *** Time consumption calculation (beginning) ***
start_time = time.time()
current_time = time.strftime("%H:%M:%S")
current_date = time.strftime("%d.%m.%Y")
print('  Today is ', current_date, ' time is ', current_time, '\n\n')

nec_file_num = len(no_of_dip)    # The number of NEC out files

Rinp = np.full((num_of_freq, array_input_num, array_input_num), 0.0)
Xinp = np.full((num_of_freq, array_input_num, array_input_num), 0.0)
IFedDipReal = np.full((num_of_freq, array_input_num, array_input_num), 0.0)
IFedDipImag = np.full((num_of_freq, array_input_num, array_input_num), 0.0) # IFedDip(frequency step, Active dipole, Passive dipole)
DistanceRP = np.full((nec_file_num, num_of_freq), 0.0)    # [file_num, FreqStep] Distance for which NEC calculates radiation pattern

Points = np.full((91, 361), 0.0 + 1j*0.0)     # [theta, phi] Matrix of E-field data for integral calculations
RSigm = np.full((num_of_freq, array_input_num, array_input_num), 0.0+1j*0.0)    #[FreqStep, Dip1, Dip2] Radiation impedance matrix
Efficiency = np.zeros((num_of_freq, array_input_num))

ETHcmplx = np.full((num_of_freq, array_input_num, 181, 361), 0.0+1j*0.0)  # [FreqStep, DipNo, theta, phi] Theta components of E-field
EPHcmplx = np.full((num_of_freq, array_input_num, 181, 361), 0.0+1j*0.0)  # [FreqStep, DipNo, theta, phi] Phi components of E-field

# *** Creating folder for results ***
result_path = 'Results'
if not os.path.exists(result_path):
    os.makedirs(result_path)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!                 Reading data from NEC output files                       !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for file_num in range(len(no_of_dip)):  # Main loop by NEC output files

    # *** Configuring the name of NEC output file ***

    # if file_num < 9:
    #     no = '0' + str(no_of_dip[file_num])
    # else:
    #     no = str(no_of_dip[file_num])

    no = str(no_of_dip[file_num])

    # file_name = path_to_data + 'UTR2_6x5-' + no + '.out'
    file_name = path_to_data + 'UTR2_SEG=1018_6x5_Ex' + no + '.out'  # UTR2_SEG=1018_6x5_Ex1
    # file_name = path_to_data + 'UTR2_6x5-' + str(file_num + 1) + '.out'
    # file_name = path_to_data + 'UTR2_6x1_EX=' + str(file_num + 1) + '.out'
    # file_name = path_to_data + 'GURT-V325_25x01_Nex=' + no + '.out.txt'
    # file_name = path_to_data + 'NEC_in_'+str(file_num+1)+'.out.txt'

    # *** Opening NECoutput file ***
    data_file = open(file_name, "r")

    print('\n **********************************************************************')
    print('\n          Reading file: ', file_name, ' \n')

    # *** Counting number of lines in the file ***
    # num_lines = sum(1 for line in data_file)
    # print ('\n  Number of lines in the file = ', num_lines, ' ')
    # data_file.seek(0)

    #  Searching the lines where loads are listed
    wires_with_loads = []
    segments_with_loads = []
    loads_values = []
    line = find_line_in_text_file(data_file, '- - - STRUCTURE IMPEDANCE LOADING - - -', 31)
    for i in range(5): 
        line = data_file.readline()  # Skip several lines
    words_in_line = 'A'  # Dummy value to make a string variable
    counter = 0
    while words_in_line[0] != 'CP' and counter < 2:
        line = data_file.readline()
        if line != '\n':
            counter = 0
            words_in_line = line.split()
            if words_in_line[4] + ' ' + words_in_line[5] == 'FIXED IMPEDANCE':
                wires_with_loads.append(int(words_in_line[0]))
                segments_with_loads.append(int(words_in_line[1]))
                loads_values.append(float(words_in_line[3]))
        else:
            counter += 1

    LoadNum = len(segments_with_loads)  # !!!!!!!!!!!! Possible error !!!!!!!!!!!!!!

    if print_or_not == 1:
        for i in range(LoadNum):
            print('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], 
                  '  Resistance = ', loads_values[i])

    # Searching number of frequencies
    data_file.seek(0)
    frequency_list = []
    num_of_frequencies = 0
    for line in data_file:
        if line.startswith('                                 - - - - - - FREQUENCY - - - - - -'):
            num_of_frequencies += 1

    print('\n  Number of frequencies analyzed = ', num_of_frequencies, '\n')

    print(' **********************************************************************')

    data_file.seek(0)  # return to the beginning of the file

    for step in range(num_of_frequencies):
        # Seek for current frequency
        line = find_line_in_text_file (data_file, 'FREQUENCY=', 36)
        words_in_line = line.split()
        frequency_list.append(float(words_in_line[1]))

        # *** Searching antenna input parameters block at the frequency ***
        line = find_line_in_text_file(data_file, '- - - ANTENNA INPUT PARAMETERS - - -', 42)

        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = find_line_in_text_file(data_file, 'IMPEDANCE', 64)
        line = data_file.readline()   # Skip one line

        line = data_file.readline()
        IFedDipReal[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1] = float(line[36 : 48])
        IFedDipImag[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1] = float(line[48 : 60])
        Rinp[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1] = float(line[60 : 72])
        Xinp[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1] = float(line[72 : 84])

        # print ('\n\n   For frequency = ', frequency_list[step], ' MHz')
        # print ('  ---------------------------- \n')
        # print ('  Self full impedance of the dipole = ', Rinp[step, no_of_dip[file_num]-1,
        # no_of_dip[file_num]-1], Xinp[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1], ' Ohm' )
        # print ('  Current in fed point of active dipole = ', IFedDipReal[step, no_of_dip[file_num]-1,
        # no_of_dip[file_num]-1],  IFedDipImag[step, no_of_dip[file_num]-1, no_of_dip[file_num]-1],' Amp. \n')

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !! Searching and reading the current in impedance load !!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # *** Seek for currents in loads ***
        line = find_line_in_text_file(data_file, '- - - CURRENTS AND LOCATION - - -', 29)
        for i in range(5): 
            line = data_file.readline()   # Skip several lines

        for j in range(LoadNum):  # loop by number of loads
            intTemp = 0
            while intTemp < wires_with_loads[j]:   # !=
                line = data_file.readline()
                intTemp = int(line[6:12])
            if segments_with_loads[j] < 1:  # if the number of segment = 1 we return to the beginning of the line
                sys.exit('ERROR!')
            if segments_with_loads[j] > 1:  # if the number of segment > 2 we skip lines
                for i in range (segments_with_loads[j] - 1): line = data_file.readline()   # Skip several lines

            DipNum = int((wires_with_loads[j]) / no_of_wires_per_dipole) + 1   # Calculation of the dipole number   !!!!
            # DipNum = int((wires_with_loads[j]) / wires_with_loads[0]) + 1   # Calculation of the dipole number   !!!!

            # *** Reading currents ***
            IFedDipReal[step, no_of_dip[file_num]-1, DipNum-1] = float(line[49: 60])
            IFedDipImag[step, no_of_dip[file_num]-1, DipNum-1] = float(line[61: 72])

            # print ('  Dip No', DipNum, ' Wire No', wires_with_loads[j], ' Segm No', segments_with_loads[j])
            # print ('  Current in the load No ', j+1, ' = ', IFedDipReal[step, no_of_dip[file_num]-1, DipNum-1],
            #        IFedDipImag[step, no_of_dip[file_num]-1, DipNum-1], ' Amp' )

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!  Reading of E values in radiation patterns  !!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        line = find_line_in_text_file(data_file, '- - - RADIATION PATTERNS - - -', 48)
        line = data_file.readline()
        DistanceRP[file_num, step] = float(line[60:74])
        # print ('\n  Distance of RP calculations = ', DistanceRP[file_num, step], ' m. \n')

        for i in range(7):
            line = data_file.readline()   # Skip several lines of table header

        # *** loop by lines to read e-field amplitudes and phases for all theta and phi ***

        for k in range(65341):  # 65341 = 181 * 361 degrees
            line = data_file.readline()
            theta1 = float(line[0: 8])
            phi1 = float(line[9: 17])
            Etheta = float(line[76: 87])
            ang1 = float(line[89: 97])
            Ephi = float(line[100: 112])
            ang2 = float(line[112: 120])

            thetaInt = (int(theta1)) + 90
            phiInt = int(phi1)
            ETHcmplx[step, no_of_dip[file_num]-1, thetaInt, phiInt] = \
                complex((Etheta * np.cos(ang1 * pi/180)), (Etheta * np.sin((ang1 * pi/180))))
            EPHcmplx[step, no_of_dip[file_num]-1, thetaInt, phiInt] = \
                complex((Ephi * np.cos(ang2 * pi/180)), (Ephi * np.sin((ang2 * pi/180))))

    data_file.close()

IFedDip = IFedDipReal + 1j * IFedDipImag
Zinp = Rinp + 1j * Xinp
del IFedDipReal, IFedDipImag, Rinp, Xinp


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!              Processing data and calculating results                     !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n    Calculations started at: ', time.strftime("%H:%M:%S"), '  \n\n')


for step in range (num_of_frequencies):     # Loop by frequencies

    # *** Mirroring of currents for cases where dipoles were not excited ***

    for i in range(len(mirror_into)):
        for j in range(array_input_num):
            IFedDip[step, mirror_into[i]-1, j] = IFedDip[step, mirror_from[i]-1, array_input_num-1-j] 

    # *** Calculations of Impedance matrices ***

    for i in range(nec_file_num):           # Loop by files
        for k in range(array_input_num):    # Loop by dipoles
            if k != (no_of_dip[i]-1):
                Zinp[step, no_of_dip[i]-1, k] = - (IFedDip[step, no_of_dip[i]-1, k] * loads_values[i] /
                                                   IFedDip[step, no_of_dip[i]-1, no_of_dip[i]-1])  # !!!-!!!

    # Mirroring the Impedances for dipoles that were not excited

    for i in range(len(mirror_into)):     # Loop by active dipoles below array diagonal
        for j in range(array_input_num):         # Loop by dipoles
            Zinp[step, mirror_into[i]-1, array_input_num - 1 - j] = Zinp[step, mirror_from[i]-1, j] 

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!  Mirroring of RP for noncalculated dipoles  !!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for i in range(len(mirror_into)):  # loop by dipoles
        for t in range(181):           # loop by theta
            for p in range(180):       # loop by phi
                ETHcmplx[step, mirror_into[i]-1, t, p+180] = - ETHcmplx[step, mirror_from[i]-1, t, p]
                EPHcmplx[step, mirror_into[i]-1, t, p+180] = - EPHcmplx[step, mirror_from[i]-1, t, p]
            for p in range(180, 361):
                ETHcmplx[step, mirror_into[i]-1, t, p-180] = - ETHcmplx[step, mirror_from[i]-1, t, p]
                EPHcmplx[step, mirror_into[i]-1, t, p-180] = - EPHcmplx[step, mirror_from[i]-1, t, p]

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!         Integrals calculation              !!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for i in range(array_input_num):       # Loop by the first dipole
        for j in range(array_input_num):   # Loop by the second dipole
            if i >= j:

                # Forming the integration elements for Theta angle
                for t in range(91):
                    for p in range(361):
                        Points[t, p] = ETHcmplx[step, i, t+90, p] * np.conj(ETHcmplx[step, j, t+90, p]) * \
                                       (np.sin(t * pi/180.) + 1j*0.0)

                # y = np.linspace(1, 1, 91)
                # x = np.linspace(1, 1, 361)
                # Integral1 = simps(simps(Points, x), y)
                Integral1 = intergration_trap_2d(0, 90, 0, 360, 90, 360, Points)

                # Forming the integration elements for Phi component
                for t in range(91):
                    for p in range(361):
                        Points[t, p] = EPHcmplx[step, i, t+90, p] * np.conj(EPHcmplx[step, j, t+90, p]) * \
                                       (np.sin(t * pi/180.) + 1j*0.0)

                # Calculating of integrals by Phi angle

                # y = np.linspace(1, 1, 91)
                # x = np.linspace(1, 1, 361)
                # Integral2 = simps(simps(Points, x), y)
                Integral2 = intergration_trap_2d(0, 90, 0, 360, 90, 360, Points)

                # Calculating the elements of radiation impedance matrix !!! Distance[1] !!!

                RSigm[step, i, j] = ((Integral1 + Integral2) * np.power(DistanceRP[0, 0], 2)) / \
                                    (120.0 * pi * (IFedDip[step, i, i] * np.conj(IFedDip[step, j, j])))
                RSigm[step, j, i] = np.conj(RSigm[step, i, j])

    # *** Calculations of self radiation efficiencies of dipols ***
    for i in range(array_input_num):
        Efficiency[step, i] = (np.real(RSigm[step, i, i]) / np.real(Zinp[step, i, i])) * 100

    current_time = time.strftime("%H:%M:%S")
    print('\n    Done for frequency = ', frequency_list[step], ' MHz at: ', current_time, ' ')

# End of loop by frequencies


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!             Displaying data and results on the screen                    !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


for step in range(num_of_frequencies):	 # Loop by frequencies for displaying results on the screen

    print('\n\n  For frequency = ', frequency_list[step], ' MHz')
    print('  --------------------------')
    print('    Zinp (1, 1) = ', Zinp[step, 0, 0])
    print('   Rsigm (1, 1) = ', RSigm[step, 0, 0])
    print('\n  Self radiation efficiency of the first dipole = ', round(Efficiency[step, 0],3), ' %')   # i -> 0


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!              Writing data and results to output files                    !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print('\n\n\n     Writing the results to output files... \n\n\n')


# *** Creating files for results ***
CurrentMatrixFile = open(result_path + "/Matrices of Currents.txt", "w")              # 24
ImpedanceMatrixFile = open(result_path + "/Matrices of Impedances.txt", "w")          # 26
RadImpMatrixFile = open(result_path + "/Matrices of Radiation Impedances.txt", "w")   # 27

for step in range(num_of_frequencies):	 # Loop by frequencies for results files writing

    # *** Matrix of complex currents in dipoles.txt ***
    CurrentMatrixFile.write('\n\n    *** Matrix of currents in fed points of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
    data = ['                                  {:2d}'.format(i+2) for i in range(array_input_num-1)]
    CurrentMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(array_input_num-1)]
    CurrentMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(IFedDip[step, k, i]), np.imag(IFedDip[step, k, i])) for i in range(array_input_num)]
        CurrentMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')

    # *** Writing impedances to file ***
    ImpedanceMatrixFile.write('\n\n    *** Matrix of impedances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
    data = ['                                  {:2d}'.format(i+2) for i in range(array_input_num-1)]
    ImpedanceMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(array_input_num-1)]
    ImpedanceMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(Zinp[step, k, i]), np.imag(Zinp[step, k, i])) for i in range(array_input_num)]
        ImpedanceMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')

    # *** Writing radiation resistances to file ***
    RadImpMatrixFile.write('\n\n    *** Matrix of radiation resistances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
    data = ['                                  {:2d}'.format(i+2) for i in range(array_input_num-1)]
    RadImpMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(array_input_num-1)]
    RadImpMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(RSigm[step, k, i]), np.imag(RSigm[step, k, i])) for i in range(array_input_num)]
        RadImpMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')

    # *** Writing currents on dipols inputs to files for MatCAD
    file_name = (result_path + '/Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()

    # *** Writing normalized currents on dipols inputs to files for MatCAD ***
    file_name = (result_path + '/Normalized Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()

    # *** Writing impedances of dipols to files for MatCAD ***
    file_name = (result_path + '/Impedances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.real(Zinp[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data) + '  \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.imag(Zinp[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data) + '  \n')
    CurrentFile.close()


    # *** Writing radiation resistances of dipols to files for MatCAD ***
    file_name = (result_path + '/Radiation resistances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.real(RSigm[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data) + '  \n')
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(np.imag(RSigm[step, i, k])) for i in range(array_input_num)]
        CurrentFile.write(' ' + ''.join(data) + '  \n')
    CurrentFile.close()

    # *** Writing radiation efficiencies of dipols to files for MatCAD ***
    file_name = (result_path + '/Radiation efficiencies of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        data = ['   {:+14.8e}'.format(Efficiency[step, i])]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()

    # *** Writing normalized Etheta RP to files for MatCAD ***
    file_name = (result_path + '/Normalized to current Etheta RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        for t in range(91):
            data = [' {:+14.8e}'.format(np.real(ETHcmplx[step, k, t+90, p] / IFedDip[step, k, k])) for p in range(361)]
            CurrentFile.write('' + ''.join(data) + ' \n')
    for k in range(array_input_num):
        for t in range(91):
            data = [' {:+14.8e}'.format(np.imag(ETHcmplx[step, k, t+90, p] / IFedDip[step, k, k])) for p in range(361)]
            CurrentFile.write('' + ''.join(data) + ' \n')
    CurrentFile.close()

    # *** Writing normalized Ephi RP to files for MatCAD ***
    file_name = (result_path + '/Normalized to current Ephi RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
    CurrentFile = open(file_name, "w")  # 6
    for k in range(array_input_num):
        for t in range(91):
            data = [' {:+14.8e}'.format(np.real(EPHcmplx[step, k, t+90, p] / IFedDip[step, k, k])) for p in range(361)]
            CurrentFile.write('' + ''.join(data) + ' \n')
    for k in range(array_input_num):
        for t in range(91):
            data = [' {:+14.8e}'.format(np.imag(EPHcmplx[step, k, t+90, p] / IFedDip[step, k, k])) for p in range(361)]
            CurrentFile.write('' + ''.join(data) + ' \n')
    CurrentFile.close()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!                            F I G U R E S                                 !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('     Making figures... \n\n')

# *** Figure of current absolute value ***
for step in range(len(frequency_list)):

    data = np.zeros((array_input_num))
    for i in range(array_input_num):
        data[i] = np.absolute(IFedDip[step, i, i])
    x_values = np.linspace(1, array_input_num, num=array_input_num)
    plt.figure()
    rc('font', weight='normal')
    plt.plot(x_values, data[:], 'bo', label='Abs(IFedDip)')
    plt.xlabel('Number of dipole')
    plt.ylabel('Current at dipole output, A')
    # plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
    plt.title('NEC modeling results, self currents - absolute values', fontsize=10)
    plt.grid(b=True, which='both', color='silver', linestyle='-')
    plt.legend(loc='lower center', fontsize=10)
    pylab.savefig(result_path + "/Dipole currents (abs) f = %5.2f MHz.png" % frequency_list[step],
                  bbox_inches='tight', dpi=160)
    plt.close('all')


# *** Figure of current absolute value matrix nondiagonal map ***
for step in range(len(frequency_list)):
    data = np.absolute(IFedDip[step, :, :])
    for i in range(array_input_num):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole','Number of dipole',
                     'Matrix of currents at inputs of antenna array at %5.2f MHz' % frequency_list[step],
                     'NEC modeling results',
                     result_path + "/Dipole currents matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step])


# *** Figure of input impedance absolute value at the main diagonal ***
for step in range(len(frequency_list)):
    data = np.zeros((array_input_num))
    for i in range(array_input_num):
        data[i] = np.absolute(Zinp[step, i, i])
    x_values = np.linspace(1, array_input_num, num=array_input_num)
    plt.figure()
    rc('font', weight='normal')
    plt.bar(x_values, data[:], bottom=np.min(data[:]) * 0.99)
    plt.xlabel('Number of dipole')
    plt.ylabel('|Z|, Ohm')
    plt.title('NEC modeling results', fontsize=7, x=0.46, y=1.005)
    plt.grid(b=True, which='both', color='0.00', linestyle='--')
    # plt.legend(loc = 'lower center', fontsize = 10)
    pylab.savefig(result_path + "/Self full impedances (abs) f = %5.2f MHz.png" % frequency_list[step],
                  bbox_inches='tight', dpi=160)
    plt.close('all')


# *** Figure of impedance absolute value matrix nondiagonal map ***
for step in range(len(frequency_list)):
    data = np.absolute(Zinp[step, :, :])
    for i in range(array_input_num):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole', 'Number of dipole',
                     'Matrix of impedances at inputs of antenna array at %5.2f MHz' % frequency_list[step],
                     'NEC modeling results, self values set to zeros to show mutual pattern',
                     result_path + "/Dipole impedances matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step])


# *** Figure of radiation resistance absolute value matrix nondiagonal map ***
for step in range(len(frequency_list)):
    data = np.absolute(RSigm[step, :, :])
    for i in range(array_input_num):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole', 'Number of dipole',
                     'Matrix of mutual radiation resistances of dipoles in antenna array at %5.2f MHz' % frequency_list[step],
                     'NEC modeling results, self values set to zeros to show mutual pattern',
                     result_path + "/Radiation resistance matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step])

central = int(array_input_num / 2) + 1

# Frequency dependence of input impedance of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list[1: len(frequency_list)], np.real(Zinp[0: len(frequency_list)-1, 0, 0]),
         color='C0', linestyle='-', linewidth='1.00', label='Re(Zinp) dipole # 1')
plt.plot(frequency_list[1: len(frequency_list)], np.imag(Zinp[0: len(frequency_list)-1, 0, 0]),
         color='C1', linestyle='-', linewidth='1.00', label='Im(Zinp) dipole # 1')
plt.plot(frequency_list[1: len(frequency_list)], np.absolute(Zinp[0: len(frequency_list)-1, 0, 0]),
         color='C4', linestyle='-', linewidth='1.00', label='Abs(Zinp) dipole # 1')
plt.plot(frequency_list[1: len(frequency_list)], np.real(Zinp[0: len(frequency_list)-1, central, central]),
         color='C0', linestyle='--', linewidth='1.00', label='Re(Zinp) central dipole # ' + str(central))
plt.plot(frequency_list[1: len(frequency_list)], np.imag(Zinp[0: len(frequency_list)-1, central, central]),
         color='C1', linestyle='--', linewidth='1.00', label='Im(Zinp) central dipole # ' + str(central))
plt.plot(frequency_list[1: len(frequency_list)], np.absolute(Zinp[0: len(frequency_list)-1, central, central]),
         color='C4', linestyle='--', linewidth='1.00', label='Abs(Zinp) central dipole # ' + str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel('R, X, |Z|, Ohm')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results', fontsize=7, x=0.46, y=1.005)
plt.grid(b=True, which='both', color='silver', linestyle='-')
plt.legend(loc='lower center', fontsize=8)
pylab.savefig(result_path + "/Self full impedances of the first dipole.png", bbox_inches='tight', dpi=160)
plt.close('all')

# Frequency dependence of radiation resistance of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list[1: len(frequency_list)], np.real(RSigm[0: len(frequency_list)-1, 0, 0]),
         color='C0', linestyle='-', linewidth='1.00', label='Rsigm dipole # 1')
plt.plot(frequency_list[1: len(frequency_list)], np.real(RSigm[0: len(frequency_list)-1, central, central]),
         color='C1', linestyle='-', linewidth='1.00', label='Rsigm central dipole # ' + str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel('Rsigm, Ohm')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results', fontsize=7, x=0.46, y=1.005)
plt.grid(b=True, which='both', color='silver', linestyle='-')
plt.legend(loc='lower center', fontsize=8)
pylab.savefig(result_path + "/Self radiation resistances of the first dipole.png", bbox_inches='tight', dpi=160)
plt.close('all')

# Frequency dependence of radiation efficiency of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list[1:], Efficiency[0: len(frequency_list)-1, 0],
         color='C0', linestyle='-', linewidth='1.00', label='Efficiency dipole # 1')
plt.plot(frequency_list[1:], Efficiency[0: len(frequency_list)-1, central],
         color='C1', linestyle='-', linewidth='1.00', label='Efficiency of central dipole # ' + str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel('Efficiency, %')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results', fontsize=7, x=0.46, y=1.005)
plt.grid(b=True, which='both', color='silver', linestyle='-')
plt.legend(loc='lower center', fontsize=8)
pylab.savefig(result_path + "/Self radiation efficiency.png", bbox_inches='tight', dpi=160)
plt.close('all')


finish_time = time.time()
print('\n\n  The program execution lasted for ', round((finish_time - start_time), 2), 'seconds (',
      round((finish_time - start_time)/60, 2), 'min. ) \n')
print('\n      *** Program ' + Software_name + ' has finished! *** \n\n\n')
