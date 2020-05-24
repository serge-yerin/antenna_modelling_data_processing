# Python3
Software_version = '2020.05.02'
Software_name = 'Array Input and Radiation Impedance Matrix Calculator for NEC4 results'

#             array_input_and_rad_impedance_matrix_calc
#    ***  Antenna Array Input and Radiation Impedance Matrices Calculator   ***
# Program reads the result files of NEC4 modeling (.nec) where each (of non mirrored) dipole
# of the antenna array is excited one by one (in each file), finds frequencies of analysis,
# excitation sources, loads and their values, currents in sources and loads, antenna patterns.
# Then the program calculates mutual impedances, radiation resistances, and store them to
# text files (.txt).
# The initial version was created for AA of 25 GURT dipoles of 1 polarization of
# incoming waves, number of segments per dipole is 325 (Multi frequencies) then it was
# extended for other possibilities, but still you need to analyze carefully the results

# *************************************************************
#                         PARAMETERS                          *
# *************************************************************
path_to_data = 'DATA/'

NoOfWiresPerDipole = 98     # Number of wires per dipole (not segments!) needed to find loads  # 434
ArrayInputNum = 30          # Array inputs number (total number of dipoles in array or inputs of dipoles)
num_of_freq = 5             # Maximal possible number of frequencies analyzed
print_or_not = 0
make_txt_or_not = 0

# Rectangular UTR-2 antenna array 6 * 5 = 30 dipoles
NoOfDip =      [1,   2,  3,  7,  8,  9, 13, 14, 15]       # Dipoles being excited
MirrorFrom =   [1,   2,  3,  7,  8,  9, 13, 14, 15]	      # Active dipoles first quadrant
MirrorInto =   [30, 29, 28, 24, 23, 22, 18, 17, 16]       # Passive dipoles fourth quadrant # Cardinal point symmetry
MirrorInto_1 = [25, 26, 27, 19, 20, 21]                   # Passive dipoles third quadrant  # Horizontal mirror
MirrorInto_2 = [ 6,  5,  4, 12, 11, 10]                   # Passive dipoles second quadrant  # Vertical mirror
# All dipoles [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
Partners_1 =  [25,26,27,28,29,30,19,20,21,22,23,24,13,14,15,16,17,18, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6]
Partners_2 =  [ 6, 5, 4, 3, 2, 1,12,11,10, 9, 8, 7,18,17,16,15,14,13,24,23,22,21,20,19,30,29,28,27,26,25]

#pi = 3.141593
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

print ('\n\n\n\n\n    ---------------------------------------------------------------------------------')
print ('    ***                                                                          ***')
print ('    *** ', Software_name, ' ***')
print ('    ***                            v.', Software_version, '                                ***')
print ('    ***       Program analyzes NEC output files with it turn excitation of       ***')
print ('    ***         array elements (one source other terminals are loaded),          ***')
print ('    ***             reads and calculates self and mutual impedances,             ***')
print ('    ***                radiation resistances of array elements.                  ***')
print ('    ***               Works with frequency sweep results of NEC4                 ***')
print ('    ***                              (c) YeS 2020                                ***')
print ('    --------------------------------------------------------------------------------\n\n\n')

# *** Time consumption calculation (beginning) ***
startTime = time.time()
currentTime = time.strftime("%H:%M:%S")
currentDate = time.strftime("%d.%m.%Y")
print ('  Today is ', currentDate, ' time is ', currentTime, '\n\n')

NECfileNum = len(NoOfDip)    # TotFileNum = The number of NEC out files

Rinp = np.full((num_of_freq, ArrayInputNum, ArrayInputNum), 0.0)
Xinp = np.full((num_of_freq, ArrayInputNum, ArrayInputNum), 0.0)
IFedDipReal = np.full((num_of_freq, ArrayInputNum, ArrayInputNum), 0.0)
IFedDipImag = np.full((num_of_freq, ArrayInputNum, ArrayInputNum), 0.0) # IFedDip(frequency step, Active dipole, Passive dipole)
DistanceRP = np.full((NECfileNum, num_of_freq), 0.0)    # [FileNum, FreqStep] Distance for which NEC calculates radiation pattern

Points = np.full((91, 361), 0.0+1j*0.0)     #[theta, phi] Matrix of E-field data for integral calculations
RSigm = np.full((num_of_freq, ArrayInputNum, ArrayInputNum), 0.0+1j*0.0)    #[FreqStep, Dip1, Dip2] Radiation impedance matrix
Efficiency = np.zeros((num_of_freq, ArrayInputNum))

ETHcmplx = np.full((num_of_freq, ArrayInputNum, 181, 361), 0.0+1j*0.0)  # [FreqStep, DipNo, theta, phi] Theta components of E-field
EPHcmplx = np.full((num_of_freq, ArrayInputNum, 181, 361), 0.0+1j*0.0)  # [FreqStep, DipNo, theta, phi] Phi components of E-field

# *** Creating folder for results ***
result_path = 'Results'
if not os.path.exists(result_path):
    os.makedirs(result_path)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!                 Reading data from NEC output files                       !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for FileNum in range (len(NoOfDip)):	# Main loop by NEC output files

    # *** Configuring the name of NEC output file ***

    if NoOfDip[FileNum] < 10:
        no = '0' + str(NoOfDip[FileNum])
    else:
        no = str(NoOfDip[FileNum])

    FileName = path_to_data + 'UTR2_6x5-' + no + '.out'
    #FileName = path_to_data + 'UTR2_6x1_EX=' + str(FileNum + 1) + '.out'
    #FileName = path_to_data + 'GURT-V325_25x01_Nex=' + no + '.out.txt'
    #FileName = path_to_data + 'NEC_in_'+str(FileNum+1)+'.out.txt'


    # *** Opening NECoutput file ***
    Data_File = open(FileName, "r")

    print (' **********************************************************************')
    print ('\n * Reading file: ', FileName, ' ')


    #  Searching the lines where loads are listed
    wires_with_loads = []
    segments_with_loads = []
    loads_values = []
    line = find_line_in_text_file(Data_File, '- - - STRUCTURE IMPEDANCE LOADING - - -', 31)
    for i in range(5): line = Data_File.readline()  # Skip several lines
    words_in_line = 'A'
    #print('0 - ',  line)
    counter = 0
    while (words_in_line[0] != 'CP' and counter < 2):
        line = Data_File.readline()
        #print('1 - ', line)
        if line != '\n':
            #print('2 - ', line)
            counter = 0
            words_in_line = line.split()
            if words_in_line[4]+' '+words_in_line[5] == 'FIXED IMPEDANCE':
                wires_with_loads.append(int(words_in_line[0]))
                segments_with_loads.append(int(words_in_line[1]))
                loads_values.append(float(words_in_line[3]))
        else:
            counter += 1

    LoadNum = len(segments_with_loads)  # !!!!!!!!!!!! Possible error !!!!!!!!!!!!!!

    if print_or_not > 0:
        for i in range (LoadNum):
            print ('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], '  Resistance = ', loads_values[i])


    # Searching number of frequencies
    Data_File.seek(0)
    frequency_list = []
    num_of_frequencies = 0
    for line in Data_File:
        if line.startswith('                                 - - - - - - FREQUENCY - - - - - -'):
            num_of_frequencies += 1

    print ('   Number of frequencies analyzed = ', num_of_frequencies, ' ')

    if print_or_not > 0: print (' **********************************************************************')

    Data_File.seek(0) # return to the beginning of the file

    for step in range (num_of_frequencies):
        # Seek for current frequency
        line = find_line_in_text_file (Data_File, 'FREQUENCY=', 36)
        words_in_line = line.split()
        frequency_list.append(float(words_in_line[1]))

        # *** Searching antenna input parameters block at the frequency ***
        line = find_line_in_text_file (Data_File, '- - - ANTENNA INPUT PARAMETERS - - -', 42)

        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = find_line_in_text_file (Data_File, 'IMPEDANCE', 64)
        line = Data_File.readline()   # Skip one line

        line = Data_File.readline()
        IFedDipReal[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[36 : 48])
        IFedDipImag[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[48 : 60])
        Rinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[60 : 72])
        Xinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[72 : 84])

        if print_or_not > 0:
            print ('\n\n   For frequency = ', frequency_list[step], ' MHz')
            print ('  ---------------------------- \n')
            print ('  Self full impedance of the dipole = ', Rinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1], Xinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1], ' Ohm' )
            print ('  Current in fed point of active dipole = ', IFedDipReal[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1],  IFedDipImag[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1],' Amp. \n')


        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !! Searching and reading the current in impedance load !!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # *** Seek for currents in loads ***
        line = find_line_in_text_file (Data_File, '- - - CURRENTS AND LOCATION - - -', 29)
        for i in range (5): line = Data_File.readline()   # Skip several lines

        for j in range (LoadNum):  # loop by number of loads
            intTemp = 0
            while (intTemp < wires_with_loads[j]):   # !=
                line = Data_File.readline()
                intTemp = int(line[6:12])
            if (segments_with_loads[j] < 1):    # if the number of segment = 1 we need to return to the begining of the line
                sys.exit('ERROR!')
            if (segments_with_loads[j] > 1):	# if the number of segment > 2 we need to skip lines
                for i in range (segments_with_loads[j] - 1): line = Data_File.readline()   # Skip several lines

            DipNum = int((wires_with_loads[j]) / NoOfWiresPerDipole) + 1   # Calculation of the dipole number   !!!!
            #DipNum = int((wires_with_loads[j]) / wires_with_loads[0]) + 1   # Calculation of the dipole number   !!!!

            # *** Reading currents ***
            IFedDipReal[step, NoOfDip[FileNum]-1, DipNum-1] = float(line[49 : 60])
            IFedDipImag[step, NoOfDip[FileNum]-1, DipNum-1] = float(line[61 : 72])

            if print_or_not > 0:
                print ('  Dip No', DipNum, ' Wire No', wires_with_loads[j], ' Segm No', segments_with_loads[j])
                print ('  Current in the load No ', j+1, ' = ', IFedDipReal[step, NoOfDip[FileNum]-1, DipNum-1], IFedDipImag[step, NoOfDip[FileNum]-1, DipNum-1], ' Amp' )


        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!  Reading of E values in radiation patterns  !!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        line = find_line_in_text_file (Data_File, '- - - RADIATION PATTERNS - - -', 48)
        line = Data_File.readline()
        DistanceRP[FileNum, step] = float(line[60 : 74])
        #print ('\n  Distance of RP calculations = ', DistanceRP[FileNum, step], ' m. \n')

        for i in range (7): line = Data_File.readline()   # Skip several lines of table header

        # *** loop by lines to read e-field amplitudes and phases for all theta and phi ***

        for k in range (65341): # 65341 = 181 * 361
            line = Data_File.readline()
            theta1 = float(line[0 : 8])
            phi1 = float(line[9 : 17])
            Etheta = float(line[76 : 87])
            ang1 = float(line[89 : 97])
            Ephi = float(line[100 : 112])
            ang2 = float(line[112 : 120])

            thetaInt = (int(theta1))+90
            phiInt = int(phi1)
            ETHcmplx[step,NoOfDip[FileNum]-1,thetaInt,phiInt] = np.complex((Etheta*np.cos(ang1*np.pi/180)),(Etheta*np.sin((ang1*np.pi/180))))
            EPHcmplx[step,NoOfDip[FileNum]-1,thetaInt,phiInt] = np.complex((Ephi*np.cos(ang2*np.pi/180)),(Ephi*np.sin((ang2*np.pi/180))))

    Data_File.close()

IFedDip = IFedDipReal + 1j * IFedDipImag
Zinp = Rinp + 1j * Xinp
del IFedDipReal, IFedDipImag, Rinp, Xinp


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!              Processing data and calculating results                     !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n    Calculations started at: ', time.strftime("%H:%M:%S"), '  \n\n')

# CURRENTS MIRRORING
# Diagonal mirror
for i in range(len(MirrorInto)):
    for j in range(ArrayInputNum):
        IFedDip[:, MirrorInto[i] - 1, j] = (IFedDip[:, MirrorFrom[i] - 1, ArrayInputNum - 1 - j])  # !!!  ArrayInputNum+1

# Quadrant 1 mirror
for i in range (len(MirrorInto_1)):
    for j in range (ArrayInputNum):
        IFedDip[:, MirrorInto_1[i]-1, j] = (IFedDip[:, MirrorFrom[i]-1, Partners_1[j]-1])

# Quadrant 2 mirror
for i in range (len(MirrorInto_2)):
    for j in range (ArrayInputNum):
        IFedDip[:, MirrorInto_2[i]-1, j] = (IFedDip[:, MirrorFrom[i]-1, Partners_2[j]-1])


# *** Calculations of Impedance matrices ***

for i in range (NECfileNum):           # Loop by files
    for k in range (ArrayInputNum):    # Loop by dipoles
        if k != (NoOfDip[i]-1):
            Zinp[:, NoOfDip[i]-1, k] = - (IFedDip[:, NoOfDip[i]-1, k] * loads_values[i] / IFedDip[:, NoOfDip[i]-1, NoOfDip[i]-1]) # !!!-!!!


# Mirroring the Impedances for dipoles that were not excited

for i in range (len(MirrorInto)):     # Loop by active dipoles below array diagonal
    for j in range (ArrayInputNum):         # Loop by dipoles
        Zinp[:, MirrorInto[i]-1, ArrayInputNum - 1 - j] = Zinp[:, MirrorFrom[i]-1, j]  # ArrayInputNum + 1

# Quadrant 1 mirror
for i in range (len(MirrorInto_1)):     # Loop by
    for j in range (ArrayInputNum):         # Loop by dipoles
        Zinp[:, MirrorInto_1[i]-1, Partners_1[j]-1] = Zinp[:, MirrorFrom[i]-1, j]

# Quadrant 2 mirror
for i in range (len(MirrorInto_2)):     # Loop by
    for j in range (ArrayInputNum):         # Loop by dipoles
        Zinp[:, MirrorInto_2[i]-1, Partners_2[j]-1] = Zinp[:, MirrorFrom[i]-1, j]


# PATTERNS MIRRORING
#'''
# Cardinal point mirror
ETHcmplx_tmp = np.zeros_like(ETHcmplx)
EPHcmplx_tmp = np.zeros_like(EPHcmplx)
for i in range (len(MirrorInto)):   # loop by dipoles
    for p in range (361):       # loop by phi
        ETHcmplx_tmp[:, MirrorInto[i]-1, :, p] = - ETHcmplx[:, MirrorFrom[i]-1, :, 360-p]  # -
        EPHcmplx_tmp[:, MirrorInto[i]-1, :, p] = - EPHcmplx[:, MirrorFrom[i]-1, :, 360-p]  # -
    for p in range(181):
        ETHcmplx[:, MirrorInto[i]-1, :,     p] = ETHcmplx_tmp[:, MirrorInto[i]-1, :, 180-p]
        ETHcmplx[:, MirrorInto[i]-1, :, 180+p] = ETHcmplx_tmp[:, MirrorInto[i]-1, :, 360-p]
        EPHcmplx[:, MirrorInto[i]-1, :,     p] = EPHcmplx_tmp[:, MirrorInto[i]-1, :, 180-p]
        EPHcmplx[:, MirrorInto[i]-1, :, 180+p] = EPHcmplx_tmp[:, MirrorInto[i]-1, :, 360-p]
del ETHcmplx_tmp, EPHcmplx_tmp

# Horizontal mirror (right to left)
for i in range (len(MirrorInto_1)):   # loop by dipoles
    for p in range (361):             # loop by phi
        ETHcmplx[:, MirrorInto_1[i]-1, :, p] =   (ETHcmplx[:, MirrorFrom[i]-1, :, 360-p])
        EPHcmplx[:, MirrorInto_1[i]-1, :, p] = - (EPHcmplx[:, MirrorFrom[i]-1, :, 360-p]) # - !!!
    #ETHcmplx[:, MirrorInto_1[i] - 1, :, :].real *= -1
    #EPHcmplx[:, MirrorInto_1[i] - 1, :, :].real *= -1



# Vertical mirror (up to down)
for i in range (len(MirrorInto_2)):   # loop by dipoles
    for p in range (181):             # loop by phi
        ETHcmplx[:, MirrorInto_2[i]-1, :,     p] = - (ETHcmplx[:, MirrorFrom[i]-1, :, 180-p]) # - !!!  np.conj() ! a.real *= factarr
        EPHcmplx[:, MirrorInto_2[i]-1, :,     p] =   (EPHcmplx[:, MirrorFrom[i]-1, :, 180-p])
        ETHcmplx[:, MirrorInto_2[i]-1, :, 180+p] = - (ETHcmplx[:, MirrorFrom[i]-1, :, 360-p]) # - !!!
        EPHcmplx[:, MirrorInto_2[i]-1, :, 180+p] =   (EPHcmplx[:, MirrorFrom[i]-1, :, 360-p])
    #ETHcmplx[:, MirrorInto_2[i] - 1, :, :].imag *= -1
    #EPHcmplx[:, MirrorInto_2[i] - 1, :, :].imag *= -1

for step in range (num_of_frequencies):     # Loop by frequencies

     # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!         Integrals calculation              !!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i in range (ArrayInputNum):       # Loop by first dipole
        for j in range (ArrayInputNum):   # Loop by second dipole
            if (i >= j):
            #if (i != j):

                # Forming the integration elements for Theta angle
                for t in range (91):
                    Points[t, :] = ETHcmplx[step, i, t+90, :] * np.conj(ETHcmplx[step, j, t+90, :]) * (np.sin(t*np.pi/180.) + 1j*0.0)

                Integral1 = intergration_trap_2d (0, 90, 0, 360, 90, 360, Points)

                # Forming the integration elements for Phi component
                for t in range (91):
                    Points[t, :] = EPHcmplx[step, i, t+90, :] * np.conj(EPHcmplx[step, j, t+90, :]) * (np.sin(t*np.pi/180.) + 1j*0.0)

                # Calculating of integrals by Phi angle
                Integral2 = intergration_trap_2d (0, 90, 0, 360, 90, 360, Points)

                # Calculating the elements of radiation impedance matrix !!! Distance[1] !!!
                RSigm[step, i, j] = ((Integral1 + Integral2) * np.power(DistanceRP[0, 0], 2)) / (120.0 * np.pi * (IFedDip[step, i, i] * np.conj(IFedDip[step, j, j])))
                RSigm[step, j, i] = np.conj(RSigm[step, i, j])

    if print_or_not > 0:
        currentTime = time.strftime("%H:%M:%S")
        print ('\n    Done for frequency = ', frequency_list[step], ' MHz at: ', currentTime, ' ')


# *** Calculations of self radiation efficiencies of dipoles ***
for i in range (ArrayInputNum):
    Efficiency[:, i] = (np.real(RSigm[:,i,i]) / np.real(Zinp[:,i,i])) * 100



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!             Displaying data and results on the screen                    !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


for step in range (num_of_frequencies):	 # Loop by frequencies for displaying results on the screen

    print ('\n\n  For frequency = ', frequency_list[step], ' MHz')
    print ('  --------------------------')
    print ('    Zinp (1, 1) = ', Zinp[step, 0, 0])
    print ('   Rsigm (1, 1) = ', RSigm[step, 0, 0])
    print ('\n  Self radiation efficiency of the first dipole = ', round(Efficiency[step, 0],3), ' %')   # i -> 0



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!              Writing data and results to output files                    !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if make_txt_or_not > 0:
    print ('\n\n\n     Writing the results to output files... \n')


    # *** Creating files for results ***
    CurrentMatrixFile = open(result_path + "/Matrices of Currents.txt", "w")              # 24
    ImpedanceMatrixFile = open(result_path + "/Matrices of Impedances.txt", "w")          # 26
    RadImpMatrixFile = open(result_path + "/Matrices of Radiation Impedances.txt", "w")   # 27

    for step in range (num_of_frequencies):	 # Loop by frequencies for results files writing


        # *** Matrix of complex currents in dipoles.txt ***
        CurrentMatrixFile.write('\n\n    *** Matrix of currents in fed points of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
        data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
        CurrentMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
        data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
        CurrentMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
        for k in range(ArrayInputNum):
            data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(IFedDip[step, k, i]), np.imag(IFedDip[step, k, i])) for i in range(ArrayInputNum)]
            CurrentMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')


        # *** Writing impedances to file ***
        ImpedanceMatrixFile.write('\n\n    *** Matrix of impedances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
        data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
        ImpedanceMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
        data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
        ImpedanceMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
        for k in range(ArrayInputNum):
            data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(Zinp[step, k, i]), np.imag(Zinp[step, k, i])) for i in range(ArrayInputNum)]
            ImpedanceMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')


        # *** Writing radiation resistances to file ***
        RadImpMatrixFile.write('\n\n    *** Matrix of radiation resistances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step])
        data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
        RadImpMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
        data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
        RadImpMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
        for k in range(ArrayInputNum):
            data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(RSigm[step, k, i]), np.imag(RSigm[step, k, i])) for i in range(ArrayInputNum)]
            RadImpMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')



        # *** Writing currents on dipols inputs to files for MatCAD
        FileName = (result_path + '/Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        CurrentFile.close()


        # *** Writing normalized currents on dipols inputs to files for MatCAD ***
        FileName = (result_path + '/Normalized Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        CurrentFile.close()


        # *** Writing impedances of dipols to files for MatCAD ***
        FileName = (result_path + '/Impedances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.real(Zinp[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.imag(Zinp[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        CurrentFile.close()


        # *** Writing radiation resistances of dipols to files for MatCAD ***
        FileName = (result_path + '/Radiation resistances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.real(RSigm[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(np.imag(RSigm[step, i, k])) for i in range(ArrayInputNum)]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        CurrentFile.close()



        # *** Writing radiation efficiencies of dipols to files for MatCAD ***
        FileName = (result_path + '/Radiation efficiencies of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            data = ['   {:+14.8e}'.format(Efficiency[step, i])]
            CurrentFile.write(' ' + ''.join(data)+ '  \n')
        CurrentFile.close()


        # *** Writing normalized Etheta RP to files for MatCAD ***
        FileName = (result_path + '/Normalized to current Etheta RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            for t in range (91):
                data = [' {:+14.8e}'.format(np.real(ETHcmplx[step,k,t+90,p] / IFedDip[step,k,k])) for p in range(361)]
                CurrentFile.write('' + ''.join(data)+ ' \n')
        for k in range (ArrayInputNum):
            for t in range (91):
                data = [' {:+14.8e}'.format(np.imag(ETHcmplx[step,k,t+90,p] / IFedDip[step,k,k])) for p in range(361)]
                CurrentFile.write('' + ''.join(data)+ ' \n')
        CurrentFile.close()



        # *** Writing normalized Ephi RP to files for MatCAD ***
        FileName = (result_path + '/Normalized to current Ephi RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step])
        CurrentFile = open(FileName, "w")  # 6
        for k in range (ArrayInputNum):
            for t in range (91):
                data = [' {:+14.8e}'.format(np.real(EPHcmplx[step,k,t+90,p] / IFedDip[step,k,k])) for p in range(361)]
                CurrentFile.write('' + ''.join(data)+ ' \n')
        for k in range (ArrayInputNum):
            for t in range (91):
                data = [' {:+14.8e}'.format(np.imag(EPHcmplx[step,k,t+90,p] / IFedDip[step,k,k])) for p in range(361)]
                CurrentFile.write('' + ''.join(data)+ ' \n')
        CurrentFile.close()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!                            F I G U R E S                                 !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n     Making figures... \n\n')

# Figures of diagonal matrix elements vs. dipole number for various frequencies

x_values = np.linspace(1, ArrayInputNum, num = ArrayInputNum) # Numbers of dipoles
data = np.zeros((ArrayInputNum))

for step in range (len(frequency_list)):

    # *** Figure of current absolute value ***
    rc('font', weight='normal', size=7)
    for i in range (ArrayInputNum):
        data[i] = np.absolute(IFedDip[step, i, i])
    plt.figure()
    plt.scatter(x_values, data[:], color = 'C0', label = r'|$I_{inp}$|')
    plt.xlabel('Number of dipole')
    plt.ylabel('Current at dipole input, A')
    plt.title('NEC4 modeling results: input currents (absolute values) f = %5.2f MHz' % frequency_list[step], fontsize=8)
    plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
    plt.legend(loc = 'lower center', fontsize = 8)
    pylab.savefig(result_path + "/Dipole currents (abs) f = %5.2f MHz.png" % frequency_list[step], bbox_inches='tight', dpi = 180)
    plt.close('all')

    # *** Figure of input impedance absolute value at the main diagonal ***
    for i in range (ArrayInputNum):
        data[i] = np.absolute(Zinp[step, i, i])
    plt.figure()
    plt.scatter(x_values, data[:], color = 'C3', label = r'|$Z_{inp}$|')
    plt.xlabel('Number of dipole')
    plt.ylabel(r'|$Z_{inp}$|, Ohm')
    plt.title('NEC4 modeling results: input impedances (absolute values) f = %5.2f MHz' % frequency_list[step], fontsize=8)
    plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
    plt.legend(loc = 'lower center', fontsize = 8)
    pylab.savefig(result_path + "/Self full impedances (abs) f = %5.2f MHz.png" % frequency_list[step], bbox_inches='tight', dpi = 180)
    plt.close('all')

    # *** Figure of self radiation resistances (at the matrix main diagonal) ***
    for i in range (ArrayInputNum):
        data[i] = np.real(RSigm[step, i, i])
    plt.figure()
    plt.scatter(x_values, data[:], color = 'C4', label = r'$R_{\Sigma}$')
    plt.xlabel('Number of dipole')
    plt.ylabel(r'$R_{\Sigma}$, Ohm')
    plt.title('NEC4 modeling results: radiation resistances f = %5.2f MHz' % frequency_list[step], fontsize=8)
    plt.legend(loc = 'lower center', fontsize = 8)
    plt.grid(b=True, which='both', color='silver', linestyle='-')
    pylab.savefig(result_path + "/Self radiation resistances f = %5.2f MHz.png" % frequency_list[step], bbox_inches='tight', dpi = 180)
    plt.close('all')




# *** Figure of current absolute value matrix nondiagonal map ***
for step in range (len(frequency_list)):
    data = np.absolute(IFedDip[step, :, :])
    for i in range (ArrayInputNum):
        data[i,i] = 0

    figure_color_map(data, 'Number of dipole','Number of dipole',
                       'Matrix of currents at inputs of antenna array at %5.2f MHz'% frequency_list[step],
                       'NEC4 modeling results',
                       result_path + "/Dipole currents matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step])


# *** Figure of impedance absolute value matrix nondiagonal map ***
for step in range (len(frequency_list)):
    data = np.absolute(Zinp[step, :, :])
    for i in range (ArrayInputNum):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole','Number of dipole',
                       'Matrix of impedances at inputs of antenna array at %5.2f MHz'% frequency_list[step],
                       'NEC4 modeling results, self values set to zeros to show mutual pattern',
                       result_path + "/Dipole impedances matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step])



# *** Figure of radiation impedance real part matrix nondiagonal map ***
for step in range (len(frequency_list)):
    data = np.zeros((ArrayInputNum, ArrayInputNum), dtype=np.float)
    data[:,:] = np.real(RSigm[step, :, :])
    for i in range (ArrayInputNum):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole','Number of dipole',
                       'Matrix of mutual radiation resistances of dipoles in antenna array at %5.2f MHz'% frequency_list[step],
                       'NEC4 modeling results, self values set to zeros to show mutual pattern',
                       result_path + "/Radiation impedance matrix nondiagonal (Re) f = %5.2f MHz.png" % frequency_list[step])



# *** Figure of radiation impedance imaginary part matrix nondiagonal map ***
for step in range (len(frequency_list)):
    data = np.zeros((ArrayInputNum, ArrayInputNum), dtype=np.float)
    data[:,:] = np.imag(RSigm[step, :, :])
    for i in range (ArrayInputNum):
        data[i, i] = 0

    figure_color_map(data, 'Number of dipole','Number of dipole',
                       'Matrix of mutual radiation reactance of dipoles in antenna array at %5.2f MHz'% frequency_list[step],
                       'NEC4 modeling results, self values set to zeros to show mutual pattern',
                       result_path + "/Radiation impedance matrix nondiagonal (Im) f = %5.2f MHz.png" % frequency_list[step])



# Calculate number od central dipole (or the nearest to phase center)
central = int(ArrayInputNum / 2) + 1

# Frequency dependence of input impedance of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list, np.real(Zinp[:, 0, 0]), color = 'C0', linestyle = '-', linewidth = '1.00', label = r'Re($Z_{inp}$) dipole # 1')
plt.plot(frequency_list, np.imag(Zinp[:, 0, 0]), color = 'C1', linestyle = '-', linewidth = '1.00', label = r'Im($Z_{inp}$) dipole # 1')
plt.plot(frequency_list, np.absolute(Zinp[:, 0, 0]), color = 'C4', linestyle = '-', linewidth = '1.00', label = r'|$Z_{inp}$| dipole # 1')
plt.plot(frequency_list, np.real(Zinp[:, central, central]), color = 'C0', linestyle = '--', linewidth = '1.00', label = r'Re($Z_{inp}$) dipole # '+str(central))
plt.plot(frequency_list, np.imag(Zinp[:, central, central]), color = 'C1', linestyle = '--', linewidth = '1.00', label = r'Im($Z_{inp}$) dipole # '+str(central))
plt.plot(frequency_list, np.absolute(Zinp[:, central, central]), color = 'C4', linestyle = '--', linewidth = '1.00', label = r'|$Z_{inp}$| dipole # '+str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel(r'$R_{inp}$, $X_{inp}$, |$Z_{inp}$|, Ohm' )  # 'R, X, |Z|, Ohm'
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC4 modeling results: input impedances vs. frequency', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.legend(loc = 'upper left', fontsize = 6)
pylab.savefig(result_path + "/Self full impedances of first dipole.png", bbox_inches='tight', dpi = 180)
plt.close('all')



# Frequency dependence of current of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list, np.real(IFedDip[:, 0, 0]), color = 'C0', linestyle = '-', linewidth = '1.00', label = r'Re($I_{inp}$) dipole # 1')  # 'Re(Iinp) dipole # 1'
plt.plot(frequency_list, np.imag(IFedDip[:, 0, 0]), color = 'C1', linestyle = '-', linewidth = '1.00', label = r'Im($I_{inp}$) dipole # 1')
plt.plot(frequency_list, np.absolute(IFedDip[:, 0, 0]), color = 'C4', linestyle = '-', linewidth = '1.00', label = r'|$I_{inp}$| dipole # 1')
plt.plot(frequency_list, np.real(IFedDip[:, central, central]), color = 'C0', linestyle = '--', linewidth = '1.00', label = r'Re($I_{inp}$) dipole # '+str(central))
plt.plot(frequency_list, np.imag(IFedDip[:, central, central]), color = 'C1', linestyle = '--', linewidth = '1.00', label = r'Im($I_{inp}$) dipole # '+str(central))
plt.plot(frequency_list, np.absolute(IFedDip[:, central, central]), color = 'C4', linestyle = '--', linewidth = '1.00', label = r'|$I_{inp}$| dipole # '+str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel(r'Im($I_{inp}$), Re($I_{inp}$), |$I_{inp}$|, A')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC4 modeling results: input currents vs. frequency', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.legend(loc = 'upper left', fontsize = 6)
pylab.savefig(result_path + "/Self current of first dipole.png", bbox_inches='tight', dpi = 180)
plt.close('all')



# Frequency dependence of radiation resistance of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list, np.real(RSigm[:, 0, 0]), color = 'C0', linestyle = '-', linewidth = '1.00', label = r'$R_{\Sigma}$ dipole # 1')
plt.plot(frequency_list, np.real(RSigm[:, central, central]), color = 'C1', linestyle = '-', linewidth = '1.00', label = r'$R_{\Sigma}$ dipole # '+str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel(r'$R_{\Sigma}$, Ohm')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC4 modeling results: radiation resistance vs. frequency', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.legend(loc = 'upper left', fontsize = 6)
pylab.savefig(result_path + "/Self radiation resistances of first dipole.png", bbox_inches='tight', dpi = 180)
plt.close('all')

# Frequency dependence of radiation efficiency of the first dipole and the central one
plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list, Efficiency[:, 0], color = 'C0', linestyle = '-', linewidth = '1.00', label = 'Efficiency dipole # 1')
plt.plot(frequency_list, Efficiency[:, central], color = 'C1', linestyle = '-', linewidth = '1.00', label = 'Efficiency dipole # '+str(central))
plt.xlabel('Frequency, MHz')
plt.ylabel('Efficiency, %')
# plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC4 modeling results: radiation efficiency vs. frequency', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.legend(loc = 'upper left', fontsize = 6)
pylab.savefig(result_path + "/Self radiation efficiency.png", bbox_inches='tight', dpi = 180)
plt.close('all')

'''
#import mpl_toolkits.mplot3d.axes3d as axes3d

theta, phi = np.linspace(0, 2 * np.pi, 361), np.linspace(0, np.pi, 91)
THETA, PHI = np.meshgrid(theta, phi)
R = np.abs(EPHcmplx[0,0,90:,:])
X = R * np.sin(PHI) * np.cos(THETA)
Y = R * np.sin(PHI) * np.sin(THETA)
Z = R * np.cos(PHI)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet', linewidth=0, antialiased=False, alpha=0.5)
pylab.savefig(result_path + "/Pattern.png", bbox_inches='tight', dpi = 180)
'''






endTime = time.time()
print('\n\n  The program execution lasted for ', round((endTime - startTime), 2), 'seconds (',
                round((endTime - startTime)/60, 2), 'min. ) \n')
print('\n      *** Program ' + Software_name + ' has finished! *** \n\n\n')
