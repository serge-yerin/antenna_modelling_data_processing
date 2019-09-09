# Python3
Software_version = '2019.09.07'
Software_name = 'Array Input Impedance Matrix Calculator'

#                 array_input_imppedance_matrix_calc
#    ***  Antenna Array Input Impedance Matrices Calculator   ***
# Program reads the result file of NEC modeling (.nec) where one dipole
# of the antenna array is excited and others are loaded, finds frequencies of
# analysis, excitation source, loads and their values, currents in source and
# loads. Then the program calculates mutual impedances, and store them
# to text files (.txt).

# The initial version was created for AA of 25 GURT dipols of 1 polarization of
# incoming waves, number of segments per dipole is 325 (Multi frequencies) then
# it was extended for other possibilities, but still you need to analyze
# carefully the results

# Unlike full version it needs only one file nec.out (not many files),
# mutual impedances are calculated based on induced currents.

#*******************************************************************************
#                           P A R A M E T E R S                                *
#*******************************************************************************
path_to_data = 'DATA/'
file_name = 'NEC_in_1.out.txt'
num_of_freq = 81            # Maximal possible number of frequencies analyzed
max_loads_num = 50          # Maximal number of loads in the file
print_or_not = 0            # 1 - print additional info, 0 - do not print
#*******************************************************************************
#                     I M P O R T   L I B R A R I E S                          *
#*******************************************************************************
import os
import time
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from package_common_modules.find_line_in_text_file import find_line_in_text_file

#*******************************************************************************
#                          M A I N   P R O G R A M                             *
#*******************************************************************************

print ('\n\n\n\n\n\n ')
print ('    ------------------------------------------------------------------')
print ('    ***                                                            ***')
print ('    ***  ', Software_name,' v.', Software_version, '  ***    (c) YeS 2019')
print ('    ***                                                            ***')
print ('    ------------------------------------------------------------------ \n\n\n')

# *** Time consumption calculation (beginning) ***
startTime = time.time()
currentTime = time.strftime("%H:%M:%S")
currentDate = time.strftime("%d.%m.%Y")
print ('  Today is ', currentDate, ' time is ', currentTime, '\n\n')

Rinp = np.full((1, num_of_freq, 5, max_loads_num), 0.0)
Xinp = np.full((1, num_of_freq, 5, max_loads_num), 0.0)
IFedDipReal = np.full((1, num_of_freq, 5, max_loads_num), 0.0)
IFedDipImag = np.full((1, num_of_freq, 5, max_loads_num), 0.0) # IFedDip(Number of files analyzed, freq step, dipole, LoadNum)

# *** Creating folder for results ***
result_path = 'Results'
if not os.path.exists(result_path):
    os.makedirs(result_path)

# *** Creating files for results ***
MutFullImpReFile = open(result_path + "/Mutual full impedances (real part).txt", "w")
MutFullImpImFile = open(result_path + "/Mutual full impedances (imag part).txt", "w")
MutFullImpModFile = open(result_path + "/Mutual full impedances (Module).txt", "w")
SelfFullImpReFile = open(result_path + "/Self full impedances (real part).txt", "w")
SelfFullImpImFile = open(result_path + "/Self full impedances (imag part).txt", "w")
SelfFullImpModFile = open(result_path + "/Self full impedances (Module).txt", "w")


TotFileNum = 1
for FileNum in range (TotFileNum):
    Active = 1

    # ***  Configuring the name of NEC output file ***

    Data_File = open(path_to_data + file_name, "r")

    print (' ')
    print ('**********************************************************************')

    print ('\n    Reading file: ', file_name + file_name, ' \n\n')


    # *** Counting number of lines in the file ***
    num_lines = sum(1 for line in Data_File)
    print ('\n  Number of lines in file = ', num_lines, '\n')
    Data_File.seek(0)


    # ***  Searching the lines where loads are listed in cards ***
    wires_with_loads = []
    segments_with_loads = []
    loads_values = []
    for line in Data_File:
        if line.startswith(' ***** DATA CARD NO.'):
            words_in_line = line.split()
            if words_in_line[5] == 'LD':
                wires_with_loads.append(int(words_in_line[7]))
                segments_with_loads.append(int(words_in_line[8]))
                loads_values.append(float(words_in_line[10]))
    LoadNum = len(segments_with_loads) # !!!!!!!!!!!! Possible error !!!!!!!!!!!!!!

    if print_or_not == 1:
        for i in range (LoadNum):
            print ('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], '  Resistance = ', loads_values[i])


    # ***  Searching and reading the reference frequency ***

    Data_File.seek(0) # return to the beginning of the file
    frequency_list = []
    for line in Data_File:
        if line.startswith('                                    FREQUENCY='):
            words_in_line = line.split()
            frequency_list.append(float(words_in_line[1]))

    if print_or_not == 1: print ('\n  Central frerquency = ', frequency_list[0], ' MHz \n')
    if print_or_not == 1:
        for i in range (1, len(frequency_list)):
            print ('  Frerquency = ', frequency_list[i], ' MHz')

    NoOfFrequencies = len(frequency_list) - 1
    print ('\n  Number of frequencies analyzed = ', NoOfFrequencies, '\n\n')

    print ('**********************************************************************')


    Data_File.seek(0) # return to the beginning of the file

    # Finding the main (central) frequency data and skip it
    line = find_line_in_text_file (Data_File, '                                    FREQUENCY=', 0, line)


    for step in range (NoOfFrequencies):
        # Seek for frequencies of analyses one by one
        line = find_line_in_text_file (Data_File, 'FREQUENCY=', 36, line)

        # *** Searching antenna input parameters block at the frequency ***
        line = find_line_in_text_file (Data_File, '- - - ANTENNA INPUT PARAMETERS - - -', 42, line)

        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = find_line_in_text_file (Data_File, 'IMPEDANCE', 64, line)


        line = Data_File.readline()   # Skip one line
        line = Data_File.readline()
        IFedDipReal[FileNum, step, 1+(Active-1)*2, 0] = float(line[36 : 48])
        IFedDipImag[FileNum, step, 1+(Active-1)*2, 0] = float(line[48 : 60])
        Rinp[FileNum, step, 1+(Active-1)*2, 0] = float(line[60 : 72])
        Xinp[FileNum, step, 1+(Active-1)*2, 0] = float(line[72 : 84])


        # Filling the mirror elements of the matrices
        #IFedDipReal[FileNum, step, 1+(Active)*2, 0] = IFedDipReal[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        #IFedDipImag[FileNum, step, 1+(Active)*2, 0] = IFedDipImag[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        #Rinp[FileNum, step, 1+(Active)*2, 0] = Rinp[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        #Xinp[FileNum, step, 1+(Active)*2, 0] = Xinp[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!


        print ('\n\n\n  For frequency = ', frequency_list[step+1], ' MHz')
        print ('  ------------------------------ \n')
        print ('  Self full impedance of the dipole = ', Rinp[FileNum, step, 1+(Active-1)*2, 0], Xinp[FileNum-1, step, 1+(Active-1)*2, 0], ' Ohm' )
        print ('  Current in fed point of active dipole = ', IFedDipReal[FileNum, step, 1+(Active-1)*2, 0],  IFedDipImag[FileNum-1, step, 1+(Active-1)*2, 0],' Amp. \n\n')

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Searching and reading the currennt in impedance load !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Seek for currents in loads
        line = find_line_in_text_file (Data_File, '- - - CURRENTS AND LOCATION - - -', 29, line)
        for i in range (6): line = Data_File.readline()   # Skip several lines


        for j in range (LoadNum):  # loop by number of loads
            intTemp = 0
            while (intTemp != wires_with_loads[j]):
                line = Data_File.readline()
                intTemp = int(line[6:12])
            if (segments_with_loads[j] < 1):   # if the number of segment = 1 we need to return to the begining of the line
                stop
            if (segments_with_loads[j] > 1):	# if the number of segment > 2 we need to skip lines
                for i in range (1, segments_with_loads[j] - 2): line = Data_File.readline()   # Skip several lines

            # Reading currents
            IFedDipReal[FileNum, step, 2+(Active-1)*2, j] = float(line[49 : 60])
            IFedDipImag[FileNum, step, 2+(Active-1)*2, j] = float(line[61 : 72])

            # Mirroring currents
            #IFedDipReal[FileNum, step, 2+(Active)*2, j] = IFedDipReal[FileNum, step, 2+(Active-1)*2, j]
            #IFedDipImag[FileNum, step, 2+(Active)*2, j] = IFedDipImag[FileNum, step, 2+(Active-1)*2, j]

            if print_or_not == 1: print ('  Current in the load No ', j+1, ' = ', IFedDipReal[FileNum, step, 2+(Active-1)*2, j], IFedDipImag[FileNum, step, 2+(Active-1)*2, j], ' Amp' )

    Data_File.close()

    print ('\n\n\n   ***********************************************************')
    print ('                  Calculations of mutual impedances ')
    print ('   ***********************************************************\n\n\n')


    IFedDip = IFedDipReal + 1j * IFedDipImag
    Zinp = Rinp + 1j * Xinp

    for step in range (NoOfFrequencies): # Loop by frequencies for calculations

        print ('  For frequency = ', frequency_list[step+1], ' MHz')
        print ('  ------------------------------')

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!!    Calculation of full mutual impedance    !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        for k in range (LoadNum):
            Zinp[FileNum, step, 2, k] = - ((IFedDip[FileNum, step, 2, k] * loads_values[k]) / IFedDip[FileNum, step, 1, 0])

        print ('    Zinp (0, 0)  = ', Zinp[FileNum, step, 2, 0] ,', Ohm \n')
        print (' ***********************************************************')



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!  Making files with calculations results  !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n     Writing the results to output files... \n\n\n')

# Mutual full impedances (real part).txt
MutFullImpReFile.write('  Frequency,MHz  ' + ''.join('    WireNo='+str(wires_with_loads[i]).zfill(5)+'   ' for i in range(LoadNum)) + '\n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.real(Zinp[FileNum, LineNum, 2, i]) for i in range(LoadNum)]
    MutFullImpReFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ '   \n')


# Mutual full impedances (imag part).txt
MutFullImpImFile.write('  Frequency,MHz  ' + ''.join('    WireNo='+str(wires_with_loads[i]).zfill(5)+'   ' for i in range(LoadNum)) + '\n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.imag(Zinp[FileNum, LineNum, 2, i]) for i in range(LoadNum)]
    MutFullImpImFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ '   \n')


# Mutual full impedances (Module).txt
MutFullImpModFile.write('  Frequency,MHz  ' + ''.join('    WireNo='+str(wires_with_loads[i]).zfill(5)+'   ' for i in range(LoadNum)) + '\n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.absolute(Zinp[FileNum, LineNum, 2, i]) for i in range(LoadNum)]
    MutFullImpModFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ '   \n')


# Self full impedances (real part).txt
SelfFullImpReFile.write('  Frequency,MHz         Re(Z11) \n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.real(Zinp[FileNum, LineNum, 1, 0])]
    SelfFullImpReFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ ' \n')

# Self full impedances (imag part).txt
SelfFullImpImFile.write('  Frequency,MHz         Im(Z11) \n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.imag(Zinp[FileNum, LineNum, 1, 0])]
    SelfFullImpImFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ ' \n')

# Self full impedances (Module).txt
SelfFullImpModFile.write('  Frequency,MHz         Abs(Z11) \n')
for LineNum in range (NoOfFrequencies):
    data = ['   %14.6e  ' % np.absolute(Zinp[FileNum, LineNum, 1, 0])]
    SelfFullImpModFile.write(('    %8.4f    ' % frequency_list[LineNum+1]) + ''.join(data)+ ' \n')

# Frequency dependance of dipole input impedance
fig = plt.figure(figsize = (10, 6))
rc('font', weight='normal')
plt.plot(frequency_list[1:len(frequency_list)], np.real(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'r', linestyle = '-', linewidth = '1.00', label = 'Re(Zinp)')
plt.plot(frequency_list[1:len(frequency_list)], np.imag(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'b', linestyle = '-', linewidth = '1.00', label = 'Im(Zinp)')
plt.plot(frequency_list[1:len(frequency_list)], np.absolute(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'violet', linestyle = '-', linewidth = '1.00', label = 'Abs(Zinp)')
plt.xlabel('Frequency, MHz')
plt.ylabel('R, X, |Z|, Ohm')
#plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.title('NEC modeling results', fontsize=7, x = 0.46, y = 1.005)
plt.legend(loc = 'lower center', fontsize = 10)
pylab.savefig("Results/Self full impedances.png", bbox_inches='tight', dpi = 160)
# plt.show()
plt.close('all')


loads = np.linspace(1, LoadNum, num = LoadNum)

# Dependence of mutual input impedance on dipole number for all frequencies
fig = plt.figure(figsize = (20, 12))
rc('font', weight='normal')
for step in range (NoOfFrequencies):
    plt.plot(loads, np.absolute(Zinp[FileNum, step, 2, 0:LoadNum]), linestyle = '-', linewidth = '1.00', label = str(frequency_list[step+1])+' MHz')
plt.xlabel('Number of load')
plt.ylabel('Abs (Z), Ohm')
#plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results. Dependence of mutual input impedance on dipole number for all frequencies', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = 'silver', linestyle = '-')
plt.legend(loc='center left', fontsize = 6, bbox_to_anchor=(1, 0.5))
pylab.savefig("Results/Mutual impedances (modules).png", bbox_inches='tight', dpi = 160)
# plt.show()
plt.close('all')


endTime = time.time()
print ('\n\n\n  The program execution lasted for ', round((endTime - startTime), 2), 'seconds (',
                                                round((endTime - startTime)/60, 2), 'min. ) \n')
print ('\n           *** Program ' + Software_name + ' has finished! *** \n\n\n')
