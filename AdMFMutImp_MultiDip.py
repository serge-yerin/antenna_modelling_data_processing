# Python3
Software_version = '2018.07.31'
#     AdMFrMutRadImp
# *** Advanced Multi Frequency Mutual&Self Radiation&Total Impedances reader and calculator  ***
# Program reads the result of one file of NEC modeling (.nec), finds frequencies of analysis,
# sources  loads and their values, currents in sources and loads, calculates mutual
# impedances and store them to text files (.txt).
# Unlike full version it does not need many files, mutual impedances are calculated based on induced currents

#*************************************************************
#                        PARAMETERS                          *
#*************************************************************
pi = 3.141593

#*************************************************************
#                   IMPORT LIBRARIES                         *
#*************************************************************
import os
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#*************************************************************
#                       FUNCTIONS                            *
#*************************************************************

# *** Find_line finds a line in .txt file which contains 'needed_text'  ***
#     needed_text - text to find
#     start_char - number of character where the needed_text starts in line
#     line - returns a text line from the file with needed_text
def Find_line (needed_text, start_char, line):
    tempChar = ''
    while tempChar != needed_text:
        line = Data_File.readline()
        tempChar = line[start_char : len(needed_text) + start_char]
    return line

#*************************************************************
#                      MAIN PROGRAM                          *
#*************************************************************

Rinp = np.full((150, 100, 5, 50), 0.0)
Xinp = np.full((150, 100, 5, 50), 0.0)
IFedDipReal = np.full((150, 100, 5, 50), 0.0)
IFedDipImag = np.full((150, 100, 5, 50), 0.0) # IFedDip(Number of files analyzed, freq step, dipole, LoadNum)



# *** Creating folder for results ***
if not os.path.exists('Results'):
    os.makedirs('Results')

# *** Creating files for results ***
MutFullImpReFile = open("Results/Mutual full impedances (real part).txt", "w")  # 9
MutFullImpImFile = open("Results/Mutual full impedances (imag part).txt", "w")  # 10
MutFullImpModFile = open("Results/Mutual full impedances (Module).txt", "w") # 17
SelfFullImpReFile = open("Results/Self full impedances (real part).txt", "w")   # 13
SelfFullImpImFile = open("Results/Self full impedances (imag part).txt", "w")   # 11
SelfFullImpModFile = open("Results/Self full impedances (Module).txt", "w")     # 19


for i in range(6): print (' ')
print ('    ----------------------------------------------------------')
print ('    ***                                                    ***')
print ('    ***        NEC output calculator v.', Software_version, '        ***')
print ('    ***           of self and mutual impedances            ***   (c) YeS 2018')
print ('    ***    of dipoles with one source and several loads    ***')
print ('    ***                                                    ***')
print ('    ----------------------------------------------------------')
for i in range(3): print (' ')



TotFileNum = 1
for FileNum in range (TotFileNum):

    Active = 1

    # ***  Configuring the name of NEC output file ***

    FileName = 'Data/NEC.out'
    Data_File = open(FileName, "r")

    for i in range (1): print (' ')
    for i in range (1): print ('**********************************************************************')

    print ('\n    Reading file: ', FileName, ' \n\n')


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

    for i in range (LoadNum):
        print ('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], '  Resistance = ', loads_values[i])


    # ***  Searching and reading the reference frequency ***

    Data_File.seek(0) # return to the beginning of the file
    frequency_list = []
    for line in Data_File:
        if line.startswith('                                    FREQUENCY='):
            words_in_line = line.split()
            frequency_list.append(float(words_in_line[1]))

    print ('\n  Central frerquency = ', frequency_list[0], ' MHz \n')
    for i in range (1, len(frequency_list)):
        print ('  Frerquency = ', frequency_list[i], ' MHz')

    NoOfFrequencies = len(frequency_list) - 1
    print ('\n  Number of frequencies analyzed = ', NoOfFrequencies, '\n\n')

    for i in range (1): print ('**********************************************************************')


    Data_File.seek(0) # return to the beginning of the file

    # Finding the main (central) frequency data and skip it
    line = Find_line ('                                    FREQUENCY=', 0, line)


    for step in range (NoOfFrequencies):
        # Seek for frequencies of analyses one by one
        line = Find_line ('FREQUENCY=', 36, line)

        # *** Searching antenna input parameters block at the frequency ***
        line = Find_line ('- - - ANTENNA INPUT PARAMETERS - - -', 42, line)

        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = Find_line ('IMPEDANCE', 64, line)


        line = Data_File.readline()   # Skip one line
        line = Data_File.readline()
        IFedDipReal[FileNum, step, 1+(Active-1)*2, 0] = float(line[36 : 48])
        IFedDipImag[FileNum, step, 1+(Active-1)*2, 0] = float(line[48 : 60])
        Rinp[FileNum, step, 1+(Active-1)*2, 0] = float(line[60 : 72])
        Xinp[FileNum, step, 1+(Active-1)*2, 0] = float(line[72 : 84])


        # Filling the mirror elements of the matrices
        IFedDipReal[FileNum, step, 1+(Active)*2, 0] = IFedDipReal[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        IFedDipImag[FileNum, step, 1+(Active)*2, 0] = IFedDipImag[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        Rinp[FileNum, step, 1+(Active)*2, 0] = Rinp[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!
        Xinp[FileNum, step, 1+(Active)*2, 0] = Xinp[FileNum, step, 1+(Active-1)*2, 0] # !!!!!!+++!!!!!!!


        print ('\n\n\n  For frequency = ', frequency_list[step+1], ' MHz')
        print ('  ------------------------------ \n')
        print ('  Self full impedance of the dipole = ', Rinp[FileNum, step, 1+(Active-1)*2, 0], Xinp[FileNum, step, 1+(Active-1)*2, 0], ' Ohm' )
        print ('  Current in fed point of active dipole = ', IFedDipReal[FileNum, step, 1+(Active-1)*2, 0],  IFedDipImag[FileNum, step, 1+(Active-1)*2, 0],' Amp. \n\n')

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Searching and reading the currennt in impedance load !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Seek for currents in loads
        line = Find_line ('- - - CURRENTS AND LOCATION - - -', 29, line)
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
            IFedDipReal[FileNum, step, 2+(Active)*2, j] = IFedDipReal[FileNum, step, 2+(Active-1)*2, j]
            IFedDipImag[FileNum, step, 2+(Active)*2, j] = IFedDipImag[FileNum, step, 2+(Active-1)*2, j]

            print ('  Current in the load No ', j+1, ' = ', IFedDipReal[FileNum, step, 2+(Active-1)*2, j], IFedDipImag[FileNum, step, 2+(Active-1)*2, j], ' Amp' )

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

#MutFullImpReFile = open("Results/Mutual full impedances (real part).txt", "w")  # 9
#MutFullImpImFile = open("Results/Mutual full impedances (imag part).txt", "w")  # 10
#MutFullImpModFile = open("Results/Mutual full impedances (Module).txt",   "w")  # 17
#SelfFullImpReFile = open("Results/Self full impedances (real part).txt",  "w")  # 13
#SelfFullImpImFile = open("Results/Self full impedances (imag part).txt",  "w")  # 11
#SelfFullImpModFile = open("Results/Self full impedances (Module).txt",    "w")  # 19


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


plt.figure()
rc('font', weight='normal')
plt.plot(frequency_list[1:len(frequency_list)], np.real(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'r', linestyle = '-', linewidth = '1.00', label = 'Re(Zinp)')
plt.plot(frequency_list[1:len(frequency_list)], np.imag(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'b', linestyle = '-', linewidth = '1.00', label = 'Im(Zinp)')
plt.plot(frequency_list[1:len(frequency_list)], np.absolute(Zinp[FileNum, 0:len(frequency_list)-1, 1, 0]), color = 'violet', linestyle = '-', linewidth = '1.00', label = 'Abs(Zinp)')
plt.xlabel('Frequency, MHz')
plt.ylabel('R, X, |Z|, Ohm')
#plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = '0.00',linestyle = '--')
plt.legend(loc = 'lower center', fontsize = 10)
pylab.savefig("Results/Self full impedances.png", bbox_inches='tight', dpi = 160)
# plt.show()
plt.close('all')


loads = np.linspace(1, LoadNum, num = LoadNum)


plt.figure()
rc('font', weight='normal')
for step in range (NoOfFrequencies):
    plt.plot(loads, np.absolute(Zinp[FileNum, step, 2, 0:LoadNum]), linestyle = '-', linewidth = '1.00', label = str(frequency_list[step+1])+' MHz')
plt.xlabel('Number of load')
plt.ylabel('Abs (Z), Ohm')
#plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
plt.title('NEC modeling results', fontsize=7, x = 0.46, y = 1.005)
plt.grid(b = True, which = 'both', color = '0.00',linestyle = '--')
plt.legend(loc = 'upper left', fontsize = 7)
pylab.savefig("Results/Mutual impedances (modules).png", bbox_inches='tight', dpi = 160)
# plt.show()
plt.close('all')
