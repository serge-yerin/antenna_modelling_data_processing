# Python3
Software_version = '2019.09.02'
Software_name = 'Array Impedance Matrix Calculator'
#                 array_imppedance_matrix_calc
#             ***  Antenna Array Impedance Matrices Calculator   ***
# Program reads the result files of NEC modeling (.nec) where each (of non mirrored) dipole
# of the antenna array is excited one by one (in each file), finds frequencies of analysis,
# excitation sources, loads and their values, currents in sources and loads, antenna patterns.
# Then the program calculates mutual impedances, radiation resistances, and store them to
# text files (.txt).

# The initial version was created for AA of 25 GURT dipols of 1 polarization of
# incoming waves, number of segments per dipole is 325 (Multi frequencies) then it was
# extended for other possibilities, but still you need to analyze carefully the results


#*************************************************************
#                        PARAMETERS                          *
#*************************************************************
pi = 3.141593

NoOfWiresPerDipole = 45
ArrayInputNum = 25 # Array inputs number (total number of dipoles in array or inputs of dipoles)

NoOfDip =    [1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13]   # Dipols being excited
MirrorFrom = [1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12]	    # Active dipols under the diagonal
MirrorInto = [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14]       # Passive dipols above the diagonal

NECfileNum = len(NoOfDip)    # TotFileNum = The number of NECout files

#*************************************************************
#                   IMPORT LIBRARIES                         *
#*************************************************************
import os
import time
#from scipy.integrate import simps
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from package_common_modules.find_line_in_text_file import find_line_in_text_file
from package_common_modules.intergration_trap_2d import intergration_trap_2d
from package_common_modules.figure_color_map import figure_color_map
#*************************************************************
#                       FUNCTIONS                            *
#*************************************************************

# *** Find_line finds a line in .txt file which contains 'needed_text'  ***
#     needed_text - text to find
#     start_char - number of character where the needed_text starts in line
#     line - returns a text line from the file with needed_text
#def Find_line (needed_text, start_char, line):
#    tempChar = ''
#    while tempChar != needed_text:
#        line = Data_File.readline()
#        tempChar = line[start_char : len(needed_text) + start_char]
#    return line

#def IntegrTrap2x (a, b, c, d, m, n, Points):
#    pi = 3.141593

    # a - the start integration limit by the first axis
    # b - the end integration limit by the first axis
    # c - the start integration limit by the second axis
    # d - the end integration limit by the second axis
    # n - number of points on 1 axis
    # m - number of points on 2 axis

#    Sum = 0.0
#    for i in range (m):
#        for k in range (n):
#            Sum = Sum + Points[i, k]

#    IntegrResult = ((((b)-(a))*pi/180.0)/m) * ((((d)-(c))*pi/180.0)/n) * ((Points[a,c] + Points[b,d]) + Sum)
#    return IntegrResult


# Plotting color map figures
#def MakeFigureColorMap(data, XLabel, YLabel, SupTitle, Title, FileName):
#    plt.figure(figsize=(10, 6))
#    rc('font', weight='normal', size = 8)
#    plt.pcolor(data, cmap = 'jet')
#    plt.axes().set_aspect('equal')
#    plt.xlabel(XLabel)
#    plt.ylabel(YLabel)
#    plt.colorbar()
#    plt.suptitle(SupTitle, fontsize=8, fontweight='bold')
#    plt.title(Title, fontsize = 8)
    # plt.subplots_adjust(left=0.2, right=0.3)
#    pylab.savefig(FileName, bbox_inches='tight', dpi = 160)
#    plt.close('all')

#*************************************************************
#                      MAIN PROGRAM                          *
#*************************************************************

print (' \n\n\n\n\n\n')
print ('    -----------------------------------------------------------------------')
print ('    ***                                                                 ***')
print ('    ***       ', Software_name,' v.', Software_version, '        ***')
print ('    ***   Program analyzes NEC output files with it turn excitation of  ***')
print ('    ***     array elements (one source other terminals are loaded),     ***')
print ('    ***         reads and calculates self and mutual impedances,        ***')
print ('    ***            radiation resistances of array elements.             ***   (c) YeS 2019')
print ('    ***           Works with frequency sweep results of NEC             ***')
print ('    ***                                                                 ***')
print ('    -----------------------------------------------------------------------\n\n\n')

# *** Time consumption calculation (beginning) ***
startTime = time.time()
currentTime = time.strftime("%H:%M:%S")
currentDate = time.strftime("%d.%m.%Y")
print ('  Today is ', currentDate, ' time is ', currentTime, '\n\n')

Rinp = np.full((100, 50, 50), 0.0)
Xinp = np.full((100, 50, 50), 0.0)
IFedDipReal = np.full((100, 50, 50), 0.0)
IFedDipImag = np.full((100, 50, 50), 0.0)   # IFedDip(frequency step, Active dipole, Passive dipole)
DistanceRP = np.full((15, 100), 0.0)        # [FileNum, FreqStep] Distance for which NEC calculates radiation pattern

Points = np.full((91, 361), 0.0+1j*0.0)     #[theta, phi] Matrix of E-field data for integral calculations
RSigm = np.full((100, 50, 50), 0.0+1j*0.0) #[FreqStep, Dip1, Dip2] Radiation impedance matrix
Efficiency = np.zeros((100, 100))

ETHcmplx = np.full((100, 50, 181, 361), 0.0+1j*0.0) # [FreqStep, DipNo, theta, phi] Theta components of E-field
EPHcmplx = np.full((100, 50, 181, 361), 0.0+1j*0.0) # [FreqStep, DipNo, theta, phi] Phi components of E-field

# *** Creating folder for results ***
if not os.path.exists('Results'):
    os.makedirs('Results')


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!                 Reading data from NEC output files                       !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for FileNum in range (len(NoOfDip)):	# Main loop by NEC output files

    # *** Configuring the name of NECoutput file ***

    if FileNum < 9:
        no = '0' + str(NoOfDip[FileNum])
    else:
        no = str(NoOfDip[FileNum])

    FileName = 'DATA/GURT-V325_25x01_Nex=' + no + '.out.txt'

    # *** Opening NECoutput file ***
    Data_File = open(FileName, "r")

    print ('\n **********************************************************************')
    print ('\n          Reading file: ', FileName)


    # *** Counting number of lines in the file ***
    num_lines = sum(1 for line in Data_File)
    print ('\n  Number of lines in the file = ', num_lines, '\n')
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

    # !!!!!!! possible error "loads_values" and others have single dimension !!!

    #for i in range (LoadNum):
    #    print ('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], '  Resistance = ', loads_values[i])


    # ***  Searching and reading the reference frequency *** (before f0)

    Data_File.seek(0) # return to the beginning of the file
    frequency_list = []
    for line in Data_File:
        if line.startswith('                                    FREQUENCY='):
            words_in_line = line.split()
            frequency_list.append(float(words_in_line[1]))

    #print ('\n  Central frerquency = ', frequency_list[0], ' MHz')
    #for i in range (1, len(frequency_list)):
    #    print ('  Frerquency = ', frequency_list[i], ' MHz')

    NoOfFrequencies = len(frequency_list) - 1
    print ('  Number of frequencies analyzed = ', NoOfFrequencies, '\n')

    for i in range (1): print (' **********************************************************************')


    Data_File.seek(0) # return to the beginning of the file

    # Finding the main (central) frequency data and skip it
    line = find_line_in_text_file (Data_File, '                                    FREQUENCY=', 0, line)


    for step in range (NoOfFrequencies):
        # Seek for current frequency
        line = find_line_in_text_file (Data_File, 'FREQUENCY=', 36, line)

        # *** Searching antenna input parameters block at the frequency ***
        line = find_line_in_text_file (Data_File, '- - - ANTENNA INPUT PARAMETERS - - -', 42, line)

        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = find_line_in_text_file (Data_File, 'IMPEDANCE', 64, line)
        line = Data_File.readline()   # Skip one line

        line = Data_File.readline()
        IFedDipReal[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[36 : 48])
        IFedDipImag[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[48 : 60])
        Rinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[60 : 72])
        Xinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1] = float(line[72 : 84])

        #print ('\n\n   For frequency = ', frequency_list[step+1], ' MHz')
        #print ('  ---------------------------- \n')
        #print ('  Self full impedance of the dipole = ', Rinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1], Xinp[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1], ' Ohm' )
        #print ('  Current in fed point of active dipole = ', IFedDipReal[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1],  IFedDipImag[step, NoOfDip[FileNum]-1, NoOfDip[FileNum]-1],' Amp. \n')


        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Searching and reading the currennt in impedance load !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # *** Seek for currents in loads ***
        line = find_line_in_text_file (Data_File, '- - - CURRENTS AND LOCATION - - -', 29, line)
        for i in range (6): line = Data_File.readline()   # Skip several lines


        for j in range (LoadNum):  # loop by number of loads
            intTemp = 0
            while (intTemp != wires_with_loads[j]):
                line = Data_File.readline()
                intTemp = int(line[6:12])
            if (segments_with_loads[j] < 1):    # if the number of segment = 1 we need to return to the begining of the line
                stop
            if (segments_with_loads[j] > 1):	# if the number of segment > 2 we need to skip lines
                for i in range (segments_with_loads[j] - 1): line = Data_File.readline()   # Skip several lines


            DipNum = int((wires_with_loads[j]) / NoOfWiresPerDipole) + 1   # Calculation of the dipole number   !!!!


            # *** Reading currents ***
            IFedDipReal[step, NoOfDip[FileNum]-1, DipNum-1] = float(line[49 : 60])
            IFedDipImag[step, NoOfDip[FileNum]-1, DipNum-1] = float(line[61 : 72])

            #print ('  Dip No', DipNum, ' Wire No', wires_with_loads[j], ' Segm No', segments_with_loads[j])
            #print ('  Current in the load No ', j+1, ' = ', IFedDipReal[step, NoOfDip[FileNum]-1, DipNum-1], IFedDipImag[step, NoOfDip[FileNum]-1, DipNum-1], ' Amp' )

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!!  Reading of E values in radiation patterns  !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        line = find_line_in_text_file (Data_File, '- - - RADIATION PATTERNS - - -', 48, line)
        line = Data_File.readline()
        DistanceRP[FileNum, step] = float(line[60 : 74])
        #print ('  Distance of RP calculations = ', DistanceRP[FileNum, step], ' m. \n')

        for i in range (6): line = Data_File.readline()   # Skip several lines of table header

        # *** loop by lines to read e-field amplitudes and phases for all theta and phi ***

        for k in range (65341):
            line = Data_File.readline()
            theta1 = float(line[0 : 8])
            phi1 = float(line[9 : 17])
            Etheta = float(line[76 : 87])
            ang1 = float(line[89 : 97])
            Ephi = float(line[100 : 112])
            ang2 = float(line[112 : 120])

            thetaInt = (int(theta1))+90
            phiInt = int(phi1)
            ETHcmplx[step,NoOfDip[FileNum]-1,thetaInt,phiInt] = np.complex((Etheta*np.cos(ang1*pi/180)),(Etheta*np.sin((ang1*pi/180))))
            EPHcmplx[step,NoOfDip[FileNum]-1,thetaInt,phiInt] = np.complex((Ephi*np.cos(ang2*pi/180)),(Ephi*np.sin((ang2*pi/180))))

    Data_File.close()

IFedDip = IFedDipReal + 1j * IFedDipImag
Zinp = Rinp + 1j * Xinp
del IFedDipReal, IFedDipImag, Rinp, Xinp

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!              Processing data and calculating results                     !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n    Performing calculations...  \n\n')


for step in range (NoOfFrequencies):     # Loop by frequencies

	# *** Mirroring of currents for cases where dipols were not excited ***

    for i in range (len(MirrorInto)):
        for j in range (ArrayInputNum):
            IFedDip[step, MirrorInto[i]-1, j] = (IFedDip[step, MirrorFrom[i]-1, ArrayInputNum-1-j]) # !!!  ArrayInputNum+1


	# *** Calculations of Impedance matrices ***

    for i in range (NECfileNum):           # Loop by files
        for k in range (ArrayInputNum):    # Loop by dipols
            if k != (NoOfDip[i]-1):
                Zinp[step, NoOfDip[i]-1, k] = - (IFedDip[step, NoOfDip[i]-1, k] * loads_values[i] / IFedDip[step, NoOfDip[i]-1, NoOfDip[i]-1]) # !!!-!!!


    # Mirroring the Impedances for dipols that were not excited

    for i in range (len(MirrorInto)):     # Loop by active dipols below array diagonal
      for j in range (ArrayInputNum):         # Loop by dipols
        Zinp[step, MirrorInto[i]-1, ArrayInputNum - 1 - j] = Zinp[step, MirrorFrom[i]-1, j]  # ArrayInputNum + 1


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!  Mirroring of RP for noncalculated dipols  !!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for i in range (len(MirrorInto)):   # loop by dipoles
        for t in range (181):           # loop by theta
            for p in range (180):       # loop by phi
                ETHcmplx[step, MirrorInto[i]-1, t, p+180] = - ETHcmplx[step, MirrorFrom[i]-1, t, p]
                EPHcmplx[step, MirrorInto[i]-1, t, p+180] = - EPHcmplx[step, MirrorFrom[i]-1, t, p]
            for p in range (180, 361):
                ETHcmplx[step, MirrorInto[i]-1, t, p-180] = - ETHcmplx[step, MirrorFrom[i]-1, t, p]
                EPHcmplx[step, MirrorInto[i]-1, t, p-180] = - EPHcmplx[step, MirrorFrom[i]-1, t, p]



    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!         Integrals calculation              !!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i in range (ArrayInputNum):       # Loop by first dipole
        for j in range (ArrayInputNum):   # Loop by second dipole
            if (i >= j):

                # Forming the integration elements for Theta angle
                for t in range (91):
                    for p in range (361):
                        Points[t, p] = ETHcmplx[step, i, t+90, p] * np.conj(ETHcmplx[step, j, t+90, p]) * (np.sin(t*pi/180.) + 1j*0.0)

                #y = np.linspace(1, 1, 91)
                #x = np.linspace(1, 1, 361)
                #Integral1 = simps(simps(Points, x), y)
                Integral1 = intergration_trap_2d (0, 90, 0, 360, 90, 360, Points)

                # Forming the integration elements for Phi component
                for t in range (91):
                    for p in range (361):
                        Points[t, p] = EPHcmplx[step, i, t+90, p] * np.conj(EPHcmplx[step, j, t+90, p]) * (np.sin(t*pi/180.) + 1j*0.0)

                # Calculating of integrals by Phi angle

                #y = np.linspace(1, 1, 91)
                #x = np.linspace(1, 1, 361)
                #Integral2 = simps(simps(Points, x), y)
                Integral2 = intergration_trap_2d (0, 90, 0, 360, 90, 360, Points)

                # Calculating the elements of radiation impedance matrix !!! Distance[1] !!!

                RSigm[step, i, j] = ((Integral1 + Integral2) * np.power(DistanceRP[0, 0], 2)) / (120.0 * pi * (IFedDip[step, i, i] * np.conj(IFedDip[step, j, j])))
                RSigm[step, j, i] = np.conj(RSigm[step, i, j])


    # *** Calculations of self radiation efficiencies of dipols ***
    for i in range (ArrayInputNum):
        Efficiency[step, i] = (np.real(RSigm[step,i,i]) / np.real(Zinp[step,i,i])) * 100

    currentTime = time.strftime("%H:%M:%S")
    print ('\n    Done for frequency = ', frequency_list[step+1], ' MHz at: ', currentTime, ' ')

# End of loop by frequencies


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!             Displaying data and results on the screen                    !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


for step in range (NoOfFrequencies):	 # Loop by frequencies for displaying results on the screen

    print ('\n\n  For frequency = ', frequency_list[step+1], ' MHz')
    print ('  --------------------------')
    print ('    Zinp (1, 1) = ', Zinp[step, 0, 0])
    print ('   Rsigm (1, 1) = ', RSigm[step, 0, 0])
    print ('\n  Self radiation efficiency of the first dipole = ', round(Efficiency[step, i],3), ' %')



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!              Writing data and results to output files                    !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n     Writing the results to output files... \n\n\n')


# *** Creating files for results ***
CurrentMatrixFile = open("Results/Matrices of Currents.txt", "w")              # 24
ImpedanceMatrixFile = open("Results/Matrices of Impedances.txt", "w")          # 26
RadImpMatrixFile = open("Results/Matrices of Radiation Impedances.txt", "w")   # 27

for step in range (NoOfFrequencies):	 # Loop by frequencies for results files writing


    # *** Matrix of complex currents in dipoles.txt ***
    CurrentMatrixFile.write('\n\n    *** Matrix of currents in fed points of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step+1])
    data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
    CurrentMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
    CurrentMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(ArrayInputNum):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(IFedDip[step, k, i]), np.imag(IFedDip[step, k, i])) for i in range(ArrayInputNum)]
        CurrentMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')


    # *** Writing impedances to file ***
    ImpedanceMatrixFile.write('\n\n    *** Matrix of impedances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step+1])
    data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
    ImpedanceMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
    ImpedanceMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(ArrayInputNum):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(Zinp[step, k, i]), np.imag(Zinp[step, k, i])) for i in range(ArrayInputNum)]
        ImpedanceMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')


    # *** Writing radiation resistances to file ***
    RadImpMatrixFile.write('\n\n    *** Matrix of radiation resistances of dipoles f = %6.3f MHz ***  \n\n' % frequency_list[step+1])
    data = ['                                  {:2d}'.format(i+2) for i in range(ArrayInputNum-1)]
    RadImpMatrixFile.write('Dip No                1' + ''.join(data)+ '   \n\n')
    data = ['          Re                Im      '.format() for i in range(ArrayInputNum-1)]
    RadImpMatrixFile.write('Dip No        Re                Im      ' + ''.join(data)+ '   \n')
    for k in range(ArrayInputNum):
        data = ['   {:+14.8e}   {:+14.8e}'.format(np.real(RSigm[step, k, i]), np.imag(RSigm[step, k, i])) for i in range(ArrayInputNum)]
        RadImpMatrixFile.write(' {:2d}'.format(k+1) + ''.join(data)+ '   \n')



    # *** Writing currents on dipols inputs to files for MatCAD
    FileName = ('Results/Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
    CurrentFile = open(FileName, "w")  # 6
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()


    # *** Writing normalized currents on dipols inputs to files for MatCAD ***
    FileName = ('Results/Normalized Currents on dipols inputs f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
    CurrentFile = open(FileName, "w")  # 6
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.real(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.imag(IFedDip[step, i, k] / IFedDip[step, i, i])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()


    # *** Writing impedances of dipols to files for MatCAD ***
    FileName = ('Results/Impedances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
    CurrentFile = open(FileName, "w")  # 6
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.real(Zinp[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.imag(Zinp[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()


    # *** Writing radiation resistances of dipols to files for MatCAD ***
    FileName = ('Results/Radiation resistances of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
    CurrentFile = open(FileName, "w")  # 6
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.real(RSigm[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(np.imag(RSigm[step, i, k])) for i in range(ArrayInputNum)]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()



    # *** Writing radiation efficiencies of dipols to files for MatCAD ***
    FileName = ('Results/Radiation efficiencies of dipols f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
    CurrentFile = open(FileName, "w")  # 6
    for k in range (ArrayInputNum):
        data = ['   {:+14.8e}'.format(Efficiency[step, i])]
        CurrentFile.write(' ' + ''.join(data)+ '  \n')
    CurrentFile.close()


    # *** Writing normalized Etheta RP to files for MatCAD ***
    FileName = ('Results/Normalized to current Etheta RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
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
    FileName = ('Results/Normalized to current Ephi RP f = %5.2f MHz for MatCAD.txt' % frequency_list[step+1])
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


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!                            F I G U R E S                                 !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# *** Figure of current absolute value ***
for step in range (len(frequency_list)-1):

    data = np.zeros((ArrayInputNum))
    for i in range (ArrayInputNum):
        data[i] = np.absolute(IFedDip[step, i, i])
    x_values = np.linspace(1, ArrayInputNum, num = ArrayInputNum)
    plt.figure()
    rc('font', weight='normal')
    plt.plot(x_values, data[:], 'bo', label = 'Abs(IFedDip)')
    plt.xlabel('Number of dipole')
    plt.ylabel('Current at dipole output, A')
    #plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
    plt.title('NEC modeling results', fontsize=8)
    plt.grid(b = True, which = 'both', color = '0.00',linestyle = '--')
    plt.legend(loc = 'lower center', fontsize = 10)
    pylab.savefig("Results/Dipole currents (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    # plt.show()
    plt.close('all')


# *** Figure of current absolute value matrix nondiagonal map ***
for step in range (len(frequency_list)-1):
    data = np.zeros((ArrayInputNum, ArrayInputNum))
    data = np.absolute(IFedDip[step, :, :])
    for i in range (ArrayInputNum):
        data[i,i] = 0

    MakeFigureColorMap(data, 'Number of dipole','Number of dipole',
                       'Matrix of currents at inputs of antenna array at %5.2f MHz'% frequency_list[step+1],
                       'NEC modeling results',
                       "Results/Dipole currents matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1])


# *** Figure of input impedance absolute value at the main diagonal ***
for step in range (len(frequency_list)-1):
    data = np.zeros((ArrayInputNum))
    for i in range (ArrayInputNum):
        data[i] = np.absolute(Zinp[step, i, i])
    x_values = np.linspace(1, ArrayInputNum, num = ArrayInputNum)
    plt.figure()
    rc('font', weight='normal')
    #plt.plot(x_values, data[:], color = 'violet', linestyle = '-', linewidth = '1.00', label = 'Abs(Zinp)')
    plt.bar(x_values, data[:], bottom = np.min(data[:])*0.99)
    plt.xlabel('Number of dipole')
    plt.ylabel('|Z|, Ohm')
    #plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
    plt.title('NEC modeling results', fontsize=7, x = 0.46, y = 1.005)
    plt.grid(b = True, which = 'both', color = '0.00',linestyle = '--')
    #plt.legend(loc = 'lower center', fontsize = 10)
    pylab.savefig("Results/Self full impedances (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    plt.close('all')


# *** Figure of impedance absolute value matrix nondiagonal map ***
for step in range (len(frequency_list)-1):
    data = np.zeros((ArrayInputNum, ArrayInputNum))
    data = np.absolute(Zinp[step, :, :])
    for i in range (ArrayInputNum):
        data[i,i] = 0

    MakeFigureColorMap(data, 'Number of dipole','Number of dipole',
                       'Matrix of impedances at inputs of antenna array at %5.2f MHz'% frequency_list[step+1],
                       'NEC modeling results',
                       "Results/Dipole impedances matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1])



# *** Figure of radiation resistance absolute value matrix nondiagonal map ***
for step in range (len(frequency_list)-1):
    data = np.zeros((ArrayInputNum, ArrayInputNum))
    data = np.absolute(RSigm[step, :, :])
    for i in range (ArrayInputNum):
        data[i,i] = 0

    MakeFigureColorMap(data, 'Number of dipole','Number of dipole',
                       'Matrix of mutual radiation resistances of dipoles in antenna array at %5.2f MHz'% frequency_list[step+1],
                       'NEC modeling results',
                       "Results/Radiation resistance matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1])


endTime = time.time()
print ('\n\n\n  The program execution lasted for ', round((endTime - startTime), 2), 'seconds (',
                                                round((endTime - startTime)/60, 2), 'min. ) \n')
print ('\n           *** Program ' + Software_name + ' has finished! *** \n\n\n')
