# Python3
Software_version = '2018.08.28'
#                   GURT_AA_ImpMatCalc
# *** GURT Antenna Array Impedance Matrices Calculator   ***
# Program reads the result files of NEC modeling (.nec) where each (of non mirrored) dipole
# of the antenna array is excited one by one (in each file), finds frequencies of analysis, 
# excitation sources, loads and their values, currents in sources and loads, antenna patterns.
# Then the program calculates mutual impedances, radiation resistances, and store them to 
# text files (.txt).

# The initial version was created for AA of 25 GURT dipols of 1 polarization of 
# incoming waves, number of segments per dipole is 325 (Multi frequencies) then it was 
# extended for other possibilities, but still you need to analyze carefully the results

#input("Press Enter to continue...") = pause command in Fortran

#*************************************************************
#                        PARAMETERS                          *
#*************************************************************
pi = 3.141593

#*************************************************************
#                   IMPORT LIBRARIES                         *
#*************************************************************
import os
import time
from scipy.integrate import simps
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pyplot import *


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

# Integration with trapezion rule
def IntegrTrap2x (a, b, c, d, m, n, Points):
    pi = 3.141593

    # a - the start integration limit by the first axis
    # b - the end integration limit by the first axis
    # c - the start integration limit by the second axis
    # d - the end integration limit by the second axis
    # n - number of points on 1 axis
    # m - number of points on 2 axis
 
    Sum = np.sum(Points)
    IntegrResult = ((((b)-(a))*pi/180.0)/m) * ((((d)-(c))*pi/180.0)/n) * ((Points[a,c] + Points[b,d]) + Sum)
    return IntegrResult

# Plotting color map figures
def MakeFigureColorMap(data, XLabel, YLabel, SupTitle, Title, FileName):
    plt.figure(figsize=(10, 6))
    rc('font', weight='normal', size = 8)
    plt.pcolor(data, cmap = 'jet')
    plt.axes().set_aspect('equal')
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.colorbar()
    plt.suptitle(SupTitle, fontsize=8, fontweight='bold')
    plt.title(Title, fontsize = 8)
    # plt.subplots_adjust(left=0.2, right=0.3)
    pylab.savefig(FileName, bbox_inches='tight', dpi = 160)
    plt.close('all')
    
    
#*************************************************************
#                      MAIN PROGRAM                          *
#*************************************************************

for i in range(6): print (' ')
print ('    ------------------------------------------------------------')
print ('    ***                                                      ***')
print ('    ***        NEC output calculator v.', Software_version, '          ***')
print ('    *** of self and mutual impedances, radiation resistances ***   (c) YeS 2018')
print ('    ***    of dipoles with one source and multiple loads     ***')
print ('    ***   of GURT antenna subarray  (MultiFreq and TwoPol)   ***')
print ('    ***                                                      ***')
print ('    ------------------------------------------------------------')
for i in range(3): print (' ')

# *** Time consumption calculation (beginning) ***
startTime = time.time()
currentTime = time.strftime("%H:%M:%S")
currentDate = time.strftime("%d.%m.%Y")
print ('  Today is ', currentDate, ' time is ', currentTime)
print (' ')


NoOfWiresPerDipole = 51.0
ArrayInputNum = 50 # ArrInpNum = Array inputs number (numder of dipols) GURT = 25 (1 pol), 50 (2 pol)
NECfileNum = 30    # The number of NECout files to analyze

NoOfDip = [1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14, 15, 19, 20, 25,
        26, 27, 28, 29, 30, 32, 33, 34, 35, 38, 39, 40, 44, 45, 50] # 30 dipols being excited
MirrorCurFrom = [2,  3,  4,  5,  8,  9,  10, 14, 15, 20] #, 27, 28, 29, 30, 33, 34, 35, 39, 40, 48] # Active dipols under the diagonal  Zerk1
MirrorCurInto = [24, 23, 22, 21, 18, 17, 16, 12, 11,  6] #, 49, 48, 47, 46, 43, 42, 41, 37, 36, 31] # Passive dipols above the diagonal Zerk2

MirrorPolFrom = [1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
MirrorPolInto = [46,41,36,31,26,47,42,37,32,27,48,43,38,33,28,49,44,39,34,29,50,45,40,35,30]


Rinp = np.full((10, 50, 50), 0.0) # [freq. step, Dip1, Dip2] Matrix of self and mutual input impedances (real part)
Xinp = np.full((10, 50, 50), 0.0) # [freq. step, Dip1, Dip2] Matrix of self and mutual input impedances (imag part)
IFedDipReal = np.full((10, 50, 50), 0.0) # IFedDip[freq. step, Active dipole, Passive dipole] current in dipole fed point
IFedDipImag = np.full((10, 50, 50), 0.0) # IFedDip[freq. step, Active dipole, Passive dipole] current in dipole fed point
DistanceRP = np.full((NECfileNum, 10), 0.0) # [FileNum, freq. step] distance of RP calculation

Points = np.full((91, 361), 0.0+1j*0.0)   # [theta, phi] Matrix of E-field data for integral calculations
RSigm = np.full((10, 50, 50), 0.0+1j*0.0) # [freq. step, Dip1, Dip2] Radiation impedance matrix
Efficiency = np.zeros((100, 100))         # Self radiation efficiency [freq. step, DipNum]

ETHcmplx = np.full((10, 50, 181, 361), 0.0+1j*0.0) # [freq. step, DipNo, theta, phi] Theta components of E-field
EPHcmplx = np.full((10, 50, 181, 361), 0.0+1j*0.0) # [freq. step, DipNo, theta, phi] Phi components of E-field


# *** Creating folder for results ***
if not os.path.exists('Results'):
    os.makedirs('Results')
    

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!                 Reading data from NEC output files                       !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for FileNum in range (NECfileNum):	# Main loop by NEC output files

    # *** Configuring the name of NECoutput file ***
    
    if FileNum < 15:
        FileName = 'Data/GURT-391_5x5_cop_N' + '{:2d}'.format(NoOfDip[FileNum]) + '.out.txt'
    if FileNum >= 15:
        FileName = 'Data/GURT-391_5x5_cross_N' + '{:2d}'.format(NoOfDip[FileNum]-25) + '.out.txt'
  
    # *** Opening NECoutput file ***
    Data_File = open(FileName, "r")

    print ('\n **********************************************************************')
    print ('\n          Reading file: ', FileName)

     
   
    # *** Counting number of lines in the file ***
    num_lines = sum(1 for line in Data_File)
    print ('\n  Number of lines in the file = ', num_lines, '\n')
    Data_File.seek(0)
    
    
    # ***  Searching the lines where loads are listed in cards *** 
    wires_with_loads = []       # numbers of wires which contain loads (load number)
    segments_with_loads = []    # numbers of segments which contain loads (load number)
    loads_values = []           # value of the load (load number)
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
    line = Find_line ('                                    FREQUENCY=', 0, line)
    

    for step in range (NoOfFrequencies):
        # Seek for current frequency
        line = Find_line ('FREQUENCY=', 36, line)
                
        # *** Searching antenna input parameters block at the frequency ***
        line = Find_line ('- - - ANTENNA INPUT PARAMETERS - - -', 42, line)
        
        # *** Searching and reading the self impedance and initial current at the frequency ***
        line = Find_line ('IMPEDANCE', 64, line)
        line = Data_File.readline()   # Skip one line
        
        line = Data_File.readline()
        
        # Read current and impedance only for first polarization (mirror later more correct)
        if FileNum < int(NECfileNum/2):
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
        line = Find_line ('- - - CURRENTS AND LOCATION - - -', 29, line)
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
            
            # Calculation of the number of passive dipole with load on the file and load number
            if FileNum < int(NECfileNum/2):
                DipNum = int(wires_with_loads[j] / NoOfWiresPerDipole) + 1   # Calculation of the dipole number   !!!!
            if FileNum >= int(NECfileNum/2):
                DipNum = int((wires_with_loads[j] - 802 + 52) / NoOfWiresPerDipole) + 26  # Calculation of the dipole number   !!!!
            
            # Calculation of the number of active dipole on the numeber of file
            if FileNum < int(NECfileNum/2):
                ActDipNo = NoOfDip[FileNum]-1
            if FileNum >= int(NECfileNum/2):
                ActDipNo = NoOfDip[FileNum-int(NECfileNum/2)]-1
            
            
            # *** Reading currents ***
            IFedDipReal[step, ActDipNo, DipNum-1] = float(line[49 : 60]) # NoOfDip[FileNum]-1
            IFedDipImag[step, ActDipNo, DipNum-1] = float(line[61 : 72]) # NoOfDip[FileNum]-1
            
            #print ('  Dip No', DipNum, ' Wire No', wires_with_loads[j], ' Segm No', segments_with_loads[j])
            #print ('  Current in the load No ', j+1, ' = ', IFedDipReal[step, NoOfDip[FileNum]-1, DipNum-1], IFedDipImag[step, NoOfDip[FileNum]-1, DipNum-1], ' Amp \n' )

            
        #input("Press Enter to continue...")
            
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!!  Reading of E values in radiation patterns  !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
        line = Find_line ('- - - RADIATION PATTERNS - - -', 48, line)
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

currentTime = time.strftime("%H:%M:%S")
print ('\n\n\n    Performing calculations started at: ', currentTime, '  \n\n')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!  Mirroring of RP for noncalculated dipols  !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#print ('   Mirroring of RPs ... ')


# *** Mirroring of RPs for dipoles of first polarization (that were not excited) ***
for i in range (len(MirrorCurInto)):
    for p in range (181):           # loop by phi
        ETHcmplx[:, MirrorCurInto[i]-1, :, p+180] = - ETHcmplx[:, MirrorCurFrom[i]-1, :, p]	
        EPHcmplx[:, MirrorCurInto[i]-1, :, p+180] = - EPHcmplx[:, MirrorCurFrom[i]-1, :, p]	
    for p in range (180, 361):
        ETHcmplx[:, MirrorCurInto[i]-1, :, p-180] = - ETHcmplx[:, MirrorCurFrom[i]-1, :, p]
        EPHcmplx[:, MirrorCurInto[i]-1, :, p-180] = - EPHcmplx[:, MirrorCurFrom[i]-1, :, p]

    
    
# *** Mirroring of RPs for dipoles of second polarization ***
for i in range (len(MirrorPolInto)):
    for p in range (91):           # loop by phi # 270
        ETHcmplx[:, MirrorPolInto[i]-1, :, p] = - ETHcmplx[:, MirrorPolFrom[i]-1, :, p+270]  #+90 
        EPHcmplx[:, MirrorPolInto[i]-1, :, p] = - EPHcmplx[:, MirrorPolFrom[i]-1, :, p+270]	
    for p in range (90, 361): # 270
        ETHcmplx[:, MirrorPolInto[i]-1, :, p] = - ETHcmplx[:, MirrorPolFrom[i]-1, :, p-90]   #-270
        EPHcmplx[:, MirrorPolInto[i]-1, :, p] = - EPHcmplx[:, MirrorPolFrom[i]-1, :, p-90]
          
          
          
currentTime = time.strftime("%H:%M:%S")  
print ('\n\n    RP were mirrored at: ', currentTime, '  \n\n')



# Filling the currents of perpendicular dipoles of the same antenna to zeros ***
for i in range (int(ArrayInputNum/2)):
    IFedDip[:, i, i + int(ArrayInputNum/2)] = 0 + 1j*0
    IFedDip[:, i + int(ArrayInputNum/2), i] = 0 + 1j*0

# *** Mirroring of currents for cases where dipols were not excited ***
for i in range (len(MirrorCurInto)):
    for j in range (int(ArrayInputNum/2)):
        IFedDip[:, MirrorCurInto[i]-1, j] = (IFedDip[:, MirrorCurFrom[i]-1, int(ArrayInputNum/2) - 1 - j]) 
        IFedDip[:, MirrorCurInto[i]-1, j + int(ArrayInputNum/2)] = (IFedDip[:, MirrorCurFrom[i]-1, ArrayInputNum - 1 - j]) 
     
# (*) *** Mirroring of excitation currents for dipols of second polarization (all were not excited) part 4 of the matrix ***
for i in range (int(ArrayInputNum/2)):
    for j in range (int(ArrayInputNum/2)):
        IFedDip[:, MirrorPolInto[i]-1,  MirrorPolInto[j]-1] = IFedDip[:, MirrorPolFrom[i]-1, MirrorPolFrom[j]-1]
    
# (*) *** Mirroring of mutual currents within both polarizations (part 3 of matrix) ***
for i in range (int(ArrayInputNum/2)):
    for j in range (int(ArrayInputNum/2)):
        IFedDip[:, i+int(ArrayInputNum/2), j] = IFedDip[:, int(ArrayInputNum/2) - 1 - i, ArrayInputNum - 1 - j]
     
print ('\n\n    Currents were mirrored at: ', currentTime, '  \n\n')





for step in range (NoOfFrequencies):     # Loop by frequencies


	# *** Calculations of Impedance matrices ***
    for i in range (ArrayInputNum):        # Loop by dipols
        for k in range (ArrayInputNum):    # Loop by dipols
            if k != i:
                if (IFedDip[step, i, i] == 0+1j*0): # Check for zero value in denominator
                    Zinp[step, i, k] = 0+1j*0
                else:
                    Zinp[step, i, k] = - (IFedDip[step, i, k] * loads_values[0] / IFedDip[step, i, i]) # !!!-!!!

    
    #print ('    Calculations of Impedance Matrix done!')
    
    
    # *** Mirroring the Impedances for dipols that were not excited (diagonal elements) ***
    
    for i in range (len(MirrorCurInto)):     # Loop by number of active dipols below array diagonal
        Zinp[step, MirrorCurInto[i]-1, MirrorCurInto[i]-1] = Zinp[step, MirrorCurFrom[i]-1, MirrorCurFrom[i]-1]  # ArrayInputNum + 1

    
    # *** Mirroring the Impedances for dipols of second polarization (main diagonal) ***
    for i in range (len(MirrorPolFrom)):     # Loop by number of dipols of second polarization
        Zinp[step, MirrorPolInto[i]-1, MirrorPolInto[i]-1] = Zinp[step, MirrorPolFrom[i]-1, MirrorPolFrom[i]-1]  # ArrayInputNum + 1
        
    #print ('   Mirroring of Impedance Matrix done!')
    
    

    
 
 
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

                # Calculating of integrals by Theta angle
                Integral1 = IntegrTrap2x (0, 90, 0, 360, 90, 360, Points)

                # Forming the integration elements for Phi component
                for t in range (91):
                    for p in range (361):
                        Points[t, p] = EPHcmplx[step, i, t+90, p] * np.conj(EPHcmplx[step, j, t+90, p]) * (np.sin(t*pi/180.) + 1j*0.0)
    
                # Calculating of integrals by Phi angle
                Integral2 = IntegrTrap2x (0, 90, 0, 360, 90, 360, Points)
            
                # Calculating the elements of radiation impedance matrix !!! Distance[1] !!!
	 
                RSigm[step, i, j] = ((Integral1 + Integral2) * np.power(DistanceRP[0, 0], 2)) / (120.0 * pi * (IFedDip[step, i, i] * np.conj(IFedDip[step, j, j])))
                RSigm[step, j, i] = np.conj(RSigm[step, i, j])

    
    # *** Calculations of self radiation efficiencies of dipols ***
    for i in range (ArrayInputNum):
        Efficiency[step, i] = (np.real(RSigm[step,i,i]) / np.real(Zinp[step,i,i])) * 100
    
    
    currentTime = time.strftime("%H:%M:%S")
    print ('\n    Done for frequency = ', frequency_list[step+1], ' MHz at: ', currentTime, ' \n')

# End of loop by frequencies


    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!             Displaying data and results on the screen                    !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


for step in range (NoOfFrequencies):	 # Loop by frequencies for displaying results on the screen

    print ('\n\n  For frequency = ', frequency_list[step+1], ' MHz') 
    print ('  --------------------------')
    print ('    Zinp (1, 1) = ', Zinp[step, 0, 0])
    print ('   Rsigm (1, 1) = ', RSigm[step, 0, 0]) 
    print ('\n  Self radiation efficiency of the first dipole = ', round(Efficiency[step, 0], 3), ' %')



    
    

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
    

    '''
    # These files are needed for check of pattern mirroring
    
    if step == 4:
        # *** Writing Etheta RP to file for check ***
        FileName = ('Results/Check Etheta RP f = %5.2f MHz.txt' % frequency_list[step+1])
        CurrentFile = open(FileName, "w")  # 6
        data = [' {:+14.8e}'.format(np.real(ETHcmplx[step,0,45+90,p])) for p in range(361)]
        CurrentFile.write('' + ''.join(data)+ ' \n')
        data = [' {:+14.8e}'.format(np.real(ETHcmplx[step,45,45+90,p])) for p in range(361)]
        CurrentFile.write('' + ''.join(data)+ ' \n')
        CurrentFile.close()
    
    
    if step == 4:
        # *** Writing Etheta RP to file for check ***
        FileName = ('Results/Check Etheta RP f = %5.2f MHz dipoles 1 and 46.txt' % frequency_list[step+1])
        CurrentFile = open(FileName, "w")  # 6
        for p in range (271):   
            data = [' {:+14.8e} {:+14.8e}'.format(np.real(ETHcmplx[step,0,45+90,p+90]), np.real(ETHcmplx[step,45,45+90,p]))] 
            CurrentFile.write('' + ''.join(data)+ ' \n')
        for p in range (270, 361):   
            data = [' {:+14.8e} {:+14.8e}'.format(np.real(ETHcmplx[step,0,45+90,p-270]), np.real(ETHcmplx[step,45,45+90,p]))] 
            CurrentFile.write('' + ''.join(data)+ ' \n')
        
        CurrentFile.close()
        
    if step == 4:
        # *** Writing Etheta RP to file for check ***
        FileName = ('Results/Check Etheta RP f = %5.2f MHz dipoles 2 and 24.txt' % frequency_list[step+1])
        CurrentFile = open(FileName, "w")  # 6
        for p in range (181):   
            data = [' {:+14.8e} {:+14.8e}'.format(np.real(ETHcmplx[step,1,45+90,p+180]), np.real(ETHcmplx[step,23,45+90,p]))] 
            CurrentFile.write('' + ''.join(data)+ ' \n')
        for p in range (181, 361):   
            data = [' {:+14.8e} {:+14.8e}'.format(np.real(ETHcmplx[step,1,45+90,p-180]), np.real(ETHcmplx[step,23,45+90,p]))] 
            CurrentFile.write('' + ''.join(data)+ ' \n')
        
        CurrentFile.close()

    '''
   
    
    
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
    
CurrentMatrixFile.close()
ImpedanceMatrixFile.close()
RadImpMatrixFile.close()










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

'''
# *** Figure of current absolute value matrix diagonal map ***
for step in range (len(frequency_list)-1):
    plt.figure()
    rc('font', weight='normal')
    plt.pcolor(np.absolute(IFedDip[step, :, :]), cmap = 'jet') # imshow
    plt.xlabel('Number of dipole')
    plt.ylabel('Number of dipole')
    #plt.suptitle(SupTitle, fontsize=10, fontweight='bold')
    plt.title('NEC modeling results', fontsize=7, x = 0.46, y = 1.005)
    #plt.grid(b = True, which = 'both', color = '0.00',linestyle = '--')
    #plt.legend(loc = 'lower center', fontsize = 10)
    pylab.savefig("Results/Dipole currents matrix diagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    # plt.show()
    plt.close('all')
'''    
    
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
    
    

    
    '''
    plt.figure(figsize=(10, 6))
    rc('font', weight='normal', size = 8)
    plt.pcolor(data, cmap = 'jet')
    plt.xlabel('Number of dipole')
    plt.ylabel('Number of dipole')
    plt.colorbar(pad=0.005)
    plt.suptitle('Matrix of currents at inputs of antenna array at %5.2f MHz'% frequency_list[step+1], fontsize=8, fontweight='bold')
    plt.title('NEC modeling results', fontsize = 8)  #, x = 0.46, y = 1.005
    pylab.savefig("Results/Dipole currents matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    plt.close('all')
    '''


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
        
    '''   
    plt.figure(figsize=(12, 6))
    rc('font', weight='normal', size = 8)
    plt.pcolor(data, cmap = 'jet') 
    plt.xlabel('Number of dipole')
    plt.ylabel('Number of dipole')
    plt.colorbar()
    plt.suptitle('Matrix of impedances at inputs of antenna array at %5.2f MHz'% frequency_list[step+1], fontsize=10, fontweight='bold')
    plt.title('NEC modeling results', fontsize=10, x = 0.46, y = 1.005)
    pylab.savefig("Results/Dipole impedances matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    plt.close('all')  
    '''
    
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
       
    '''
    plt.figure(figsize=(12, 6))
    rc('font', weight='normal', size = 8)
    plt.pcolor(data, cmap = 'jet') # imshow
    plt.xlabel('Number of dipole')
    plt.ylabel('Number of dipole')
    plt.colorbar()
    plt.axes().set_aspect('equal') # !!!
    plt.suptitle('Matrix of mutual radiation resistances of dipoles in antenna array at %5.2f MHz'% frequency_list[step+1], fontsize=10, fontweight='bold')
    plt.title('NEC modeling results', fontsize=10, x = 0.46, y = 1.005)
    pylab.savefig("Results/Radiation resistance matrix nondiagonal (abs) f = %5.2f MHz.png" % frequency_list[step+1], bbox_inches='tight', dpi = 160)
    plt.close('all')   
    '''
    
    
    
    
# *** Time consumption calculation (ending) ***    
endTime = time.time()    # Time of calculations      

print (' ')
print ('  The program execution lasted for ', round((endTime - startTime),2), 'seconds')
for i in range (0,2) : print (' ')
print ('    *** Program GURT_AA_ImpMatCalc reader has finished! ***')
for i in range (0,3) : print (' ')

