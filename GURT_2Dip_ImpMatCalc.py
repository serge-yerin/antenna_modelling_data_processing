# Python3
Software_version = '2018.08.16'
#                   GURT_2Dip_ImpMatCalc
# *** GURT 2 Dipole Impedance Matrices Calculator   ***
# Program reads the result files of NEC modeling (.nec) where each of 2 dipoles
# are excited one by one (in each file), finds frequencies of analysis, 
# excitation sources, loads and their values, currents in sources and loads, antenna patterns.
# Then the program calculates mutual impedances, radiation resistances, and store them to 
# text files (.txt).

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
#from scipy.integrate import simps
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
    
def IntegrTrap2x (a, b, c, d, m, n, Points):
    pi = 3.141593

    # a - the start integration limit by the first axis
    # b - the end integration limit by the first axis
    # c - the start integration limit by the second axis
    # d - the end integration limit by the second axis
    # n - number of points on 1 axis
    # m - number of points on 2 axis
    '''
    Sum = 0.0
    for i in range (m):
        for k in range (n):
            Sum = Sum + Points[i, k]
    '''
    Sum = np.sum(Points)
    
    IntegrResult = ((((b)-(a))*pi/180.0)/m) * ((((d)-(c))*pi/180.0)/n) * ((Points[a,c] + Points[b,d]) + Sum)
    return IntegrResult


#*************************************************************
#                      MAIN PROGRAM                          *
#*************************************************************

for i in range(6): print (' ')
print ('    ------------------------------------------------------------')
print ('    ***                                                      ***')
print ('    ***        NEC output calculator v.', Software_version, '          ***')
print ('    *** of self and mutual impedances, radiation resistances ***   (c) YeS 2018')
print ('    ***    of 2 dipoles with 1 source and multiple loads     ***')
print ('    ***            of GURT dipoles  (MultiFreq )             ***')
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
ArrayInputNum = 2 # ArrInpNum = Array inputs number
NECfileNum = 2    # TotFileNum = The number of NECout files

Rinp = np.full((10, ArrayInputNum, ArrayInputNum), 0.0)
Xinp = np.full((10, ArrayInputNum, ArrayInputNum), 0.0)
IFedDipReal = np.full((10, ArrayInputNum, ArrayInputNum), 0.0)
IFedDipImag = np.full((10, ArrayInputNum, ArrayInputNum), 0.0)  # IFedDip(frequency step, Active dipole, Passive dipole)
DistanceRP = np.full((15, 10), 0.0) # [FileNum, FreqStep]

Points = np.full((91, 361), 0.0+1j*0.0)     #[theta, phi] Matrix of E-field data for integral calculations
RSigm = np.full((10, ArrayInputNum, ArrayInputNum), 0.0+1j*0.0) #[FreqStep, Dip1, Dip2] Radiation impedance matrix
Efficiency = np.zeros((100, 100))

ETHcmplx = np.full((10, ArrayInputNum, 181, 361), 0.0+1j*0.0) # [FreqStep, DipNo, theta, phi] Theta components of E-field
EPHcmplx = np.full((10, ArrayInputNum, 181, 361), 0.0+1j*0.0) # [FreqStep, DipNo, theta, phi] Phi components of E-field


# *** Creating folder for results ***
if not os.path.exists('Results'):
    os.makedirs('Results')
    

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!                 Reading data from NEC output files                       !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for FileNum in range (NECfileNum):	# Main loop by NEC output files

    # *** Configuring the name of NECoutput file ***
    FileName = 'Data/GURT-391_1x2_cop_N ' + str(FileNum+1) + '.out'
  
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
    
    
    for i in range (LoadNum):
        print ('  Wire No ', wires_with_loads[i], '  Segm No', segments_with_loads[i], '  Resistance = ', loads_values[i])
    
        
    # ***  Searching and reading the reference frequency *** (before f0)
    
    Data_File.seek(0) # return to the beginning of the file
    frequency_list = []
    for line in Data_File:
        if line.startswith('                                    FREQUENCY='):
            words_in_line = line.split()
            frequency_list.append(float(words_in_line[1]))
            
    #print ('\n  Central frerquency = ', frequency_list[0], ' MHz')
    for i in range (1, len(frequency_list)):
        print ('  Frerquency = ', frequency_list[i], ' MHz')
 
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
        IFedDipReal[step, FileNum, FileNum] = float(line[36 : 48])
        IFedDipImag[step, FileNum, FileNum] = float(line[48 : 60])
        Rinp[step, FileNum, FileNum] = float(line[60 : 72])
        Xinp[step, FileNum, FileNum] = float(line[72 : 84])
      
        print ('\n\n   For frequency = ', frequency_list[step+1], ' MHz')
        print ('  ---------------------------- \n')
        print ('  Self full impedance of the dipole = ', Rinp[step, FileNum, FileNum], Xinp[step, FileNum, FileNum], ' Ohm' )
        print ('  Current in fed point of active dipole = ', IFedDipReal[step, FileNum, FileNum],  IFedDipImag[step, FileNum, FileNum],' Amp. \n')
        
        
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
            
            if FileNum == 0: k = 1
            if FileNum == 1: k = 0
           
            # *** Reading currents ***
            IFedDipReal[step, FileNum, k] = float(line[49 : 60])
            IFedDipImag[step, FileNum, k] = float(line[61 : 72])
            
            print ('  Dip No', FileNum+1, ' Wire No', wires_with_loads[j], ' Segm No', segments_with_loads[j])
            print ('  Current in the load No ', j+1, ' = ', IFedDipReal[step, FileNum, k], IFedDipImag[step, FileNum, k], ' Amp' )

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
            ETHcmplx[step,FileNum,thetaInt,phiInt] = np.complex((Etheta*np.cos(ang1*pi/180)),(Etheta*np.sin((ang1*pi/180)))) 
            EPHcmplx[step,FileNum,thetaInt,phiInt] = np.complex((Ephi*np.cos(ang2*pi/180)),(Ephi*np.sin((ang2*pi/180))))

    Data_File.close()

IFedDip = IFedDipReal + 1j * IFedDipImag
Zinp = Rinp + 1j * Xinp
del IFedDipReal, IFedDipImag, Rinp, Xinp
    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!              Processing data and calculating results                     !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print ('\n\n\n    Performing calculations...  \n\n')


for step in range (NoOfFrequencies):     # Loop by frequencies

	# *** Calculations of Impedance matrices ***
	
    for i in range (ArrayInputNum):        # Loop by dipoles
        for k in range (ArrayInputNum):    # Loop by dipoles
            if k != i:  
                Zinp[step, i, k] = - (IFedDip[step, i, k] * loads_values[0] / IFedDip[step, i, i]) # !!!-!!!

    #print ('    Calculations of Impedance Matrix done!')

    
 
 
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
    print ('    Zinp (2, 2) = ', Zinp[step, 1, 1])
    print ('    Zinp (1, 2) = ', Zinp[step, 0, 1])
    print ('    Zinp (2, 1) = ', Zinp[step, 1, 0], '\n')
    
    print ('   Rsigm (1, 1) = ', RSigm[step, 0, 0])
    print ('   Rsigm (2, 2) = ', RSigm[step, 1, 1])
    print ('   Rsigm (1, 2) = ', RSigm[step, 0, 1])
    print ('   Rsigm (2, 1) = ', RSigm[step, 1, 0])    
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
    
    
   
    
# *** Time consumption calculation (ending) ***    
endTime = time.time()    # Time of calculations      

print (' ')
print ('  The program execution lasted for ', round((endTime - startTime),2), 'seconds')
for i in range (0,2) : print (' ')
print ('    *** Program GURT_AA_ImpMatCalc reader has finished! ***')
for i in range (0,3) : print (' ')

