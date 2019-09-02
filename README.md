# "antenna_modelling_data_processing" 
Processing of data from the results of NEC engine antenna modeling. For calculations usually 4NEC2X software was used. 
This soft was available from: http://home.ict.nl/~arivoors/. Now this link does not work.

##        Array Impedance Matrix Calculator (array_imppedance_matrix_calc)
                 
Program reads the result files of NEC modeling (.nec) where each (of non mirrored) dipole of the antenna array 
is excited one by one (one in each file), finds frequencies of analysis,  excitation sources, loads and their values, 
currents in sources and loads, antenna patterns.
Then the program calculates mutual impedances, radiation resistances, and store them to text files (.txt).
The initial version was created for AA of 25 GURT dipols of 1 polarization of incoming waves, 
number of segments per dipole is 325 then it was  extended for other possibilities, 
but still you need to analyze carefully the results.

To save some time of modelling the symmetry of antenna arrays is used. It means that for some of array elements you
do not need to model their excitation with NEC because their radiation pattern is the same for mirror elements
except for you need to "mirror" it in respect to symmetry point. It is implemented in the code. But you have to
specify which number of dipoles were modeled with excitation, which of them can be "mirrored" and which of non
modeled dipoles will have "mirrored" patterns.
