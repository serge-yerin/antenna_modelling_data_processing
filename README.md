# "antenna_modelling_data_processing" repository
Processing of data from the results of NEC engine antenna modeling. For
calculations usually 4NEC2X software was used.
This soft was available from: http://home.ict.nl/~arivoors/.
Now this link does not work.
This code was used for modelling of Giant Urainian Radio Telescope active
antenna and subarray. The mathematical model of the active antenna and
the subarray are presented in papers:
https://ieeexplore.ieee.org/abstract/document/7987810
Tokarsky, P.L., Konovalenko, A.A. and Yerin, S.N., 2017. Sensitivity of an
active antenna array element for the low-frequency radio telescope GURT.
IEEE Transactions on Antennas and Propagation, 65(9), pp.4636-4644.
https://ieeexplore.ieee.org/abstract/document/8782153
Tokarsky, P.L., Konovalenko, A.A., Yerin, S.N. and Bubnov, I.N., 2019.
An Active Antenna Subarray for the Low-Frequency Radio Telescope GURT–Part I:
Design and Theoretical Model. IEEE Transactions on Antennas and Propagation.
https://ieeexplore.ieee.org/document/8809883
Tokarsky, P.L., Konovalenko, A.A., Yerin, S.N. and Bubnov, I.N., 2019.
An Active Antenna Subarray for the Low-Frequency Radio Telescope GURT–Part II:
Numerical Analysis and Experiment. IEEE Transactions on Antennas and Propagation.


##        Array Impedance Matrix Calculator (array_imppedance_matrix_calc.py)

Program reads the result files of NEC modeling (.nec) where each (of non
mirrored) dipole of the antenna array is excited one by one (one in each
file), finds frequencies of analysis,  excitation sources, loads and their
values, currents in sources and loads, antenna patterns.
Then the program calculates mutual impedances, radiation resistances,
and store them to text files (.txt).
The initial version was created for AA of 25 GURT dipols of 1 polarization
of incoming waves, number of segments per dipole is 325 then it was  extended
for other possibilities, but still you need to analyze carefully the results.

To save some time of modelling the symmetry of antenna arrays is used. It
means that for some of array elements you do not need to model their
excitation with NEC because their radiation pattern is the same for mirror
elements except for you need to "mirror" it in respect to symmetry point.
It is implemented in the code. But you have to specify which number of
dipoles were modeled with excitation, which of them can be "mirrored" and
which of non modeled dipoles will have "mirrored" patterns.

Program was translated from FOTRAN in a simple way, so do not expect it to be
clear, efficient and fast. it just works.


##        Array Impedance Matrix Calculator (array_imppedance_matrix_calc_Ort.py)
is a previous version of the code adjusted to model antenna array of 25
crossed dipoles (50 dipoles in total). Needs check.
