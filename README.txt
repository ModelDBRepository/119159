This is the readme for the models associated with the paper

Talathi SS, Hwang DU, Ditto WL (2008) Spike timing dependent
plasticity promotes synchrony of inhibitory networks in the presence
of heterogeneity. J Comput Neurosci 25:262-81

Usage:

The Main Code is ./C++HH/Codes/InhibLearning/Test_Slow_Syn.cc

The accompanying make file Make_Inhib_Slow.f will allow to compile the code
to generate the executable run_slow_compute.

Example runs are provided in the perlscript TestRun.pl which generates the
data files, used to produce the two sample figures Learning.pdf and
NoLearning.pdf

Compile the make file on terminal through commannd

$ make -f Make_Inhib_Slow.f

It will generate the executable run_slow_compute

Run the perl script TestRun.pl.

It will generate two data files BiDirecLearn.dat and BiDirecNoLearn.dat

The voltage trace of the two neurons is the 2nd and 3rd column of
these dat files.

The Figures Learning.pdf and NoLearning.pdf represent the plotted
Figures from the data files above.

-Sachin S Talathi
