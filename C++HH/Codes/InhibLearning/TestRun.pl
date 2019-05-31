#!/usr/bin/perl
#NoLearning Run
$cmd="./run_slow_compute 1 5000 2.7 2.5 .1 .15 5 2 0 0 BiDirecNoLearn.dat .1 .1 0";
system($cmd);
#Learning Run
$cmd="./run_slow_compute 1 5000 2.7 2.5 .1 .15 5 2 0 1 BiDirecLearn.dat .1 .1 0";
system($cmd);