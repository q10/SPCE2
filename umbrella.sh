#!/bin/bash

mkdir -p tmp
make clean && make o3 && make mbar

MAIN_N="A"
N1=$MAIN_N"1"
N2=$MAIN_N"2"
N3=$MAIN_N"3"
N4=$MAIN_N"4"
N5=$MAIN_N"5"
N6=$MAIN_N"6"

qsub -v WINDOW_LOWER_BOUND=0.0,WINDOW_UPPER_BOUND=5.0,SIM_NAME=$N1 -N $N1 qsub.sh
qsub -v WINDOW_LOWER_BOUND=4.5,WINDOW_UPPER_BOUND=7.0,SIM_NAME=$N2 -N $N2 qsub.sh
qsub -v WINDOW_LOWER_BOUND=6.5,WINDOW_UPPER_BOUND=9.0,SIM_NAME=$N3 -N $N3 qsub.sh
qsub -v WINDOW_LOWER_BOUND=8.5,WINDOW_UPPER_BOUND=11.0,SIM_NAME=$N4 -N $N4 qsub.sh
qsub -v WINDOW_LOWER_BOUND=10.5,WINDOW_UPPER_BOUND=13.0,SIM_NAME=$N5 -N $N5 qsub.sh
qsub -v WINDOW_LOWER_BOUND=12.5,WINDOW_UPPER_BOUND=15.0,SIM_NAME=$N6 -N $N6 qsub.sh

while [ $(ls tmp/*$MAIN_N.jobdone | wc -l) != '6' ]; do
    sleep 3600
done

# renames file extensions of results for mbar
rename .ipair_dist .txt *.ipair_dist

NWINDOWS=5        # number of windows (files)
WINDOWSPACING=2.0 # spacing between the start of one window and the start of the next
WINDOWWIDTH=2.5  # width of a single window
MINIMUMH=5.0     # distance at which the lowest windows begins
HISTOGRAMFILE=umbrella_rdf # name of output file
NITERATIONS=400   # number of iterations to build the RDF
DATAFILE=window   # name of the input files (files start as %s1.txt, %s2.txt, etc)

./mbar $NWINDOWS $WINDOWSPACING $WINDOWWIDTH $MINIMUMH $HISTOGRAMFILE $NITERATIONS $DATAFILE
