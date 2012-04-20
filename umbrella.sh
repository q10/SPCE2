#!/bin/bash

MAIN_N="A"
NWINDOWS=5
WINDOWSPACING=1.0
WINDOWWIDTH=1.5
MINIMUMH=4.0

# Make results folder if nonexistent, and check that it's empty
mkdir -p results
if [ "$(ls -A results)" ]; then
    echo "There are other files in results/, possibly from older simulation runs.  Please remove them or move them out and rerun script."
    exit 2
fi

# Build
make clean && make o3 && make mbar

# Submit the jobs
lower=$MINIMUMH
for i in $(seq 1 $NWINDOWS); do
    upper=$(echo "scale=10;$lower + $WINDOWWIDTH" | bc)
    qsub -v WINDOW_LOWER_BOUND=$lower,WINDOW_UPPER_BOUND=$upper,SIM_NAME=$MAIN_N$i -N $MAIN_N$i qsub.sh
    lower=$(echo "scale=10;$lower + $WINDOWSPACING" | bc)
    sleep 10 # so that jobs can be more randomly distributed in the grid
done


# Wait for umbrella windows to finish running on the grid
while [ $(ls results/*.jobdone | wc -l) -lt $NWINDOWS ]; do
    sleep 10
done

# Move mbar into results folder and rename file extensions of results for mbar
cp mbar results/mbar && cd results 
rename .ipair_dist .txt *.ipair_dist


# Run mbar
#NWINDOWS        # number of windows (files)
#WINDOWSPACING=2.0 # spacing between the start of one window and the start of the next
#WINDOWWIDTH=2.5  # width of a single window
#MINIMUMH=5.0     # distance at which the lowest windows begins
HISTOGRAMFILE=$MAIN_N"_umbrella_rdf" # name of output file
NITERATIONS=1000   # number of iterations to build the RDF
#DATAFILE=$MAIN_N   # name of the input files (files start as %s1.txt, %s2.txt, etc)
./mbar $NWINDOWS $WINDOWSPACING $WINDOWWIDTH $MINIMUMH $HISTOGRAMFILE $NITERATIONS $MAIN_N
rm mbar


# Plot output RDF on gnuplot
cat << EOF | gnuplot
set term gif large size 1600, 1200
set output "$HISTOGRAMFILE.gif"
plot "$HISTOGRAMFILE" using 1:(\$2/\$1**2) with lp
EOF

exit 0
