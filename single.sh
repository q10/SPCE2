#!/bin/bash

MAIN_N=SPCEX
EMAIL=bm@berkeley.edu

# Make results folder if nonexistent, and check that it's empty
mkdir -p results
if [ "$(ls -A results)" ]; then
    echo "There are other files in results/, possibly from older simulation runs.  Please remove them or move them out and rerun script."
    exit 2
fi

# Build
make clean && make o3
qsub -v SIM_NAME=$MAIN_N,NUM_ION_PAIRS=1 -N $MAIN_N qsub.sh

# Wait for simulation to finish running on the grid
while [ $(ls results/*.jobdone | wc -l) -lt "1" ]; do
    sleep 900
done

# Plot output RDF on gnuplot
cd results && cat << EOF | gnuplot
set term gif large size 1600, 1200
set output "OO-RDF.gif"
plot "$MAIN_N.out" using 1:2 with lp title 'O-O RDF'
set output "AO-RDF.gif"
plot "$MAIN_N.out" using 1:3 with lp title 'A-O RDF'
set output "CO-RDF.gif"
plot "$MAIN_N.out" using 1:4 with lp title 'C-O RDF'
set output "AC-RDF.gif"
plot "$MAIN_N.out" using 1:5 with lp title 'A-C RDF'
set output "ALL-RDF.gif"
plot "$MAIN_N.out" using 1:2 with lp title 'O-O RDF', "$MAIN_N.out" using 1:3 with lp title 'A-O RDF', "$MAIN_N.out" using 1:4 with lp title 'C-O RDF', "$MAIN_N.out" using 1:5 with lp title 'A-C RDF'
EOF

echo "" | mail -s "Simulation \"$MAIN_N\" has finished" "$EMAIL"

exit 0
