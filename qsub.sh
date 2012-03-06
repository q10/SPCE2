# Declare a name for this job
# -N "SPCE_"$SIM_NAME

# Request that email is sent when job has started, ended, aborted, or suspended
#$ -m beas
#$ -M bm@berkeley.edu

# Run program and save logfiles in the current working directory
#$ -cwd

# Run the program
if [ "x" == "x$WINDOW_LOWER_BOUND" ] || [ "x" == "x$WINDOW_UPPER_BOUND" ] || [ "x" == "x$SIM_NAME" ]; then
    echo "window bounds or name not set"
    ./SPCE
else
    echo "WINDOW_LOWER_BOUND: $WINDOW_LOWER_BOUND"
    echo "WINDOW_UPPER_BOUND: $WINDOW_UPPER_BOUND"
    ./SPCE $WINDOW_LOWER_BOUND $WINDOW_UPPER_BOUND $SIM_NAME
fi

exit 0
