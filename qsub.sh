# Declare a name for this job
# -N "SPCE_"$SIM_NAME

# Request that email is sent when job has started, ended, aborted, or suspended
#$ -m beas
#$ -M bm@berkeley.edu

# Run program and save logfiles in the current working directory
#$ -cwd


final_path=$HOME/SPCE2
output_path=/scratch/bm/$JOB_ID
[[ ! -e $output_path ]]  && mkdir -p $output_path
[[ ! -e $final_path ]] && mkdir -p $final_path

echo This job is being run on `hostname`

cp SPCE $output_path/
cd $output_path


# Run the program
if [ "x" == "x$WINDOW_LOWER_BOUND" ] || [ "x" == "x$WINDOW_UPPER_BOUND" ] || [ "x" == "x$SIM_NAME" ]; then
    echo "window bounds or name not set"
    ./SPCE
else
    echo "WINDOW_LOWER_BOUND: $WINDOW_LOWER_BOUND"
    echo "WINDOW_UPPER_BOUND: $WINDOW_UPPER_BOUND"
    ./SPCE $WINDOW_LOWER_BOUND $WINDOW_UPPER_BOUND $SIM_NAME
fi

# copy the tons of data back to local disk
rm -rf SPCE
cd $final_path && mkdir -p results
mv $output_path/* results/

# signal master script that it is done
touch results/$SIM_NAME.jobdone

exit 0
