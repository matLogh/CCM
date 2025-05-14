#!/bin/bash

# Define arrays
run_list=(1003 1008 1010 1011 1012 1013 1014 1015 1016 1017 1 2 3 5 6 7 8 1018 1019 1020 1021 1023)
crystal_list=("00A" "00B" "00C" "01A" "01C" "02A" "02B" "02C" "04A" "04B" "04C" "05B" "05C" "06A" "06B" "06C" "07A" "07B" "08A" "08B" "09A" "09B" "09C" "10A" "10B" "10C" "11A" "11B" "11C" "14A" "14B" "14C")

# Array to store failed runs and crystals
declare -a failed_runs_and_crystals

# Loop over the runs
for run in "${run_list[@]}"; do
    echo "Processing run number: $run"

    # Run the command and capture its output and exit status
    output=$(matTimeEvo_cores --allcrys --run "$run" --Tbinning 5 2>&1)
    exit_status=$?

    # Check if the command was killed or failed
    if [[ $output == *"Killed"* || $exit_status -eq 137 ]]; then
        echo "Run $run failed: Process was killed (likely due to insufficient memory)."
        for crystal in "${crystal_list[@]}"; do
            echo "Processing crystal: $crystal for run $run"
            crystal_output=$(matTimeEvo_cores --crys "$crystal" --run "$run" --Tbinning 5 2>&1)
            crystal_exit_status=$?

            if [[ $crystal_output == *"Killed"* || $crystal_exit_status -eq 137 ]]; then
                echo "Crystal $crystal for run $run failed: Process was killed."
                failed_runs_and_crystals+=("$run:$crystal")
            elif [[ $crystal_exit_status -ne 0 ]]; then
                echo "Crystal $crystal for run $run failed: Command exited with status $crystal_exit_status."
                failed_runs_and_crystals+=("$run:$crystal")
            else
                echo "Crystal $crystal for run $run completed successfully."
            fi
        done
    elif [[ $exit_status -ne 0 ]]; then
        echo "Run $run failed: Command exited with status $exit_status."
    else
        echo "Run $run completed successfully."
    fi
done

# Print all failed runs and crystals
if [[ ${#failed_runs_and_crystals[@]} -gt 0 ]]; then
    echo "The following runs and crystals failed:"
    for entry in "${failed_runs_and_crystals[@]}"; do
        echo "Run and Crystal: $entry"
    done
else
    echo "All runs and crystals completed successfully."
fi