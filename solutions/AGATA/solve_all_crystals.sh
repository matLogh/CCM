#!/bin/bash

# Define arrays
run=7
source="66Ga"
ref_time_start=3980
ref_time_end=3990
# crystal_list=("00A" "00B" "00C" "01A" "01C" "02A" "02B" "02C" "04A" "04B" "04C" "05B" "05C" "06A" "06B" "06C" "07A" "07B" "08A" "08B" "09A" "09B" "09C" "10A" "10B" "10C" "11A" "11B" "11C" "14A" "14B" "14C")
crystal_list=("00A" "00B" "00C" "01A" "01C" "02A" "02B" "02C" "04A" "04B" "04C" "05B" "05C" "06A" "06B" "06C" "07A"       "08A" "08B" "09A" "09B" "09C" "10A" "10B" "10C" "11A" "11B" "11C" "14A" "14B" "14C")

echo "Run number: $run"
echo "Source: $source"
echo "Reference time start: $ref_time_start"
echo "Reference time end: $ref_time_end"
echo "-----------------------------------------"
# Loop over the runs
for cry in "${crystal_list[@]}"; do
    echo "Processing run number: $run and crystal: $cry"

    ./solveTimeEvo --dir /home/mbalogh/perf/timeEvo --run $run --ref_time $ref_time_start $ref_time_end --ROIsource $source --crystal $cry #--fit_peak 2751.835 2720. 2780.
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