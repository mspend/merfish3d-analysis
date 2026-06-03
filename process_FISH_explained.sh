# This line is called the shebang. 
# Tells the operating system to run the script using bash.
# env finds the correct bash executable in your PATH.
#!/usr/bin/env bash

# Enables safer bash behavior:
# -e: exit immediately if a command fails.
# -u: error if you use an undefined variable.
# -o pipefail: if any command in a pipeline fails, the pipeline is considered failed.
set -Euo pipefail

# Usage: bash process_FISH.sh /path/to/data
# $1 = first command-line argument.
path_to_data="$1"

# Change this to true if you want to decode
is_MERFISH=false

# Checks whether the variable is empty.
# [[ ... ]] is bash's conditional syntax.
# -z means "string length is zero".
# Equivalent English: If the user did not provide a data path, print a usage message and quit.
# >&2 sends the message to stderr (error output) instead of normal output.
# exit 1 ends the script with a non-zero exit code, indicating failure.
if [[ -z "$path_to_data" ]]; then
  echo "Please provide the /path/to/data" >&2
  exit 1
fi

# Use these if you want to pass extra options
datastore_opts="--codebook-path /data/smfish/codebook.csv --bit-order-path /data/smfish/bit_order.csv"
preprocess_opts=""
global_register_opts=""
segment_opts=""
decode_opts=""

# conda run executes a command inside a conda environment.
# conda run is capturing stdout/stderr by default. 
# --no-capture-output (also spelled --live-stream) to “not capture stdout/stderr.”
# -n merfish3d_052926 specifies the environment name.
# "Run this command as if the merfish3d_052926 environment were activated."
# conda activate fails unless your shell has been initialized for conda
# -l = start a login shell.
# -c = execute the following command string.
# Start a new bash shell, load normal shell startup files, run qi2lab-preprocess /data, then exit.
echo "Starting conversion: $path_to_data"
conda run --no-capture-output -n merfish3d_052926  bash -lc "qi2lab-datastore $path_to_data" 
echo "Finished conversion for $path_to_data"

echo "Starting registration"
conda run -n --no-capture-output merfish3d_052926  bash -lc "qi2lab-preprocess $path_to_data"
echo "Finished registration for $path_to_data"


echo "Starting global registration"
conda run -n --no-capture-output merfish3d_052926 bash -lc "qi2lab-globalregister $path_to_data"
echo "Finished global registration for $path_to_data"
 

echo "Starting segmentation"
conda run -n --no-capture-output merfish3d_052926 bash -lc \ "qi2lab-segment $path_to_data"

if [[ "$is_MERFISH" == true ]]; then
    echo "Starting decoding"
    conda run -n --no-capture-output merfish3d_052926 bash -lc \ "qi2lab-decode $path_to_data"
    echo "Finished decoding for $path_to_data"
fi

echo "All done."

