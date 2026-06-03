#!/usr/bin/env bash
set -Euo pipefail

# Usage: bash process_FISH.sh /path/to/data
path_to_data="$1"

# Change this to true if you want to decode
is_MERFISH=false

# If no path provided, raise a warning message
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

echo "Starting conversion: $path_to_data"
conda run --live-stream -n merfish3d_052926  bash -lc "qi2lab-datastore $path_to_data $datastore_opts" 
echo "Finished conversion for $path_to_data"

echo "Starting registration"
conda run -n --live-stream merfish3d_052926  bash -lc "qi2lab-preprocess $path_to_data $preprocess_opts"
echo "Finished registration for $path_to_data"


echo "Starting global registration"
conda run -n --live-stream merfish3d_052926 bash -lc "qi2lab-globalregister $path_to_data $global_register_opts"
echo "Finished global registration for $path_to_data"
 

echo "Starting segmentation"
conda run -n --live-stream merfish3d_052926 bash -lc \ "qi2lab-segment $path_to_data $segment_opts"
echo "Finished segmentation for $path_to_data"

if [[ "$is_MERFISH" == true ]]; then
    echo "Starting decoding"
    conda run -n --live-stream merfish3d_052926 bash -lc \ "qi2lab-decode $path_to_data $decode_opts"
    echo "Finished decoding for $path_to_data"
fi

echo "All done."

