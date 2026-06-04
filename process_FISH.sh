#!/usr/bin/env bash
set -Euo pipefail

# If no path provided, raise a warning message
if [[ $# -ne 1 ]]; then
  echo "Please provide the /path/to/data" >&2
  exit 1
fi

# Usage: bash process_FISH.sh /path/to/data
path_to_data="$1"

# Change this to true if you want to decode
is_MERFISH=false

run_step() {
  local step="$1"
  shift

  echo "Starting $step"
  if "$@"; then
    echo "Finished $step for $path_to_data"
  else
    echo "$step failed for $path_to_data" >&2
    exit 1
  fi
}

run_step "conversion" conda run --live-stream -n merfish3d_052926 bash -lc "qi2lab-datastore \"$path_to_data\" --codebook-path /data/smfish/codebook.csv --bit-order-path /data/smfish/bit_order.csv"
run_step "registration" conda run --live-stream -n merfish3d_052926 bash -lc "qi2lab-preprocess \"$path_to_data\""
run_step "global registration" conda run --live-stream -n merfish3d-stitcher bash -lc "qi2lab-globalregister \"$path_to_data\""
run_step "segmentation" conda run --live-stream -n merfish3d_052926 bash -lc "qi2lab-segment \"$path_to_data\""

if [[ "$is_MERFISH" == true ]]; then
    run_step "decoding" conda run --live-stream -n merfish3d_052926 bash -lc "qi2lab-decode \"$path_to_data\""
fi

echo "All done."