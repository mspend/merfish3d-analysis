#!/usr/bin/env bash
# Run sim-* pipeline for each subdirectory of a given root.
# Usage: ./run_sim_pipeline.sh /absolute/or/relative/root/path


##### NOTE:
# Everything has moved to a command line interface. I use a bash script to analyze all of the simulated data in one go, happy to share that if helpful.

# This only does one simulation "type" and you'll need to pass some parameters or edit the script for the --smFISH option







set -Euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <root_dir>" >&2
  exit 1
fi

ROOT_DIR=$1

# --- Conda activation (non-interactive safe) ---
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
else
  # Fallback: try the common install path; adjust if your conda lives elsewhere.
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  else
    echo "Could not locate conda. Ensure conda is installed and on PATH." >&2
    exit 1
  fi
fi

# Activate your env (hard-coded as requested)
conda activate merfish3d

# ---- Hard-coded options (edit if you want to pass extra flags) ----
# Example placeholders:
# CONVERT_OPTS="--some-flag foo"
# DATASTORE_OPTS="--threads 8"
# PREPROCESS_OPTS=""
# DECODE_OPTS="--radius 3.0"
# F1SCORE_OPTS="--threshold 0.5"

CONVERT_OPTS=""
DATASTORE_OPTS=""
PREPROCESS_OPTS=""
DECODE_OPTS=""
F1SCORE_OPTS=""

echo "Root: $ROOT_DIR"
echo "Finding subdirectories (one level deep)…"

# Iterate immediate subdirectories (depth 1). Change -maxdepth to recurse.
while IFS= read -r -d '' DIR; do
  # Skip if not a directory (shouldn't happen with -type d, but just in case)
  [[ -d "$DIR" ]] || continue

  echo
  echo "=== Processing: $DIR ==="

  # Define derived paths
  ACQ_DIR="$DIR/sim_acquisition"

  # Run the pipeline; keep going on errors but report them
  if ! (
    echo "+ sim-convert \"$DIR\" $CONVERT_OPTS"
    sim-convert "$DIR" $CONVERT_OPTS

    echo "+ sim-datastore \"$ACQ_DIR\" $DATASTORE_OPTS"
    sim-datastore "$ACQ_DIR" $DATASTORE_OPTS

    echo "+ sim-preprocess \"$ACQ_DIR\" $PREPROCESS_OPTS"
    sim-preprocess "$ACQ_DIR" $PREPROCESS_OPTS

    echo "+ sim-decode \"$ACQ_DIR\" $DECODE_OPTS"
    sim-decode "$ACQ_DIR" $DECODE_OPTS

    echo "+ sim-f1score \"$DIR\" $F1SCORE_OPTS"
    sim-f1score "$DIR" $F1SCORE_OPTS
  ); then
    echo "!!! Failed on: $DIR (continuing) !!!" >&2
  else
    echo "✓ Completed: $DIR"
  fi
done < <(find "$ROOT_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

echo
echo "All done."