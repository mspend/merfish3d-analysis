# Usage: bash process_FISH.sh /path/to/data


path_to_data = $1

is_MERFISH = False


conda activate merfish3d_052926 

echo "Starting processing: $path_to_data"
# qi2lab-datastore {path_to_data} --codebook-path /data/smfish/codebook.csv --bit-order-path /data/smfish/bit_order.csv 
echo "Finished conversion to datastore for $path_to_data"

echo "Starting registration"
# qi2lab-preprocess {path_to_data}
echo "Finished registration for $path_to_data"

conda deactivate
conda activate merfish3d-stitcher 

echo "Starting global registration"
# qi2lab-globalregister {path_to_data}
echo "Finished global registration for $path_to_data"

conda deactivate
conda activate merfish3d_052926 

echo "Starting segmentation"
# qi2lab-segment {path_to_data}

if is_MERFISH:
    echo "Starting decoding"
    # qi2lab-decode {path_to_data}
    echo "Finished decoding for $path_to_data"

echo "All done."

