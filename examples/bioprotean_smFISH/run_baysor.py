from merfish3danalysis.qi2labDataStore import qi2labDataStore
from merfish3danalysis.PixelDecoder import PixelDecoder
from pathlib import Path

root_path = Path(r"/data/smFISH/02202025_Bartelle_control_smFISH_TqIB")

# initialize datastore
datastore_path = root_path / Path(r"qi2labdatastore")
datastore = qi2labDataStore(datastore_path)
merfish_bits = datastore.num_bits

# initialize decodor class
decoder = PixelDecoder(
    datastore=datastore, 
    use_mask=False, 
    merfish_bits=merfish_bits, 
    verbose=1
)



# Prep for Baysor
decoder._reformat_barcodes_for_baysor

# Run Baysor segmentation
datastore.save_spots_prepped_for_baysor
datastore._reformat_barcodes_for_baysor
datastore._run_baysor
datastore.reformat_baysor_3D_oultines
datastore.load_global_baysor_filtered_spots
datastore.reprocess_and_save_filtered_spots_with_baysor_outlines



if run_baysor:
    datastore.run_baysor()
    datastore.save_mtx()


# do I need to set the Baysor path and Baysor options?


# if not (baysor_spots_path.exists()):
#                 raise FileNotFoundError("Baysor filtered decoded spots missing.")