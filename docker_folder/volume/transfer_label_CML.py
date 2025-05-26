#conda environment: cellxgene
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import argparse
import spatialdata as sd
from spatialdata_io import xenium
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import anndata as ad
import pandas as pd
#import decoupler as dc
import os                        # For filesystem operations



parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-mtc", "--metacell_path", help="Path to metacell folder that has the file name 'metacell_output.h5ad', but not the file")
parser.add_argument("-cn", "--name_column_cell_type", default="cell_type", help="Name of the collumn to pass the information, defoult: cell_type")
parser.add_argument("-sg", "--base_dir_segger", help="Path to segger folder that has the segment_res folder in it")

args = vars(parser.parse_args())
 
# Set up parameters
metacell_path = args["metacell_path"]
base_dir_segger = args["base_dir_segger"]
name_column_cell_type = args["name_column_cell_type"]

metacell_anndata = os.path.join(metacell_path, "metacell_output.h5ad")
segger_output_path = os.path.join(base_dir_segger, "segment_res/")
from pathlib import Path
segger_output_path = Path(segger_output_path)
subfolders = [f for f in segger_output_path.iterdir() if f.is_dir()]

if subfolders:
    unknown_folder = subfolders[0]  # Assuming it's the first folder
    segger_output_path = unknown_folder / "segger_adata.h5ad"


#input
# segger output anndata file
segger_output = ad.read_h5ad(segger_output_path)
#metacell
metacell = ad.read_h5ad(metacell_anndata)

#for metacell
from segger.validation.utils import annotate_query_with_reference

segger_output = annotate_query_with_reference(reference_adata=metacell, query_adata=segger_output, transfer_column=name_column_cell_type)
spatial = segger_output.obs[["cell_centroid_x", "cell_centroid_y"]]
spatial = spatial.to_numpy()
segger_output.obsm["spatial"] = spatial

#saving object
saving_transfer = os.path.join(base_dir_segger, "segger_output.h5ad")

segger_output.write_h5ad(saving_transfer)
